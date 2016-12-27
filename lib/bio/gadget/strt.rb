require 'fileutils'
require 'open3'
require 'parallel'

require 'bio/gadget/strt/prepare_transcriptome.rb'

module Bio
  class Gadget
    class Strt < Bio::Gadget

      OPT_GENOME = [ :genome, { :default => 'hg38',
                                :desc => 'Genome assembly' } ]
      
      OPT_UMI_LENGTH = [ :umi_length, { :banner => 'NT',
                                        :default => 6,
                                        :desc => 'Length of UMI',
                                        :type => :numeric } ]

      # strt:alignment

      desc 'alignment INDEX FQGZDIR BAMDIR', 'Align reads to reference'
      long_desc <<-DESC
Align STRT reads (*.fq.gz files at FQGZDIR) to a reference (INDEX/ref.*.ht2). The results will be at BAMDIR/*.bam.
DESC

      method_option *OPT_PARALLEL

      def alignment(index0, indir0, outdir0)
        
        index = File.expand_path(index0)
        outdir = File.expand_path(outdir0)
        
        Dir.glob("#{File.expand_path(indir0)}/*.fq.gz").each do |fqgz|
          
          STDERR.puts "#{`date`.strip}: Align #{fqgz}..."

          tmpbam = get_temporary_path('strt.alignment', 'bam')
          bam = "#{outdir}/#{File.basename(fqgz, '.fq.gz')}.bam"
          
          pipeline("hisat2 --rna-strandness F --dta-cufflinks -p #{options.parallel} -x #{index}/ref -U #{fqgz}",
                   "samtools view -S -b -@ #{options.parallel} - > #{tmpbam}")
          system "samtools sort -f -@ #{options.parallel} #{tmpbam} #{bam}" or exit $?.exitstatus
          system "samtools index #{bam}" or exit $?.exitstatus
          
        end
        
      end
      
      # strt:build_index

      desc 'build_index DIR', 'Build index for alignment'
      long_desc <<-DESC
Build index for alignment of STRT reads, from the speficied GENOME, TRANSCRIPTOME and VARIATION, at DIR.
DESC

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GENOME
      method_option *OPT_GREP_PREFIX
      method_option *OPT_PARALLEL
      
      def build_index(dir0)

        dir = File.expand_path(dir0)
        FileUtils.mkdir_p(dir)

        STDERR.puts "#{`date`.strip}: Preparing data files..."
        
        Parallel.map(
          ["strt prepare_genome#{coreutils_prefix_option(options)}#{genome_option(options)} #{dir}",
           "strt prepare_variation#{coreutils_prefix_option(options)}#{genome_option(options)} --download=only #{dir}",
           "strt prepare_transcriptome #{options.genome}#{coreutils_prefix_option(options)}#{grep_prefix_option(options)} --download=only #{dir}",
           "strt prepare_spikein#{coreutils_prefix_option(options)} #{dir}",
           "strt prepare_ribosome#{coreutils_prefix_option(options)}#{genome_option(options)} #{dir}"], in_threads: options.parallel) do |cmd|
          system cmd or exit $?.exitstatus
        end

        system "unpigz -c #{dir}/genome.fa.gz #{dir}/spikein.fa.gz #{dir}/ribosome.fa.gz > #{dir}/ref.fa"
        system "samtools faidx #{dir}/ref.fa"

        Parallel.map(["strt prepare_transcriptome #{options.genome}#{coreutils_prefix_option(options)}#{grep_prefix_option(options)} --download=no #{dir}",
                      "strt prepare_variation#{coreutils_prefix_option(options)}#{genome_option(options)} --download=no #{dir}"], in_threads: options.parallel) do |cmd|
          STDERR.puts cmd
          system cmd or exit $?.exitstatus
        end
        
        STDERR.puts "#{`date`.strip}: Building index..."
        
        system "hisat2-build -f -p #{options.parallel} --snp #{dir}/variation.snp --haplotype #{dir}/variation.haplotype --ss #{dir}/transcriptome.splice_sites --exon #{dir}/transcriptome.exons #{dir}/ref.fa #{dir}/ref"
        
      end

      # strt:prepare_genome

      desc 'prepare_genome DIR', 'Prepare genome data'
      long_desc <<-DESC
Prepare data files of the specified GENOME at DIR.
DESC

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_DOWNLOAD
      method_option *OPT_GENOME

      def prepare_genome(dir0)

        dir = File.expand_path(dir0)
        tgz = "#{dir}/#{options.genome}.chromFa.tar.gz"
        ucsc = "rsync://hgdownload.cse.ucsc.edu/goldenPath/#{options.genome}/bigZips"

        if options.download != 'no'
          if options.genome == 'hg38'
            rsync_file("#{ucsc}/#{options.genome}.chromFa.tar.gz", tgz)
          else
            rsync_file("#{ucsc}/chromFa.tar.gz", tgz)
          end
        end
        pipeline("unpigz -c #{tgz}",
                 "#{options.coreutils_prefix}tar -xOf - --exclude \"*_*\"",
                 "gawk 'BEGIN{p=\"\"} /^>/{ print p $1 } !/^>/{ printf $1; p=\"\\n\" } END{ print }'",
                 "#{fold_command(options)} -w 50",
                 "pigz -c > #{dir}/genome.fa.gz"
                ) if options.download != 'only'
        
      end
      
      # strt:prepare_reads

      desc 'prepare_reads BASE MAP FQGZ ...', 'Prepare STRT reads'
      long_desc <<-DESC
Prepare STRT reads from raw sequence files before alignment. After demultiplexing, it performs
(i) exclusion of redundant reads,
(ii) exclusion of noncanonical reads, which does not begin with template switching primer,
(iii) trimming from low-quality base,
(iv) trimming from sequence similar to HiSeq universal primer,
and (v) trimming of the template switching primer.

Mandatory paramters are (1) BASE; basename for demulplexed and gzipped fastq files, (2) MAP; filename of comma-separated table between barcode and well, and (3) FQGZs; comma-separated filenames of raw sequences; each file is gzipped fastq. When MAP contains 'CAAAGT,A2' and BASE is '~/test', reads having CAAAGT-like barcode are in '~/test.A2.fq.gz' file after the preprocesses.
DESC

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_PARALLEL
      method_option *OPT_UMI_LENGTH

      method_option :maximum_memory,
                    default: 50,
                    desc: 'Maximum memory usage in percent rate',
                    type: :numeric
      
      method_option :minimum_length,
                    banner: 'NT',
                    default: 24,
                    desc: 'Minimum length after the preprocess',
                    type: :numeric
      
      method_option :reads,
                    desc: 'Number of raw reads for the preprocess',
                    type: :numeric

      # method_option :maximum_distance,
      #               default: 1,
      #               desc: 'Maximum distance between barcode and sequence',
      #               type: :numeric

      def prepare_reads(base, map, fqgz0, *fqgzs0)

        fqgzs = [fqgz0] + fqgzs0

        bcs = Hash.new
        open(map, 'r').each do |line|
          bc, well = line.rstrip.split(',')
          bcs[bc] = well
        end
        
        bcl = bcs.keys.map!{|key| key.length}.sort.uniq[0]

        tso_pattern = '.'*options.umi_length + '.'*bcl + 'GG'

        #
        
        STDERR.puts "#{`date`.strip}: Demultiplexing each raw sequence files..."
        
        fqgz2csv0 = Hash.new
        fqgz2csv1 = Hash.new
        fqgz2base = Hash.new
        fqgzs.each do |fqgz|
          fqgz2csv0[fqgz] = get_temporary_path('strt.preprocess', 'csv', false)
          fqgz2csv1[fqgz] = get_temporary_path('strt.preprocess', 'csv', false)
          fqgz2base[fqgz] = get_temporary_path('strt.preprocess', 'base', false)
        end

        Parallel.map(fqgz2csv0.keys, in_processes: options.parallel) do |fqgz|
          cmds = [
            "unpigz -c #{fqgz}",
            "#{fq1l_convert_command(options)}",
            "#{fq1l_count_command(options)} #{fqgz2csv0[fqgz]}",
            "fq1l match_5end#{grep_prefix_option(options)} #{tso_pattern}",
            "#{fq1l_count_command(options)} #{fqgz2csv1[fqgz]}",
            "fq1l annotate_index --first-cycle=#{options.umi_length+1} --last-cycle=#{options.umi_length+bcl}",
            "fq1l annotate_umi --first-cycle=1 --last-cycle=#{options.umi_length}",
            "fq1l sort_index#{coreutils_prefix_option(options)}#{parallel_option(options)} --buffer-size=#{(options.maximum_memory/(fqgz2csv0.keys.size+1)).to_i}%",
            "fq1l demultiplex #{fqgz2base[fqgz]} #{map}"
          ]
          cmds.insert(2, "#{head_command(options)} -n #{options.reads}") unless options.reads.nil?
          stats = Open3.pipeline(*cmds)
          stats.each_index do |i|
            raise "Fail at process #{i}; #{stats[i]}; #{cmds[i]}" unless stats[i].success? || (stats[i].signaled? && stats[i].termsig == 13)
          end
        end

        system "fq1l sum_counts #{fqgz2csv0.values.join(' ')} > #{base}.count.step1.csv"
        unlink_files(fqgz2csv0.values)
        
        system "fq1l sum_counts #{fqgz2csv1.values.join(' ')} > #{base}.count.step2.csv"
        unlink_files(fqgz2csv1.values)

        #
        
        (bcs.values + ['NA']).each do |well|

          STDERR.puts "#{`date`.strip}: Finishing well #{well}..."
          
          tmpfqgzs = fqgz2base.values.map {|base| "#{base}.#{well}.fq.gz"}
          csvs = Array.new(6) {|i| "#{base}.#{well}.count.step#{i+3}.csv"}
          
          pipeline("unpigz -c #{tmpfqgzs.join(' ')}",
                   "#{fq1l_convert_command(options)}",
                   "#{fq1l_count_command(options)} #{csvs[0]}",
                   "#{fq1l_sort_command(options)} --buffer-size=#{(options.maximum_memory/2).to_i}%",
                   "fq1l exclude_duplicate",
                   "#{fq1l_count_command(options)} #{csvs[1]}",
                   "fq1l trim_3end_quality",
                   "#{fq1l_count_command(options)} #{csvs[2]}",
                   "fq1l trim_3end_primer#{coreutils_prefix_option(options)}#{grep_prefix_option(options)}#{parallel_option(options)}",
                   "#{fq1l_count_command(options)} #{csvs[3]}",
                   "#{fq1l_sort_command(options)} --buffer-size=#{(options.maximum_memory/2).to_i}%",
                   "fq1l exclude_degenerate",
                   "#{fq1l_count_command(options)} #{csvs[4]}",
                   "fq1l trim_5end --minimum-length=#{options.minimum_length} #{tso_pattern}+",
                   "#{fq1l_count_command(options)} #{csvs[5]}",
                   "fq1l restore#{coreutils_prefix_option(options)}",
                   "pigz -c > #{base}.#{well}.fq.gz")
          
          unlink_files(tmpfqgzs)
          
        end
                             
      end

      # strt:prepare_ribosome

      desc 'prepare_ribosome DIR', 'Prepare ribosome data'
      long_desc <<-DESC
Prepare ribosome data files for the specified GENOME at DIR.
DESC

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_DOWNLOAD
      method_option *OPT_GENOME

      def prepare_ribosome(dir0)

        dir = File.expand_path(dir0)

        if options.genome[0..1] == 'hg'
          if options.download != 'no'
            download_file("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=555853&strand=1&rettype=fasta&retmode=text", "#{dir}/U13369.fa")
            system "pigz #{dir}/U13369.fa" or exit $?.exitstatus
          end
          pipeline("unpigz -c #{dir}/U13369.fa.gz",
                   "gawk '/^>/{ print \">RIBO_U13369.1\" } !/^>/{ printf $1 } END{ print }'",
                   "#{fold_command(options)} -w 50",
                   "pigz -c > #{dir}/ribosome.fa.gz"
                  ) if options.download != 'only'
        else
          pipeline("echo", "pigz -c > #{dir}/ribosome.fa.gz")
        end
        
      end

      # strt:prepare_spikein

      desc 'prepare_spikein DIR', 'Prepare spikein data'
      long_desc <<-DESC
Prepare spikein data files at DIR.
DESC

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_DOWNLOAD

      def prepare_spikein(dir0)
        
        dir = File.expand_path(dir0)
        zip = "#{dir}/ERCC92.zip"
        
        download_file("https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip", zip) if options.download != 'no'
        pipeline("unzip -cq #{zip} ERCC92.fa",
                 "gawk 'BEGIN{p=\"\"} /^>/{ print p \">RNA_SPIKE_\" substr($1, 2); printf(\"AATTC\" ($1 == \">ERCC-00130\" ? \"GAGCTC\" : \"\") ) } /^[ACGT]/{ printf $1; p=\"\\n\" } END{ print }'",
                 "#{fold_command(options)} -w 50",
                 "pigz -c > #{dir}/spikein.fa.gz") if options.download != 'only'
        
      end
      
      # strt:prepare_transcriptome

      register(Bio::Gadget::StrtPrepareTranscriptome,
               'prepare_transcriptome',
               'prepare_transcriptome GENOME',
               'Prepare transcriptome data')
      
      # strt:prepare_variation

      desc 'prepare_variation DIR', 'Prepare variation data'
      long_desc <<-DESC
Prepare genome variation data files for the specified GENOME dir based on common variations in dbSNP BUILD, at DIR.
DESC
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_DOWNLOAD
      method_option *OPT_GENOME

      method_option :dbsnp,
                    banner: 'BUILD',
                    default: 146,
                    desc: 'Build number of dbSNP',
                    type: :numeric
      
      def prepare_variation(dir0)
        
        dir = File.expand_path(dir0)
        snp = "#{dir}/#{options.genome}.snp#{options.dbsnp}Common.txt.gz"

        rsync_file("rsync://hgdownload.soe.ucsc.edu/goldenPath/#{options.genome}/database/snp#{options.dbsnp}Common.txt.gz", snp) if options.download != 'no'
        pipeline("unpigz -c #{dir}/genome.fa.gz",
                 "hisat2_extract_snps_haplotypes_UCSC.py - #{snp} #{dir}/variation"
                ) if options.download != 'only'
        
      end

      #

      no_commands do

        def download_option(options)
          " --download=#{options.download}"
        end

        def genome_option(options)
          " --genome=#{options.genome}"
        end

        def rsync_file(remote, local)
          system "rsync -a #{remote} #{local}" or exit $?.exitstatus
        end
        
      end
      
    end
  end
end

require 'bio/gadget/strt/bam2bed5p.rb'
require 'bio/gadget/strt/count.rb'
require 'bio/gadget/strt/depth.rb'
require 'bio/gadget/strt/qcSmp.rb'
