require 'csv'
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

      desc 'alignment REFDIR SEQDIR MAPDIR', 'Align reads to reference'
      long_desc <<-DESC
Align STRT reads (*.fq.gz files at SEQDIR) to a reference (REFDIR/ref.*.ht2). The alignments will be at MAPDIR/*.bam, and per-base 5'-end counts will be at MAPDIR/*.bed.gz.
DESC

      method_option *OPT_BUFFER_SIZE
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_PARALLEL

      def alignment(refdir, seqdir, mapdir)
        
        Dir.glob("#{File.expand_path(seqdir)}/*.fq.gz").each do |fqgz|
          base = File.basename(fqgz, '.fq.gz')
          STDERR.puts "#{`date`.strip}: Align #{base}..."
          bam = "#{mapdir}/#{base}.bam"
          pipeline(
            "hisat2 --rna-strandness F --dta-cufflinks -p #{options.parallel} -x #{refdir}/ref -U #{fqgz}",
            "#{grep_command} -v -E 'NH:i:([2-9][0-9]*|1[0-9]+)'",
            "samtools sort -@ #{options.parallel} -o #{bam}")
          sh "samtools index #{bam}"
        end

        STDERR.puts "#{`date`.strip}: Count from all alignments."
        Parallel.map(Dir.glob("#{File.expand_path(mapdir)}/*.bam"),
                     in_threads: options.parallel) do |bam|
          pipeline(
            "strt count_per_base#{buffer_size_option}#{coreutils_prefix_option(options)}#{parallel_option(options)} #{bam}",
            "pigz -c > #{mapdir}/#{File.basename(bam, '.bam')}.bed.gz")
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

      # strt:call_allele

      desc 'call_allele CSV BAMDIR REFDIR', 'Call allele frequency'
      long_desc <<-DESC
Call allele frequencies of multiple samples specified in a design CSV, based on alignment files at BAMDIR, and reference sequence 'ref.fa' and the index at REFDIR.
DESC

      method_option *OPT_GENOME

      def call_allele(csv, bamdir, refdir)

        design = CSV.table(csv)
        bams = get_temporary_path('strt.call_allele', 'bams')
        fp = open(bams, 'w')
        design[:base].each {|bam| fp.puts "#{bamdir}/#{bam}.bam" }
        fp.close
        csvdir = File.dirname(csv)
        bcf = "#{csvdir}/strt-call_allele.bcf"

        pipeline("samtools mpileup -u -t AD,ADF,ADR,DP -f #{refdir}/ref.fa -b #{bams}",
                 "bcftools call --multiallelic-caller --variants-only --output-type u",
                 "bcftools filter -s LowQual -e '%QUAL<20 || MIN(FORMAT/DP)<20' --output-type b > #{bcf}")
        pipeline("bcftools view #{bcf}",
                 "table_annovar.pl - #{refdir} --buildver #{options.genome} --outfile #{csvdir}/strt-call_allele --remove --protocol refGene --operation g --nastring . --vcfinput")

      end

      # strt:check_samples

      desc 'check_samples CSV SEQDIR BAMDIR BEDDIR REFDIR', 'Check samples'
long_desc <<-DESC
Check samples in a design CSV by counting.'
DESC

      method_option *OPT_BUFFER_SIZE
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_PARALLEL

      def check_samples(csv, seqdir, bamdir, beddir, refdir)

        count_commands = ["#{cut_command(options)} -f 5",
                          "ruby -e 'n=0; while gets; n+=$_.to_i; end; puts n'"]

        samples = CSV.read(csv, {
                             headers: true,
                             converters: :numeric
                           })
        bases = samples["BASE"]

        samples["TOTAL_READS"] =
          Parallel.map(bases, in_threads: options.parallel) do |base|
          stat = CSV.table("#{seqdir}/#{base}.count.step8.csv")
          n = 0
          stat[:reads].each {|i| n += i }
          n
        end

        samples["MAPPED_READS"] =
          Parallel.map(bases, in_threads: options.parallel) do |base|
          pipeline_readline("unpigz -c #{beddir}/#{base}.bed.gz",
                            *count_commands).to_i
        end

        tmp = Array.new
        samples.each do |row|
          tmp << row["MAPPED_READS"].to_f / row["TOTAL_READS"]
        end
        samples["MAPPED_RATE"] = tmp

        samples["RIBOSOME_READS"] = 
          Parallel.map(bases, in_threads: options.parallel) do |base|
          pipeline_readline("bedtools intersect -nonamecheck -u -s -a #{beddir}/#{base}.bed.gz -b #{refdir}/ribosome.bed.gz",
                            *count_commands).to_i
        end

        samples["SPIKEIN_READS"] =
          Parallel.map(bases, in_threads: options.parallel) do |base|
          pipeline_readline("bedtools intersect -nonamecheck -u -s -a #{beddir}/#{base}.bed.gz -b #{refdir}/spikein_whole.bed.gz",
                            *count_commands).to_i
        end

        samples["SPIKEIN_5END_READS"] =
          Parallel.map(bases, in_threads: options.parallel) do |base|
          pipeline_readline("bedtools intersect -nonamecheck -u -s -a #{beddir}/#{base}.bed.gz -b #{refdir}/spikein_5end.bed.gz",
                            *count_commands).to_i
        end

        tmp = Array.new
        samples.each do |row|
          tmp << row["SPIKEIN_5END_READS"].to_f / row["SPIKEIN_READS"]
        end
        samples["SPIKEIN_5END_RATE"] = tmp

        tmp = Array.new
        samples.each do |row|
          tmp << (row["MAPPED_READS"] - row["RIBOSOME_READS"] - row["SPIKEIN_READS"]) / row["SPIKEIN_5END_READS"].to_f
        end
        samples["RELATIVE_POLYA_RNAS"] = tmp

        samples["CODING_READS"] =
          Parallel.map(bases, in_threads: options.parallel) do |base|
          pipeline_readline("bedtools intersect -nonamecheck -u -s -a #{beddir}/#{base}.bed.gz -b #{refdir}/transcriptome.coding_whole.bed.gz",
                            *count_commands).to_i
        end

        samples["CODING_5END_READS"] =
          Parallel.map(bases, in_threads: options.parallel) do |base|
          pipeline_readline("bedtools intersect -nonamecheck -u -s -a #{beddir}/#{base}.bed.gz -b #{refdir}/transcriptome.coding_5end.bed.gz",
                            *count_commands).to_i
        end

        tmp = Array.new
        samples.each do |row|
          tmp << row["CODING_5END_READS"].to_f / row["CODING_READS"]
        end
        samples["CODING_5END_RATE"] = tmp

        tmp = Array.new
        samples.each do |row|
          tmp << row["CODING_5END_READS"].to_f / row["SPIKEIN_5END_READS"]
        end
        samples["RELATIVE_MRNAS"] = tmp

        puts samples

      end
      
      # strt:count_per_base

      desc 'count_per_base BAM',  'Count reads per base'
      long_desc <<-DESC
Count reads per base , based on an alignment BAM.
DESC

      method_option *OPT_BUFFER_SIZE
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_PARALLEL

      def count_per_base(bam)

        pipeline(
          "bedtools bamtobed -i #{bam}",
          "ruby -F'\t' -anle 'puts [$F[0], $F[5]==\"+\" ? $F[1] : $F[2].to_i-1, $F[5]==\"+\" ? $F[1].to_i+1 : $F[2], $F[5]].join(\"\t\")'",
          "#{sort_command(options)} -t '\t' -k 1,1 -k 2,2n",
          "#{uniq_command(options)} -c",
          "ruby -anle 'puts ($F[1..3]+[\"#{File.basename(bam, '.bam')}\", $F[0], $F[4]]).join(\"\t\")'")
        
      end

      # strt:count_per_region

      desc 'count_per_region COUNT REGION [REGION ...]',
           'Count reads per region'
      long_desc <<-DESC
Count reads per region. Read counts, which are summerized by a BED-format COUNT, within regions, which are defined by a BED-format REGION, are summed by the region names.
DESC
      
      method_option *OPT_COREUTILS_PREFIX

      def count_per_region(count, region0, *regions0)

        pipeline(
          "bedtools intersect -nonamecheck -s -wa -wb -a #{count} -b #{([region0]+regions0).join(' ')}",
          "#{cut_command(options)} -f 5,11",
          "ruby -F'\t' -e 'n2c={}; while gets; c,n=$_.strip.split /\\t/; n2c[n]=(n2c.key?(n) ? n2c[n] : 0)+c.to_i; end; puts \"name,count\"; n2c.each {|n,c| puts \"\#{n},\#{c}\"}'")
        
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

        def pipeline_readline(*cmds)
          fp, ths = Open3.pipeline_r(*cmds)
          line = fp.gets.strip
          fp.close
          ths[-1].join
          line
        end

        def rsync_file(remote, local)
          system "rsync -a #{remote} #{local}" or exit $?.exitstatus
        end
        
      end
      
    end
  end
end

require 'bio/gadget/strt/count.rb'
require 'bio/gadget/strt/depth.rb'
