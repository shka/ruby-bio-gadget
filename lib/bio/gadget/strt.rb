require 'open3'
require 'parallel'

module Bio
  class Gadget
    class Strt < Bio::Gadget

      OPT_UMI_LENGTH = [ :umi_length, { :banner => 'NT',
                                        :default => 6,
                                        :desc => 'Length of UMI',
                                        :type => :numeric } ]
      
      # strt:preprocess

      desc 'preprocess BASE MAP FQGZs', 'Preprocess of raw sequences for STRT'
      long_desc <<-DESC
Process raw sequences for STRT before alignment. After demultiplexing, it performs
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

      def preprocess(base, map, fqgzs0)

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
        fqgzs0.split(/,/).each do |fqgz|
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
          
          cmds = [
            "unpigz -c #{tmpfqgzs.join(' ')}",
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
            "pigz -c > #{base}.#{well}.fq.gz"
          ]
          stats = Open3.pipeline(*cmds)
          stats.each_index do |i|
            raise "Fail at process #{i}; #{stats[i]}; #{cmds[i]}" unless stats[i].success?
          end
          
          unlink_files(tmpfqgzs)
          
        end
                             
      end
      
    end
  end
end

require 'bio/gadget/strt/bam2bed5p.rb'
require 'bio/gadget/strt/count.rb'
require 'bio/gadget/strt/depth.rb'
require 'bio/gadget/strt/prepBed.rb'
require 'bio/gadget/strt/prepGtf.rb'
require 'bio/gadget/strt/qcSmp.rb'
