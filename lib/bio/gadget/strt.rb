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
Process raw sequences for STRT before alignment. It performs
(i) exclusion of redundant reads,
(ii) exclusion of noncanonical reads, which does not begin with template switching primer,
(iii) trimming from low-quality base,
(iv) trimming from sequence similar to HiSeq universal primer,
(v) trimming of the template switching primer,
and (vi) demultiplexing based on barcode.

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

        tmps = Array.new

        bcs = Hash.new
        open(map, 'r').each do |line|
          bc, well = line.rstrip.split(',')
          bcs[bc] = well
        end
        
        bcl = bcs.keys.map!{|key| key.length}.sort.uniq[0]

        tso_pattern = '.'*options.umi_length + '.'*bcl + 'GG'

        begin

          tmps << fifo1 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count1 = get_temporary_path('strt.preprocess', 'txt')
          pid1 = spawn_wcl(options, fifo1, count1)
          tmps << fifo2 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count2 = get_temporary_path('strt.preprocess', 'txt')
          pid2 = spawn_wcl(options, fifo2, count2)
          tmps << fifo3 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count3 = get_temporary_path('strt.preprocess', 'txt')
          pid3 = spawn_wcl(options, fifo3, count3)
          tmps << fifo4 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count4 = get_temporary_path('strt.preprocess', 'txt')
          pid4 = spawn_wcl(options, fifo4, count4)
          tmps << fifo5 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count5 = get_temporary_path('strt.preprocess', 'txt')
          pid5 = spawn_wcl(options, fifo5, count5)
          tmps << fifo6 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count6 = get_temporary_path('strt.preprocess', 'txt')
          pid6 = spawn_wcl(options, fifo6, count6)
          tmps << fifo7 = get_fifo('strt.preprocess', 'fq1l')
          tmps << count7 = get_temporary_path('strt.preprocess', 'txt')
          pid7 = spawn_wcl(options, fifo7, count7)
          
          cmds = [
            "gunzip -c #{fqgzs0.split(/,/).join(' ')}",
            "fq1l convert#{coreutils_prefix_option(options)}",
            "#{tee_command(options)} #{fifo1}",
            "fq1l match_5end#{grep_prefix_option(options)} #{tso_pattern}",
            "#{tee_command(options)} #{fifo2}",
            "fq1l sort#{coreutils_prefix_option(options)}#{parallel_option(options)} --buffer-size=#{(options.maximum_memory/2).to_i}%",
            "fq1l exclude_duplicate",
            "#{tee_command(options)} #{fifo3}",
            "fq1l trim_3end_quality",
            "#{tee_command(options)} #{fifo4}",
            "fq1l trim_3end_primer#{coreutils_prefix_option(options)}#{grep_prefix_option(options)}#{parallel_option(options)}",
            "#{tee_command(options)} #{fifo5}",
            "fq1l sort#{coreutils_prefix_option(options)}#{parallel_option(options)} --buffer-size=#{(options.maximum_memory/2).to_i}%",
            "fq1l exclude_degenerate",
            "#{tee_command(options)} #{fifo6}",
            "fq1l annotate_index --first-cycle=#{options.umi_length+1} --last-cycle=#{options.umi_length+bcl}",
            "fq1l annotate_umi --first-cycle=1 --last-cycle=#{options.umi_length}",
            "fq1l trim_5end --minimum-length=#{options.minimum_length} #{tso_pattern}+",
            "#{tee_command(options)} #{fifo7}",
            "fq1l sort_index#{coreutils_prefix_option(options)}#{parallel_option(options)} --buffer-size=#{(options.maximum_memory/2).to_i}%",
            "fq1l demultiplex #{base} #{map}"
          ]
          cmds.insert(2, "#{head_command(options)} -n #{options.reads}") unless options.reads.nil?
          stats = Open3.pipeline(*cmds)
          stats.each_index do |i|
            raise "Fail at process #{i}; #{stats[i]}; #{cmds[i]}" unless stats[i].success? || (stats[i].signaled? && stats[i].termsig == 13)
          end
          
          Process.waitpid(pid1)
          reads1 = `#{cat_command(options)} #{count1}`.rstrip.to_i
          STDERR.puts "#{reads1} raw reads."
          
          Process.waitpid(pid2)
          reads2 = `#{cat_command(options)} #{count2}`.rstrip.to_i
          STDERR.puts "#{reads2} canonical reads."
          
          Process.waitpid(pid3)
          reads3 = `#{cat_command(options)} #{count3}`.rstrip.to_i
          STDOUT.puts "#{reads3} nonredundant reads, before trimming."
          
          Process.waitpid(pid4)
          reads4 = `#{cat_command(options)} #{count4}`.rstrip.to_i
          STDOUT.puts "#{reads4} qualified reads after quality trimming."
          
          Process.waitpid(pid5)
          reads5 = `#{cat_command(options)} #{count5}`.rstrip.to_i
          STDOUT.puts "#{reads5} qualified reads after primer trimming."

          Process.waitpid(pid6)
          reads6 = `#{cat_command(options)} #{count6}`.rstrip.to_i
          STDOUT.puts "#{reads6} nonredundant reads, after quality/primer trimming."
          Process.waitpid(pid7)
          reads7 = `#{cat_command(options)} #{count7}`.rstrip.to_i
          STDOUT.puts "#{reads7} qualified reads, after UMI/barcode trimming."

        ensure
          unlink_files(tmps)
        end
        
      end
      
      no_commands do

        def spawn_wcl(options, infile, outfile)
          Kernel.spawn("#{wc_command(options)} -l #{infile} > #{outfile}")
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
