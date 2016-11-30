require 'open3'
require 'bio/gadget/fq1l/bm'
require 'bio/gadget/fq1l/dmp'
require 'bio/gadget/fq1l/mt5'
require 'bio/gadget/fq1l/rst'
require 'bio/gadget/fq1l/to'

module Bio
  class Gadget
    class Fq1l < Bio::Gadget
      
      OPT_INVERT_MATCH = [
        :invert_match, {
          :desc => 'The sense of matching',
          :type => :boolean
        }
      ]

      OPT_MINIMUM_LENGTH = [
        :minimum_length, {
          :banner => 'NT',
          :default => 40,
          :desc => 'Minimum length after trimming',
          :type => :numeric
        }
      ]

      # fq1l:convert
      
      desc 'convert', '(Filter) Convert fastq from 4 lines/read to 1 line/read'

      method_option *OPT_COREUTILS_PREFIX
      
      def convert
        exit unless STDIN.wait
        exec "#{options.coreutils_prefix}paste - - - -"
      end

      # fq1l:exclude_degenerate

      desc 'exclude_degenerate', '(Filter) Exclude degenerated reads in the order'

      def exclude_degenerate
        exit unless STDIN.wait
        BioGadget.nr_deg()
      end
      
      # fq1l:exclude_duplicate

      desc 'exclude_duplicate', '(Filter) Exclude duplicated reads in the order'

      def exclude_duplicate
        exit unless STDIN.wait
        BioGadget.nr_std()
      end

      # fq1l:index_to_id

      desc 'index_to_id FIRST LAST', '(Filter) Append index to sequence identifier'
      def index_to_id(first, last)
        exit unless STDIN.wait
        BioGadget.i2i(first.to_i, last.to_i)
      end
      
      # fq1l:match_3end

      desc 'match_3end PATTERN', '(Filter) Select sequences that match the 3\'-end with a given PATTERN'

      method_option *OPT_INVERT_MATCH
      method_option *OPT_GREP_PREFIX

      def match_3end(pattern)
        exec "#{grep_command(options)} #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t[^\\t]*#{pattern}\\t'"
      end
      
      # fq1l:match_5end

      desc 'match_5end PATTERN', '(Filter) Select sequences that match the 5\'-end with a given PATTERN'

      method_option *OPT_INVERT_MATCH
      method_option *OPT_GREP_PREFIX

      def match_5end(pattern)
        exit unless STDIN.wait
        exec "#{grep_command(options)} #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t#{pattern}'"
      end
      
      # fq1l:sort

      desc 'sort', '(Filter) Sort by sequences and the qualities in descending order'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL

      def sort
        exit unless STDIN.wait
        exec "#{sort_command(options)} -t '\t' -r -k2,4"
      end

      # fq1l:trim_3end

      desc 'trim_3end SEQUENCE', '(Filter) Trim 3\'-end that match with a given SEQUENCE'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_MINIMUM_LENGTH
      
      method_option :trimmed,
                    banner: 'FILE',
                    desc: 'FILE for trimmed reads; STDOUT if not speficied',
                    type: :string

      def trim_3end(sequence)
        exit unless STDIN.wait
        gPrefix = options.key?(:grep_prefix) ? " --grep-prefix=#{options.grep_prefix}" : ''
        fifo = get_fifo('fq1l.trim_3end', 'fq1l', false)
        begin
          tmpfile = options.key?(:trimmed) ? File.expand_path(options.trimmed) : get_temporary_path('fq1l.trim_3end', 'fq1l', false)
          begin 
            pid = Process.fork do
              BioGadget.t3("fq1l match_3end#{gPrefix} #{sequence} < #{fifo}", sequence.length, options.minimum_length, tmpfile)
            end
            commands = ["#{tee_command(options)} #{fifo}",
                        "fq1l match_3end#{gPrefix} #{sequence} --invert-match"]
            stats = Open3.pipeline(*commands)
            Process.waitpid(pid)
            stats.each_index do |i|
              raise "Fail at process #{i}; #{stats[i]}; #{commands[i]}" unless stats[i].success?
            end
          ensure
            unless options.key?(:trimmed)
              system "#{cat_command(options)} #{tmpfile}"
              File.unlink(tmpfile) if File.exist?(tmpfile)
            end
          end
        ensure
          File.unlink(fifo) if File.exist?(fifo)
        end
      end

      # fq1l:trim_3end_length

      desc 'trim_3end_length', '(Filter) Trim 3\'-end by a specific length'

      method_option *OPT_MINIMUM_LENGTH

      method_option :trimming_length,
                    default: 1,
                    desc: 'Length of the trimming',
                    type: :numeric

      def trim_3end_length
        exit unless STDIN.wait
        BioGadget.t3(nil, options.trimming_length, options.minimum_length, nil)
      end
      
      # fq1l:trim_3end_primer

      desc 'trim_3end_primer', '(Filter) Trim 3\'-end that match with a given primer'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_MINIMUM_LENGTH
      method_option *OPT_PARALLEL

      method_option :primer,
                    default: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG',
                    desc: 'Primer sequence that be used for trimming',
                    type: :string

      def trim_3end_primer

        opt_minimum_length = "--minimum-length=#{options.minimum_length}"
        primer = options.primer
        
        fragments = Hash.new
        max = primer.length-1
        tmp = Hash.new
        for i in 0..max do
          for j in i..max do
            fragment = primer[i..j]
            unless tmp.key?(fragment)
              l = fragment.length
              fragments[l] = Array.new unless fragments.key?(l)
              fragments[l] << fragment
              tmp[fragment] = true
            end
          end
        end
        
        exit unless STDIN.wait

        tmpfiles = Array.new
        commands = Array.new
        pids = Array.new
        begin 
          fragments.keys.sort.reverse.each do |length|
            if 4**length == fragments[length].size
              commands << "fq1l trim_3end_length --trimming-length=#{length} #{opt_minimum_length}"
              break
            else
              fragments[length].sort.reverse.each do |fragment|
                tmpfiles << tmpfile = get_temporary_path('fq1l.trim_3end_primer', 'fq1l', false)
                commands << "fq1l trim_3end#{' --coreutils-prefix='+options.coreutils_prefix if options.key?(:coreutils_prefix)}#{' --grep-prefix='+options.grep_prefix if options.key?(:grep_prefix)} #{opt_minimum_length} --trimmed=#{tmpfile} #{fragment}"
              end
            end
          end
          stats = pipeline(options.parallel*2, *commands)
          stats.each_index do |i|
            raise "Fail at process #{i}; #{stats[i]}; #{commands[i]}" unless stats[i].success?
          end
          system "#{cat_command(options)} #{tmpfiles.join(' ')}"
        ensure
          tmpfiles.each do |tmpfile|
            File.unlink(tmpfile) if FileTest.exist?(tmpfile)
          end
        end
        
      end

      # fq1l:trim_3end_quality

      desc 'trim_3end_quality', '(Filter) Trim 3\'-end from a low quality base'

      method_option *OPT_MINIMUM_LENGTH
      
      method_option :low_qualities,
                    banner: 'CHARACTERS',
                    default: '!"#',
                    desc: 'Low quality characters',
                    type: :string
      
      def trim_3end_quality
        BioGadget.t3q(options.low_qualities, options.minimum_length)
      end
      
      # fq1l:umi_to_id

      desc 'umi_to_id FIRST LAST', '(Filter) Append UMI to sequence identifier'
      def umi_to_id(first, last)
        exit unless STDIN.wait
        BioGadget.u2i(first.to_i, last.to_i)
      end
      
      #

      no_commands do
        
        def read_barcodes(map)
          bcs = Hash.new
          open(map, 'r').each do |line|
            bc, well = line.rstrip.split(',')
            bcs[bc] = well
          end
          return bcs
        end

      end

    end
  end
end
