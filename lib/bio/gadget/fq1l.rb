require 'open3'
require 'bio/gadget/fq1l/bm'
require 'bio/gadget/fq1l/dmp'
require 'bio/gadget/fq1l/mt5'
require 'bio/gadget/fq1l/nr'
require 'bio/gadget/fq1l/pt3'
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

      # fq1l:exclude_redundant

      desc 'exclude_duplicate', '(Filter) Exclude duplicates in the order'

      def exclude_duplicate
        exit unless STDIN.wait
        BioGadget.nr_std()
      end

      # fq1l:match_3end

      desc 'match_3end PATTERN', '(Filter) Select sequences that match the 3\'-end with a given PATTERN'

      method_option *OPT_INVERT_MATCH
      method_option *OPT_GREP_PREFIX

      def match_3end(pattern, *files)
        if files.length == 0
          exit unless STDIN.wait
        end
        exec "#{grep_command(options)} #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t[^\\t]+#{pattern}\\t'#{files.length > 0 ? ' '+files.join(' ') : ''}"
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

      desc 'trim_3nd PATTERN', '(Filter) Trim sequences that match the 3\'-end with a given PATTERN'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_MINIMUM_LENGTH
      
      # method_option :primer,
      #               default: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG',
      #               desc: 'Primer sequence that be used for trimming',
      #               type: :string

      method_option :trimmed,
                    banner: 'FILE',
                    desc: 'FILE for trimmed reads; STDOUT if not speficied',
                    type: :string

      def trim_3end(pattern)
        exit unless STDIN.wait
        gPrefix = options.key?(:grep_prefix) ? " --grep-prefix=#{options.grep_prefix}" : ''
        fifo = get_fifo('fq1l.trim_3end', 'fq1l', false)
        begin
          tmpfile = options.key?(:trimmed) ? File.expand_path(options.trimmed) : get_temporary_path('fq1l.trim_3end', 'fq1l', false)
          begin 
            pid = Process.fork do
              BioGadget.t3("fq1l match_3end#{gPrefix} #{pattern} #{fifo}", pattern.length, options.minimum_length, tmpfile)
            end
            stats = Open3.pipeline(
              "#{tee_command(options)} #{fifo}",
              "fq1l match_3end#{gPrefix} --invert-match #{pattern}")
            Process.waitpid(pid)
            stats.each_index {|i| raise "Fail at process #{i}; #{stats[i]}" unless stats[i].success? }
          ensure
            unless options.key?(:trimmed)
              system "cat #{tmpfile}"
              File.unlink(tmpfile)
            end
          end
        ensure
          File.unlink(fifo)
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
