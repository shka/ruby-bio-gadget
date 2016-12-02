require 'damerau-levenshtein'
require 'io/wait'
require 'open3'

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

      # fq1l:annotate_index

      desc 'annotate_index', 'Annotate sequence identifier by index sequence at the specified region'

      method_option :first_cycle,
                    default: 7,
                    desc: 'The first cycle of index',
                    type: :numeric
      
      method_option :last_cycle,
                    default: 12,
                    desc: 'The last cycle of index',
                    type: :numeric
      
      def annotate_index
        exit unless STDIN.wait
        BioGadget.i2i(options.first_cycle, options.last_cycle)
      end
      
      # fq1l:annotate_umi

      desc 'annotate_umi', 'Annotate sequence identifier by UMI sequence at the specified region'

      method_option :first_cycle,
                    default: 1,
                    desc: 'The first cycle of UMI',
                    type: :numeric
      
      method_option :last_cycle,
                    default: 6,
                    desc: 'The last cycle of UMI',
                    type: :numeric
      
      def annotate_umi
        exit unless STDIN.wait
        BioGadget.u2i(options.first_cycle, options.last_cycle)
      end
      
      # fq1l:convert
      
      desc 'convert', 'Convert fastq from 4 lines/read to 1 line/read for this utility'

      method_option *OPT_COREUTILS_PREFIX
      
      def convert
        exit unless STDIN.wait
        exec "#{options.coreutils_prefix}paste - - - -"
      end

      # fq1l:demultiplex

      desc 'demultiplex MAP BASE', 'Demultiplex based on a barcode MAP, and restore sequence files with BASE names'

      method_option :maximum_distance,
                    default: 1,
                    desc: 'Maximum distance between barcode and sequence',
                    type: :numeric

      def demultiplex(map, base)
        
        bcs = Hash.new
        open(map, 'r').each do |line|
          bc, well = line.rstrip.split(',')
          bcs[bc] = IO.popen("pigz -c > #{base}.#{well}.fq.gz", 'w:BINARY')
          bcs[bc].sync = false
        end
        na = IO.popen("pigz -c > #{base}.NA.fq.gz", 'w:BINARY')
        na.sync = false
        bcl = bcs.keys.map!{|key| key.length}.sort.uniq[0]

        dl = DamerauLevenshtein

        exit unless STDIN.wait

        STDERR.puts na.sync
        STDERR.puts STDOUT.sync
        
        fp = na
        pbc = nil
        STDIN.set_encoding('BINARY').each do |line|
          acc, seq, sep, qual = line.rstrip.split(/\t/)
          bc = acc[-bcl, bcl]
          if bc != pbc
            mindist = options.maximum_distance+1
            minbc = nil
            bcs.each_key do |key|
              dist = dl.distance(key, bc, 0, options.maximum_distance)
              if dist < mindist
                mindist = dist
                minbc = key
              end
              break if dist == 0
            end
            fp = mindist <= options.maximum_distance ? bcs[minbc] : na
            pbc = bc
          end
          fp.puts "#{acc}\n#{seq}\n#{sep}\n#{qual}"
        end

        bcs.each_value do |fp|
          fp.close
        end
        na.close
        
      end
      
      # fq1l:exclude_degenerate

      desc 'exclude_degenerate', 'Exclude degenerated reads in the order'

      def exclude_degenerate
        exit unless STDIN.wait
        BioGadget.nr_deg()
      end
      
      # fq1l:exclude_duplicate

      desc 'exclude_duplicate', 'Exclude duplicated reads in the order'

      def exclude_duplicate
        exit unless STDIN.wait
        BioGadget.nr_std()
      end

      # fq1l:match_3end

      desc 'match_3end PATTERN', 'Select sequences that match the 3\'-end with a given PATTERN'

      method_option *OPT_INVERT_MATCH
      method_option *OPT_GREP_PREFIX

      def match_3end(pattern)
        exit unless STDIN.wait
        system "#{grep_command(options)} #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t[^\\t]*#{pattern}\\t'"
        exit $?.to_i == 0 || $?.to_i == 1 ? 0 : $?.to_i
      end
      
      # fq1l:match_5end

      desc 'match_5end PATTERN', 'Select sequences that match the 5\'-end with a given PATTERN'

      method_option *OPT_INVERT_MATCH
      method_option *OPT_GREP_PREFIX

      def match_5end(pattern)
        exit unless STDIN.wait
        system "#{grep_command(options)} #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t#{pattern}'"
        exit $?.to_i == 0 || $?.to_i == 1 ? 0 : $?.to_i
      end

      # fq1l:restore

      desc 'restore', 'Convert fastq from 1 line/read to 4 lines/read'

      method_option *OPT_COREUTILS_PREFIX
      
      def restore
        exit unless STDIN.wait
        exec "#{options.coreutils_prefix}tr \"\\t\" \"\\n\""
      end
      
      # fq1l:sort

      desc 'sort', 'Sort by sequence and the quality in descending order'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL

      def sort
        exit unless STDIN.wait
        exec "#{sort_command(options)} -t '\t' -r -k2,4"
      end

      # fq1l:sort_index

      desc 'sort_index', 'Sort by index'
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL
      
      def sort_index
        exit unless STDIN.wait
        exec "#{sort_command(options)} -k2"
      end
      
      # fq1l:thin_out

      desc 'thin_out DRAW SKIP', 'Thin out the sequences'

      def to(draw, skip)
        exit unless STDIN.wait
        BioGadget.to(draw.to_i, skip.to_i)
      end
      
      # fq1l:trim_3end

      desc 'trim_3end SEQUENCE', 'Trim 3\'-end that match with a given SEQUENCE'

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
            system "#{cat_command(options)} #{tmpfile}" unless options.key?(:trimmed)
          end
        ensure
          File.unlink(fifo) if File.exist?(fifo)
          File.unlink(tmpfile) if File.exist?(tmpfile) && !options.key?(:trimmed)
        end
      end

      # fq1l:trim_3end_length

      desc 'trim_3end_length', 'Trim 3\'-end by a specific length'

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

      desc 'trim_3end_primer', 'Trim 3\'-end that match with a given primer'

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

      desc 'trim_3end_quality', 'Trim 3\'-end from a low quality base'

      method_option *OPT_MINIMUM_LENGTH
      
      method_option :low_qualities,
                    banner: 'CHARACTERS',
                    default: '!"#',
                    desc: 'Low quality characters',
                    type: :string
      
      def trim_3end_quality
        BioGadget.t3q(options.low_qualities, options.minimum_length)
      end
      
      # fq1l:trim_5end

      desc 'trim_5end PATTERN', 'Trim 5\'-end that match with a given PATTERN'

      method_option :minimum_length,
                    banner: 'NT',
                    default: 26,
                    desc: 'Minimum length after trimming',
                    type: :numeric
                      
      
      def trim_5end(pattern)
        exit unless STDIN.wait
        BioGadget.t5(pattern, options.minimum_length)
      end

    end
  end
end
