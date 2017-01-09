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

      # fq1l:count

      desc 'count [CSV]', 'Count sequences by the length'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_PARALLEL
      
      def count(csv = nil)
        exit unless STDIN.wait
        if csv.nil?
          puts "length,reads"
          pipeline("#{cut_command} -f 2",
                   "ruby -nle 'puts $_.length'",
                   "#{sort_command} -n",
                   "#{uniq_command(options)} -c",
                   "ruby -anle 'puts $F.reverse.join(\",\")'")
        else
          fifo = get_fifo('fq1l.count', 'fq1l')
          pid = Kernel.spawn("fq1l count#{coreutils_prefix_option} < #{fifo} > #{csv}")
          system "#{tee_command(options)} #{fifo}"
          Process.waitpid(pid)
        end
      end
      
      # fq1l:demultiplex

      desc 'demultiplex BASE MAP', 'Demultiplex based on a barcode MAP, and restore sequence files with BASE names'

      method_option :maximum_distance,
                    default: 1,
                    desc: 'Maximum distance between barcode and sequence',
                    type: :numeric

      def demultiplex(base, map)
        
        dl = DamerauLevenshtein

        exit unless STDIN.wait

        bc2fq = Hash.new
        open(map, 'r').each do |line|
          bc, well = line.rstrip.split(',')
          bc2fq[bc] = fq = "#{base}.#{well}.fq"
          File.unlink(fq) if File.exist?(fq)
        end
        na = "#{base}.NA.fq"
        File.unlink(na) if File.exist?(na)
        
        bcl = bc2fq.keys.map!{|key| key.length}.sort.uniq[0]

        fp = nil
        pbc = nil
        STDIN.set_encoding('BINARY').each do |line|
          acc, seq, sep, qual = line.rstrip.split(/\t/)
          bc = acc[-bcl, bcl]
          if bc != pbc
            mindist = options.maximum_distance+1
            minbc = nil
            bc2fq.each_key do |key|
              dist = dl.distance(key, bc, 0, options.maximum_distance)
              if dist < mindist
                mindist = dist
                minbc = key
              end
              break if dist == 0
            end
            fp.close unless fp.nil?
            fp = open(mindist <= options.maximum_distance ? bc2fq[minbc] : na, 'a')
            pbc = bc
          end
          fp.puts "#{acc}\n#{seq}\n#{sep}\n#{qual}"
        end
        fp.close unless fp.nil?
        
        bc2fq.each_value {|fq| system "pigz #{fq}" if File.exist?(fq) }
        system "pigz #{na}" if File.exist?(na)
        
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
        # PCRE was faster than BRE and ERE in GNU grep 2.25
        system "#{grep_command}#{options.invert_match ? ' -v' : ''} -P -e '^[^\\t]+\\t[^\\t]*#{pattern}\\t'"
        exit $?.to_i == 0 || $?.to_i == 1 ? 0 : $?.to_i
      end
      
      # fq1l:match_5end

      desc 'match_5end PATTERN', 'Select sequences that match the 5\'-end with a given PATTERN'

      method_option *OPT_INVERT_MATCH
      method_option *OPT_GREP_PREFIX

      def match_5end(pattern)
        exit unless STDIN.wait
        # PCRE was faster than BRE and ERE in GNU grep 2.25
        system "#{grep_command}#{options.invert_match ? ' -v' : ''} -P -e '^[^\\t]+\\t#{pattern}'"
        exit $?.to_i == 0 || $?.to_i == 1 ? 0 : $?.to_i
      end

      # fq1l:restore

      desc 'restore', 'Convert fastq from 1 line/read to 4 lines/read'

      method_option *OPT_COREUTILS_PREFIX
      
      def restore
        exit unless STDIN.wait
        exec "#{options.coreutils_prefix}tr \"\\t\" \"\\n\""
      end

      # fq1l:slice

      desc 'slice Nth SLICE', 'Slice the sequences'

      def slice(nth, slice)
        exit unless STDIN.wait
        BioGadget.slice(nth.to_i, slice.to_i)
      end
      
      # fq1l:sort

      desc 'sort [FQ1Ls]', 'Sort by sequence and the quality in descending order'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL

      def sort(*fq1ls)
        if fq1ls.size == 0
          exit unless STDIN.wait
          exec "#{sort_command} -t '\t' -r -k2,4"
        else
          exec "#{sort_command} -t '\t' -r -k2,4 -m #{fq1ls.join(' ')}"
        end
      end

      # fq1l:sort_index

      desc 'sort_index', 'Sort by index'
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL
      
      def sort_index
        exit unless STDIN.wait
        exec "#{sort_command} -k2"
      end

      # fq1l:sum_counts

      desc 'sum_counts CSV ...', 'Sum counts of sequences by the length'

      def sum_counts(*csvs)
        length2count = Hash.new
        csvs.each do |csv|
          open(csv).each do |line|
            l, c = line.rstrip.split(/,/)
            next if l == 'length'
            length = l.to_i
            length2count[length] = 0 unless length2count.key?(length)
            length2count[length] += c.to_i
          end
        end
        puts "length,count"
        length2count.keys.sort.each do |length|
          puts "#{length},#{length2count[length]}"
        end
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
        # exit unless STDIN.wait
        gPrefix = options.key?(:grep_prefix) ? " --grep-prefix=#{options.grep_prefix}" : ''
        fifo = get_fifo('fq1l.trim_3end', 'fq1l', false)
        begin
          tmpfile = options.key?(:trimmed) ? File.expand_path(options.trimmed) : get_temporary_path('fq1l.trim_3end', 'fq1l', false)
          begin 
            pid = Process.fork do
              BioGadget.t3("fq1l match_3end#{gPrefix} #{sequence} < #{fifo}", sequence.length, options.minimum_length, tmpfile)
            end
            pipeline("#{tee_command(options)} #{fifo}",
                     "fq1l match_3end#{gPrefix} #{sequence} --invert-match")
          ensure
            system "#{cat_command} #{tmpfile}" unless options.key?(:trimmed)
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

      method_option :primers,
                    default: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG,CTCGTATGCCGTCTTCTGCTTG',
                    desc: 'Comma-separated primer sequences that be used for trimming',
                    type: :string

      def trim_3end_primer

        opt_minimum_length = "--minimum-length=#{options.minimum_length}"
        primers = options.primers.split(',')
        
        fragments = Hash.new
        tmp = Hash.new
        primers.each do |primer|
          max = primer.length-1
          for i in 0..max do
            fragment = primer[0..i]
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

        fragments.keys.sort.reverse.each do |length|
          if 4**length == fragments[length].size
            commands << "fq1l trim_3end_length --trimming-length=#{length} #{opt_minimum_length}"
            break
          else
            fragments[length].sort.reverse.each do |fragment|
              tmpfiles << tmpfile = get_temporary_path("fq1l.trim_3end_primer.#{fragment}", 'fq1l', false)
              commands << "fq1l trim_3end#{' --coreutils-prefix='+options.coreutils_prefix if options.key?(:coreutils_prefix)}#{' --grep-prefix='+options.grep_prefix if options.key?(:grep_prefix)} #{opt_minimum_length} --trimmed=#{tmpfile} #{fragment}"
            end
          end
        end
        stats = Open3.pipeline(*commands)
        stats.each_index do |i|
          unless stats[i].success?
            unlink_files(tmpfiles)
            raise "Fail at process #{i}; #{stats[i]}; #{commands[i]}" 
          end
        end
        system "#{cat_command} #{tmpfiles.join(' ')}"
        unlink_files(tmpfiles)
        
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
                    default: 24,
                    desc: 'Minimum length after trimming',
                    type: :numeric
                      
      
      def trim_5end(pattern)
        exit unless STDIN.wait
        BioGadget.t5(pattern, options.minimum_length)
      end

      # #

      # no_commands do
        
      #   def pipeline(parallel, *commands)
      #     stats = Array.new
      #     tmpin = nil
      #     tmpout = nil
      #     begin
      #       while commands.size > 0
      #         cmds = commands.shift(parallel)
      #         cmds[0] = cmds[0] + " < #{tmpin}" unless tmpin.nil?
      #         if commands.size > 0
      #           tmpout = get_temporary_path('pipeline', 'tmp', false)
      #           cmds[-1] = cmds[-1] + " > #{tmpout}"
      #         end
      #         tmpstats = Open3.pipeline(*cmds)
      #         stats.concat(tmpstats)
      #         tmpstats.each {|tmpstat| commands = nil unless tmpstat.success? }
      #         break if commands.nil?
      #         File.unlink(tmpin) unless tmpin.nil?
      #         tmpin = tmpout
      #       end
      #     ensure
      #       File.unlink(tmpin) if !tmpin.nil? && File.exist?(tmpin)
      #       File.unlink(tmpout) if !tmpout.nil? && File.exist?(tmpout)
      #     end
      #     stats
      #   end
        
      # end
      
    end
  end
end
