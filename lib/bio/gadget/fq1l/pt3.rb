require 'io/wait'

module Bio
  class Gadget
    class Fq1l < Bio::Gadget
      
      desc 'pt3', 'Primer trimming, 3\'-end'

      method_option :primer,
                    default: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG',
                    desc: 'Primer sequence that be used for trimming',
                    type: :string
      
      method_option :minimum_length,
                    banner: 'NT',
                    default: 40,
                    desc: 'Minimum length after trimming',
                    type: :numeric
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX

      def pt3
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
        #
        exit unless STDIN.wait
        #
        prefix = options.coreutils_prefix
        main = ''
        fifos = Array.new
        tmpfiles = Array.new
        mlen = options.minimum_length
        fragments.keys.sort.reverse.each do |length|
          fifo = mkfifo('fq1l.pt3', 'fq1l', false)
          tmpfiles << tmpfile = getTmpname('fq1l.pt3', 'fq1l', false)
          if 4**length == fragments[length].size
            main += "#{prefix}cat > #{fifo}"
            Kernel.fork do
              BioGadget.pt3(length, "#{prefix}cat #{fifo}", tmpfile, mlen)
            end
          else
            patterns = (fragments[length].map { |fragment| "-e '#{fragment}\t+'" }).join ' '
            main += "#{teeCommand(options)} #{fifo} | #{grepCommand(options)} -v #{patterns} | "
            Kernel.fork do
              BioGadget.pt3(length, "#{grepCommand(options)} #{patterns} #{fifo}", tmpfile, mlen)
            end
          end
          at_exit {
            begin
              File.unlink(fifo) if FileTest.exist?(fifo)
            rescue Errno::ENOENT
            end
          }
          break if 4**length == fragments[length].size
        end
        tmpfiles.map { |tmpfile|
          at_exit {
            begin
              File.unlink(tmpfile) if FileTest.exist?(tmpfile)
            rescue Errno::ENOENT
            end
          }
        }
        system main
        Process.waitall
        system "#{prefix}cat #{tmpfiles.join(' ')}"
      end
      
    end
  end
end
