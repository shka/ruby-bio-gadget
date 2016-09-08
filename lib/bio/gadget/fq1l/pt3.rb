module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'pt3', 'Primer trimming, 3\'-end'

      method_option :primer,
                    aliases: '-p',
                    type: :string,
                    default: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG',
                    desc: 'Primer sequence that be used for trimming'
      
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
        prefix = options.prefix_coreutils
        main = ''
        fifos = Array.new
        tmpfiles = Array.new
        fragments.keys.sort.reverse.each do |length|
          fifo = Bio::Gadgets.mkfifo('fq1l.pt3', 'fq1l', false)
          tmpfiles << tmpfile = Bio::Gadgets.getTmpname('fq1l.pt3', 'fq1l', false)
          if 4**length == fragments[length].size
            main += "#{prefix}cat > #{fifo}"
            Kernel.fork do
              BioGadget.pt3(length, "#{prefix}cat #{fifo}", tmpfile)
            end
          else
            patterns = (fragments[length].map { |fragment| "-e '#{fragment}\t+'" }).join ' '
            main += "#{prefix}tee #{fifo} | #{prefix}grep -v #{patterns} | "
            Kernel.fork do
              BioGadget.pt3(length, "#{prefix}grep #{patterns} #{fifo}", tmpfile)
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
        system "#{prefix}cat #{tmpfiles.join(' ')}"
        Process.waitall
      end
      
    end
  end
end
