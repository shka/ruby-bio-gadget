module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'at3 ADAPTOR', 'Adaptor trimming, 3\'-end'

      def at3(adaptor)
        fragments = Hash.new
        max = adaptor.length-1
        tmp = Hash.new
        for i in 0..max do
          for j in i..max do
            fragment = adaptor[i..j]
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
        tmpfiles = Array.new
        fragments.keys.sort.reverse.each do |length|
          fifo = Bio::Gadgets.mkfifo('fq1l.at3', 'fq1l')
          tmpfiles << tmpfile = Bio::Gadgets.getTmpname('fq1l.at3', 'fq1l', false)
          if 4**length == fragments[length].size
            main += "#{prefix}cat > #{fifo}"
            Kernel.fork do
              BioGadget.trim3(length, "#{prefix}cat #{fifo}", tmpfile)
            end
          else
            patterns = (fragments[length].map { |fragment| "-e '#{fragment}\t+'" }).join ' '
            main += "#{prefix}tee #{fifo} | #{prefix}grep -v #{patterns} | "
            Kernel.fork do
              BioGadget.trim3(length, "#{prefix}grep #{patterns} #{fifo}", tmpfile)
            end
          end
        end
        tmpfiles.map { |tmpfile| at_exit { File.unlink(tmpfile) if FileTest.exist?(tmpfile) } }
        system main
        system "cat #{tmpfiles.join(' ')}"
        Process.waitall
      end
      
    end
  end
end
