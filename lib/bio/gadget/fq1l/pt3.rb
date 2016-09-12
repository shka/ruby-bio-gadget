module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'pt3', 'Primer trimming, 3\'-end'

      method_option :primer,
                    aliases: '-p',
                    default: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG',
                    desc: 'Primer sequence that be used for trimming',
                    type: :string
      
      method_option :minimum_length,
                    aliases: '-l',
                    banner: 'NT',
                    default: 40,
                    desc: 'Minimum length after trimming',
                    type: :numeric
      
      method_option :prefix_coreutils,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU coreutils',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : ''

      method_option :prefix_grep,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU grep',
                    default: system('which ggrep >/dev/null 2>&1') ? 'g' : ''

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
        gprefix = options.prefix_grep
        main = ''
        fifos = Array.new
        tmpfiles = Array.new
        mlen = options.minimum_length
        fragments.keys.sort.reverse.each do |length|
          fifo = Bio::Gadgets.mkfifo('fq1l.pt3', 'fq1l', false)
          tmpfiles << tmpfile = Bio::Gadgets.getTmpname('fq1l.pt3', 'fq1l', false)
          if 4**length == fragments[length].size
            main += "#{prefix}cat > #{fifo}"
            Kernel.fork do
              BioGadget.pt3(length, "#{prefix}cat #{fifo}", tmpfile, mlen)
            end
          else
            patterns = (fragments[length].map { |fragment| "-e '#{fragment}\t+'" }).join ' '
            main += "#{prefix}tee #{fifo} | #{gprefix}grep -v #{patterns} | "
            Kernel.fork do
              BioGadget.pt3(length, "#{gprefix}grep #{patterns} #{fifo}", tmpfile, mlen)
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
