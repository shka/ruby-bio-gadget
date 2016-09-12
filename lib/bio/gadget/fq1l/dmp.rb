module Bio
  module Gadget
    class Fq1l < Thor

      desc 'dmp MAP BASE', 'Dempltiplex (and restore); BASE is a basename of the demultiplexed files'

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

      def dmp(map, base)
        bcs = readBarcodeMap(map)
        #
        prefix = options.prefix_coreutils
        gprefix = options.prefix_grep
        cmds = ''
        bcs.values.each do |well|
          fifo = Bio::Gadgets.mkfifo('fq1l.dmp', 'fq1l')
          cmds += "#{prefix}tee #{fifo} | "
          Process.fork do
            exec "#{gprefix}grep -P '^[^\t]+ #{well}\t' #{fifo} | #{prefix}tr \"\\t\" \"\\n\" | pigz -c > #{base}.#{well}.fq.gz"
          end
        end
        #
        exec "#{cmds} #{gprefix}grep -P '^[^\t]+ undef\t'"
      end
      
    end
  end
end
