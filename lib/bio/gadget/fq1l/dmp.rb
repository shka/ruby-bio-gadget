module Bio
  module Gadget
    class Fq1l < Thor

      desc 'dmp MAP', 'Dempltiplex (and restore)'

      method_option :prefix,
                    aliases: '-p',
                    desc: 'prefix for output files',
                    type: :string
      
      def dmp(map)
        bcs = readBarcodeMap(map)
        #
        prefix = options.prefix_coreutils
        cmds = ''
        bcs.values.each do |well|
          fifo = Bio::Gadgets.mkfifo('fq1l.dmp', 'fq1l')
          cmds += "#{prefix}tee #{fifo} | "
          Process.fork do
            exec "#{prefix}grep -P '^[^\t]+ #{well}\t' #{fifo} | #{prefix}tr \"\\t\" \"\\n\" | pigz -c > #{options.prefix}.#{well}.fq.gz"
          end
        end
        #
        exec "#{cmds} #{prefix}grep -P '^[^\t]+ undef\t'"
      end
      
    end
  end
end
