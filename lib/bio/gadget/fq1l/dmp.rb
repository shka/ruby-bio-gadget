require 'io/wait'
require 'parallel'

module Bio
  module Gadget
    class Fq1l < Thor

      desc 'dmp MAP BASE', 'Dempltiplex (and restore); BASE is a basename of the demultiplexed files'

      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX
      method_option *Bio::Gadgets::OPT_PARALLEL

      method_option :prefix_grep,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU grep',
                    default: system('which ggrep >/dev/null 2>&1') ? 'g' : ''

      def dmp(map, base)
        bcs = readBarcodeMap(map)
        tmpfile = Bio::Gadgets.getTmpname('fq1l.dmp', 'fq1l.gz')
        p = options.parallel
        system "pigz -p #{p} -c > #{tmpfile}"
        #
        prefix = options.coreutils_prefix
        gprefix = options.prefix_grep
        Parallel.each(bcs.values, in_threads: p) do |well|
          system "unpigz -c #{tmpfile} | #{gprefix}grep -P '^[^\t]+ #{well}\t' | #{prefix}tr \"\\t\" \"\\n\" | pigz -p #{p} -c > #{base}.#{well}.fq.gz"
        end
        system "unpigz -c #{tmpfile} | #{gprefix}grep -P '^[^\t]+ undef\t'"
      end
      
    end
  end
end
