require 'io/wait'
require 'parallel'

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

      method_option :parallel,
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
      
      def dmp(map, base)
        bcs = readBarcodeMap(map)
        tmpfile = Bio::Gadgets.getTmpname('fq1l.dmp', 'fq1l.gz')
        p = options.parallel
        system "pigz -p #{p} -c > #{tmpfile}"
        #
        prefix = options.prefix_coreutils
        gprefix = options.prefix_grep
        Parallel.each(bcs.values, in_threads: p) do |well|
          system "unpigz -c #{tmpfile} | #{gprefix}grep -P '^[^\t]+ #{well}\t' | #{prefix}tr \"\\t\" \"\\n\" | pigz -p #{p} -c > #{base}.#{well}.fq.gz"
        end
        system "unpigz -c #{tmpfile} | #{gprefix}grep -P '^[^\t]+ undef\t'"
      end
      
    end
  end
end
