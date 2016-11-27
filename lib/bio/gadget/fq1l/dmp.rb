require 'io/wait'
require 'parallel'

module Bio
  class Gadget
    class Fq1l < Bio::Gadget

      desc 'dmp MAP BASE', 'Dempltiplex (and restore); BASE is a basename of the demultiplexed files'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_PARALLEL

      def dmp(map, base)
        bcs = readBarcodeMap(map)
        tmpfile = getTmpname('fq1l.dmp', 'fq1l.gz')
        p = options.parallel
        system "pigz -p #{p} -c > #{tmpfile}"
        #
        prefix = options.coreutils_prefix
        Parallel.each(bcs.values, in_threads: p) do |well|
          system "unpigz -c #{tmpfile} | #{grepCommand(options)} -P '^[^\t]+ #{well}\t' | #{prefix}tr \"\\t\" \"\\n\" | pigz -p #{p} -c > #{base}.#{well}.fq.gz"
        end
        system "unpigz -c #{tmpfile} | #{grepCommand(options)} -P '^[^\t]+ undef\t'"
      end
      
    end
  end
end
