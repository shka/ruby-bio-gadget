module Bio
  class Gadget
    class Fq1l < Bio::Gadget

      desc 'rst', 'Restore fastq from 1 line/read to 4 lines/read'
      
      method_option *OPT_COREUTILS_PREFIX

      def rst
        exec "#{options.coreutils_prefix}tr \"\\t\" \"\\n\""
      end
      
    end
  end
end
