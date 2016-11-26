module Bio
  module Gadget
    class Fq1l < Thor

      desc 'rst', 'Restore fastq from 1 line/read to 4 lines/read'
      
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX

      def rst
        exec "#{options.coreutils_prefix}tr \"\\t\" \"\\n\""
      end
      
    end
  end
end
