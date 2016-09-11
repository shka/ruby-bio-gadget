module Bio
  module Gadget
    class Fq1l < Thor

      desc 'rst', 'Restore fastq from 1 line/read to 4 lines/read'
      
      def rst
        exec "#{options.prefix_coreutils}tr \"\\t\" \"\\n\""
      end
      
    end
  end
end
