module Bio
  module Gadget
    class Fq1l < Thor

      desc 'restore', 'Restore fastq from 1 line/read to 4 lines/read'
      
      def restore
        exec "#{options.prefix_coreutils}cut -f 1-4 | #{options.prefix_coreutils}tr \"\\t\" \"\\n\""
      end
      
    end
  end
end
