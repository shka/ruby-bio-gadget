module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'cnv', 'Convert fastq from 4 lines/read to 1 line/read'
      
      def cnv
        exec "#{options.prefix_coreutils}paste - - - -"
      end
      
    end
  end
end
