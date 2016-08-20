module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'convert', 'Convert fastq from 4 lines/read to 1 line/read'
      
      def convert
        exec "#{options.prefix_coreutils}paste - - - -"
      end
      
    end
  end
end
