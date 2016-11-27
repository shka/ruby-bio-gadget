require 'bio/gadget/fq1l/bm'
require 'bio/gadget/fq1l/dmp'
require 'bio/gadget/fq1l/m5'
require 'bio/gadget/fq1l/mt5'
require 'bio/gadget/fq1l/nr'
require 'bio/gadget/fq1l/pt3'
require 'bio/gadget/fq1l/qt3'
require 'bio/gadget/fq1l/rst'
require 'bio/gadget/fq1l/to'

module Bio
  module Gadget
    class Fq1l < Thor

      # fq1l:convert
      
      desc 'convert', '(filter) convert fastq from 4 lines/read to 1 line/read'

      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX
      
      def convert
        exit unless STDIN.wait
        exec "#{options.coreutils_prefix}paste - - - -"
      end

      # fq1l:exclude_redundant

      desc 'exclude_duplicate', '(filter) exclude duplicates in the order'

      def exclude_duplicate
        exit unless STDIN.wait
        BioGadget.nr_std()
      end
      
      # fq1l:sort

      desc 'sort', '(filter) sort by sequences and the qualities in descending order'

      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX
      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_PARALLEL

      def sort
        exit unless STDIN.wait
        exec "#{Bio::Gadgets.sortCommand(options)} -t '\t' -r -k2,4"
      end
      
      #

      no_commands do
        
        def readBarcodeMap(map)
          bcs = Hash.new
          open(map, 'r').each do |line|
            bc, well = line.rstrip.split(',')
            bcs[bc] = well
          end
          return bcs
        end

      end

    end
  end
end
