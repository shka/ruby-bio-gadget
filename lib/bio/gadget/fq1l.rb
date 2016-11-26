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

      method_option *Bio::Gadgets::OPT_PREFIX_COREUTILS
      
      def convert
        exit unless STDIN.wait
        exec "#{options.prefix_coreutils}paste - - - -"
      end

      # fq1l:sort

      desc 'sort', '(filter) sort by sequences and the qualities'

      method_option *Bio::Gadgets::OPT_PREFIX_COREUTILS
      method_option *Bio::Gadgets::OPT_PARALLEL
      method_option *Bio::Gadgets::OPT_BUFFER_SIZE

      def sort
        exit unless STDIN.wait
        exec "#{options.prefix_coreutils}sort -t '\t' --parallel=#{options.parallel} -r -k2,4 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}"
      end
      
      #

      no_commands do
        
        def readBarcodeMap(map)
          bcs = Hash.new
          open(map, 'r').each do |line|
            bc, well = line.rstrip.split(',')
            bcs[bc] = well
          end
          bcs
        end
        
      end

    end
  end
end
