require 'bio/gadget/fq1l/bm'
require 'bio/gadget/fq1l/cnv'
require 'bio/gadget/fq1l/dmp'
require 'bio/gadget/fq1l/m5'
require 'bio/gadget/fq1l/mt5'
require 'bio/gadget/fq1l/nr'
require 'bio/gadget/fq1l/pt3'
require 'bio/gadget/fq1l/qt3'
require 'bio/gadget/fq1l/rst'

module Bio
  module Gadget
    class Fq1l < Thor
      
      class_option :prefix_coreutils,
                   type: :string,
                   banner: 'PREFIX',
                   desc: 'A prefix character for GNU coreutils',
                   default: system('which gnproc >/dev/null 2>&1') ? 'g' : ''

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
