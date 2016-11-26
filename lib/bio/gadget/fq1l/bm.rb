require 'io/wait'
require 'damerau-levenshtein'

module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'bm MAP', 'Barcode match; a csv file MAP w/ barcode sequence and wellname'

      method_option :begin,
                    default: 7,
                    desc: 'Start position of barcode',
                    type: :numeric
      
      method_option :end,
                    default: 12,
                    desc: 'End position of barcode',
                    type: :numeric

      method_option :maximum_distance,
                    default: 1,
                    desc: 'Maximum distance between barcode and sequence',
                    type: :numeric

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX
      method_option *Bio::Gadgets::OPT_PARALLEL
                    
      def bm(map)
        bcs = readBarcodeMap(map)
        #
        exit unless STDIN.wait
        #
        pbc = ''
        pwell = ''
        bidx = options.begin-1
        eidx = options.end-1
        bcr = bidx..eidx
        dl = DamerauLevenshtein
        open("| #{options.coreutils_prefix}sort -t '\t' --parallel=#{options.parallel} -k2.#{options.begin},2.#{options.end} #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}", 'r').each do |line|
          acc, seq, sep, qual = line.rstrip.split(/\t/)
          bc = seq[bcr]
          if bc != pbc
            mindist = bc.length
            minbc = ''
            bcs.keys.each do |k|
              dist = dl.distance(k, bc)
              if dist < mindist
                mindist = dist
                minbc = k
              end
              break if dist == 0
            end
            pbc = bc
            pwell = mindist <= options.maximum_distance ? bcs[minbc] : 'undef'
          end
          puts ["#{acc} #{pwell}", seq, sep, qual].join("\t")
        end
      end
      
    end
  end
end
