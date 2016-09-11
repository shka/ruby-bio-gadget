require 'damerau-levenshtein'

module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'bm MAP', 'Barcode match; a csv file MAP w/ barcode sequence and wellname'

      method_option :begin,
                    aliases: '-b',
                    default: 7,
                    type: :numeric,
                    desc: 'Start position of barcode'
      
      method_option :end,
                    aliases: '-e',
                    default: 12,
                    type: :numeric,
                    desc: 'End position of barcode'

      method_option :maximum_distance,
                    aliases: '-d',
                    default: 1,
                    type: :numeric,
                    desc: 'Maximum distance between barcode and sequence'

      method_option :buffer_size,
                    banner: 'SIZE',
                    aliases: '-S',
                    desc: 'Use SIZE for main memory buffer'
      
      def bm(map)
        bcs = readBarcodeMap(map)
        #
        pbc = ''
        pwell = ''
        bidx = options.begin-1
        eidx = options.end-1
        bcr = bidx..eidx
        dl = DamerauLevenshtein
        prefix = options.prefix_coreutils
        open("| #{prefix}sort -t '\t' --parallel=`#{prefix}nproc` -k2.#{options.begin},2.#{options.end} #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}", 'r').each do |line|
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
