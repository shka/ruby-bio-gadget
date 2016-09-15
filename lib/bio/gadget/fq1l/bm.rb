require 'io/wait'
require 'damerau-levenshtein'

module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'bm MAP', 'Barcode match; a csv file MAP w/ barcode sequence and wellname'

      method_option :begin,
                    aliases: '-b',
                    default: 7,
                    desc: 'Start position of barcode',
                    type: :numeric
      
      method_option :end,
                    aliases: '-e',
                    default: 12,
                    desc: 'End position of barcode',
                    type: :numeric

      method_option :maximum_distance,
                    aliases: '-d',
                    default: 1,
                    desc: 'Maximum distance between barcode and sequence',
                    type: :numeric

      method_option :buffer_size,
                    banner: 'SIZE',
                    aliases: '-S',
                    desc: 'Use SIZE for main memory buffer',
                    type: :string
      
      method_option :prefix_coreutils,
                    banner: 'PREFIX',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : '',
                    desc: 'A prefix character for GNU coreutils',
                    type: :string

      method_option :parallel,
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
                    
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
        open("| #{options.prefix_coreutils}sort -t '\t' --parallel=#{options.parallel} -k2.#{options.begin},2.#{options.end} #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}", 'r').each do |line|
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
