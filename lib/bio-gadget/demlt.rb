require 'bio-faster'
require 'levenshtein'
require 'mkfifo'
require 'parallel'

module Bio
  class Gadget < Thor

    namespace :bio

    desc 'demlt BARCODE [FASTQ]', 'demultiplex fastq from STDIN by barcodes'
    option 'output-dir', :aliases => '-o', :type => :string, :default => '.'
    option 'umi-length', :aliases => '-u', :type => :numeric, :default => 4
    option 'cdna-length', :aliases => '-c', :type => :numeric, :default => 37
    def demlt(bcfile, fastq=:stdin)

      ofs = options['umi-length']
      len = options['cdna-length']

      wells = Array.new
      bcs = Array.new
      bclens = Array.new
      open(bcfile).each do |line|
        cols = line.rstrip.split
        wells.push(cols[0])
        bcs.push(cols[1])
        bclens.push(cols[1].length)
      end

      bclens.uniq!
      if bclens.size != 1
        raise 'Inconsistent barcode sequence lengths'
      end
      bclen = bclens[0]

      procs = Parallel.processor_count

      fifo1paths = Array.new
      procs.times { |i|
        fifo1path = mytemppath('fifo1-')
        File.mkfifo(fifo1path)
        fifo1paths.push(fifo1path)
      }
      pid = Kernel.fork {
        fifo1s = Array.new
        fifo1paths.each { |fifo1path| fifo1s.push(open(fifo1path, 'w')) }
        total = 0
        Bio::Faster.new(fastq).each_record(:quality => :raw) do |vals|
          fifo1 = fifo1s[total % procs]
          fifo1.puts(vals.join("\t"))
          total += 1
        end
        fifo1s.each { |fifo1| fifo1.close }
        Kernel.exit!
      }

      fifo2paths = Array.new
      procs.times { |i|
        fifo2path = mytemppath('fifo2-')
        File.mkfifo(fifo2path)
        fifo2paths.push(fifo2path)
        pid = Kernel.fork {
          open(fifo2path, 'w') { |fifo2|
            open(fifo1paths[i], 'r').each { |line|
              seqid, seq, qvs = line.rstrip.split(/\t/)
              tmpdists = Hash.new
              bcs.each_index { |bcidx|
                tmpdists[bcidx] = Levenshtein.distance(bcs[bcidx], seq[ofs, bclen])
              }
              dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
              bc = dists[0][1] < 2 && dists[0][1] < dists[1][1] ? dists[0][0] : -1
              fifo2.puts("#{bc}\t#{seqid}\t#{seq}\t#{qvs}")
            }
          }
          Kernel.exit!
        }
      }

      tmpwells = wells + ['other']

      fifo3paths = Array.new
      tmpwells.each_index { |i|
        fifo3path = mytemppath('fifo3-')
        File.mkfifo(fifo3path)
        fifo3paths.push(fifo3path)
      }
      pid = Kernel.fork {
        fifo2s = Array.new
        fifo2paths.each { |fifo2path| fifo2s.push(open(fifo2path, 'r')) }
        fifo2done = Hash.new
        fifo3s = Array.new
        fifo3paths.each { |fifo3path| fifo3s.push(open(fifo3path, 'w')) }
        fifo2s.cycle { |fifo2|
          unless fifo2done.key?(fifo2)
            line = fifo2.gets
            if line.nil?
              fifo2done[fifo2] = ''
            else
              bcs, seqid, seq, qvs = line.rstrip.split(/\t/)
              fifo3 = fifo3s[bcs.to_i]
              fifo3.puts([seqid, seq, qvs].join("\t"))
            end
          end
          if fifo2done.size == fifo2s.size
            break
          end
        }
        fifo2s.each { |fifo2| fifo2.close }
        fifo3s.each { |fifo3| fifo3.close }
        Kernel.exit!
      }

      tmpwells.each_index { |i|
        well = tmpwells[i]
        outpath = "#{options['output-dir']}/#{well}.fq.xz"
        pid = Kernel.fork {
          left = ofs+bclen
          right = ofs+bclen+len-1
          preprocess = ofs > 0 ? <<"DEDUPandFORMAT"
ruby -F'\\t' -anle 'f1=$F[1][0..#{right}];f2=$F[2][0..#{right}];puts([f1+f2, $F[0], f2, f1].join("\\t"))' #{fifo3paths[i]} \\
| sort -k 1 -r | cut -f 2- | uniq -f 2 \\
| ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[2][#{left}..-1], "+", $F[1][#{left}..-1]].join("\\n"))' \\
DEDUPandFORMAT
          : <<"FORMAT"
ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[1][#{left}..#{right}], "+", $F[2][#{left}..#{right}].rstrip].join("\\n"))' #{fifo3paths[i]} \\
FORMAT
          exec preprocess+"| xz -z -c -e > #{outpath}"
        }
      }

      Process.waitall

    end

  end
end
