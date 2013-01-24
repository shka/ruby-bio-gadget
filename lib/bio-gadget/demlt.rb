require 'bio-faster'
require 'levenshtein'
require 'mkfifo'
require 'parallel'

module Bio
  class Gadget < Thor
    namespace :bio

    desc 'demlt BC POS LEN', 'demultiplex fastq (via STDIN) by barcodes'
    option 'output-dir', :aliases => '-o', :type => :string, :default  => '.'
    def demlt(bcfile, tmpofs, tmplen)

      ofs = tmpofs.to_i
      len = tmplen.to_i

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
      procs.times { |i| fifo1paths.push(mytempfile('fifo1-')) }
      pid = Kernel.fork {
        fifo1s = Array.new
        fifo1paths.each { |fifo1path| fifo1s.push(open(fifo1path, 'w+')) }
        total = 0
        Bio::Faster.new(:stdin).each_record(:quality => :raw) do |vals|
          fifo1 = fifo1s[total % procs]
          fifo1.puts(vals.join("\t"))
          fifo1.flush
          total += 1
        end
        fifo1s.each { |fifo1| fifo1.puts('*'); fifo1.close }
        Kernel.exit!
      }

      fifo1paths.each { |fifo1path|
        until File.exist?(fifo1path)
          sleep 1
        end
      }

      fifo2paths = Array.new
      procs.times { |i|
        fifo2path = mytempfile('fifo2-')
        fifo2paths.push(fifo2path)
        pid = Kernel.fork {
          open(fifo2path, 'w+') { |fifo2|
            fifo1 = open(fifo1paths[i], 'r+')
            while true
              line = fifo1.gets
              if line.nil?
                sleep 1
              elsif line == "*\n"
                break
              else
                seqid, seq, qvs = line.rstrip.split(/\t/)
                tmpdists = Hash.new
                bcs.each_index { |bcidx|
                  tmpdists[bcidx] = Levenshtein.distance(bcs[bcidx], seq[ofs, bclen])
                }
                dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
                bc = dists[0][1] < 2 && dists[0][1] < dists[1][1] ? dists[0][0] : -1
                fifo2.puts("#{bc}\t#{seqid}\t#{seq}\t#{qvs}")
                fifo2.flush
              end
            end
            fifo1.close
            fifo2.puts('*')
          }
          Kernel.exit!
        }
      }

      fifo2paths.each { |fifo2path|
        until File.exist?(fifo2path)
          sleep 1
        end
      }

      tmpwells = wells + ['other']

      fifo3paths = Array.new
      tmpwells.each_index { |i| fifo3paths.push(mytempfile('fifo3-')) }
      pid = Kernel.fork {
        fifo2s = Array.new
        fifo2paths.each { |fifo2path| fifo2s.push(open(fifo2path, 'r+')) }
        fifo2done = Hash.new
        fifo3s = Array.new
        fifo3paths.each { |fifo3path| fifo3s.push(open(fifo3path, 'w+')) }
        fifo2s.cycle { |fifo2|
          unless fifo2done.key?(fifo2)
            line = fifo2.gets
            if line.nil?
              sleep 1
            elsif line == "*\n"
              # puts("#{fifo2} eof.")
              fifo2done[fifo2] = ''
            else
              bcs, seqid, seq, qvs = line.rstrip.split(/\t/)
              fifo3 = fifo3s[bcs.to_i]
              fifo3.puts([seqid, seq, qvs].join("\t"))
              fifo3.flush
            end
          end
          if fifo2done.size == fifo2s.size
            break
          end
        }
        fifo2s.each { |fifo2| fifo2.close }
        fifo3s.each { |fifo3| fifo3.puts('*'); fifo3.close }
        Kernel.exit!
      }

      fifo3paths.each { |fifo3path|
        until File.exist?(fifo3path)
          sleep 1
        end
      }

      tmpwells.each_index { |i|
        well = tmpwells[i]
        outpath = "#{options['output-dir']}/#{well}.fq.xz"
        pid = Kernel.fork {
          left = ofs+bclen
          right = ofs+bclen+len-1
          preprocess = ofs > 0 ? <<"DEDUPandFORMAT"
| ruby -F'\\t' -anle 'f1=$F[1][0..#{right}];f2=$F[2][0..#{right}];puts([f1+f2, $F[0], f2, f1].join("\\t"))' \\
| sort -k 1 -r | cut -f 2- | uniq -f 2 \\
| ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[2][#{left}..-1], "+", $F[1][#{left}..-1]].join("\\n"))' \\
DEDUPandFORMAT
          : <<"FORMAT"
| ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[1][#{left}..#{right}], "+", $F[2][#{left}..#{right}].rstrip].join("\\n"))' \\
FORMAT
          preprocess += "| xz -z -c -e > #{outpath}"
          open(preprocess, 'w') { |fp|
            fifo3 = open(fifo3paths[i], 'r+')
            while true
              line = fifo3.gets
              if line.nil?
                sleep 1
              elsif line == "*\n"
                break
              else
                fp.puts(line)
                fp.flush
              end
            end
            fifo3.close
          }
          Kernel.exit!
        }
      }

      Process.waitall

    end

  end
end
