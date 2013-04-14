require 'bio-faster'
require 'levenshtein'
require 'mkfifo'
require 'parallel'

module Bio
  class Gadget < Thor

    namespace :bio

    desc 'demlt BARCODE [FASTQ]', "Demultiplex fastq from STDIN by barcodes.\n\n"
    option 'output-dir', :aliases => '-o', :type => :string, :default => '.'
    option 'umi-length', :aliases => '-u', :type => :numeric, :default => 4, :desc => '0 is no umi, means no PCR-amplicon reduction.'
    option 'cdna-length', :aliases => '-c', :type => :numeric, :default => 37, :desc => 'Trimming length before PCA-amplicon reduction. -1 is no trimming by length.'
    option 'g-trimming', :aliases => '-g', :type => :boolean, :default => false, :desc => "Trimming of 5'-end poly-G. Length of the trimmed Gs attached after the read name."
    option 'q-trimming', :aliases => '-q', :type => :string, :default => '~', :desc => "Quality threshold - nucleotides with lower quality will be trimmed, from the end of the sequence. '~' is no trimming by quality, because this is the maximum quality base character."
    option 'min-length', :aliases => '-l', :type => :numeric, :default => 0, :desc => 'Length threshold - sequences shorter than this after trimming will be filtered out. 0 is no filtering.'
    def demlt(bcfile, fastq=:stdin)

      ofs = options['umi-length']
      clen = options['cdna-length']
      gtrim = options['g-trimming']
      qtrim = options['q-trimming']
      mlen = options['min-length']

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
          right = clen > -1 ? -1 : ofs+bclen+clen-1
          preprocess = ofs > 0 ? <<"DEDUPandFORMAT"
ruby -F'\\t' -anle 'f1=$F[1][0..#{right}];f2=$F[2][0..#{right}];puts([f1+f2, $F[0], f2, f1].join("\\t"))' #{fifo3paths[i]} \\
| sort -k 1 -r | cut -f 2- | uniq -f 2 \\
| ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[2][#{left}..-1], "+", $F[1][#{left}..-1]].join("\\n"))' \\
DEDUPandFORMAT
          : <<"FORMAT"
ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[1][#{left}..#{right}], "+", $F[2][#{left}..#{right}].rstrip].join("\\n"))' #{fifo3paths[i]} \\
FORMAT

          preprocess += '| ruby -e \'require "bio-faster";Bio::Faster.new(:stdin).each_record(:quality=>:raw){|v|s=v[1].gsub(/^G+/,"");l=v[1].length-s.length;puts("@#{v[0]}|-G#{l}\\n#{s}\\n+\\n#{v[2][l,s.length]}") if s.length>0}\'' if gtrim

          if qtrim != '~' || mlen > 0
            preprocess += '| ruby -e \'require "bio-faster";Bio::Faster.new(:stdin).each_record(:quality=>:raw){|v|m=v[2].length-1;0.upto(m){|i|if v[2][i]<"'+qtrim+'" then m=i-1;break;end};puts("@#{v[0]}\n#{v[1][0..m]}\n+\n#{v[2][0..m]}") if m+1>='+mlen.to_s+'}\''
          end

          exec preprocess+"| xz -z -c -e > #{outpath}"
        }
      }

      Process.waitall

    end

  end
end
