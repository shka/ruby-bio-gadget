require 'bio-faster'
require 'levenshtein'
require 'thread'

module Bio
  class Gadget < Thor
    namespace :bio

    desc 'demlt BC POS', 'demultiplex fastq (via STDIN) by barcodes'
    def demlt(bcfile, tmpofs)

      ofs = tmpofs.to_i

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
        raise 'Inconsistent barcode sequences'
      end
      bclen = bclens[0]

      ts = Array.new
      qs = Array.new
      (wells + ['other']).each { |well|
        q = Queue.new
        t = Thread.new(well, q) do |well|
          fp = open("| gzip -c > #{well}.fq.gz", 'w')
          while vals = q.shift
            if vals == ""
              break
            else
              fp.puts(vals)
            end
          end
          fp.close()
        end
        qs.push(q)
        ts.push(t)
      }

      reads = Array.new(bcs.size+1, 0)
      tmpdist = Hash.new
      Bio::Faster.new(:stdin).each_record(:quality => :raw) do |seqid, seq, qvs|
        seqbc = seq[ofs, bclen]
        bcs.each_index do |i|
          tmpdist[i] = Levenshtein.distance(bcs[i], seqbc)
        end
        dists = tmpdist.sort { |a, b| a[1] <=> b[1] }
        if dists[0][1] < dists[1][1] && dists[0][1] < 2
          idx = dists[0][0]
          qs[idx].push(">#{seqid}\n#{seq}\n+\n#{qvs}")
          reads[idx] = reads[idx]+1
        else
          qs[-1].push(">#{seqid}\n#{seq}\n+\n#{qvs}")
          reads[-1] = reads[-1]+1
        end
      end

      qs.each { |q| q.push('') }
      ts.each { |t| t.join }

      total = 0
      bcs.each_index { |i|
        r = reads[i]
        puts "#{bcs[i]}\t#{r}\t#{wells[i]}.fq.gz"
        total = total+r
      }
      puts "Other\t#{reads[-1]}\tother.fq.gz"
      puts '===='
      puts "Total\t#{total+reads[-1]}"

    end
  end
end
