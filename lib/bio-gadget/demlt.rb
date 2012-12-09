require 'bio-faster'
require 'levenshtein'
require 'parallel'
require 'thread'

module Bio
  class Gadget < Thor
    namespace :bio

    desc 'demlt BC POS', 'demultiplex fastq (via STDIN) by barcodes'
    option 'output-dir', :type => :string, :default  => '.', :aliases => '-o'
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
        q = SizedQueue.new(100000)
        t = Thread.new(well, q) do |well, q|
          tc = Thread.current
          tc[:file] = "#{options['output-dir']}/#{well}.fq.xz"
          fp = open("| xz -z -c -e > #{tc[:file]}", 'w')
          tc[:read] = 0
          while vals = q.shift
            if vals == ""
              break
            else
              fp.puts(vals)
              tc[:read] = tc[:read] + 1
            end
          end
          fp.close()
        end
        qs.push(q)
        ts.push(t)
      }

      rq = Queue.new
      Thread.new(rq) {
        Bio::Faster.new(:stdin).each_record(:quality => :raw) do |seqid, seq, qvs|
          rq.push([seqid, seq, qvs])
        end
        rq.push('')
      }

      seqs = Array.new
      while vals = rq.shift
        if vals != ""
          seqs.push(vals)
        end
        if vals == "" || seqs.size == 100000 * Parallel.processor_count
          parallel_Levenshtein(seqs, bcs, ofs, bclen, qs)
          seqs = Array.new
        end
        if vals == ""
          qs.each { |q| q.push('') }
          break
        end
      end
      ts.each { |t| t.join }

      total = 0
      bcs.each_index { |i|
        t = ts[i]
        r = t[:read]
        puts "#{bcs[i]}\t#{r}\t#{t[:file]}"
        total = total+r
      }
      t = ts[-1]
      r = t[:read]
      puts "Other\t#{r}\t#{t[:file]}"
      puts '===='
      puts "Total\t#{total+r}"

    end

    protected

    def parallel_Levenshtein(seqs, bcs, ofs, bclen, qs)

      tmpdists = Parallel.map_with_index(bcs, :in_processes => Parallel.processor_count) do |bc, bcidx|
        tmpdist = Array.new
        seqs.each_index do |seqidx|
          seqbc = seqs[seqidx][1][ofs, bclen]
          tmpdist.push(Levenshtein.distance(bc, seqbc))
        end
        tmpdist
      end

      tmpdist = Hash.new
      seqs.each_index do |seqidx|
        seqid, seq, qvs = seqs[seqidx]
        bcs.each_index do |bcidx|
          tmpdist[bcidx] = tmpdists[bcidx][seqidx]
        end
        dists = tmpdist.sort { |a, b| a[1] <=> b[1] }
        if dists[0][1] < dists[1][1] && dists[0][1] < 2
          idx = dists[0][0]
          qs[idx].push("@#{seqid}\n#{seq}\n+\n#{qvs}")
        else
          qs[-1].push("@#{seqid}\n#{seq}\n+\n#{qvs}")
        end
      end

    end

  end
end
