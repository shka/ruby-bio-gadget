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
        q = Queue.new
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

      Bio::Faster.new(:stdin).each_record(:quality => :raw) do |seqid, seq, qvs|
        tmpdists = Hash.new
        bcs.each_index { |bcidx|
          tmpdists[bcidx] = Levenshtein.distance(bcs[bcidx], seq[ofs, bclen])
        }
        dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
        if dists[0][1] < 2 && dists[0][1] < dists[1][1]
          qs[dists[0][0]].push("@#{seqid}\n#{seq}\n+\n#{qvs}")
        else
          qs[-1].push("@#{seqid}\n#{seq}\n+\n#{qvs}")
        end
        Thread.pass
      end

      qs.each { |q| q.push('') }
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

  end
end
