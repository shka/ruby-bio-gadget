module Bio
  class Gadget < Thor

    desc 'peak WIG1,WIG2,... [GTF]', <<DESC
Find peak within each exon from (gzipped) variableStep wigs by a majority vote. It will read from a standard input if no GTF option.
DESC
    def peak(wigs, gtf="/dev/stdin")

      nchrpos2val = Hash.new
      wigs.split(/,/).each { |wig|
        n = wigs.split(/,/).index(wig)
        nchrpos2val[n] = Hash.new unless nchrpos2val.key?(n)
        chr = ''
        fp = open("| zcat #{wig} | grep -v ^track")
        fp.each { |line|
          cols = line.rstrip.split(/\s+/)
          if cols[0] == 'variableStep'
            chr = cols[1].match(/chrom=(\S+)$/).to_a[1]
            nchrpos2val[n][chr] = Hash.new unless nchrpos2val[n].key?(chr)
          else
            nchrpos2val[n][chr][cols[0].to_i] = cols[1].to_f
          end
        }
        fp.close
      }

      open("| grep exon #{gtf}").each { |line|
        cols = line.rstrip.split(/\t/)
        oid = cols[8].match(/oId \"([^\"]+)/).to_a[1]
        oid = "#{cols[0]}.#{oid}" if cols[0] =~ /^RNA_SPIKE_/
        exn = cols[8].match(/exon_number \"(\d+)\"/).to_a[1]
        chr = cols[0]
        str = cols[6]
        start = cols[3].to_i
        stop = cols[4].to_i
        peak = ''
        #
        poss = Hash.new
        nchrpos2val.each { |n, chrpos2val|
          if chrpos2val.key?(chr)
            pos2val = chrpos2val[chr]
            tmppos2val = Hash.new
            pos2val.each { |pos, val|
              tmppos2val[pos] = val if start <= pos && pos <= stop
            }
            if tmppos2val.size > 0
              tmpposs = tmppos2val.keys.sort { |a, b|
                tmppos2val[b] == tmppos2val[a] ? (str == '+' ? (a <=> b) : (b <=> a)) : (tmppos2val[b] <=> tmppos2val[a])
              }
              tmppos = tmpposs[0]
              # puts "#{n} | #{chr}:#{start}-#{stop} #{str} | #{tmpposs}"
              poss[tmppos] = poss.key?(tmppos) ? poss[tmppos]+1 : 1
            end
          end
        }
        if poss.size > 0
          peaks = poss.keys.sort { |a, b|
            poss[b] == poss[a] ? (str == '+' ? (a <=> b) : (b <=> a)) : (poss[b] <=> poss[a])
          }
          peak = peaks[0]
        end
        #
        puts [oid, exn, peak].join("\t")
      }

    end

  end
end
