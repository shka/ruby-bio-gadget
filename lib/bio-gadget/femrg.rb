module Bio
  class Gadget < Thor

    desc 'femrg [GTF...]', 'Extract and merge overlapping first exons.'

    def femrg(*gtf_files)

      exon1s = Hash.new
      fp = open("| grep -h '\texon\t' #{gtf_files.nil? ? '-' : gtf_files.join(' ')} | cut -f 1,4,5,7,9")
      fp.each { |line|
        chr, sstart, sstop, str, attr = line.rstrip.split(/\t/)
        exon1s[chr] = Hash.new unless exon1s.key?(chr)
        exon1s[chr][str] = Hash.new unless exon1s[chr].key?(str)
        id = attr.match(/transcript_id \"([^\"]+)\"/).to_a[1]
        en = attr.match(/exon_number \"(\d+)\"/).to_a[1].to_i
        if (!exon1s[chr][str].key?(id) ||
            (str == "+" && en < exon1s[chr][str][id][2]) ||
            (str == "-" && exon1s[chr][str][id][2] < en))
          exon1s[chr][str][id] = [sstart.to_i, sstop.to_i, en]
        end
      }
      fp.close

      idx = 0
      exon1s.each { |chr, exon1schr|
        exon1schr.each { |str, exon1schrstr|
          ids = exon1schrstr.keys.sort { |a, b|
            exon1schrstr[a][0] <=> exon1schrstr[b][0]
          }
          #
          clusters = Array.new
          members = [id = ids.shift]
          start, stop, en = exon1schrstr[id]
          while ids.length > 0
            nid = ids.shift
            nstart, nstop, nen = exon1schrstr[nid]
            if stop < nstart
              clusters.push([members, start, stop])
              members = [nid]
              start = nstart
              stop = nstop
            else
              members.push(nid)
              start, stop = [start, stop, nstart, nstop].sort.values_at(0, 3)
            end
          end
          clusters.push([members, start, stop])
          #
          clusters.each { |members, start, stop|
            attr = "gene_id \"FE#{idx}\"; transcript_id \"FE#{idx}\"; "
            puts [chr, 'bio-gadget:femrg', 'transcript', start, stop, 1000, str, '.', attr].join("\t")
            puts [chr, 'bio-gadget:femrg', 'exon', start, stop, 1000, str, '.', "#{attr}exon_number \"1\"; member_ids \"#{members.join('|')}\""].join("\t")
            idx += 1
          }
        }
      }

    end

  end
end
