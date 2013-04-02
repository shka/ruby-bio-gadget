module Bio
  class Gadget < Thor

    desc 'ucscann ISOFORMSGZ XREFGZ [KNOWNGENE]', 'create knownGene.gtf'
    def ucscann(isoforms, xref, knowngene='/dev/stdin')

      acc2id = Hash.new
      open("| gunzip -c #{isoforms}").each { |line|
        id, acc = line.rstrip.split(/\t/)
        acc2id[acc] = id
      }

      acc2sym = Hash.new
      open("| gunzip -c #{xref} | cut -f 1,5").each { |line|
        acc, sym = line.rstrip.split(/\t/)
        acc2sym[acc] = sym
      }

      open("| sort -k1,1 -k4,4n", 'w') { |fp|
        open(knowngene).each { |line|
          cols = line.rstrip.split(/\t/)
          acc, chr, str = cols.values_at(0, 1, 2)
          cs = cols[5].to_i
          ce = cols[6].to_i
          lefts = cols[8].split(/,/)
          rights = cols[9].split(/,/)
          prop = "gene_id \"#{acc2id[acc]}\"; transcript_id \"#{acc}\""
          prop += "; gene_name \"#{acc2sym[acc]}\"" if acc2sym.has_key?(acc)
          lefts.each_index { |i|
            next if lefts[i].nil?
            l = lefts[i].to_i
            r = rights[i].to_i
            fp.puts [chr, 'knownGene', 'exon', l+1, r, 0, str, '.', prop].join("\t")
          }
          if cs != ce
            fp.puts [chr, 'knownGene', 'start_codon', (str == '+' ? cs+1 : ce-2), (str == '+' ? cs+3 : ce), 0, str, '.', prop].join("\t")
            lefts.each_index { |i|
              next if lefts[i].nil?
              l = lefts[i].to_i
              next if ce-1 < l
              r = rights[i].to_i
              next if r-1 < cs
              fp.puts [chr, 'knownGene', 'CDS', (cs < l ? l : cs), (r < ce ? r : ce), 0, str, '.', prop].join("\t")
            }
          end
        }
      }

    end

    desc 'ensann GENENAMEGZ [ENSGENE]', 'create ensGene.gtf'
    def ensann(genename, ensgene='/dev/stdin')

      acc2sym = Hash.new
      open("| gunzip -c #{genename}").each { |line|
        acc, sym = line.rstrip.split(/\t/)
        acc2sym[acc] = sym
      }

      open("| sort -k1,1 -k4,4n", 'w') { |fp|
        open(ensgene).each { |line|
          cols = line.rstrip.split(/\t/)
          acc, chr, str, id = cols.values_at(1, 2, 3, 12)
          cs = cols[6].to_i
          ce = cols[7].to_i
          lefts = cols[9].split(/,/)
          rights = cols[10].split(/,/)
          prop = "gene_id \"#{id}\"; transcript_id \"#{acc}\""
          prop += "; gene_name \"#{acc2sym[acc]}\"" if acc2sym.has_key?(acc)
          lefts.each_index { |i|
            next if lefts[i].nil?
            l = lefts[i].to_i
            r = rights[i].to_i
            fp.puts [chr, 'ensGene', 'exon', l+1, r, 0, str, '.', prop].join("\t")
          }
          if cs != ce
            fp.puts [chr, 'ensGene', 'start_codon', (str == '+' ? cs+1 : ce-2), (str == '+' ? cs+3 : ce), 0, str, '.', prop].join("\t")
            lefts.each_index { |i|
              next if lefts[i].nil?
              l = lefts[i].to_i
              next if ce-1 < l
              r = rights[i].to_i
              next if r-1 < cs
              fp.puts [chr, 'ensGene', 'CDS', (cs < l ? l : cs), (r < ce ? r : ce), 0, str, '.', prop].join("\t")
            }
          end
        }
      }

    end

  end
end
