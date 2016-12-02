module Bio
  class Gadget
    class Strt < Bio::Gadget

      desc 'prepBed FAI GTF [ANNBASE]',
           'Preprocess annotation files to create bed files at ANNBASE for STRT'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      
      def prepBed(fai, gtf, base='./')

        cPrefix = options.coreutils_prefix
        gPrefix = options.grep_prefix
        reSym = /gene_symbol "([^"]+)/
        reSep = /\t/
        sortBed = "#{cPrefix}sort -k 1,1 -k 2,2n -k 3,3n -k 4,4"
        
        isCoding = Hash.new
        open("| #{gPrefix}grep '\tCDS\t' #{gtf} | #{cPrefix}cut -f 9").each do |line|
          isCoding[reSym.match(line).to_a[1]] = ''
        end

        # create bed for spike, whole
        fp = open("| #{sortBed} > #{base}spike_whole.bed", 'w')
        open("| #{gPrefix}grep ^RNA_SPIKE_ #{fai} | #{cPrefix}cut -f 1,2").each do |line|
          acc, len = line.rstrip.split(reSep)
          fp.puts [acc, 0, len, acc, 0, '+'].join("\t")
        end
        fp.close
        
        # create bed for spike, 5'-end
        fp = open("| #{sortBed} > #{base}spike_5end.bed", 'w')
        open("| #{gPrefix}grep ^RNA_SPIKE_ #{fai} | #{cPrefix}cut -f 1").each do |line|
          acc = line.rstrip
          fp.puts [acc, 0, 50, acc, 0, '+'].join("\t")
        end
        fp.close

        # create bed for coding, whole
        fp = open("| #{sortBed} -u > #{base}coding_whole.bed", 'w')
        open("| #{gPrefix}grep -E '\t(promoter|exon)\t' #{gtf}").each do |line|
          cols = line.rstrip.split(reSep)
          sym = reSym.match(cols[8]).to_a[1]
          if isCoding.key?(sym)
            fp.puts [cols[0], cols[3].to_i-1, cols[4], sym, 0, cols[6]].join("\t")
          end
        end
        fp.close
        
        # create bed for coding, 5end
        fp = open("| #{sortBed} -u > #{base}coding_5end.bed", 'w')
        open("| #{gPrefix}grep -E '\t(promoter|5UTR)\t' #{gtf}").each do |line|
          cols = line.rstrip.split(reSep)
          sym = reSym.match(cols[8]).to_a[1]
          if isCoding.key?(sym)
            fp.puts [cols[0], cols[3].to_i-1, cols[4], sym, 0, cols[6]].join("\t")
          end
        end
        fp.close
        
      end

    end
  end
end
