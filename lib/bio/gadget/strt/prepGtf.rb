module Bio
  module Gadget
    class STRT < Thor

      desc 'prepGtf xRef genePred', 'Preprocess annotation and create gtf from gzipped tables by UCSC'

      method_option :promoter_size,
                    banner: 'NT',
                    default: 500,
                    desc: 'Length of promoter',
                    type: :numeric
      
      method_option :prefix_coreutils,
                    banner: 'PREFIX',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : '',
                    desc: 'A prefix character for GNU coreutils',
                    type: :string

      def prepGtf(xRef, genePred)

        cPrefix = options.prefix_coreutils
        pLength = options.promoter_size

        tid2sym = Hash.new
        open("| gunzip -c #{xRef} | #{cPrefix}cut -f 1,5").each do |line|
          tid, sym = line.rstrip.split(/\t/)
          tid2sym[tid] = sym
        end

        tmpfile = Bio::Gadgets::getTmpname('strt.prepGtf', 'gtf')
        system "gunzip -c #{genePred} | #{cPrefix}cut -f 1-10 | genePredToGtf -utr file stdin #{tmpfile}"

        reSep = /\t/
        reTid = /transcript_id "([^"]+)/
        tid2exon1 = Hash.new
        tmpfile2 = Bio::Gadgets::getTmpname('strt.prepGtf', 'gtf')
        fp = open(tmpfile2, 'w')
        open(tmpfile).each do |line|
          cols = line.rstrip.split(reSep)
          tid = reTid.match(cols[8]).to_a[1]
          newAttr = tid2sym.key?(tid) ? "#{cols[8]} gene_symbol \"#{tid2sym[tid]}\"; " : cols[8]
          if cols[2] == 'transcript'
            fp.puts [cols[0], cols[1], 'promoter',
                     cols[6] == '+' ? cols[3].to_i-pLength : cols[4].to_i+1,
                     cols[6] == '+' ? cols[3].to_i-1 : cols[4].to_i+pLength,
                     cols[5], cols[6], cols[7], newAttr].join("\t")
          elsif cols[2] == 'exon'
            if (cols[6] == '+' && !tid2exon1.key?(tid)) || cols[6] == '-'
              tid2exon1[tid] = (cols[0..1]+['1stExon']+cols[3..7]+[newAttr]).join("\t")
            end
          end
          fp.puts (cols[0..7] + [newAttr]).join("\t")
        end
        tid2exon1.each_value do |line|
          fp.puts line
        end
        fp.close

        system "#{cPrefix}sort -t '\t' -k 1,1 -k 4,4n -k 5,5n -k 3,3 -k 9,9 #{tmpfile2}"
      end
      
    end
  end
end
