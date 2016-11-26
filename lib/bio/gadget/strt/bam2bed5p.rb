module Bio
  module Gadget
    class STRT < Thor

      desc 'bam2bed5p BAM', "Convert bam to bed of 5'-end of alignments"

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX
      method_option *Bio::Gadgets::OPT_PARALLEL

      method_option :minimum_quality,
                    default: 14,
                    desc: 'Minimum mapping quality',
                    type: :numeric
      
      def bam2bed5p(bam)

        reSep = /\t/
        
        fp = open("| #{Bio::Gadgets.sortCommand(options)} -t '\t' -k 1,1 -k 2,2n -k 3,3n -k 4,4", 'w')
        open("| samtools view #{bam} -b -q #{options.minimum_quality} | bedtools bamtobed -i stdin").each do |line|
          cols = line.rstrip.split(reSep)
          fp.puts ([cols[0],
                    cols[5] == '+' ? cols[1] : cols[2],
                    cols[5] == '+' ? cols[1].to_i+1 : cols[2].to_i+1] +
                   cols[3..5]).join("\t")
        end
        fp.close
        
      end
      
    end
  end
end
