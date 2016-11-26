module Bio
  module Gadget
    class STRT < Thor

      desc 'bam2bed5p BAM', "Convert bam to bed of 5'-end of alignments"

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX

      method_option :minimum_quality,
                    default: 14,
                    desc: 'Minimum mapping quality',
                    type: :numeric
      
      method_option :parallel,
                    aliases: '-p',
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
      
      def bam2bed5p(bam)

        bSize = options.key?('buffer_size') ? '--buffer-size='+options.buffer_size : ''
        parallel = "--parallel=#{options.parallel}"
        reSep = /\t/
        
        fp = open("| #{options.coreutils_prefix}sort #{bSize} #{parallel} -k 1,1 -k 2,2n -k 3,3n -k 4,4", 'w')
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
