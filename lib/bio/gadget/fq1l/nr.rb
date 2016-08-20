module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'
      
      def nr
        pseq = ''
        prefix = options.prefix_coreutils
        open("| #{prefix}sort --parallel=`#{prefix}nproc` -r -k2,4").each do |line|
          acc, seq, tmp, qual = line.split /\t/
          if pseq != seq
            puts line.strip
            pseq = seq
          end
        end
      end
      
    end
  end
end

