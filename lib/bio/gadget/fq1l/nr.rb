module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'
      
      method_option :buffer_size,
                    banner: 'SIZE',
                    aliases: '-S',
                    desc: 'Use SIZE for main memory buffer'
      
      def nr
        pseq = ''
        prefix = options.prefix_coreutils
        command = <<CMD
| #{prefix}sort --parallel=`#{prefix}nproc` -r -k2,4 \
    #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}
CMD
        open(command).each do |line|
          acc, seq, tmp, qual = line.split /\t/
          if pseq != seq
            puts line.rstrip
            pseq = seq
          end
        end
      end
      
    end
  end
end

