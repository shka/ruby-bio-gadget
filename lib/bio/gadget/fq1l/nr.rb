module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'

      method_option :buffer_size,
                    banner: 'SIZE',
                    aliases: '-S',
                    desc: 'Use SIZE for main memory buffer'

      method_option :degenerated_mode,
                    aliases: '-d',
                    type: :boolean,
                    default: false,
                    desc: 'Exclude redundant and shorter sequece'
      
      def nr
        pseq = ''
        prefix = options.prefix_coreutils
        if options.degenerated_mode
          command = <<CMD
| #{prefix}sort -t '\t' --parallel=`#{prefix}nproc` -r -k2,4 \
    #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}
CMD
          open(command).each do |line|
            acc, seq, tmp, qual = line.split /\t/
            if pseq != seq
              puts line.rstrip
              pseq = seq
            end
          end
        else
          command = <<CMD
| #{prefix}sort -t '\t' --parallel=`#{prefix}nproc` -r -k2,2 \
    #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}
CMD
          pseql = 2**(0.size * 8 -2) -1
          open(command).each do |line|
            acc, seq, tmp, qual = line.split /\t/
            if seq.length >= pseq.length || !Regexp.new("^#{seq}").match(pseq).nil?
              puts line.rstrip
              pseq = seq
              pseql = pseq.length
            end
          end
        end
      end
      
    end
  end
end

