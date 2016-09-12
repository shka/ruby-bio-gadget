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
      
      method_option :prefix_coreutils,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU coreutils',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : ''

      def nr
        pseq = ''
        prefix = options.prefix_coreutils
        if options.degenerated_mode
          BioGadget.nr_deg("#{prefix}sort -t '\t' --parallel=`#{prefix}nproc` -r -k2,2 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}")
        else
          BioGadget.nr_std("#{prefix}sort -t '\t' --parallel=`#{prefix}nproc` -r -k2,4 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}")
        end
      end
      
    end
  end
end

