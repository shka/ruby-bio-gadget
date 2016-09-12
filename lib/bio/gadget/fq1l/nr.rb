module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'

      method_option :buffer_size,
                    aliases: '-S',
                    banner: 'SIZE',
                    desc: 'Use SIZE for main memory buffer',
                    type: :string

      method_option :degenerated_mode,
                    aliases: '-d',
                    default: false,
                    desc: 'Exclude redundant and shorter sequece',
                    type: :boolean
      
      method_option :prefix_coreutils,
                    banner: 'PREFIX',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : '',
                    desc: 'A prefix character for GNU coreutils',
                    type: :string

      method_option :parallel,
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
                    
      def nr
        pseq = ''
        if options.degenerated_mode
          BioGadget.nr_deg("#{options.prefix_coreutils}sort -t '\t' --parallel=#{options.parallel} -r -k2,2 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}")
        else
          BioGadget.nr_std("#{options.prefix_coreutils}sort -t '\t' --parallel=#{options.parallel} -r -k2,4 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}")
        end
      end
      
    end
  end
end

