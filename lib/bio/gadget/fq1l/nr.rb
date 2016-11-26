require 'io/wait'

module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX

      method_option :degenerated_mode,
                    default: false,
                    desc: 'Exclude redundant and shorter sequece',
                    type: :boolean
      
      method_option :parallel,
                    aliases: '-p',
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
                    
      def nr
        exit unless STDIN.wait
        #
        pseq = ''
        if options.degenerated_mode
          BioGadget.nr_deg("#{options.coreutils_prefix}sort -t '\t' --parallel=#{options.parallel} -r -k2,2 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}")
        else
          BioGadget.nr_std("#{options.coreutils_prefix}sort -t '\t' --parallel=#{options.parallel} -r -k2,4 #{options.key?('buffer_size') ? '-S '+options.buffer_size : ''}")
        end
      end
      
    end
  end
end

