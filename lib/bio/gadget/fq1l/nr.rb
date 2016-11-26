require 'io/wait'

module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX
      method_option *Bio::Gadgets::OPT_PARALLEL

      method_option :degenerated_mode,
                    default: false,
                    desc: 'Exclude redundant and shorter sequece',
                    type: :boolean
      
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

