require 'io/wait'

module Bio
  class Gadget
    class Fq1l < Bio::Gadget
      
      desc 'nr', 'Extract non-redundant best-quality sequences, esp. for UMI'

      method_option *OPT_BUFFER_SIZE
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_PARALLEL

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

