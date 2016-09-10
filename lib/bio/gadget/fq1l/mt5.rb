module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'mt5 PATTERN', 'Match trimming, 5\'-end'

      method_option :minimum_length,
                    alias: '-l',
                    banner: 'NT',
                    default: 25,
                    desc: 'Minimum length after trimming',
                    type: :numeric
      
      def mt5(pattern)
        BioGadget.mt5(pattern, options.minimum_length)
      end
      
    end
  end
end
