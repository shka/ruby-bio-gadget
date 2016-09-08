module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'qt3', 'Quality trimming, 3\'-end'
      
      method_option :low_qualities,
                    aliases: '-q',
                    banner: 'CHARACTERS',
                    default: '!"#',
                    desc: 'Low quality characters',
                    type: :string
      
      method_option :minimum_length,
                    aliases: '-l',
                    banner: 'NT',
                    default: 40,
                    desc: 'Minimum length after trimming',
                    type: :numeric
      
      def qt3
        BioGadget.qt3(options.low_qualities, options.minimum_length)
      end
      
    end
  end
end
