module Bio
  class Gadget
    class Fq1l < Bio::Gadget
      
      desc 'qt3', 'Quality trimming, 3\'-end'
      
      method_option :low_qualities,
                    banner: 'CHARACTERS',
                    default: '!"#',
                    desc: 'Low quality characters',
                    type: :string
      
      method_option :minimum_length,
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
