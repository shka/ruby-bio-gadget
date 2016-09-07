module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'qt3', 'Quality trimming, 3\'-end'

      method_option :low_qualities,
                    :type => :string,
                    :default => '!"#',
                    :banner => 'CHARACTERS',
                    :desc => 'Low quality characters'
      
      def qt3
        BioGadget.qt3(options.low_qualities)
      end
      
    end
  end
end
