module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'm5 PATTERN',
           'Select sequences that match the 5\'-end with a given PATTERN'
      
      def m5(pattern)
        exec "#{options.prefix_coreutils}grep -P -e '^[^\\t]+\\t#{pattern}'"
      end
      
    end
  end
end
