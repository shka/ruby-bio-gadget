module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'm5 PATTERN',
           'Select sequences that match the 5\'-end with a given PATTERN'

      method_option :invert_match,
                    type: :boolean,
                    aliases: '-v',
                    desc: 'Invert the sense of matching, to select non-matching lines'
      
      def m5(pattern)
        exec "#{options.prefix_coreutils}grep #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t#{pattern}'"
      end
      
    end
  end
end
