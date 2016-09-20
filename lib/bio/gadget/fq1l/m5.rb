module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'm5 PATTERN',
           'Select sequences that match the 5\'-end with a given PATTERN'

      method_option :invert_match,
                    type: :boolean,
                    desc: 'Invert the sense of matching, to select non-matching lines'
      
      method_option :prefix_grep,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU grep',
                    default: system('which ggrep >/dev/null 2>&1') ? 'g' : ''

      def m5(pattern)
        exec "#{options.prefix_grep}grep #{options.invert_match ? '-v' : ''} -P -e '^[^\\t]+\\t#{pattern}'"
      end
      
    end
  end
end
