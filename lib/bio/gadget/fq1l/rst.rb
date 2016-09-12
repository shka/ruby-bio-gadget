module Bio
  module Gadget
    class Fq1l < Thor

      desc 'rst', 'Restore fastq from 1 line/read to 4 lines/read'
      
      method_option :prefix_coreutils,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU coreutils',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : ''

      def rst
        exec "#{options.prefix_coreutils}tr \"\\t\" \"\\n\""
      end
      
    end
  end
end
