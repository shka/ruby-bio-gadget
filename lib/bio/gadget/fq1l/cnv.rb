module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'cnv', 'Convert fastq from 4 lines/read to 1 line/read'
      
      method_option :prefix_coreutils,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU coreutils',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : ''

      def cnv
        exec "#{options.prefix_coreutils}paste - - - -"
      end
      
    end
  end
end
