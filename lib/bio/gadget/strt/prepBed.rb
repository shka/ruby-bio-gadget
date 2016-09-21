module Bio
  module Gadget
    class STRT < Thor

      desc 'prepBed FAI [BASE]',
           'Preprocess annotation files to create bed files at BASE for STRT'

      method_option :prefix_coreutils,
                    banner: 'PREFIX',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : '',
                    desc: 'A prefix character for GNU coreutils',
                    type: :string
      
      def prepBed(fai, base='./')

        cPrefix = options.prefix_coreutils
        reSep = /\t/
        
        # create bed for spike, whole
        fp = open("| #{cPrefix}sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 > #{base}spike_whole.bed", 'w')
        open("| grep ^RNA_SPIKE_ #{fai} | #{cPrefix}cut -f 1,2").each do |line|
          acc, len = line.rstrip.split(reSep)
          fp.puts [acc, 0, len, acc, 0, '+'].join("\t")
        end
        fp.close
        
        # create bed for spike, 5'-end
        fp = open("| #{cPrefix}sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 > #{base}spike_5end.bed", 'w')
        open("| grep ^RNA_SPIKE_ #{fai} | #{cPrefix}cut -f 1").each do |line|
          acc = line.rstrip
          fp.puts [acc, 0, 50, acc, 0, '+'].join("\t")
        end
        fp.close
        
      end

    end
  end
end
