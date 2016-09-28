module Bio
  module Gadget
    class STRT < Thor

      desc 'depth FQGZ [FQGZ ...]',
           'Count nonredundant reads according to the sequencing depths'

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_PARALLEL
      method_option *Bio::Gadgets::OPT_PREFIX_COREUTILS
      method_option *Bio::Gadgets::OPT_PREFIX_GREP

      method_option *OPT_LENGTH_BARCODE
      method_option *OPT_LENGTH_GAP
      method_option *OPT_LENGTH_MINIMUM
      method_option *OPT_LENGTH_UMI
      method_option *OPT_LOW_QUALITIES

      def depth(fqgz, *fqgzs0)
        
        bSize,
        cPfx,
        gPfx,
        par,
        bLen,
        gLen,
        mLen,
        pLen,
        uLen,
        match,
        cPfx0 = Bio::Gadget::STRT.configure_prepSeq(options)

        fqgzs = [fqgz] + fqgzs0
        reads = Hash.new

        1.upto(12).each do |draw|
          fifo = Bio::Gadgets.mkfifo('strt.depth', 'fq1l')
          fp0 = open("| #{cPfx0}wc -l #{fifo}")
          fp1 = open(<<CMD
| LC_ALL=C unpigz -c #{fqgzs.join(' ')} \
| fq1l cnv #{cPfx} \
| fq1l to #{draw} #{12-draw} \
| #{cPfx0}tee #{fifo} \
| fq1l nr #{bSize} #{cPfx} #{par} \
| fq1l m5 #{gPfx} #{match} \
| fq1l m5 #{gPfx} --invert-match '[^\\t]*N' \
| fq1l qt3 --low-qualities='#{options.low_qualities}' --minimum-length=#{pLen} \
| fq1l pt3 --primer=AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG --minimum-length=#{pLen} #{cPfx} #{gPfx} \
| fq1l nr #{bSize} --degenerated-mode #{cPfx} #{par} \
| #{cPfx0}wc -l
CMD
                    )
          raw = fp0.gets.strip.split(/\s+/)[0]
          fp0.close
          nr = fp1.gets.strip
          fp1.close
          puts [raw, nr].join(',')
        end
        
      end
      
    end
  end
end