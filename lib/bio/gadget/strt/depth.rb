require 'parallel'

module Bio
  class Gadget
    class STRT < Bio::Gadget

      desc 'depth FQGZ [FQGZ ...]',
           'Count nonredundant reads according to the sequencing depths'

      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX

      method_option *OPT_LENGTH_BARCODE
      method_option *OPT_LENGTH_GAP
      method_option *OPT_LENGTH_MINIMUM
      method_option *OPT_LENGTH_UMI
      method_option *OPT_LOW_QUALITIES

      method_option :tss,
                    default: false,
                    desc: 'Check number of TSSs, instead of STRT reads',
                    type: :boolean

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
        tmpfiles = Array.new(fqgzs.length) do |i|
          getTmpname('strt.depth', 'fq1l')
        end
        tsscmd =
          options.tss ? "fq1l mt5 --minimum-length=#{mLen} #{match}+ | #{cPfx0}cut -f 2 | #{sortCommand(options)} -u |" : ''
        indexes = Array.new(fqgzs.length) { |i| i }
        Parallel.each(indexes, in_threads: options.parallel) do |i|
          system "gunzip -c #{fqgzs[i]} | fq1l convert #{cPfx} > #{tmpfiles[i]}"
        end
        
        1.upto(12).each do |draw|
          fifo = mkfifo('strt.depth', 'fq1l')
          fp0 = open("| #{cPfx0}wc -l #{fifo}")
          fp1 = open(<<CMD
| LC_ALL=C cat #{tmpfiles.join(' ')} \
| fq1l to #{draw} #{12-draw} \
| #{cPfx0}tee #{fifo} \
| fq1l nr #{bSize} #{cPfx} #{par} \
| fq1l m5 #{gPfx} #{match} \
| fq1l m5 #{gPfx} --invert-match '[^\\t]*N' \
| fq1l qt3 --low-qualities='#{options.low_qualities}' --minimum-length=#{pLen} \
| fq1l pt3 --primer=AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG --minimum-length=#{pLen} #{cPfx} #{gPfx} \
| fq1l nr #{bSize} --degenerated-mode #{cPfx} #{par} \
| #{tsscmd} #{cPfx0}wc -l
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
