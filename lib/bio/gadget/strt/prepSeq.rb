module Bio
  module Gadget
    class STRT < Thor

      desc 'prepSeq MAP BASE FQGZ [FQGZ ...]', 'Preprocess of FQGZs with demuliplex based on MAP and BASE, and also the quality reports'

      method_option *Bio::Gadgets::OPT_BUFFER_SIZE
      method_option *Bio::Gadgets::OPT_PARALLEL
      method_option *Bio::Gadgets::OPT_PREFIX_COREUTILS
      method_option *Bio::Gadgets::OPT_PREFIX_GREP

      method_option *OPT_LENGTH_BARCODE
      method_option *OPT_LENGTH_GAP
      method_option *OPT_LENGTH_MINIMUM
      method_option *OPT_LENGTH_UMI
      method_option *OPT_LOW_QUALITIES
      
      method_option :maximum_distance,
                    default: 1,
                    desc: 'Maximum distance between barcode and sequence',
                    type: :numeric

      method_option :maximum_reads,
                    aliases: '-n',
                    desc: 'Maximum number of raw reads for (debug of) the preprocess',
                    type: :numeric
      
      def prepSeq(map, base, fqgz, *fqgzs)

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
        cPfx0 = self.configure_prepSeq(options)

        fifos = Array.new(6) { Bio::Gadgets.mkfifo('strt.prepSeq', 'fq1l') }
        fifos.each_index do |i|
          Process.fork do
            exec "#{cPfx0}cut -f 2 #{fifos[i]}| ruby -nle 'puts $_.length' | #{cPfx0}sort -n #{par} | #{cPfx0}uniq -c | ruby -nle \"puts \\$_.lstrip.tr(' ',',')\" > #{base}.stat#{i}.csv"
          end
        end
        
        cmd = <<CMD
LC_ALL=C unpigz -c #{fqgz} #{fqgzs.join(' ')} \
| fq1l cnv #{cPfx} #{options.key?('maximum_reads') ? '| '+cPfx0+'head -n '+options.maximum_reads.to_s : ''} \
| #{cPfx0}tee #{fifos[0]} \
| fq1l nr #{bSize} #{cPfx} #{par} \
| #{cPfx0}tee #{fifos[1]} \
| fq1l m5 #{gPfx} #{match} \
| fq1l m5 #{gPfx} --invert-match '[^\\t]*N' \
| #{cPfx0}tee #{fifos[2]} \
| fq1l qt3 --low-qualities='#{options.low_qualities}' --minimum-length=#{pLen} \
| #{cPfx0}tee #{fifos[3]} \
| fq1l pt3 --primer=AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG --minimum-length=#{pLen} #{cPfx} #{gPfx} \
| #{cPfx0}tee #{fifos[4]} \
| fq1l nr #{bSize} --degenerated-mode #{cPfx} #{par} \
| #{cPfx0}tee #{fifos[5]} \
| fq1l bm --begin=#{uLen+1} --end=#{uLen+bLen} --maximum-distance=#{options.maximum_distance} #{bSize} #{cPfx} #{par} #{map} \
| fq1l mt5 --minimum-length=#{mLen} #{match}+ \
| fq1l dmp #{cPfx} #{gPfx} #{par} #{map} #{base} \
| pigz --processes #{options.parallel} -c > #{base}.undef.fq1l.gz
CMD
        system cmd
        Process.waitall
        
      end
      
    end
  end
end
