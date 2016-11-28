module Bio
  class Gadget
    class STRT < Bio::Gadget

      desc 'prepSeq MAP BASE FQGZ [FQGZ ...]', 'Preprocess of FQGZs with demuliplex based on MAP and BASE, and also the quality reports'

      method_option *OPT_BUFFER_SIZE
      method_option *OPT_PARALLEL
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX

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
      
      def prepSeq(map, base, fqgz, *fqgzs0)

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

        fqgzs = [fqgz] + fqgzs0
        tmpfiles = Array.new(fqgzs.length) do |i|
          get_temporary_path('strt.depth', 'fq1l')
        end
        indexes = Array.new(fqgzs.length) { |i| i }
        Parallel.each(indexes, in_threads: options.parallel) do |i|
          system "gunzip -c #{fqgzs[i]} | fq1l convert #{cPfx} > #{tmpfiles[i]}"
        end
        
        fifos = Array.new(6) { get_fifo('strt.prepSeq', 'fq1l') }
        fifos.each_index do |i|
          Process.fork do
            exec "#{cPfx0}cut -f 2 #{fifos[i]}| ruby -nle 'puts $_.length' | #{cPfx0}sort -n #{par} | #{cPfx0}uniq -c | ruby -nle \"puts \\$_.lstrip.tr(' ',',')\" > #{base}.stat#{i}.csv"
          end
        end
        
        cmd = <<CMD
LC_ALL=C cat #{tmpfiles.join(' ')} #{options.key?('maximum_reads') ? '| '+cPfx0+'head -n '+options.maximum_reads.to_s : ''} \
| #{tee_command(options)} #{fifos[0]} \
| fq1l nr #{bSize} #{cPfx} #{par} \
| #{tee_command(options)} #{fifos[1]} \
| fq1l m5 #{gPfx} #{match} \
| fq1l m5 #{gPfx} --invert-match '[^\\t]*N' \
| #{tee_command(options)} #{fifos[2]} \
| fq1l qt3 --low-qualities='#{options.low_qualities}' --minimum-length=#{pLen} \
| #{tee_command(options)} #{fifos[3]} \
| fq1l pt3 --primer=AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG --minimum-length=#{pLen} #{cPfx} #{gPfx} \
| #{tee_command(options)} #{fifos[4]} \
| fq1l nr #{bSize} --degenerated-mode #{cPfx} #{par} \
| #{tee_command(options)} #{fifos[5]} \
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
