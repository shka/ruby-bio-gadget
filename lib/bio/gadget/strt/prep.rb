module Bio
  module Gadget
    class STRT < Thor

      desc 'prep MAP BASE FQGZ [FQGZ ...]', 'Preprocess of FQGZs with demuliplex based on MAP and BASE, and also the quality reports'

      method_option :buffer_size,
                    aliases: '-S',
                    banner: 'SIZE',
                    desc: 'Use SIZE for main memory buffer',
                    type: :string

      method_option :length_barcode,
                    banner: 'NT',
                    default: 6,
                    desc: 'Length of barcode',
                    type: :numeric
      
      method_option :length_gap,
                    banner: 'NT',
                    default: 3,
                    desc: 'Length of gap (polyG)',
                    type: :numeric
      
      method_option :length_minimum,
                    banner: 'NT',
                    default: 25,
                    desc: 'Minimum length after preprocess',
                    type: :numeric
      
      method_option :length_umi,
                    banner: 'NT',
                    default: 6,
                    desc: 'Length of UMI',
                    type: :numeric

      method_option :low_qualities,
                    banner: 'CHARACTERS',
                    default: '!"#',
                    desc: 'Low quality characters',
                    type: :string
      
      method_option :maximum_distance,
                    default: 1,
                    desc: 'Maximum distance between barcode and sequence',
                    type: :numeric

      method_option :maximum_reads,
                    aliases: '-n',
                    desc: 'Maximum number of raw reads for (debug of) the preprocess',
                    type: :numeric
      
      method_option :parallel,
                    alias: '-p',
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
      
      method_option :prefix_coreutils,
                    banner: 'PREFIX',
                    default: system('which gnproc >/dev/null 2>&1') ? 'g' : '',
                    desc: 'A prefix character for GNU coreutils',
                    type: :string

      method_option :prefix_grep,
                    type: :string,
                    banner: 'PREFIX',
                    desc: 'A prefix character for GNU grep',
                    default: system('which ggrep >/dev/null 2>&1') ? 'g' : ''

      def prep(map, base, fqgz, *fqgzs)
        bSize = options.key?('buffer_size') ? '--buffer-size='+options.buffer_size : ''
        parallel = "--parallel=#{options.parallel}"
        cPrefix0 = options.prefix_coreutils
        cPrefix = options.prefix_coreutils == '' ? '' : "--prefix-coreutils=#{options.prefix_coreutils}"
        gPrefix = options.prefix_grep == '' ? '' : "--prefix-grep=#{options.prefix_grep}"
        uLength = options.length_umi
        bLength = options.length_barcode
        gLength = options.length_gap
        mLength = options.length_minimum
        mLength0 = mLength + uLength + bLength + gLength
        match = "#{'.' * uLength}#{'.' * bLength}#{'G' * (gLength-1)}"
        fifos = Array.new(6) { Bio::Gadgets.mkfifo('strt', 'fq1l') }
        fifos.each_index do |i|
          Process.fork do
            exec "#{cPrefix0}cut -f 2 #{fifos[i]}| ruby -nle 'puts $_.length' | #{cPrefix0}sort -n #{parallel} | #{cPrefix0}uniq -c | ruby -nle \"puts \\$_.lstrip.tr(' ',',')\" > #{base}.stat#{i}.csv"
          end
        end
        cmd = <<CMD
LC_ALL=C \
gunzip -c #{fqgz} #{fqgzs.join(' ')} \
| fq1l cnv #{cPrefix} \
 #{options.key?('maximum_reads') ? '| '+cPrefix0+'head -n '+options.maximum_reads.to_s : ''} \
| #{cPrefix0}tee #{fifos[0]} \
| fq1l nr #{bSize} #{cPrefix} #{parallel} \
| #{cPrefix0}tee #{fifos[1]} \
| fq1l m5 #{gPrefix} #{match}\
| fq1l m5 #{gPrefix} --invert-match '[^\\t]*N' \
| #{cPrefix0}tee #{fifos[2]} \
| fq1l qt3 --low-qualities='#{options.low_qualities}' --minimum-length=#{mLength0} \
| #{cPrefix0}tee #{fifos[3]} \
| fq1l pt3 --primer=AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG --minimum-length=#{mLength0} #{cPrefix} #{gPrefix} \
| #{cPrefix0}tee #{fifos[4]} \
| fq1l nr #{bSize} --degenerated-mode #{cPrefix} #{parallel} \
| #{cPrefix0}tee #{fifos[5]} \
| fq1l bm --begin=#{uLength+1} --end=#{uLength+bLength} --maximum-distance=#{options.maximum_distance} #{bSize} #{cPrefix} #{parallel} #{map} \
| fq1l mt5 --minimum-length=#{mLength} #{match}+ \
| fq1l dmp #{cPrefix} #{gPrefix} #{parallel} #{map} #{base} \
| pigz --processes #{options.parallel} -c > #{base}.undef.fq1l.gz
CMD
        STDERR.puts cmd
        system cmd
        Process.waitall
      end
      
    end
  end
end

  # time \
  #   unpigz -c \
  #     /proj/b2014069/nobackup/private/STRTprep3.NATHALIEB/src/run42/lane2_NoIndex_L002_R1_001.fastq.gz \
  #     /proj/b2014069/nobackup/private/STRTprep3.NATHALIEB/src/run42/lane3_NoIndex_L003_R1_001.fastq.gz \
  #     /proj/b2014069/nobackup/private/STRTprep3.NATHALIEB/src/run42/lane4_NoIndex_L004_R1_001.fastq.gz \
  # | bundle exec exe/fq1l cnv \
  # | bundle exec exe/fq1l nr -S 12% --parallel=\$SLURM_NNODES \
  # | bundle exec exe/fq1l m5 ............GG \
  # | bundle exec exe/fq1l m5 -v '[^\t]*N' \
  # | bundle exec exe/fq1l qt3 \
  # | bundle exec exe/fq1l pt3 \
  # | bundle exec exe/fq1l nr --degenerated-mode -S 12% --parallel=\$SLURM_NNODES \
  # | bundle exec exe/fq1l bm -S 12% --parallel=\$SLURM_NNODES barcodes.csv \
  # | bundle exec exe/fq1l mt5 ............GG+ \
  # | bundle exec exe/fq1l dmp barcodes.csv ../seq/NATHALIE2B --parallel=\$SLURM_NNODES \
  # | pigz --processes \$SLURM_NNODES -c > ../seq/NATHALIE2B.undef.fq1l.gz
  # EOF
