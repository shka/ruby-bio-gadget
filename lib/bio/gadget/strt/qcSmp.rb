require 'parallel'

module Bio
  module Gadget
    class STRT < Thor

      desc 'qcSmp BAM 5pBED NAME [ANNBASE]', 'Measure for sample quality check'

      method_option :parallel,
                    aliases: '-p',
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
      
      def qcSmp(bam, bed, name, base='./')

        cPrefix = options.prefix_coreutils

        bedtool = 'bedtools intersect -s -nonamecheck -a stdin -b'
        count = "| #{cPrefix}cut -f 4 | #{cPrefix}sort -u | #{cPrefix}wc -l | ruby -nle 'puts $_.lstrip'"
        grep = "#{options.prefix_grep}grep"
        
        commands =
          [ '',
            "| #{grep} ^RIBO_",
            "| #{bedtool} #{base}spike_whole.bed",
            "| #{bedtool} #{base}spike_5end.bed",
            "| #{grep} -v ^RIBO_ | #{grep} -v ^RNA_SPIKE_",
            "| #{bedtool} #{base}coding_whole.bed",
            "| #{bedtool} #{base}coding_5end.bed" ].map { |cmd| "gunzip -c #{bed} #{cmd} #{count}" }
        commands.unshift "samtools view #{bam} | #{cPrefix}cut -f 1 | #{cPrefix}sort -u | #{cPrefix}wc -l | ruby -nle 'puts $_.lstrip'"
        
        tmpfiles = Array.new(commands.length) { Bio::Gadgets.getTmpname('strt.qcSmp', 'txt') }

        Parallel.map(0.upto(commands.length-1),
                     in_threads: options.parallel) do |i|
          system "#{commands[i]} > #{tmpfiles[i]}"
        end

        fp = open("| paste #{tmpfiles.join(' ')}")
        values = fp.gets.rstrip.split(/\t/).map { |n| n.to_f }
        fp.close
        
        puts (
          [ name ] + values +
          [ values[1] / values[0], # mapping rate
            values[2] / values[1], # ribosomal rate
            values[3] / values[1], # spike rate
            values[5] / values[1], # polyA+ RNA rate
            values[5] / values[3], # relative polyA+ RNA amount
            values[6] / values[1], # mRNA rate
            values[4] / values[3], # 5'-end capture rate among spike-ins
            values[7] / values[6]  # 5'-end capture rate among mRNAs
          ]
        ).join("\t")
        
      end
      
    end
  end
end
