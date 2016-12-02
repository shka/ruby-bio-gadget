require 'parallel'

module Bio
  class Gadget
    class Strt < Bio::Gadget

      desc 'qcSmp BAM 5pBED NAME [ANNBASE]', 'Measure for sample quality check'

      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_GREP_PREFIX
      method_option *OPT_PARALLEL

      def qcSmp(bam, bed, name, base='./')

        cPrefix = options.coreutils_prefix

        bedtool = 'bedtools intersect -s -nonamecheck -a stdin -b'
        count = "| #{cPrefix}cut -f 4 | #{cPrefix}sort -u | #{cPrefix}wc -l | ruby -nle 'puts $_.lstrip'"
        grep = "#{options.grep_prefix}grep"
        
        commands =
          [ '',
            "| #{grep} ^RIBO_",
            "| #{bedtool} #{base}spike_whole.bed",
            "| #{bedtool} #{base}spike_5end.bed",
            "| #{grep} -v ^RIBO_ | #{grep} -v ^RNA_SPIKE_",
            "| #{bedtool} #{base}coding_whole.bed",
            "| #{bedtool} #{base}coding_5end.bed" ].map { |cmd| "gunzip -c #{bed} #{cmd} #{count}" }
        commands.unshift "samtools view #{bam} | #{cPrefix}cut -f 1 | #{cPrefix}sort -u | #{cPrefix}wc -l | ruby -nle 'puts $_.lstrip'"
        
        tmpfiles = Array.new(commands.length) { get_temporary_path('strt.qcSmp', 'txt') }

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
        ).join(',')
        
      end
      
    end
  end
end
