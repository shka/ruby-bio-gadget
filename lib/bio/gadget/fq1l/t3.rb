require 'mkfifo'
require 'tempfile'

module Bio
  module Gadget
    class Fq1l < Thor

      desc 't3 TRIMMED SEQ1 [SEQ2 ...]', 'Trim sequences that match the 3\'-end with SEQs'
      
      def t3(trimmed, *seqs)
        len = seqs[0].length
        seqs.each do |seq|
          abort "SEQs must be a same length." if seq.length != len
        end
        patterns = (seqs.map { |seq| "-e '#{seq}\t'" }).join ' '
        #
        fifo = Dir::Tmpname.create(['rbg.fq1l.t3.', '.fq1l']) {  }
        File.mkfifo(fifo)
        prefix = options.prefix_coreutils 
        pid = Kernel.fork do
          exec "#{prefix}tee #{fifo} | #{prefix}grep -v #{patterns}"
        end
        #
        range = 0..(-len-1)
        out = open(trimmed, 'w')
        open("| #{prefix}grep #{patterns} #{fifo}").each do |line|
          acc, raw, tmp, qual = line.rstrip.split /\t/
          out.puts [acc, raw[range], tmp, qual[range]].join("\t")
        end
        out.close
        Process.wait(pid)
      ensure
        File.unlink(fifo)
      end
      
    end
  end
end
