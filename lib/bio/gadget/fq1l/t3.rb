require 'mkfifo'
require 'tempfile'

module Bio
  module Gadget
    class Fq1l < Thor

      ## Change to 't3 TRIMMED SEQ1 SEQ2 ...', which internally supports 'grep -e SEQ1 -e SEQ2 ...'
      desc 't3 SEQ TRIMMED', 'Trim sequences that match the 3\'-end with SEQ'
      
      def t3(seq, trimmed)
        fifo = Dir::Tmpname.create(['rbg.fq1l.t3.', '.fq1l']) {  }
        File.mkfifo(fifo)
        prefix = options.prefix_coreutils 
        pid = Kernel.fork do
          exec "#{prefix}tee #{fifo} | #{prefix}grep -v '#{seq}\t+'"
        end
        range = 0..(-seq.length-1)
        out = open(trimmed, 'w')
        open("| #{prefix}grep '#{seq}\t+' #{fifo}").each do |line|
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
