require 'mkfifo'
require 'tempfile'

module Bio
  module Gadget
    class Fq1l < Thor
      
      desc 'at3 ADAPTOR', 'Adaptor trimming, 3\'-end'

      def at3(adaptor)
        seqs = Hash.new
        max = adaptor.length-1
        tmp = Hash.new
        for i in 0..max do
          for j in i..max do
            seq = adaptor[i..j]
            unless tmp.key?(seq)
              len = seq.length
              seqs[len] = Array.new unless seqs.key?(len)
              seqs[len] << seq
              tmp[seq] = true
            end
          end
        end
        #
        keys = seqs.keys.sort
        tmps = Array.new(keys.length) do
          Dir::Tmpname.create(['rbg.fq1l.at3.', '.fq1l']) {  }
        end
        pipes = Array.new(keys.length-1) { IO.pipe }
        pipes = [nil, pipes.flatten.reverse, nil].flatten
        keys.each_index do |i|
          len = keys[i]
          tmp = tmps[i]
          r, w = pipes.shift 2
          Kernel.fork do
            STDIN.reopen r if r
            STDOUT.reopen w if w
            t3(tmp, seqs[len])
          end
          r.close if r
          w.close if w
        end
        Process.waitall
        #
        exec "cat #{tmps.join(' ')}"
      ensure
        tmps.map { |tmp| File.unlink(tmp) } unless tmps.nil?
      end
      
      no_commands do
        
        def t3(trimmed, seqs)
          len = seqs[0].length
          seqs.each do |seq|
            abort "SEQs must be a same length." if seq.length != len
          end
          patterns = (seqs.map { |seq| "-e '#{seq}\t+'" }).join ' '
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
          File.unlink(fifo) unless fifo.nil?
        end
        
      end
      
    end
  end
end
