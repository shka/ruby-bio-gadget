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
        keys = Array.new
        seqs.keys.sort.reverse.each do |len|
          keys << len
          break if 4**len == seqs[len].length
        end
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
            t3(len, seqs[len], tmp)
          end
          r.close if r
          w.close if w
        end
        Process.waitall
        #
        exec "unpigz -c #{tmps.join(' ')}"
      ensure
        tmps.map { |tmp| File.unlink(tmp) } unless tmps.nil?
      end
      
      no_commands do
        
        def t3(len, seqs, trimmed)
          fifo = Dir::Tmpname.create(['rbg.fq1l.t3.', '.fq1l']) {  }
          File.mkfifo(fifo)
          prefix = options.prefix_coreutils
          if 4**len == seqs.length
            pCmds = "#{prefix}cat > #{fifo}"
            cCmds = "#{prefix}cat #{fifo}"
          else
            patterns = (seqs.map { |seq| "-e '#{seq}\t+'" }).join ' '
            pCmds = "#{prefix}tee #{fifo} | #{prefix}grep -v #{patterns}"
            cCmds = "#{prefix}grep #{patterns} #{fifo}"
          end
          pid = Kernel.fork do
            exec pCmds
          end
          BioGadget.trim3(len, cCmds, trimmed)
          Process.wait(pid)
          File.unlink(fifo) unless fifo.nil?
        end
        
      end
      
    end
  end
end
