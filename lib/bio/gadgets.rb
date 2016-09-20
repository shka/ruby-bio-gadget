require 'mkfifo'
require 'tempfile'
require 'thor'
require 'bio/gadget/fq1l'
require 'bio/gadget/strt'

module Bio
  class Gadgets < Thor
    
    def self.getTmpname(prefix, suffix, cleanup=true)
      tmpname = Dir::Tmpname.create(["rbg.#{prefix}.", ".#{suffix}"]) {  }
      if cleanup
        at_exit {
          begin
            File.unlink(tmpname) if FileTest.exist?(tmpname)
          rescue Errno::ENOENT
          end
        }
      end
      tmpname
    end

    def self.mkfifo(prefix, suffix, cleanup=true)
      fifo = self.getTmpname("#{prefix}.fifo", suffix, cleanup)
      File.mkfifo(fifo)
      fifo
    end
    
  end
end

require 'bio/gadget/bio_gadget'
