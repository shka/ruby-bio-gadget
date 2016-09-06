require 'mkfifo'
require 'tempfile'
require 'thor'
require 'bio/gadget/fq1l'

module Bio
  class Gadgets < Thor
    
    VERSION = '0.5.0'
    
    def self.getTmpname(prefix, suffix, cleanup=true)
      tmpname = Dir::Tmpname.create(["rbg.#{prefix}.", ".#{suffix}"]) {  }
      if cleanup
        at_exit { File.unlink(tmpname) if FileTest.exist?(tmpname) || FileTest.pipe?(tmpname) }
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
