require 'mkfifo'
require 'tempfile'
require 'thor'

module Bio
  class Gadget < Thor

    OPT_BUFFER_SIZE = [
      :buffer_size, {
        :aliases => '-S',
        :banner => 'SIZE',
        :desc => 'Use SIZE for main memory buffer',
        :type => :string
      }
    ]
    
    OPT_PARALLEL = [
      :parallel, {
        :aliases => '-p',
        :banner => 'N',
        :default => (
          system('which gnproc >/dev/null 2>&1') ?
            `gnproc`.to_i :
            (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2)
        ),
        :desc => 'Change the number of sorts run concurrently to N',
        :type => :numeric
      }
    ]

    OPT_COREUTILS_PREFIX = [
      :coreutils_prefix, {
        :banner => 'PREFIX',
        :default => system('which gnproc >/dev/null 2>&1') ? 'g' : '',
        :desc => 'A prefix character for GNU coreutils',
        :type => :string
      }
    ]

    OPT_GREP_PREFIX = [
      :grep_prefix, {
        :banner => 'PREFIX',
        :default => system('which ggrep >/dev/null 2>&1') ? 'g' : '',
        :desc => 'A prefix character for GNU grep',
        :type => :string
      }
    ]
    
    #

    no_commands do
      
      def get_temporary_path(prefix, suffix, cleanup=true)
        tmpname = Dir::Tmpname.create(["rbg.#{prefix}.", ".#{suffix}"]) {  }
        if cleanup
          at_exit { File.unlink(tmpname) if FileTest.exist?(tmpname) }
        end
        tmpname
      end
      
      def get_fifo(prefix, suffix, cleanup=true)
        fifo = get_temporary_path("#{prefix}.fifo", suffix, cleanup)
        File.mkfifo(fifo)
        fifo
      end
      
      def grep_command(options)
        "#{options.grep_prefix}grep"
      end
      
      def sort_command(options)
        "#{options.coreutils_prefix}sort --parallel=#{options.parallel}#{options.key?('buffer_size') ? ' --buffer-size='+options.buffer_size+' ' : ''}"
      end
      
      def tee_command(options)
        "#{options.coreutils_prefix}tee"
      end

    end
    
  end
end

require 'bio/gadgets'
require 'bio/gadget/fq1l'
require 'bio/gadget/strt'
require 'bio/gadget/bio_gadget'
