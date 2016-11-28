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
      
      def cat_command(options)
        "#{options.coreutils_prefix}cat"
      end
      
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
      
      def pipeline(parallel, *commands)
        stats = Array.new
        tmpfiles = Array.new
        begin
          while commands.size > 0
            cmds = commands.shift(parallel)
            tmpin = tmpfiles[0]
            cmds[0] = cmds[0] + " < #{tmpin}" unless tmpin.nil?
            tmpfiles << tmpout = get_temporary_path('pipeline', 'tmp', false)
            cmds[-1] = cmds[-1] + " > #{tmpout}" if commands.size > 0
            tmpstats = Open3.pipeline(*cmds)
            stats.concat(tmpstats)
            tmpstats.each do |tmpstat|
              commands = nil unless tmpstat.success?
            end
            unless commands.nil?
              File.unlink(tmpin) unless tmpin.nil?
              tmpfiles.shift if tmpfiles.size > 1
            else
              break
            end
          end
        ensure
          tmpfiles.each do |tmpfile|
            File.unlink(tmpfile) if File.exist?(tmpfile)
          end
        end
        stats
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
