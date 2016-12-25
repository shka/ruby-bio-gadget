require 'mkfifo'
require 'open3'
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
    
    OPT_DOWNLOAD = [ :download, {
                       :banner => 'BEHAVIOR',
                       :default => 'yes',
                       :desc => 'Download and process, no download or only',
                       :enum => ['yes', 'no', 'only'] } ]
    
    OPT_PARALLEL = [
      :parallel, {
        :banner => 'N',
        :default => (
          system('which gnproc >/dev/null 2>&1') ?
            `gnproc`.to_i :
            (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2)
        ),
        :desc => 'Change the number of sorts run concurrently',
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

      def coreutils_prefix_option(options)
         options.key?(:coreutils_prefix) ? " --coreutils-prefix=#{options.coreutils_prefix}" : ''
      end
      
      def cut_command(options)
        "#{options.coreutils_prefix}cut"
      end

      def download_file(url, path)
        system "curl -R -f -s -S -o #{path} '#{url}'" or exit $?.exitstatus
      end

      def fold_command(options)
        "#{options.coreutils_prefix}fold"
      end
      
      def fq1l_convert_command(options)
        "fq1l convert#{coreutils_prefix_option(options)}"
      end
      
      def fq1l_count_command(options)
        "fq1l count#{coreutils_prefix_option(options)}#{parallel_option(options)}"
      end

      def fq1l_sort_command(options)
        "fq1l sort#{coreutils_prefix_option(options)}#{parallel_option(options)}"
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

      def grep_prefix_option(options)
         options.key?(:grep_prefix) ? " --grep-prefix=#{options.grep_prefix}" : ''
      end
      
      def head_command(options)
        "#{options.coreutils_prefix}head"
      end

      def parallel_option(options)
         options.key?(:parallel) ? " --parallel=#{options.parallel}" : ''
      end
      
      def pipeline(*cmds)
        stats = Open3.pipeline(*cmds)
        stats.each_index do |i|
          raise "Fail at process #{i}; #{stats[i]}; #{cmds[i]}" unless stats[i].success?
        end
      end
      
      def sort_command(options)
        "#{options.coreutils_prefix}sort#{options.key?(:parallel) ? ' --parallel='+options.parallel.to_s : ''}#{options.key?(:buffer_size) ? ' --buffer-size='+options.buffer_size+' ' : ''} --compress-program=pigz"
      end
      
      def tee_command(options)
        "#{options.coreutils_prefix}tee"
      end

      def uniq_command(options)
        "#{options.coreutils_prefix}uniq"
      end

      def unlink_files(files)
        files.each do |file|
          File.unlink(file) if File.exist?(file)
        end
      end
      
      def wc_command(options)
        "#{options.coreutils_prefix}wc"
      end

    end
    
  end
end

require 'bio/gadgets'
require 'bio/gadget/fq1l'
require 'bio/gadget/strt'
require 'bio/gadget/bio_gadget'
