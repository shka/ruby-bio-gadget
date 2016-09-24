require 'bio'
require 'mkfifo'
require 'tempfile'
require 'thor'

module Bio
  class Gadgets < Thor
    
    def self.getTmpname(prefix, suffix, cleanup=true)
      tmpname = Dir::Tmpname.create(["rbg.#{prefix}.", ".#{suffix}"]) {  }
      if cleanup
        at_exit { File.unlink(tmpname) if FileTest.exist?(tmpname) }
      end
      tmpname
    end

    def self.mkfifo(prefix, suffix, cleanup=true)
      fifo = self.getTmpname("#{prefix}.fifo", suffix, cleanup)
      File.mkfifo(fifo)
      fifo
    end

    #

    option_buffer_size = [
      :buffer_size, {
        :aliases => '-S',
        :banner => 'SIZE',
        :desc => 'Use SIZE for main memory buffer',
        :type => :string
      }
    ]
    
    option_parallel = [
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

    option_prefix_coreutils = [
      :prefix_coreutils, {
        :banner => 'PREFIX',
        :default => system('which gnproc >/dev/null 2>&1') ? 'g' : '',
        :desc => 'A prefix character for GNU coreutils',
        :type => :string
      }
    ]

    #
    
    desc 'find PATTERN [NAME]',
         'Find fragments matching with regexp PATTERN from FASTA-format STDIN'

    method_option *option_buffer_size
    method_option *option_parallel
    method_option *option_prefix_coreutils

    method_option :ignore_case,
                  default: true,
                  desc: 'Fold lower case to upper case characters',
                  type: :boolean
    
    def find(pattern, name0 = '')

      bSize = options.key?('buffer_size') ? '--buffer-size='+options.buffer_size : ''
      cPrefix = options.prefix_coreutils
      re = Regexp.new("(#{pattern})", options.ignore_case)
      name = name0 == '' ? pattern : name0
      
      pids = Array.new
      tmpfiles = Array.new
      ff = Bio::FlatFile.open(Bio::FastaFormat, STDIN)
      ff.each do |entry|
        tmpfiles << tmpfile = Bio::Gadgets.getTmpname('find', 'bed', false)
        acc = entry.entry_id
        seq = entry.seq
        pids << Process.fork do
          fp = open("| #{cPrefix}sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 #{bSize} > #{tmpfile}", 'w')
          #
          pos = 0
          match = re.match(seq, pos)
          while !match.nil?
            fp.puts [acc, match.begin(1), match.end(1), name, '0', '+'].join("\t")
            pos = match.begin(1)+1
            match = re.match(seq, pos)
          end
          #
          pos = 0
          seq = seq.reverse.tr('acgtACGT', 'tgcaTGCA')
          len = seq.length
          match = re.match(seq, pos)
          while !match.nil?
            fp.puts [acc, len-match.end(1), len-match.begin(1), name, '0', '-'].join("\t")
            pos = match.begin(1)+1
            match = re.match(seq, pos)
          end
          #
          fp.close
        end
        while pids.length == options.parallel
          pids.delete(Process.wait)
        end
      end
      ff.close
      while pids.length > 0
        pids.delete(Process.wait)
      end
      
      system "#{cPrefix}sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 --merge #{bSize} #{tmpfiles.join(' ')}"
      tmpfiles.each do |tmpfile|
        File.unlink(tmpfile) if FileTest.exist?(tmpfile)
      end

    end
  end
end

require 'bio/gadget/fq1l'
require 'bio/gadget/strt'
require 'bio/gadget/bio_gadget'
