require 'bio-faster'
require 'parallel'

module Bio
  class Gadget < Thor
    namespace :bio

    desc 'dedup', 'deduplicate fastq (via STDIN)'
    def dedup

      p1in, p1out = IO.pipe

      fork {
        p1in.close
        $stdout.reopen(p1out)
        open("| sort -k 1 -r -S #{sprintf('%2d', 100/(Parallel.processor_count+1))}% -T $TMPDIR | cut -f 2- | uniq -f 2", 'w') { |fp|
          Bio::Faster.new(:stdin).each_record(:quality => :raw) do |seqid, seq, qvs|
            fp.puts "#{seq}#{qvs}\t#{seqid}\t#{qvs}\t#{seq}"
          end
        }
      }

      p1out.close

      p1in.each_line { |line|
        seqid, qvs, seq = line.rstrip.split
        puts "@#{seqid}\n#{seq}\n+\n#{qvs}"
      }

    end

  end
end
