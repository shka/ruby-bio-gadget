require 'bio'
require 'parallel'

module Bio
  class Gadgets < Bio::Gadget
    
    desc 'find PATTERN [NAME]',
         'Find fragments matching with regexp PATTERN from FASTA-format STDIN'
    
    method_option *OPT_BUFFER_SIZE
    method_option *OPT_PARALLEL
    method_option *OPT_COREUTILS_PREFIX

    method_option :ignore_case,
                  default: true,
                  desc: 'Fold lower case to upper case characters',
                  type: :boolean
    
    def find(pattern, name0 = '')

      bSize = options.key?('buffer_size') ? '--buffer-size='+options.buffer_size : ''
      cPrefix = options.coreutils_prefix
      re = Regexp.new("(#{pattern})", options.ignore_case)
      name = name0 == '' ? pattern : name0
      
      pids = Array.new
      tmpfiles = Array.new
      ff = Bio::FlatFile.open(Bio::FastaFormat, STDIN)
      ff.each do |entry|
        tmpfiles << tmpfile = get_temporary_path('find', 'bed', false)
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
      unlink_files(tmpfiles)
      
    end

    #

    desc 'gap BEDgz1 BEDgz2',
         "Calculate gap distances from 5'-end of fragments 1 to 3'-end of fragments 2"

    method_option *OPT_PARALLEL
    method_option *OPT_COREUTILS_PREFIX

    method_option :minimum_gap,
                  default: -10000,
                  desc: 'Minimum gap distans to be reported',
                  type: :numeric

    method_option :maximum_gap,
                  default: 2500,
                  desc: 'Maximum gap distans to be reported',
                  type: :numeric

    def gap(bedgz1, bedgz2)

      cPrefix = options.coreutils_prefix

      chrs = Hash.new
      open("| unpigz -c #{bedgz1} #{bedgz2} | #{cPrefix}cut -f 1 | #{cPrefix}uniq | #{cPrefix}sort -u").each do |line|
        chrs[line.rstrip] = ''
      end

      max = options.maximum_gap
      min = options.minimum_gap
      reSep = /\t/
      tmpfiles = Hash.new
      chrs.keys.each do |chr|
        tmpfiles[chr] = get_temporary_path('gap', 'csv', false)
      end
      Parallel.each(chrs.keys, in_processes: options.parallel) do |chr|
        bed2 = Array.new
        open("| gunzip -c #{bedgz2} | grep '^#{chr}\t'").each do |line|
          chr0, *cols = line.rstrip.split(reSep)
          cols[0] = cols[0].to_i
          cols[1] = cols[1].to_i
          bed2 << cols
        end
        fp = open(tmpfiles[chr], 'w')
        open("| gunzip -c #{bedgz1} | grep '^#{chr}\t'").each do |line|
          chr0, start, stop, name, score, str = line.rstrip.split(reSep)
          if str == '+'
            bed2.each do |bed|
              dist = bed[1] - start.to_i + 1
              fp.puts [name, bed[2], dist].join(',') if min <= dist && dist < max
            end
          else
            bed2.each do |bed|
              dist = stop.to_i - bed[0] + 1
              fp.puts [name, bed[2], dist].join(',') if min <= dist && dist < max
            end
          end
        end
        fp.close
      end
      system "cat #{tmpfiles.values.join(' ')}"
      unlink_files(tmpfiles)
      
    end
    
  end
end
