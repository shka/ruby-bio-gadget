module Bio
  class Gadget < Thor

    namespace :bio

    desc 'qvstat QUAL', 'statistics of quality values in *.qual file'
    def qvstat(qualfile)
      stat = Hash.new
      myopen(qualfile) { |fp|
        fp.each { |line|
          next if /^[\#\>]/ =~ line
          qvs = line.rstrip.split
          qvs.each_index { |i|
            qv = qvs[i]
            stat[qv] = Array.new unless stat.key?(qv)
            if stat[qv][i].nil?
              stat[qv][i] = 1
            else
              stat[qv][i] += 1
            end
          }
        }
      }
      statfile = qualfile.sub(/.qual(.gz)?$/, '.qvstat')
      open(statfile, 'w') { |out|
        qvs = stat.keys.sort { |a, b| a.to_i <=> b.to_i }
        qvs.each { |qv|
          out.puts "#{qv} #{stat[qv].join(' ')}"
        }
      }
    end

    private

    def myopen(file, &block)
      # how to write?
      f = (/\|/ !~ file && /\.gz$/ =~ file) ? "| gunzip -c #{file}" : file
      unless block.nil?
        o = open(f); block.call(o); o.close
      else
        open(f)
      end
    end

  end
end
