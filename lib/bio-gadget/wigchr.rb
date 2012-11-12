module Bio
  class Gadget < Thor

    namespace :bio

    desc 'wigchr WIG CHR', 'extract wiggle track on specified chromosome'
    def wigchr(wigfile, chr)
      target = false
      myopen(wigfile) { |fp|
        fp.each { |line|
          if (/^(fixed|variable)Step/ =~ line)
            if (/chrom=#{chr}\s/ =~ line)
              target = true
              puts line
            else
              target = false
            end
          elsif (/^\d/ =~ line)
            puts line if target
          else
            puts line
          end
        }
      }
    end

  end
end
