require 'parallel'
require 'pathname'

module Bio
  class Gadget < Thor

    namespace :bio

    desc 'fqxz', 'automatic (re)compression of *.fq(.gz|.bz2) files'
    def fqxz
      Parallel.map(Pathname.glob('*.fq{.gz,.bz2,}')) { |fqfilename|
        xzfilename = fqfilename.sub(/\.fq(\.(gz|bz2))*$/, '.fq.xz')
        if !xzfilename.exist?
          case fqfilename.extname
          when '.gz'
            decompressor = 'gunzip -c'
          when '.bz2'
            decompressor = 'bunzip2 -c'
          else
            decompressor = 'cat'
          end
          puts "compressing #{xzfilename}..."
          system "#{decompressor} #{fqfilename} | xz -z -e -c > #{xzfilename} 2> #{xzfilename}.log"
          system "xz -t #{xzfilename} >> #{xzfilename}.log 2>&1"
        end
      }
    end

  end
end
