require 'parallel'
require 'pathname'

module Bio
  class Gadget < Thor

    namespace :bio

    desc 'fqlzma', 'automatic (re)compression of *.fq(.gz|.bz2) files'
    def fqlzma
      Parallel.map(Pathname.glob('*.fq{.gz,.bz2,}')) { |fqfilename|
        lzmafilename = fqfilename.sub(/\.fq(\.(gz|bz2))*$/, '.fq.lzma')
        if !lzmafilename.exist?
          case fqfilename.extname
          when '.gz'
            decompressor = 'gunzip -c'
          when '.bz2'
            decompressor = 'bunzip2 -c'
          else
            decompressor = 'cat'
          end
          puts "compressing #{lzmafilename}..."
          system "#{decompressor} #{fqfilename} | lzma -c > #{lzmafilename} 2> #{lzmafilename}.log"
          system "lzma -t #{lzmafilename} >> #{lzmafilename}.log 2>&1"
        end
      }
    end

  end
end
