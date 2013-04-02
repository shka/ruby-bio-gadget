require 'bio-gadget/version'
require 'bio-gadget/dedup'
require 'bio-gadget/demlt'
require 'bio-gadget/fqxz'
require 'bio-gadget/qvstat'
require 'bio-gadget/rgt2mtx'
require 'bio-gadget/wig5p'
require 'bio-gadget/wigchr'

require 'tempfile'

module Bio
  class Gadget < Thor

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

    @@mytemppaths = Array.new

    def mytemppath(basename, tmpdir = Dir::tmpdir)
      fp = Tempfile.open(basename, tmpdir)
      path = fp.path
      @@mytemppaths.push(path)
      fp.close!
      path
    end

    END {
      @@mytemppaths.each { |path| File.unlink(path) if File.exist?(path) }
    }

  end
end
