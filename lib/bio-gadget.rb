require 'bio-gadget/version'
require 'bio-gadget/dedup'
require 'bio-gadget/demlt'
require 'bio-gadget/fqxz'
require 'bio-gadget/qvstat'
require 'bio-gadget/wigchr'

require 'tempfile'

module Bio
  class Gadget < Thor

    namespace :bio

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

    def mytempfile(basename, tmpdir = Dir::tmpdir)
      fp = Tempfile.open(basename, tmpdir)
      path = fp.path
      fp.close #!
      path
    end

  end
end
