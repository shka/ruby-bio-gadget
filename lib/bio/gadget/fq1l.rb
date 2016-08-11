require "bio/gadget/fq1l/convert"
require "thor"

module Bio
  class Gadget < Thor

    register(Fq1l, "fq1l", "fq1l [COMMAND]", "Tools for oneline-fastq")
    
  end
end
