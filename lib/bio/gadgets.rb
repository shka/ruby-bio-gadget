require 'thor'
require 'bio/gadget/fq1l'

module Bio
  class Gadgets < Thor

    VERSION = '0.5.0'

    register(Bio::Gadget::Fq1l, 'fq1l', 'fq1l [COMMAND]', 'Tools for oneline-fastq')

  end
end
