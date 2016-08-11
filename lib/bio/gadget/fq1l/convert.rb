class Fq1l < Thor

  desc "convert", "convert fastq from 4 lines/read to 1 line/read"

  def convert
    exec "paste - - - -"
  end
  
end
