class Fq1l < Thor

  desc "restore", "restore fastq from 1 line/read to 4 lines/read"

  def restore
    exec "cut -f 1-4 | tr \"\\t\" \"\\n\""
  end
  
end
