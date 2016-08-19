class Fq1l < Thor

  desc "nr", "extract non-redundant best-quality sequences, esp. for UMI"

  def nr
    pseq = ''
    open("| gsort --parallel=`gnproc` -r -k2,4").each do |line|
      acc, seq, tmp, qual = line.split /\t/
      if pseq != seq
        puts line.strip
        pseq = seq
      end
    end
  end
  
end
