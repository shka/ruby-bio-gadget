module Bio
  class Gadget < Thor

    desc 'rgt2mtx [RGT]', <<DESC
Convert cuffdiff read group tracking file into tab-separated matrix. If no given name of tracking file, it reads from standard input.
DESC
    def rgt2mtx(rgt="/dev/stdin")

      id = nil
      header = true
      raws = Hash.new
      open("| tail -n +2 #{rgt} | sort -k 1").each { |line|
        cols = line.rstrip.split(/\t/)
        id = cols[0] if id.nil?
        if id != cols[0]
          if header
            puts (['tracking_id'] + raws.keys.sort).join("\t")
            header = false
          end
          tmp = [id]
          raws.keys.sort.each { |k| tmp.push(raws[k]) }
          puts tmp.join("\t")
          id = cols[0]
        end
        raws["#{cols[1]}|#{cols[2]}"] = cols[3]
      }
      tmp = [id]
      raws.keys.sort.each { |k| tmp.push(raws[k]) }
      puts tmp.join("\t")
    end

  end
end
