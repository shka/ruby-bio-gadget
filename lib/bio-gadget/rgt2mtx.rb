module Bio
  class Gadget < Thor

    desc 'rgt2mtx [RGT]', <<DESC
Convert cuffdiff read group tracking file (*.read_group_tracking) into tab-separated matrix. If no given name of tracking file, it reads from standard input.
DESC
    option 'gtf', :aliases => '-g', :type => :string, :desc => 'GTF to revert old transcript_id (oId) renamed by cuffcompare. Moreover, "RNA_SPIKE_*" chromosome names will be inserted for the transcripts aligned on the forward strand of spike-in sequences.'
    option 'sample', :aliases => '-s', :type => :string, :desc => 'Mapping from condition/replicate to sample ID for the column names of output matrix. Tab-separated text with three columns of condition, replicate and the sample ID.'
    def rgt2mtx(rgt="/dev/stdin")

      tid2oid = Hash.new
      unless options['gtf'].nil?
        open(options['gtf']).each { |line|
          cols = line.rstrip.split(/\t/)
          tid = cols[8].match(/transcript_id \"([^\"]+)/).to_a[1]
          oid = cols[8].match(/oId \"([^\"]+)/).to_a[1]
          oid = "#{cols[0]}.#{oid}" if cols[0] =~ /^RNA_SPIKE_/
          tid2oid[tid] = oid
        }
      end

      cr2sid = Hash.new
      unless options['sample'].nil?
        open(options['sample']).each { |line|
          c, r, sid = line.rstrip.split(/\t/)
          cr2sid["#{c}|#{r}"] = sid
        }
      end

      id = nil
      header = true
      raws = Hash.new
      open("| tail -n +2 #{rgt} | sort -k 1").each { |line|
        cols = line.rstrip.split(/\t/)
        id = cols[0] if id.nil?
        if id != cols[0]
          if header
            if options['sample'].nil?
              puts ([''] + raws.keys.sort).join("\t")
            else
              tmp = ['']
              raws.keys.sort.each { |cr| tmp.push(cr2sid[cr]) }
              puts tmp.join("\t")
            end
            header = false
          end
          tmp = [tid2oid.key?(id) ? tid2oid[id] : id]
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
