module Bio
  class Gadget < Thor

    desc 'wig5p BAM', <<DESC
Convert bam-format alignments into wig-format table. Count is at 5'-end of alignments. Moreover the counts are normalized when there are alignments to references named 'RNA_SPIKE_*'. This procedure requires samtools and bamToBed in BEDTools.
DESC
    option 'reverse', :aliases => '-r', :type => :boolean, :default => false
    option 'name', :aliases => '-n', :type => :string, :default => ' '
    def wig5p(bam)

      str = options['reverse'] ? '-' : '+'

      tmp = 0.0
      open("| samtools view #{bam} | grep RNA_SPIKE_ | grep 'XS:A:+' | cut -f 1 | sort -u").each { |line| tmp += 1.0 }
      abort "No spike-in reads." if tmp == 0.0
      spike = 2234.8/tmp

      acc2dups = Hash.new
      open("| samtools view #{bam} | cut -f 1 | sort | uniq -c").each { |line|
        dups, acc = line.strip.split(/\s/)
        acc2dups[acc] = dups.to_f
      }

      cnts = Hash.new
      tmpbam = mytemppath('.bam')
      abort unless system("samtools view -h #{bam} | grep -E 'XS:A:\\#{str}|^@SQ' | samtools view -S -b - > #{tmpbam}")

      open("| bamToBed -i #{tmpbam} | cut -f 1,#{options['reverse'] ? 3 : 2},4").each { |line|
        chr, poss, acc = line.rstrip.split(/\t/)
        cnts[chr] = Hash.new if !cnts.key?(chr)
        pos = poss.to_i
        if !cnts[chr].key?(pos)
          cnts[chr][pos] = 1.0/acc2dups[acc]
        else
          cnts[chr][pos] = cnts[chr][pos]+1.0/acc2dups[acc]
        end
      }

      puts "track type=wiggle_0 name=\"#{options['name']}\" description=\" \" alwaysZero=on visivility=full maxHeightPixels=128:64:16 color=#{options['reverse'] ? '0,128,255' : '255,128,0'}"

      offset = options['reverse'] ? 0 : 1
      signal = options['reverse'] ? -spike : spike
      cnts.each { |chr, posvals|
        puts "variableStep chrom=#{chr}"
        posvals.keys.sort.each { |pos| puts [pos+offset, signal*posvals[pos]].join("\t") }
      }

    end

  end
end
