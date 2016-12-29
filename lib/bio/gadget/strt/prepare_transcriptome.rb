require 'open3'

module Bio
  class Gadget
    class StrtPrepareTranscriptome < Bio::Gadget

      package_name :prepare_transcriptome

      #
      
      desc 'hg38 DIR', 'GRCh38/hg38 - human'
      long_desc <<-DESC
Prepare transcriptome data files based on GENCODE gene annotation RELEASE for GRCh38/h38 at DIR, where it has 'ref.fa.fai' genome index file.
DESC
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_DOWNLOAD
      method_option *OPT_GREP_PREFIX

      method_option :gencode,
                    banner: 'RELEASE',
                    default: 25,
                    desc: 'Release number of GENCODE',
                    type: :numeric
      
      def hg38(dir0)

        dir = File.expand_path(dir0)
        gtf = "#{dir}/hg38.gencode.v#{options.gencode}.annotation.gtf.gz"
        
        if options.download != 'no'
          download_file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v#{options.gencode}.annotation.gtf.gz",
                        gtf)
        end
        
        if options.download != 'only'
          
          pipeline("unpigz -c #{gtf}",
                   "hisat2_extract_splice_sites.py - > #{dir}/transcriptome.splice_sites")
          
          pipeline("unpigz -c #{gtf}",
                   "hisat2_extract_exons.py - > #{dir}/transcriptome.exons")
          
          fp_ribosome = open_bed_w("#{dir}/ribosome.bed")
          fp_whole = open_bed_w("#{dir}/spikein_whole.bed")
          fp_5end = open_bed_w("#{dir}/spikein_5end.bed")
          open("#{dir}/ref.fa.fai").each do |line|
            acc, len, *tmp = line.rstrip.split
            fp_ribosome.puts [acc, 0, len, acc, 0, '+'].join("\t") if acc =~ /^RIBO_/
            if acc =~ /^RNA_SPIKE_/
              fp_whole.puts [acc, 0, len, acc, 0, '+'].join("\t")
              fp_5end.puts [acc, 0, 50, acc, 0, '+'].join("\t")
            end
          end
          fp_ribosome.close
          fp_whole.close
          fp_5end.close

          atgs = Hash.new
          regex_transcript_id = /transcript_id "([^"]+)"/
          regex_gene_id = /gene_id "([^"]+)"/
          regex_gene_name = /gene_name "([^"]+)"/
          regex_exon_number = /exon_number (\d+)/
          Open3.pipeline_r(
            "unpigz -c #{gtf}",
            "#{grep_command(options)} '\tstart_codon\t'") do |fp, threads|
            fp.each do |line|
              cols = line.rstrip.split /\t/
              atgs[regex_transcript_id.match(cols[8]).to_a[1]] =
                cols[cols[6] == '+' ? 3 : 4].to_i
            end
            fp.close
          end
          
          bed_coding_exon = "#{dir}/transcriptome.coding_exon.bed"
          bed_coding_5utr = "#{dir}/transcriptome.coding_5utr.bed"
          bed_coding_promoter = "#{dir}/transcriptome.coding_promoter.bed"
          
          fp_coding_gene = open_bed_w("#{dir}/transcriptome.coding_gene.bed")
          fp_coding_exon = open_bed_w(bed_coding_exon)
          fp_coding_5utr = open_bed_w(bed_coding_5utr)
          fp_coding_promoter = open_bed_w(bed_coding_promoter)
          fp_other_gene = open_bed_w("#{dir}/transcriptome.other_gene.bed")
          fp_other_exon = open_bed_w("#{dir}/transcriptome.other_exon.bed")
          fp_other_1st_exon = open_bed_w("#{dir}/transcriptome.other_1st_exon.bed")
          fp_other_promoter = open_bed_w("#{dir}/transcriptome.other_promoter.bed")
          Open3.pipeline_r(
            "unpigz -c #{gtf}",
            "#{grep_command(options)} -E '\t(exon|transcript)\t'") do |fp, threads|
            fp.each do |line|
              cols = line.rstrip.split /\t/
              ann = cols[8]
              transcript_id = regex_transcript_id.match(ann).to_a[1]
              gene_id = regex_gene_id.match(ann).to_a[1]
              gene_name = regex_gene_name.match(ann).to_a[1]
              exon_number = regex_exon_number.match(ann).to_a[1].to_i
              chr = cols[0]
              left = cols[3].to_i
              right = cols[4].to_i
              str = cols[6]
              acc = "#{gene_name}|#{gene_id}|#{transcript_id}"
              exon = [chr, left-1, right, acc, 0, str].join("\t")
              if cols[2] == 'transcript'
                if atgs.key?(transcript_id)
                  fp_coding_gene.puts exon
                else
                  fp_other_gene.puts exon
                end
                next
              end
              if atgs.key?(transcript_id)
                fp_coding_exon.puts exon
                atg = atgs[transcript_id]
                if (str == '+' && right < atg) || (str == '-' && atg < left)
                  fp_coding_5utr.puts exon
                elsif (str == '+' && left < atg && atg <= right)
                  fp_coding_5utr.puts [chr, left-1, atg-1, acc, 0, '+'].join("\t")
                elsif (str == '-' && left <= atg && atg < right)
                  fp_coding_5utr.puts [chr, atg, right, acc, 0, '-'].join("\t")
                end
                if exon_number == 1
                  if str == '+'
                    fp_coding_promoter.puts [chr, left-500, left-1, acc, 0, '+'].join("\t")
                  else
                    fp_coding_promoter.puts [chr, right, right+499, acc, 0, '-'].join("\t")
                  end
                end
              else
                fp_other_exon.puts exon
                if exon_number == 1
                  fp_other_1st_exon.puts exon
                  if str == '+'
                    fp_other_promoter.puts [chr, left-500, left-1, acc, 0, '+'].join("\t")
                  else
                    fp_other_promoter.puts [chr, right, right+499, acc, 0, '-'].join("\t")
                  end
                end
              end
            end
            fp.close
          end
          fp_coding_gene.close
          fp_coding_exon.close
          fp_coding_5utr.close
          fp_coding_promoter.close
          fp_other_gene.close
          fp_other_exon.close
          fp_other_1st_exon.close
          fp_other_promoter.close

          merge_bed_by_gene(options,
                            "#{dir}/transcriptome.coding_5end.bed",
                            "Coding5end",
                            "5UTR and the proximal upstream of protein coding genes based on GENCODE version #{options.gencode} for GRCh38/hg38",
                            bed_coding_promoter, bed_coding_5utr)
          merge_bed_by_gene(options,
                            "#{dir}/transcriptome.coding_whole.bed",
                            "CodingWhole",
                            "Exon and the proximal upstream of protein coding genes based on GENCODE version #{options.gencode} for GRCh38/hg38",
                            bed_coding_promoter, bed_coding_exon)

          system "pigz #{dir}/*.bed" or exit $?.exitstatus

          genePred = "#{dir}/hg38_refGene.txt"
          pipelin("unpigz -c #{gtf}",
                  "#{grep_command(options)} 'tag \"basic\"'",
                  "gtfToGenePred -geneNameAsName2 -genePredExt stdin #{genePred}")
          system "retrieve_seq_from_fasta.pl --format refGene --seqfile #{dir}/ref.fa --outfile #{dir}/hg38_refGeneMrna.fa #{genePred}" or exit $?.exitstatus

        end
        
      end

      no_commands do

        def open_bed_w(bed)
          open("| bedtools sort -i stdin > #{bed}", 'w')
        end

        def merge_bed_by_gene(options, outbed, tname, tdesc, *inbeds)
          tmpbed = get_temporary_path('strt.prepare_transcriptome', 'bed')
          Open3.pipeline_r(
            "#{cat_command(options)} #{inbeds.join(' ')}",
            "#{sort_command(options)} -t '|' -k 2") do |infp, inths|

            pregacc = ''
            outfp = nil
            outths = nil
            infp.each do |line|
              chr, left, right, name, *cols = line.rstrip.split /\t/
              sym, gacc, tacc = name.split /\|/
              if pregacc != gacc
                unless outfp.nil?
                  outfp.close
                  outths[1].join
                end
                outfp, outths = Open3.pipeline_w(
                         "bedtools sort -i stdin",
                         "bedtools merge -s -c 4 -o distinct >> #{tmpbed}")
                pregacc = gacc
              end
              outfp.puts ([chr, left, right, "#{sym}|#{gacc}"]+cols).join("\t")
            end
            unless outfp.nil?
              outfp.close
              outths[1].join
            end
          end
          system "echo 'track name=#{tname} description=\"#{tdesc}\" visibility=3 colorByStrand=\"38,139,210 203,75,22\"' > #{outbed}"
          pipeline("ruby -anle 'puts ($F.values_at(0, 1, 2) + [$F[4], 0, $F[3]]).join(\"\t\")' < #{tmpbed}",
                   "bedtools sort -i stdin >> #{outbed}")
        end
      end
      
    end
  end
end
