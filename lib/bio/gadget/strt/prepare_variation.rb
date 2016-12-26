require 'open3'

module Bio
  class Gadget
    class StrtPrepareVariation < Bio::Gadget

      package_name :prepare_variation

      #

      desc 'hg38 DIR', 'GRCh38/hg38 - human'
      long_desc <<-DESC
Prepare variation data files based on dbSNP BUILD for GRCh38/h38 at DIR, where it has 'ref.fa.fai' genome index file.
DESC
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_DOWNLOAD
      method_option *OPT_GREP_PREFIX

      method_option :dbsnp,
                    banner: 'BUILD',
                    default: 149,
                    desc: 'Build number of dbSNP',
                    type: :numeric
      
      def hg38(dir0)

        dir = File.expand_path(dir0)
        vcf = "#{dir}/hg38.snp#{options.dbsnp}Common.vcf.gz"

        if options.download != 'no'
          download_file("ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b#{options.dbsnp}_GRCh38p7/VCF/common_all_20161122.vcf.gz", vcf)
        end
        
        if options.download != 'only'
        
          outfp = open("#{dir}/variation.snp", 'w')
          bedfp, bedthreads = Open3.pipeline_w("bedtools sort -i stdin",
                                               "pigz -c > #{dir}/variation.snp.bed.gz")
          Open3.pipeline_r("unpigz -c #{vcf}",
                           "#{grep_command(options)} -v '^#'") do |infp, inthreads|
            infp.each do |line|
              chr0, pos0, name0, ref, alt0, *tmp = line.rstrip.split /\t/
              chr = "chr#{chr0}"
              pos = pos0.to_i-1
              alts = alt0.split /,/
              refl = ref.length
              alts.each_index do |i|
                name = alts.length > 1 ? "#{name0}.#{i}" : name0
                alt = alts[i]
                altl = alt.length
                if refl == 1 && altl == 1 && ref != alt
                    outfp.puts [name, 'single', chr, pos, alt].join("\t")
                    bedfp.puts [chr, pos, pos0, name].join("\t")
                elsif refl < altl && !/^#{ref}/.match(alt).nil?
                  posr = pos+refl
                  outfp.puts [name, 'insertion', chr, posr, alt[refl..-1]].join("\t")
                  bedfp.puts [chr, posr-1, posr, name].join("\t")
                elsif refl > altl && !/^#{alt}/.match(ref).nil?
                  posr = pos+altl
                  outfp.puts [name, 'deletion', chr, posr, refl-altl].join("\t")
                  bedfp.puts [chr, posr-1, posr, name].join("\t")
                end
              end
            end
            infp.close
          end
          outfp.close
          bedfp.close
          bedthreads[1].join

        end

      end

    end
  end
end
