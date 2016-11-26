require 'parallel'

module Bio
  module Gadget
    class STRT < Thor

      desc 'count SMP BASE BED [BED ...]',
           "Count 5'-ends at BASE in each region defined by BEDs"
      
      method_option :parallel,
                    aliases: '-p',
                    banner: 'N',
                    default: system('which gnproc >/dev/null 2>&1') ? `gnproc`.to_i : (system('which nproc >/dev/null 2>&1') ? `nproc`.to_i : 2),
                    desc: 'Change the number of sorts run concurrently to N',
                    type: :numeric
      
      method_option *Bio::Gadgets::OPT_COREUTILS_PREFIX

      def count(smp, base, bed0, *beds)

        cPrefix = options.coreutils_prefix
        
        smps = Hash.new
        fp = open(smp)
        header = fp.gets.rstrip.split(',')
        idxName = header.index('NAME')
        idxBeds = header.index('5pBEDs')
        fp.each do |line|
          cols = line.rstrip.split(',')
          smps[cols[idxName]] = cols[idxBeds].split(';')
        end
        fp.close

        tmpfile = Bio::Gadgets.getTmpname('strt.count', 'bed')
        system "#{cPrefix}sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 #{bed0} #{beds.join(' ')} > #{tmpfile}"

        counts = Hash.new
        mutex = Mutex.new
        Parallel.map(smps.keys, in_threads: options.parallel) do |name|
          bed5ps = smps[name].map { |bed| "#{base}#{bed}.5p.bed.gz" }
          open("| bedtools intersect -nonamecheck -wa -wb -s -sorted -a #{tmpfile} -b #{bed5ps.join(' ')} | #{cPrefix}cut -f 4,10 | #{cPrefix}sort -u | #{cPrefix}cut -f 1 | #{cPrefix}uniq -c").each do |line|
            cnt, id = line.strip.split(' ')
            mutex.synchronize do
              counts[id] = Hash.new unless counts.key?(id)
              counts[id][name] = cnt.to_i
            end
          end
        end

        names = smps.keys.sort
        puts (['ID'] + names.map { |name| "R|#{name}" }).join(',')
        counts.each do |id, name2count|
          puts ([id] + names.map { |name| name2count.key?(name) ? name2count[name] : 0 }).join(',')
        end
      end
      
    end
  end
end
