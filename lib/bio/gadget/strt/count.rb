require 'parallel'

module Bio
  class Gadget
    class STRT < Bio::Gadget

      desc 'count SMP BASE BED [BED ...]',
           "Count 5'-ends at BASE in each region defined by BEDs"
      
      method_option *OPT_COREUTILS_PREFIX
      method_option *OPT_PARALLEL

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

        tmpfile = get_temporary_path('strt.count', 'bed')
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
