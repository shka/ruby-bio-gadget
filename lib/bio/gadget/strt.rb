module Bio
  module Gadget
    class STRT < Thor

      OPT_LENGTH_BARCODE = [ :length_barcode, { :banner => 'NT',
                                                :default => 6,
                                                :desc => 'Length of barcode',
                                                :type => :numeric } ]

      OPT_LENGTH_GAP = [ :length_gap, { :banner => 'NT',
                                        :default => 3,
                                        :desc => 'Length of gap (polyG)',
                                        :type => :numeric } ]

      OPT_LENGTH_MINIMUM = [ :length_minimum,
                             { :banner => 'NT',
                               :default => 25,
                               :desc => 'Minimum length after preprocess',
                               :type => :numeric } ]

      OPT_LENGTH_UMI = [ :length_umi, { :banner => 'NT',
                                        :default => 6,
                                        :desc => 'Length of UMI',
                                        :type => :numeric } ]
      
      OPT_LOW_QUALITIES = [ :low_qualities, { :banner => 'CHARACTERS',
                                              :default => '!"#',
                                              :desc => 'Low quality characters',
                                              :type => :string } ]
      
      def self.configure_prepSeq(options)
        uLength = options.length_umi
        bLength = options.length_barcode 
        gLength = options.length_gap
        mLength = options.length_minimum
        return [ options.key?('buffer_size') ?
                   '--buffer-size='+options.buffer_size : '',
                 options.coreutils_prefix == '' ?
                   '' : "--prefix-coreutils=#{options.coreutils_prefix}",
                 options.prefix_grep == '' ?
                   '' : "--prefix-grep=#{options.prefix_grep}",
                 "--parallel=#{options.parallel}",
                 bLength,
                 gLength,
                 mLength,
                 mLength + uLength + bLength + gLength,
                 uLength,
                 "#{'.' * uLength}#{'.' * bLength}#{'G' * (gLength-1)}",
                 options.coreutils_prefix ]
      end

    end
  end
end

require 'bio/gadget/strt/bam2bed5p.rb'
require 'bio/gadget/strt/count.rb'
require 'bio/gadget/strt/depth.rb'
require 'bio/gadget/strt/prepBed.rb'
require 'bio/gadget/strt/prepGtf.rb'
require 'bio/gadget/strt/prepSeq.rb'
require 'bio/gadget/strt/qcSmp.rb'
