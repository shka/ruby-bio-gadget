module Bio
  class Gadget
    class Fq1l < Bio::Gadget

      desc 'to DRAW SKIP', 'Thin the sequences out'

      def to(draw, skip)
        BioGadget.to(draw.to_i, skip.to_i)
      end

    end
  end
end
