require 'test_helper'

class Bio::GadgetTest < Minitest::Test
  def test_that_it_has_a_version_number
    refute_nil ::Bio::Gadget::VERSION
  end

  def test_it_does_something_useful
    assert false
  end
end
