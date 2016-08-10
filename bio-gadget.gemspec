require File.expand_path('../lib/bio-gadget/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name          = "bio-gadget"
  gem.version       = Bio::Gadget::VERSION
  gem.licesnses     = ['MIT']
  gem.summary       = gem.description
  gem.description   = %q{Gadgets for bioinformatics}
  gem.authors       = ["Shintaro Katayama"]
  gem.email         = ["shintaro.katayama@gmail.com"]
  gem.homepage      = "https://github.com/shka/ruby-bio-gadget"

  gem.files         = `git ls-files`.split($\)
  gem.executables   = gem.files.grep(%r{^bin/}).map{ |f| File.basename(f) }
  gem.test_files    = gem.files.grep(%r{^(test|spec|features)/})
  gem.require_paths = ["lib"]

  gem.add_dependency 'thor'
end
