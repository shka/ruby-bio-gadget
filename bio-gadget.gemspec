require File.expand_path('../lib/bio-gadget/version', __FILE__)

Gem::Specification.new do |gem|
  gem.authors       = ["Shintaro Katayama"]
  gem.email         = ["shintaro.katayama@gmail.com"]
  gem.description   = %q{Gadgets for bioinformatics}
  gem.summary       = gem.description
  gem.homepage      = "https://github.com/shka/ruby-bio-gadget"

  gem.files         = `git ls-files`.split($\)
  gem.executables   = gem.files.grep(%r{^bin/}).map{ |f| File.basename(f) }
  gem.test_files    = gem.files.grep(%r{^(test|spec|features)/})
  gem.name          = "bio-gadget"
  gem.require_paths = ["lib"]
  gem.version       = Bio::Gadget::VERSION

  gem.add_dependency 'thor'
  gem.add_dependency 'parallel'
  gem.add_dependency 'levenshtein-ffi'
  gem.add_dependency 'bio-faster'
  gem.add_dependency 'mkfifo'
end
