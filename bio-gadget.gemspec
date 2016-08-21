# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'bio/gadgets'

Gem::Specification.new do |gem|
  gem.name          = "bio-gadget"
  gem.version       = Bio::Gadgets::VERSION
  gem.licenses      = ['MIT']
  gem.summary       = gem.description
  gem.description   = %q{Gadgets for bioinformatics}
  gem.authors       = ["Shintaro Katayama"]
  gem.email         = ["shintaro.katayama@gmail.com"]
  gem.homepage      = "https://github.com/shka/ruby-bio-gadget"

  gem.files         = `git ls-files`.split($\)
  gem.bindir        = "exe"                                           
  gem.executables   = gem.files.grep(%r{^exe/}).map{ |f| File.basename(f) }
  gem.test_files    = gem.files.grep(%r{^(test|spec|features)/})
  gem.require_paths = ["lib"]

  gem.add_development_dependency "bundler", "~> 1.12"
  gem.add_development_dependency "rake", "~> 10.0"
  gem.add_development_dependency "minitest", "~> 5.0"

  gem.add_dependency 'mkfifo'
  gem.add_dependency 'thor'
end
