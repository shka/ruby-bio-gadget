Gem::Specification.new do |gem|
  gem.name          = 'bio-gadget'
  gem.version       = '0.5.0'
  gem.licenses      = ['MIT']
  gem.summary       = gem.description
  gem.description   = %q{Gadgets for bioinformatics}
  gem.authors       = ['Shintaro Katayama']
  gem.email         = ['shintaro.katayama@gmail.com']
  gem.homepage      = 'https://github.com/shka/ruby-bio-gadget'
  gem.extensions    = %w[ext/bio_gadget/extconf.rb]

  gem.files         = `git ls-files`.split($\)
  gem.bindir        = 'exe'                                           
  gem.executables   = gem.files.grep(%r{^exe/}).map{ |f| File.basename(f) }
  gem.test_files    = gem.files.grep(%r{^(test|spec|features)/})
  gem.require_paths = ['lib']

  gem.add_development_dependency 'bundler', '~> 1.12'
  gem.add_development_dependency 'rake', '~> 10.0'
  gem.add_development_dependency 'rake-compiler'
  gem.add_development_dependency 'minitest', '~> 5.0'

  gem.add_dependency 'bio'
  gem.add_dependency 'damerau-levenshtein'
  gem.add_dependency 'mkfifo'
  gem.add_dependency 'parallel'
  gem.add_dependency 'thor' , '~> 0.19.3'
end
