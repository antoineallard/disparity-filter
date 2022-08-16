from distutils.core import setup
setup(
  name = 'backbone_network',         # How you named your package folder (MyLib)
  packages = ['backbone_network'],   # Chose the same as "name"
  version = '1.1.3',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Extracts the multiscale backbone of complex weighted networks',   # Give a short description about your library
  author = 'Malcolm van Raalte',                   # Type in your name
  author_email = 'malcolm@van.raalte.ca',      # Type in your E-Mail
  url = 'https://github.com/malcolmvr/backbone_network',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/malcolmvr/backbone_network/archive/v1.1.2.tar.gz',    # I explain this later on
  keywords = ['network', 'backbone', 'multiscale', 'serrano et al. 2009'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
    'networkx',
    'numpy',
    'scipy',
  ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
  ],
)