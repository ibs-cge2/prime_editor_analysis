from setuptools import setup, find_packages

setup(name='pea',
      version='0.4.8',
      packages=find_packages(include=['pea']),
      package_data={'pea':['log.ini']},
      install_requires=['absl-py',
                        'pandas', 'fastparquet',
                        'biopython',
                        'termcolor',
                        'pyyaml',
                       ]
)

