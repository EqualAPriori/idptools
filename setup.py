from setuptools import setup
setup(name = 'idptools',
      version = '0.1dev',
      description = 'idp tools',
      packages = ['idptools'],
      author = 'Kevin Shen',
      author_email = 'kevin.shen@ucsb.edu',
      entry_points = {          # this here is the magic that binds your function into a callable script
          'console_scripts': ['idptools=idptools.cmdln:dispatch',
              'pyleap=idptools.tleap:cmdln'],
      }
)
