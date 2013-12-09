try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name = 'pypartitions',
      version= '0.1dev',
      description = 'Efficient Sampling Algorithims for Integer Partitioning',
      author = "Ken Locey and Dan McGlinn",
      url = 'https://github.com/klocey/partitions',
      packages = ['pypartitions', 'metrics'],
      package_data = {'pypartitions': ['testfiles/*.txt']},
      license = 'GNU General Public License',
)