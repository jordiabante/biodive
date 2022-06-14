import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.0.0'
PACKAGE_NAME = 'pybiodive'
AUTHOR = 'Jordi Abante'
AUTHOR_EMAIL = 'jordiabante@protonmail.com'
URL = 'https://github.com/jordiabante/biodive'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Package to detect diversity generating mechanisms from sequencing data alone.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'biopython',
      'scipy',
      'sklearn',
      'matplotlib',
      'tqdm',
      'torch',
      'pyro-ppl',
      'mmh3',
      'editdistance',
      'statsmodels'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      keywords=['python', 'crispr', 'dgr'],
      packages=find_packages()
)
