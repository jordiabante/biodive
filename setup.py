import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.1'
PACKAGE_NAME = 'biodive'
AUTHOR = 'Jordi Abante'
AUTHOR_EMAIL = 'jordiabante@protonmail.com'
URL = 'https://github.com/jordiabante/biodive'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Discovery of k-mer sequences associated with high rates of sequence diversification.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'scipy',
      'biopython',
      'matplotlib',
      'statsmodels',
      'scikit-learn',
      'editdistance'
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
      keywords=['python', 'horizontal gene transfer', 'mobile genetic elements', 'crispr', 'diversity-generating mechanisms'],
      packages=find_packages()
)
