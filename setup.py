import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'biodive'
AUTHOR = 'Jordi Abante'
AUTHOR_EMAIL = 'jordiabante@protonmail.com'
URL = 'https://github.com/jordiabante/biodive'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'Discovery of k-mer sequences associated with high rates of sequence diversification.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy==1.19.5',
      'scipy==1.5.4',
      'biopython==1.79',
      'matplotlib==3.3.4',
      'statsmodels==0.12.2',
      'scikit-learn==0.24.2',
      'editdistance==0.6.0'
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
