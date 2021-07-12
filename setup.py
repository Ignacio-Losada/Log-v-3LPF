from setuptools import setup, find_packages
from os.path import dirname, realpath

def _read_requirements_file():
    """Return the elements in requirements.txt."""
    req_file_path = '%s/requirements.txt' % dirname(realpath("__file__"))
    print(req_file_path)
    with open(req_file_path) as f:
        return [line.strip() for line in f]

setup(
    name='pysoda',
    version="0.1.0",
    packages=find_packages(),
    license='MIT',
    description='SoDa is an irradiance-based synthetic Solar Data generation tool to generate realistic sub-minute solar photovoltaic (PV) power time series',
    install_requires=_read_requirements_file(),
    url='https://github.com/Ignacio-Losada/SoDa',
    author='Ignacio Losada Carreno',
    author_email='ilosadac@asu.edu',
    scripts=['soda/__init__.py',]
)
