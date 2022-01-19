from setuptools import setup, find_packages
from os.path import dirname, realpath

def _read_requirements_file():
    """Return the elements in requirements.txt."""
    req_file_path = '%s/requirements.txt' % dirname(realpath("__file__"))
    print(req_file_path)
    with open(req_file_path) as f:
        return [line.strip() for line in f]

setup(
    name='logv3lpf',
    version="0.1.0",
    packages=find_packages(),
    license='MIT',
    description='Log(v) 3LPF: A linear power flow solver for three-phase unbalanced distribution systems',
    install_requires=_read_requirements_file(),
    url='https://github.com/Ignacio-Losada/Log-v-3LPF',
    author='Ignacio Losada Carreno',
    author_email='il244@cornell.edu',
    scripts=['logv3lpf/__init__.py',]
)
