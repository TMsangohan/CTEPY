# load installer class to use for running post install scripts
from distutils.command.install import install as _install
from os import path
import setuptools

# path for README.rst
here = path.abspath(path.dirname(__file__))

# adding post install commands
# can also run a python script from file
class install(_install):
    def run(self):
        _install.run(self)
	from subprocess import call
        call("cp ./tfs/*.tfs .", shell=True)	
        call("cp ./tfs/*.csv .", shell=True)

# To use a consistent encoding
from codecs import open

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# call the extension in order to compile the fortran files
# on the target machine and make them available in python
from numpy.distutils.core import Extension

# external fortran routines
ext = [Extension(name = 'addparticles', sources =['./src/addparticles.f95']),
       Extension(name = 'blowup', sources = ['./src/blowup.f95']),
       Extension(name = 'collimation', sources = ['./src/collimation.f95']),
       Extension(name = 'collision', sources = ['./src/collision.f95']),
       Extension(name = 'fmohl', sources = ['./src/fmohl.f95']),
       Extension(name = 'getemittance', sources = ['./src/getemittance.f95']),
       Extension(name = 'getgauss', sources = ['./src/getgauss.f95']),
       Extension(name = 'ham', sources = ['./src/ham.f95']),
       Extension(name = 'ibsbane', sources = ['./src/ibsbane.f95']),
       Extension(name = 'ibsinterpolat', sources = ['./src/ibsinterpolat.f95']),
       Extension(name = 'ibslong', sources = ['./src/ibslong.f95']),
       Extension(name = 'ibsmodpiwilattice', sources = ['./src/ibsmodpiwilattice.f95']),
       Extension(name = 'ibsnagaitsev', sources = ['./src/ibsnagaitsev.f95']),
       Extension(name = 'ibspiwilattice', sources = ['./src/ibspiwilattice.f95']),
       Extension(name = 'ibspiwismooth', sources = ['./src/ibspiwismooth.f95']),
       Extension(name = 'IO', sources = ['./src/IO.f95']),
       Extension(name = 'keeptransvprof', sources = ['./src/keeptransvprof.f95']),
       Extension(name = 'keepturn', sources = ['./src/keepturn.f95']),
       Extension(name = 'longmatch', sources = ['./src/longmatch.f95']),
       Extension(name = 'rds', sources = ['./src/rds.f95']),
       Extension(name = 'rfupdate', sources = ['./src/rfupdate.f95']),
       Extension(name = 'raddamp',sources = ['./src/raddamp.f95']),
       Extension(name = 'sptranschrom',sources = ['./src/sptranschrom.f95'])]

# running the setup from distutils
# $ python setup.py sdist OR
# $ python setup.py bdist OR
# $ python setup.py build
# for installing after extracting the tarball
# $ python setup.py install
if __name__ =="__main__" :
  from numpy.distutils.core import setup
  setup(name    = 'CTEPY',
        version = '1.0.2.dev',
        description = 'Collider Time Evolution Python FORTRAN version',
        long_description = long_description,
        author   = 'T. Mertens, R. Bruce, M. Blasciwisc,M. Schaumann',
        author_email = 'tom.mertens@cern.ch',
        packages = ['src','ctepy'],
        license  = 'MIT',
        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        # Indicate who your project is intended for
        'Intended Audience :: Accelerator Physicists',
        'Topic :: Simulation Tools :: CTE',
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
        #'Programming Language :: Python :: 3',
        #'Programming Language :: Python :: 3.3',
        #'Programming Language :: Python :: 3.4',
        #'Programming Language :: Python :: 3.5',
        ],
        # List run-time dependencies here.  These will be installed by pip when
        # your project is installed. For an analysis of "install_requires" vs pip's
        # requirements files see:
        # https://packaging.python.org/en/latest/requirements.html
#        install_requires=['numpy','pandas','scipy'],
#        only works with setuptools but then extension breaks
        # If there are data files included in your packages that need to be
        # installed, specify them here.  If using Python 2.6 or less, then these
        # have to be included in MANIFEST.in as well.
        package_data={
        "": ['./tfs/lhcb1-PbPb-6500.tfs','./tfs/lhcb2-PbPb-6500.tfs',
             './tfs/lhcb1-pPb-6500.tfs','./tfs/lhcb2-pPb-6500.tfs',
             './tfs/lhcb1-pp-6500.tfs','./tfs/lhcb2-pp-6500.tfs',
             './tfs/lhcb1-PbPb-HLLHC-7000.tfs','./tfs/lhcb2-PbPb-HLLHC-7000.tfs',
             './tfs/CTEPYinputdf.csv','./tfs/b1bunches.csv','./tfs/b2bunches.csv','/mad/madinCTE.madx','./tfs/ctefsinput.csv'
             ],
        },
#	data_files=[('.',['./tfs/lhcb1-PbPb-6500.tfs','./tfs/lhcb2-PbPb-6500.tfs',
#                './tfs/lhcb1-pPb-6500.tfs','./tfs/lhcb2-pPb-6500.tfs',
#                './tfs/lhcb1-pp-6500.tfs','./tfs/lhcb2-pp-6500.tfs',
#                './tfs/lhcb1-PbPb-HLLHC-7000.tfs','./tfs/lhcb2-PbPb-HLLHC-7000.tfs','./scripts/ctepyfs.py'])],
        ext_modules = ext,
        use_2to3=True,
        scripts=['./scripts/ctepyfs.py','./scripts/ctepylxplusscript.sh','./scripts/submitlxplus.sh','./src/bin/CTEPYv5.py','./mad/madinCTE.madx'],
        cmdclass={'install': install}
        )

