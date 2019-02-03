
from numpy.distutils.core import Extension

ext = [Extension(name = 'raddamp',sources = ['raddamp.f95']),
       Extension(name = 'ibslong',sources = ['ibslong.f95'])]

if __name__ =="__main__" :
  from numpy.distutils.core import setup
  setup(name='CTEPY',
        version='1.0',
        author='R. Bruce, M. Blasciwisc, T. Mertens, M. Schaumann',
       ext_modules = ext)
