
from numpy.distutils.core import Extension

ext = Extension(name = 'raddamp',
                sources = ['raddamp.f95'])

if __name__ =="__main__" :
  from numpy.distutils.core import setup
  setup(name='f2py_example',
       ext_modules = [ext])
