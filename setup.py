from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['spkmeansmodule.c'])
setup(name='mykmeanssp',
     version='2.3',
     description='Python wrapper for custom C extension',
     ext_modules=[module])