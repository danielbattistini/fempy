from setuptools import setup

setup(
   name='fempy',
   version='0.0.0',
   description='Package for femtoscopic studies',
   author='Daniel Battistini',
   author_email='daniel.battistini@cern.ch',
   packages=['fempy'],  #same as name
   install_requires=['pyyaml'], #external packages as dependencies
)
