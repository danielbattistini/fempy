from setuptools import setup, find_packages

setup(
   name='fempy',
   version='0.0.0',
   description='Package for femtoscopic studies',
   author='Daniel Battistini',
   author_email='daniel.battistini@cern.ch',
   packages=find_packages(),  # To load all the submodules automatically
   install_requires=['pyyaml'], # External packages as dependencies
)
