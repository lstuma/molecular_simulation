from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

c_ext = Extension("fastsim", sources=["_fastsim.c", "fastsim.c"])

setup(
    name='fastsim',
    version='1.0',
    description='Calculating atom Van-Der-Waals interactions.',
    ext_modules=[c_ext]
)