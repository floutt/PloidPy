import setuptools
from Cython.Build import cythonize
from pysam import get_include as pysam_get_include
from distutils.extension import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()

extensions = [
    Extension(
       'PloidPy.process_bam', ['PloidPy/process_bam.pyx'],
       include_dirs=pysam_get_include()),
]

setuptools.setup(
    name="PloidPy-oaolayin",
    version="0.1.0",
    author="Oluwatosin Olayinka",
    author_email="oaolayin@live.unc.edu",
    description="Discrete mixture model based ploidy inference tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/floutt/PloidPy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['scripts/PloidPy'],
    python_requires='>=3.6',
    # ext_modules=cythonize('PloidPy/process_bam.pyx',
    #                      include_path=pysam_get_include())
    ext_modules=cythonize(extensions)
)
