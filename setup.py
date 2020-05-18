import setuptools
from Cython.Build import cythonize
from distutils.extension import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()

extensions = [
    Extension(
       'PloidPy.process_bam', ['PloidPy/process_bam.pyx',
                               'PloidPy/parse_bam.c'],
        libraries=["m"],
    )
]

setuptools.setup(
    name="PloidPy",
    version="1.1.0",
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
    install_requires=[
        'numpy',
        'statsmodels',
        'matplotlib',
        'seaborn'
    ],
    scripts=['scripts/PloidPy'],
    python_requires='>=3.6',
    headers=["PloidPy/parse_bam.h"],
    ext_modules=cythonize(extensions)
)
