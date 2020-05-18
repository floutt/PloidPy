import setuptools
import subprocess
import os
from Cython.Build import cythonize
from distutils.extension import Extension


def is_library_installed(lib):
    o = subprocess.run(['whereis %s' % lib], stdout=subprocess.PIPE,
                       shell=True)
    if o.returncode != 0:
        raise Exception("""whereis command not found. Unable to determine if
                        library %s is present""" % lib)
    return not ((o.stdout[len(lib)+1:] == b"\n") or
                (o.stdout[len(lib)+1:] == b""))


def install_htslib():
    print("Compiling C htslib library")
    old_dir = os.getcwd()
    os.chdir("htslib")  # change into htslib directory

    def run_cmd_raise_if_error(cmd):
        o = subprocess.run([cmd], shell=True)
        if o.returncode != 0:
            raise Exception("An error occured when running '%s' command" % cmd)
    # run compilation commands
    run_cmd_raise_if_error("autoheader")
    run_cmd_raise_if_error("autoconf")
    run_cmd_raise_if_error("./configure")
    run_cmd_raise_if_error("make")
    os.chdir(old_dir)
    print("\tsuccess!")


source_files = ['PloidPy/process_bam.pyx', 'PloidPy/parse_bam.c']
libraries = ["m"]
include_dirs = []
headers = ["PloidPy/parse_bam.h"]


if is_library_installed("libhts"):
    print("""htslib library found in system, using system version in
          installation.""")
    libraries.append("hts")
else:
    print("htslib library not found in system.")
    install_htslib()
    include_dirs.append("htslib")
    libraries.append("htslib/hts")


with open("README.md", "r") as fh:
    long_description = fh.read()

extensions = [
    Extension('PloidPy.process_bam', source_files, libraries=libraries,
              include_dirs=include_dirs)
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
    headers=headers,
    ext_modules=cythonize(extensions)
)
