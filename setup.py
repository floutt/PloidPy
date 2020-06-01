import setuptools
import subprocess
import glob
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


def run_configure():
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
    include_dirs.append("htslib")
    run_configure()
    libraries += ['z', 'bz2', 'lzma', 'curl', 'crypt']
    c_files = glob.glob('htslib/*.c') + glob.glob('htslib/cram/*.c')
    elim = ['irods', 'plugin']
    source_files += [x for x in c_files if not any(e in x for e in elim)]
    source_files = list(filter(lambda x: not x.endswith(
        ('htsfile.c', 'tabix.c', 'bgzip.c')), source_files))

with open("README.md", "r") as fh:
    long_description = fh.read()

cmpl_args = ["-Wno-sign-compare", "-Wno-unused-function",
             "-Wno-strict-prototypes", "-Wno-unused-result",
             "-Wno-discarded-qualifiers"]

extensions = [
    Extension('PloidPy.process_bam', source_files, libraries=libraries,
              include_dirs=include_dirs, extra_compile_args=cmpl_args)
]

setuptools.setup(
    name="PloidPy",
    version="1.2.1",
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
