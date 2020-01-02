import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

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
)
