import sys
try:
    import setuptools
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages


install_requires = ['numpy', 'astropy>=1.2, !=1.3.3', 'scipy']
if sys.version_info < (3, 2):  # add backported configparser module if py < 3.2
    install_requires.append('configparser')


extras = {
    'test': ['pytest', 'mpmath'],
}

setup(
    # package description and version
    name="line_classifier",
    version="0.0,0",
    author="Daniel Farrow",
    author_email="dfarrow@mpe.mpg.de",
    description="Tools to classify emission lines detected in HETDEX",

    # list of packages and data
    packages=find_packages(),
    include_package_data=True,
    # don't zip when installing
    zip_safe=False,

    # dependences
    install_requires=install_requires,
    extras_require=extras,

    classifiers=["Development Status :: 3 - Alpha",
                 "Environment :: Console",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License (GPL)",
                 "Operating System :: Unix",
                 "Programming Language :: Python :: 2.7",
                 "Topic :: Scientific/Engineering :: Astronomy",
                 "Topic :: Utilities",
                 ]
)
