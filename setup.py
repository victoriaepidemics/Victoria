#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="victoriaepi",
    version="0.0.1",
    author="Marcos Capistran, Antonio Capella, Andres Christen",
    author_email="",
    description="""Victoria epidemics ODE model, and parameter inference using
            Bayesian inference (MCMC with the twalk)""",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Operating System :: POSIX',
        'Operating System :: MacOS :: MacOS X'
    ],
    install_requires=[
        "scipy>=1.5.0",
        "numpy>=1.18.5",
        "matplotlib>=3.2.2",
        "dacite>=1.5.1",
        "pandas>=1.0.5",
        "hjson>=3.0.2",
        "dataclasses"],
    requires=[
        "scipy(>=1.5.0)",
        "numpy(>=1.18.5)",
        "matplotlib(>=3.2.2)",
        "dacite(>=1.5.1)",
        "pandas(>=1.0.5)",
        "hjson(>=3.0.2)",
        "dataclasses"],
    python_requires='>=3.7',
)
# python3 setup.py clean sdist bdist_wheel
# pipreqs --ignore calls,models  victoria/
