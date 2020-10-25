#!/usr/bin/env python3
# coding: utf-8

from setuptools import setup, find_packages

setup(
    name = "biologi",
    version = "0.1",
    url = "https://github.com/RandolphLiu/BioLoGI",
    author = "Yujia Liu",
    author_email = "rainl199922@gmail.com",
    license = "GPLv3",
    description = "A computational tool for automated modeling of" +\
        "logic gates in syntheic circuits and human transcriptome",
    install_requires = [
        "numpy",
        "scipy",
        "sympy",
        "pymc3",
        "matplotlib",
        "dnaplotlib == 1.0",
        "python-libsbml",
    ],
    entry_points = {
        "console_scripts": [
            "biologi=biologi.command:main",
        ]
    },
    packages = find_packages(),
)

