#!/usr/bin/env python
from setuptools import setup, find_packages
from cellphonedb import version

with open("README.md", "rt", encoding="utf-8") as fh:
    long_description_readme_md = fh.read()

with open("requirements.txt", "rt", encoding="utf-8") as fh:
    install_requirements_txt = [line.strip() for line in fh.readlines()]

setup(
    name='CellphoneDB',
    author='TeichLab/VentoLab',
    author_email='contact@cellphonedb.org',
    version=version,
    description='Inferring cell-cell communication',
    long_description=long_description_readme_md,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    maintainer='VentoLab',
    maintainer_email='contact@cellphonedb.org',
    url='https://cellphonedb.org',
    license='MIT',
    exclude_package_data={'': []},
    entry_points={},
    install_requires=install_requirements_txt,
)