import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.md').read()
setup(
    name='edicer',
    version='1.0.1',
    description = "A script to simulate all possible Dicer sequence fragments",
    long_description=readme,
    license = "GPLv2",
    author='Finlay Maguire',
    author_email='finlaymaguire@gmail.com',
    url='https://github.com/fmaguire/eDicer',
    packages=[
        'edicer',
    ],
    package_dir={'edicer': 'edicer'},
    include_package_data=True,
    install_requires=[
    ],
    zip_safe=False,
    keywords='edicer',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
