from setuptools import setup, find_packages

setup(
    name='ataqc',
    version='0.1',

    description='ATAqC - quality control for ATAC-seq',

    url='https://github.com/kundajelab/ataqc',

    author='Chuan-Sheng Foo, Daniel Kim',
    author_email='csfoo@cs.stanford.edu',

    license='BSD-3',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Operating System :: POSIX :: Linux',
    ],

    keywords='ATAC-seq bioinformatics QC genomics',

    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    install_requires=['numpy', 'scipy', 'pandas', 'matplotlib',
                      'pysam', 'pybedtools', 'metaseq', 'jinja2' ],

    include_package_data=True,
)
