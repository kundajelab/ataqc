from setuptools import setup, find_packages

setup(
    name='ataqc',
    version='0.2',

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

    packages=find_packages(include=['ataqc'],
                           exclude=['contrib', 'docs', 'tests*']),
    package_dir={'ataqc': 'ataqc'},

    scripts=['ataqc/ataqc'],

    install_requires=['numpy >= 1.10.2', 
                      'scipy', 
                      'pandas', 
                      'matplotlib >= 1.5.1',
                      'pysam >= 0.8.2.1', 
                      'pybedtools >= 0.6.9', 
                      'metaseq == 0.5.6', 
                      'jinja2'],

    include_package_data=True
)
