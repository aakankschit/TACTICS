from setuptools import setup, find_packages

setup(
    name='PRISMS',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'polars',
        'pandas',
        'rdkit',
        'tqdm',
        'openeye-toolkits'




    ],
    entry_points={
        'console_scripts': [
            # Add command-line scripts here if needed
        ],
    },
    author='AK Nandkeolyar',
    author_email='anandkeo@uci.edu',
    description='A package for probabilistic virtual screening of chemical libraries',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/library_enumeration',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)