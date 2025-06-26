from setuptools import setup, find_packages

setup(
    name='luxetenebrae',
    url='https://github.com/AstroMusers/luxetenebrae',
    author='Aysu Ece Saricaoglu',
    author_email='a.saricaoglu@wustl.edu',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        'numpy',
        'matplotlib',
        'seaborn',
        'astropy',
        'pandas',
        'datetime',
        'h5py'
    ],
    version='0.1',
    license='MIT',
    description='An example of a python package from pre-existing code',
    # long_description=open('README.txt').read(),
)
