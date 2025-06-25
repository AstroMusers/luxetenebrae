from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='luxetenebrae',
    url='https://github.com/AstroMusers/luxetenebrae',
    author='Aysu Ece Saricaoglu',
    author_email='a.saricaoglu@wustl.edu',
    # Needed to actually package something
    packages=['luxetenebrae'],
    # Needed for dependencies
    install_requires=[
        'numpy',
        'matplotlib',
        'seaborn',
        'astropy',
        'pandas',
        'datetime',
        'h5py',
        'git+https://github.com/TeamCOMPAS/COMPAS.git'
    ],
    # *strongly* suggested for sharing
    version='0.1',
    # The license can be anything you like
    license='MIT',
    description='An example of a python package from pre-existing code',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)
