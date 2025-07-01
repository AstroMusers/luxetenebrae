from setuptools import setup, find_packages
import os

# Ensure __init__.py exists in scripts and subfolders for importability
working_dir = os.path.dirname(os.path.abspath(__file__))
print("Working directory:", working_dir)
init_paths = [
    'scripts',
    'scripts/detailed_evolution',
]
for path in init_paths:
    init_file = os.path.join(working_dir, path, '__init__.py')
    if not os.path.exists(init_file):
        with open(init_file, 'w'):
            pass

setup(
    name='luxetenebrae',
    url='https://github.com/AstroMusers/luxetenebrae',
    author='Aysu Ece Saricaoglu',
    author_email='a.saricaoglu@wustl.edu',
    packages=find_packages(where="."),
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
