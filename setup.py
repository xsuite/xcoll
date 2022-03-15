from setuptools import setup, find_packages, Extension

#######################################
# Prepare list of compiled extensions #
#######################################

extensions = []


#########
# Setup #
#########

setup(
    name='xcoll',
    version='0.0.0',
    description='Xsuite collimation package',
    url='https://github.com/xsuite/xcoll',
    author='Giovanni Iadarola, A. Abramov, F. Van der Veken',
    packages=find_packages(),
    ext_modules = extensions,
    include_package_data=True,
    install_requires=[
        'numpy>=1.0',
        'scipy',
        'xobjects',
        'xpart',
        'xdeps',
        'xtrack'
        ]
    )
