from setuptools import setup,find_packages, Extension
import os
from shutil import copyfile,copytree,rmtree
from platform import system
import pathlib
from Cython.Distutils import build_ext 
from Cython.Build import cythonize


# Copy C source
try:
    src_path = pathlib.Path(os.environ["PWD"], "../../../daqp")
except:
    src_path = []
csrc_dir = pathlib.Path('./csrc')
if src_path and os.path.exists(src_path) and not os.path.exists(csrc_dir):
    os.mkdir(csrc_dir)
    copytree(os.path.join(src_path,'src'),os.path.join(csrc_dir,'src'))
    copytree(os.path.join(src_path,'include'),os.path.join(csrc_dir,'include'))
    copytree(os.path.join(src_path,'codegen'),os.path.join(csrc_dir,'codegen'))
    copyfile(os.path.join(src_path,'CMakeLists.txt'),os.path.join(csrc_dir,'CMakeLists.txt'))
    copyfile(os.path.join(src_path,'LICENSE'),pathlib.Path('./LICENSE'))
else:
    print("Could not find daqp directory")


cython_ext = Extension('daqp',
        sources = ['daqp.pyx','daqp.pxd', 
            'csrc/src/api.c', 
            'csrc/src/auxiliary.c', 
            'csrc/src/bnb.c', 
            'csrc/src/daqp.c', 
            'csrc/src/daqp_prox.c', 
            'csrc/src/factorization.c', 
            'csrc/src/utils.c', 
            ],
        library_dirs=['.'],
        extra_compile_args=["-O3"],
        include_dirs=['csrc/include'])

setup(name='daqp',
        version='0.4.0',
        description='DAQP: A dual active-set QP solver',
        url='http://github.com/darnstrom/daqp',
        author='Daniel Arnström',
        author_email='daniel.arnstrom@liu.se',
        license='MIT',
        long_description=open('README.md','r').read(),
        long_description_content_type='text/markdown',
        ext_modules=cythonize(cython_ext),
        cmdclass={'build_ext': build_ext},
        package_data = {'':["daqp.pyx","daqp.pxd"]},
        include_package_data = True,
        zip_safe=False)

# Cleanup C-source
if src_path and os.path.exists(src_path):
    rmtree(csrc_dir)
    os.remove(pathlib.Path('./LICENSE'))
