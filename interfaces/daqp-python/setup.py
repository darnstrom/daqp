from setuptools import setup,find_packages, Extension
import os
from shutil import copyfile,copytree,rmtree
from platform import system
import pathlib
from Cython.Distutils import build_ext 
from Cython.Build import cythonize


# Directory containing this setup.py file (works regardless of where pip is invoked from)
setup_dir = pathlib.Path(__file__).parent.resolve()

# Copy C source
src_path = (setup_dir / "../..").resolve()
csrc_dir = setup_dir / 'csrc'
# Check that src_path actually contains the daqp source tree (not just that it exists)
daqp_src_exists = (src_path / 'src').exists() and (src_path / 'include').exists()
if not daqp_src_exists:
    if not csrc_dir.exists():
        print("Could not find daqp directory")
elif not csrc_dir.exists():
    os.mkdir(csrc_dir)
    copytree(str(src_path / 'src'), str(csrc_dir / 'src'))
    copytree(str(src_path / 'include'), str(csrc_dir / 'include'))
    copytree(str(src_path / 'codegen'), str(csrc_dir / 'codegen'))
    copyfile(str(src_path / 'CMakeLists.txt'), str(csrc_dir / 'CMakeLists.txt'))
    copyfile(str(src_path / 'LICENSE'), str(setup_dir / 'LICENSE'))


cython_ext = Extension('daqp',
        sources = [str(setup_dir / 'daqp.pyx'),
            str(csrc_dir / 'src/api.c'),
            str(csrc_dir / 'src/avi.c'),
            str(csrc_dir / 'src/auxiliary.c'),
            str(csrc_dir / 'src/bnb.c'),
            str(csrc_dir / 'src/daqp.c'),
            str(csrc_dir / 'src/daqp_prox.c'),
            str(csrc_dir / 'src/factorization.c'),
            str(csrc_dir / 'src/hierarchical.c'),
            str(csrc_dir / 'src/utils.c'),
            ],
        library_dirs=[str(setup_dir)],
        extra_compile_args=["-O3", "-DPROFILING", "-fassociative-math",
                            "-fno-signed-zeros", "-fno-trapping-math"],
        include_dirs=[str(csrc_dir / 'include')])

setup(name='daqp',
        version='0.8.5',
        description='DAQP: A dual active-set QP solver',
        url='http://github.com/darnstrom/daqp',
        author='Daniel Arnström',
        author_email='daniel.arnstrom@gmail.com',
        license='MIT',
        long_description=open(str(setup_dir / 'README.md'),'r').read(),
        long_description_content_type='text/markdown',
        ext_modules=cythonize(cython_ext),
        cmdclass={'build_ext': build_ext},
        zip_safe=False)

# Cleanup C-source
if daqp_src_exists:
    rmtree(str(csrc_dir))
    os.remove(str(setup_dir / 'LICENSE'))
