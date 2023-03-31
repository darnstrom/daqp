from setuptools import setup,find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
import os
from shutil import copyfile,copytree,rmtree
from platform import system
import pathlib


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


class CMakeExtension(Extension):

    def __init__(self, name):
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()
        cmake_path = os.path.join(pathlib.Path().absolute(),'csrc')


        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config
        ]

        if system() == 'Windows':
            cmake_args += ['-G', 'MinGW Makefiles']
        else:  # Unix 
            cmake_args += ['-G', 'Unix Makefiles']

        build_args = [
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cmake_path)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(str(cwd))

        if system() == 'Windows':
            copyfile(os.path.join(build_temp,'libdaqp.dll'),
                     os.path.join(str(extdir.parent.absolute()),'libdaqp.dll'))


setup(name='daqp',
        version='0.0',
        description='DAQP: A dual active-set QP solver',
        url='http://github.com/darnstrom/daqp',
        author='Daniel Arnström',
        author_email='daniel.arnstrom@liu.se',
        license='MIT',
        packages=find_packages(
            where='src',
            include=['daqp']),
        package_dir={"": "src"},
        long_description=open('README.md','r').read(),
        ext_modules=[CMakeExtension('daqp/daqp')],
        cmdclass={
            'build_ext': build_ext,
            },
        zip_safe=False)

# Cleanup C-source
if src_path and os.path.exists(src_path):
    rmtree(csrc_dir)
    os.remove(pathlib.Path('./LICENSE'))
