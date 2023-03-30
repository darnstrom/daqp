from setuptools import setup,find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
import os
from platform import system
import pathlib



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
        cmake_path = pathlib.Path().absolute().parents[1]

        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute())
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
        ext_modules=[CMakeExtension('daqp/c')],
        cmdclass={
            'build_ext': build_ext,
            },
        zip_safe=False)
