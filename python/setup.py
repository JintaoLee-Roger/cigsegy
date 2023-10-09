from setuptools import setup, find_packages

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):

        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False
except ImportError:
    bdist_wheel = None

setup(
    name='cigsegy',
    version='1.1.7',
    description=
    'A tool for segy-format file reading and segy-format creating from binary file',
    author='roger',
    url='https://github.com/JintaoLee-Roger/cigsegy',
    license='MIT',
    cmdclass={'bdist_wheel': bdist_wheel},
    packages=find_packages(),
    package_data={'': ['*.so', '*.dylib', '*.pyd', '*.pyi']})
