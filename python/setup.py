from setuptools import setup, find_packages

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):

        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False
except ImportError:
    bdist_wheel = None

package_name = "cigsegy"
version = "1.1.7"

setup(
    name=package_name,
    version=version,
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    author='Jintao Li',
    url='https://github.com/JintaoLee-Roger/cigsegy',
    license='MIT',
    install_requires=['numpy'],
    python_requires=">=3.6",
    cmdclass={'bdist_wheel': bdist_wheel},
    packages=find_packages(),
    package_data={'': ['*.so', '*.dylib', '*.pyd', '*.pyi']})
