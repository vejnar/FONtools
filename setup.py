import os
import shutil
from distutils.core import setup
from distutils.command.install_scripts import install_scripts

class remove_extension_install_scripts(install_scripts):
    def run(self):
        install_scripts.run(self)
        for script in self.get_outputs():
            if script.endswith('.py'):
                shutil.move(script, script[:-3])

setup(
    name='fontools',
    version='0.1',
    packages=['fontools'],
    author='Charles E. Vejnar',
    url='https://github.com/vejnar/fontools',
    license='Mozilla Public License 2.0 (MPL 2.0)',
    long_description=open('README.md').read(),
    install_requires=['pyfaidx', 'pyfnutils'],
    cmdclass = {'install_scripts': remove_extension_install_scripts},
    scripts=[os.path.join('script', p) for p in os.listdir('script')]
)
