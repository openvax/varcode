# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import os

from setuptools import setup, find_packages
import versioneer

readme_filename = "README.md"
current_directory = os.path.dirname(__file__)
readme_path = os.path.join(current_directory, readme_filename)

readme = ""
try:
    with open(readme_path, 'r') as f:
        readme = f.read()
except Exception as e:
    print(e)
    print("Failed to open %s" % readme_path)

try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except Exception as e:
    print(e)
    print("Failed to convert %s from Markdown to reStructuredText" % readme_filename)


if __name__ == '__main__':
    cmdclass = versioneer.get_cmdclass()
    print("Command class = %s" % cmdclass)

    version = versioneer.get_version()
    print("Version: %s" % version)
    setup(
        name='varcode',
        packages=find_packages(),
        package_data={'varcode.cli': ['logging.conf']},
        version=version,
        cmdclass=cmdclass,
        description="Variant annotation in Python",
        long_description=readme,
        url="https://github.com/hammerlab/varcode",
        author="Alex Rubinsteyn",
        author_email="alex rubinsteyn at gmail's fine email service",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'numpy>=1.7, <2.0',
            'pandas>=0.15',
            'pyensembl>=1.0.3',
            'biopython>=1.64',
            'pyvcf>=0.6.7',
            'memoized_property>=1.0.2',
            'serializable>=0.0.8',
            'sercol>=0.0.2',
        ],
        entry_points={
            'console_scripts': [
                'varcode-variants = varcode.cli.variants_script:main'
            ]
        })
