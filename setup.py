# Copyright (c) 2014-2018. Mount Sinai School of Medicine
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
import re

from setuptools import setup, find_packages

readme_filename = "README.md"
current_directory = os.path.dirname(__file__)
readme_path = os.path.join(current_directory, readme_filename)

try:
    with open(readme_path, 'r') as f:
        readme_markdown = f.read()
except Exception as e:
    readme_markdown = ""
    print(e)
    print("Failed to open %s" % readme_path)

# convert README to restructured text format required by PyPI
try:
    import pypandoc
    readme_restructured = pypandoc.convert(readme_markdown, to='rst', format='md')
except Exception as e:
    readme_restructured = readme_markdown
    print(e)
    print("Failed to convert %s from Markdown to reStructuredText" % readme_filename)

# Determine version number
with open('varcode/__init__.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)
    print("Version: %s" % version)

if __name__ == '__main__':
    setup(
        name='varcode',
        packages=find_packages(),
        package_data={'varcode.cli': ['logging.conf']},
        version=version,
        description="Variant annotation in Python",
        long_description=readme_restructured,
        url="https://github.com/openvax/varcode",
        author="Alex Rubinsteyn",
        author_email="alex.rubinsteyn@mssm.edu",
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
            'sercol>=0.1.0',
        ],
        entry_points={
            'console_scripts': [
                'varcode-variants = varcode.cli.variants_script:main'
            ]
        })
