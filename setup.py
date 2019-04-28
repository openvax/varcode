# Copyright (c) 2014-2019. Mount Sinai School of Medicine
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
        long_description=readme_markdown,
        long_description_content_type='text/markdown',
        url="https://github.com/openvax/varcode",
        author="Alex Rubinsteyn",
        author_email="alex@openvax.org",
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
            'pyensembl>=1.7.4',
            'biopython>=1.64',
            'pyvcf>=0.6.7',
            'memoized_property>=1.0.2',
            'serializable>=0.1.1',
            'sercol>=0.1.4',
        ],
        entry_points={
            'console_scripts': [
                'varcode-genes = varcode.cli.genes_script:main',
                'varcode = varcode.cli.effects_script:main',
            ]
        })
