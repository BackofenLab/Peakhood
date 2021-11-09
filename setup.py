#!/usr/bin/env python3

from setuptools import setup


"""Setup peakhood"""

setup(
    name='peakhood',
    version='0.4',
    description='Individual site context extraction for CLIP-Seq peak regions',
    long_description=open('README.md').read(),
    author='Michael Uhl',
    author_email='uhlm@informatik.uni-freiburg.de',
    url='https://github.com/BackofenLab/Peakhood',
    license='MIT',
    scripts=['bin/peakhood', 'bin/gtf_add_transcript_biotype_info.py', 'bin/bed_generate_unique_ids.py'],
    packages=['peakhood'],
    package_data={'peakhood': ['content/*']},
    zip_safe=False,
)

