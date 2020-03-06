import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

package_data = {
    '': ['../data/*', 'omega_routines/*'],
    }

setuptools.setup(
    name='omegapy',
    version='1.2',
    author='Aurélien Stcherbinine',
    author_email='aurelien.stcherbinine@ias.u-psud.fr',
    description='Python tools for OMEGA/MEx observations analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://git.ias.u-psud.fr/astcherb1/omegapy',
    packages=setuptools.find_packages(),
    package_data=package_data,
    python_requires='>=3.7',
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy'
        ]
    )
