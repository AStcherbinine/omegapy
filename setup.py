import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

package_data = {
    '': ['OMEGA_dataref/*', 'res_findcube/*'],
    }

setuptools.setup(
    name='omegapy',
    version='2.2.10',
    author='Aurélien Stcherbinine',
    author_email='aurelien@stcherbinine.net',
    description='Python tools for OMEGA/MEx observations analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://astcherbinine.github.io/omegapy',
    project_urls={
        'Source' : 'https://github.com/AStcherbinine/omegapy',
    },
    packages=setuptools.find_packages(),
    package_data=package_data,
    python_requires='>=3.7',
    setup_requires=['wheel'],
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy'
        ]
    )
