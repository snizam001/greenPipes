from setuptools import find_packages,setup

requirements = open("requirements.txt").read().strip().split("\n")

exec(open('greenPipe/__init__.py').read())

setup(
	name='greenPipe',
	packages=find_packages(),
	entry_points={
		"console_scripts": ['greenPipe = greenPipe.run:main']
		},
	version=__version__,
	description="pipeline for the analysis of greenCUTRUN datasets and integrated analysis with RNAseq and MassSpectrometry",
	author="Sheikh Nizamuddin",
	license="DKFZ",
	include_package_data=True,
	package_data={"greenPipe": ["data/*"]},
#	scripts=['rscripts/*'],
#	package_rscript={"greenCUTRUN": ["rscripts/*"]},
#	package_dir={"":'greenCUTRUN'},
	install_requires=requirements,
	url='https://github.com/snizam001/greenPipe'
) 
