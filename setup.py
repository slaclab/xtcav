from setuptools import setup

setup(name='xtcav',
	version='0.1',
	description='Updated XTCAV analysis code',
	packages=['xtcav'],
	scripts=['bin/xtcavDisp', 'bin/xtcavDark', 'bin/xtcavLasingOff', 'bin/xtcavLasingOn'])

