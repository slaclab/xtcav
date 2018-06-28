from setuptools import setup

setup(name='xtcav2',
      version='0.1',
      description='Updated XTCAV analysis code',
      packages=['xtcav2'],
      package_dir = {'xtcav2': 'xtcav'},
      scripts=['bin/xtcavDisplay', 'bin/xtcavDark', 'bin/xtcavLasingOff', 'bin/xtcavLasingOn'])

