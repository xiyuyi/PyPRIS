from setuptools import setup

setup(name='PyPRIS',
      version='1',
      description='performing 3D super-resolution reconstruction using progressive '
                  'refinement method for sparse recovery',
      url='[to be added]',
      author='Xiyu Yi, Xingjia Wang',
      author_email='xiyu.yi@gmail.com, xingjia.wang9805@gmail.com',
      license='MIT',
      packages=['PyPRIS'],
      install_requires=[
          'torch', 'numpy',
      ],
      zip_safe=False)

