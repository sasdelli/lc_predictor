from setuptools import setup

setup(name='spfit',
      version='0.1',
      description='A spectral series fitter for SNIa',
      url='',
      author='Michele Sasdelli',
      author_email='sasdelli@mpa-garching.mpg.de',
      license='GPL3',
      packages=['lc_predictor'],
      install_requires=[
          'scipy',
          'scikit-learn>=0.10'
      ],
      dependency_links=['http://github.com/sbailey/empca/raw/master/empca.py'],
      scripts=['lc_predictor/bin/fit_lc.py'],
      include_package_data=True,
      package_dir= {'spfit': 'lc_predictor',
                    'data': 'lc_predictor/trained_data', 
                    'input_examples': 'lc_predictor/input_data'},
      package_data= {'trained_data': ['lc_predictor/trained_data/*']},
      zip_safe=False)
