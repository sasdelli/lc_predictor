==========================================
lc_predictor
==========================================

## Installation

We recommend to run this package inside a dedicated anaconda/miniconda python
environment. However, it has also been successfully tested on Ubuntu 14.04.

### Anaconda/Miniconda

We assume that you have anacondan or miniconda installed on your system. You
can download it from

Use the requirements file, ``lcp.yml`` (or ``lcp_OSX.yml`` if you are on Mac Os),
distributed with the repository to create a new environment

    conda-env create -f lcp.yml

This will in set up a anaconda/miniconda python environment ``LCP``, which can
be used by invoking

    source activate LCP

Now install ``empca.py`` (see below)

### Ubuntu 14.04 instructions

Follow this procedure to get install all dependencies on a fresh Ubuntu 14.04:

* Add the universe repositories for scipy (if they are not already enabled):
```cmd
sudo add-apt-repository universe
sudo apt-get update
```

* Install necessary packages:
```cmd
sudo apt-get install python-dev
sudo apt-get install python-numpy python-scipy python-matplotlib
```

* On new Ubuntu systems you can install sklearn from the default repository:
```cmd
sudo apt-get install python-sklearn
```

Finally, install ``empca.py`` (see below)

### Installing empca.py

The lc_predictor module requires the ``empca.py`` (Bailey 2012). Download it
and put it in a location where python will find it, e.g.:
```cmd
wget https://github.com/sbailey/empca/raw/master/empca.py \
    -O $HOME/.local/lib/python2.7/site-packages/empca.py
```

If you are using anaconda/miniconda, you could also move it to the environment
directory. In the case of miniconda (assuming a standard installation) this
would be
```cmd
$HOME/miniconda2/envs/LCP/lib/python2.7/site-packages/
```

### Install and test the package

* Go in the directory where you extracted the code:
```cmd
cd lc_predictor/
```

* If your are using anaconda/miniconda, activate the environment
```
source activate LCP
```

* Install the package globally:
```cmd
python setup.py install
```

* OR Install the package locally:
```cmd
python setup.py install --user
```

* Test it!
```cmd
cd lc_predictor/bin/
python fit_lc.py
```

If everything worked, a number of ``pdf`` files should have been created in the
``out_dir`` folder.
