==========================================
lc_predictor
==========================================

## Instalation

Steps to install the package on a fresh Ubuntu 14.04:

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

* download and install empca.py (Bailey 2012):
```cmd
wget https://github.com/sbailey/empca/raw/master/empca.py \
    -O $HOME/.local/lib/python2.7/site-packages/empca.py
```

* Go in the directory where you extracted the code:
```cmd
cd lc_predictor/
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

* Outputs

It creates some output pdfs in out_dir.

For example:

```cmd
evince out_dir/spectra.py &
```
