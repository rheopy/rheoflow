# rheoflow

Lbrary for non-Newtonian flow calculators

# Try it in the cloud
[![Binder](http://mybinder.org/badge_logo.svg)](http://beta.mybinder.org/v2/gh/rheopy/rheoflow/master)
[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/whharttiv/rheoflow/blob/master/notebooks/index.ipynb)

# Installation

* Currently since the repository is not public download the zip file of this branch and from any folder (e.g. download folder) run:

```
pip install downloaded_zip_file.zip
```

* if the repository is set to public within the intranet

* Requires git installed

* to use ssh protocol looks like gitlab.pg.com requires to set:

run from shell
```
git config --global http.sslVerify false
```

* install the library wit pip:
```
pip install git+https://github.com/rheopy/rheoflow.git
```

* if you change your mind:
```
pip uninstall rheoflow
```

* to upgrade (will uninstall current version and install newest from master branch)

```
pip install git+https://github.com/rheopy/rheoflow.git --upgrade
```

## Notebooks

* friction factor documentation.ipynb

This notebook is a start at documenting models and numerical methods used in friction_factor_property.py.  This will develop into training material.

* friction factor example template.ipynb


Quick notebook for calculations

* friction factor newtonian and powerlaw.ipynb


Notebook showing how to implement Newtonian and power-law viscosity models and calculate turbulent friction factors and such

* friction factor class demonstration carreau.ipynb


Notebook with Carreau viscosity model example for pipe flow friction factors

* friction factor class demonstration herschel_bulkley.ipynb


Notebook with Herschel-Bulkey viscosity model example for pipe flow friction factors

* laminar_pipe_demonstration.ipynb

* laminar_pipe_dashboard.ipynb

* laminar_pipe_deomstration_powerlaw.ipynb

* von_karmen_pipe_flow.ipynb

* von_karmen_pipe_flow_nn.ipynb




