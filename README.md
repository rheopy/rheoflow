# rheoflow

Non-Newtonian flow calculators

Short term usage
Be sure laminar.py and friction_factor.py are in current working directory.  laminar.py is for laminar non-Newtonian pipe flow.  friction_factor_property.py is for non-Newtonian turbulent pipe flow.

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


## Phase I

Get laminar.py and friction_factor_property.py working and in MVP status

## Phase II

Create notebooks for documentation and training



