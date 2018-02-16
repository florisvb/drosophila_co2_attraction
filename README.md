# CO2 attraction in Drosophila
This repository contains the software associated with the paper "Drosophila have distinct activity-gated pathways that mediate attraction and aversion to CO2". For a pre-print, see: https://www.biorxiv.org/content/early/2017/12/03/227991

The data (4 TB) will be available upon reasonable request upon formal publication.

This readme assumes working knowledge of Ubuntu and python. This code is not actively maintained. It worked on 2018-01-16 using up-to-date versions of the required software below.

Code and data are licensed under a [Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) License](https://creativecommons.org/licenses/by-nc-sa/4.0/ "CC BY-NC-SA 4.0").

## What you need to run our analysis
* Ubuntu (we used Ubuntu 12-16)
* Python (2.7)
* ROS (Robot Operating System, Kinetic): http://wiki.ros.org/kinetic/Installation/Ubuntu
* apt-get repositories: git python-pip python-scipy python-h5py python-progressbar python-sympy python-networkx
* Manual downloads: 
  * https://github.com/taynaud/python-louvain/
  * http://www.pyqtgraph.org/
* pip installs: pandas (0.19), statsmodels
* My packages:
  * FigureFirst: https://github.com/FlyRanch/figurefirst
  * FlyPlotLib: https://github.com/florisvb/FlyPlotLib
  * DataFit: https://github.com/florisvb/DataFit
  * FlyStat: https://github.com/florisvb/FlyStat
* Inkscape

You may wish to do all of this in a virtual environment.

## The data

If you would like to re-run our analysis, please contact the authors (florisvb@gmail.com). We will make the 4 TB data available through a resilio-sync readonly key.

Once you have the data, you will need to follow the following instructions for making the data accessible to the analysis (below).

## Making the data automatically accessible to the analysis
We ran our analysis on several different computers, so to keep track of everything, we created a python package that points to the data and figure template locations. In order to run our analysis, you will need to add your machine and local paths to this repository. 

In `co2_paper_locations/co2_paper_locations`, create a duplicate of `data_locations_analysiscavetech_organized.py`, e.g. `data_locations_yourname.py`. Edit the file so that the paths correspond to the data locations on your machine. Next, you will need to create an environmental variable called `co2_paper_locations` (e.g. type `export co2_paper_locations co2_paper_yourname` in any terminal window in which you plan to run our code, or add that to your .bashrc). Add an `elif` statement to the `__init__.py` file in `co2_paper_locations` that matches, for example, `co2_paper_yourname` to  `data_locations_yourname.py` as in the other if and elif statements.

In `co2_paper_locations/co2_paper_locations`, edit the file `figure_template_locations.py`, so that the paths match your system.

Install the package (from `co2_paper_locations` type `python setup.py install`). 

## Installing our analysis

In addition to co2_paper_locations, you need to install the following included python packages:
* Orchard
* multicat_analysis

## Running the analysis

In each "figure" folder there is a make_figureX.py file. Run this file (`python ./make_figureX.py`) to rerun the analysis and update the associated svg figure files in that directory. You can use this to trace backwards our analysis.

