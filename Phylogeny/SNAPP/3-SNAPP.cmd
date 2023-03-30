#!/bin/bash

beast -overwrite -threads 19 highcov_n2_5.xml

ruby ~/bin/add_theta_to_log.rb -l ./highcov_n2_4.log -t highcov_n2_4.trees -g 1.5 -o ./highcov_n2_4_popsize.log

# Check the outputed trees with tensitree

# Create a consensus tree with treeannotator
