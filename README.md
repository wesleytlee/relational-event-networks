# Supplementary material for paper "Inferring social structure from continuous-time interaction data"

This repository contains code and figures for the two applications featured in our paper "[Inferring social structure from continuous-time interaction data] (https://arxiv.org/abs/1609.02629)." Please visit https://wesleytlee.github.io/relational-event-networks/ for more details.

## Authors

* [Wesley Lee] (https://www.stat.washington.edu/people/wtlee/)
* [Bailey K. Fosdick] (http://www.stat.colostate.edu/~bailey/)
* [Tyler H. McCormick] (http://www.stat.washington.edu/~tylermc/)

## Directory Structure

Each folder, which corresponds to a single application, consists of the following files:
* process_data.R - Cleans and processes raw data (not included)
* students.stan OR swallow.stan - Stan file with corresponding model
* stan_sampling.R - Samples parameters and network using appropriate .stan file
* create_figs.R - Creates figures and animations using processed data and results from Stan
* /figures - A subfolder containing various figures and animations used in the paper
