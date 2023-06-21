# astro_scripts_uibk

[![pipeline status](https://git.uibk.ac.at/csap5791/astro_scripts_uibk/badges/master/pipeline.svg)](https://git.uibk.ac.at/csap5791/astro_scripts_uibk/commits/master)
[![coverage report](https://git.uibk.ac.at/csap5791/astro_scripts_uibk/badges/master/coverage.svg)](https://git.uibk.ac.at/csap5791/astro_scripts_uibk/commits/master)

## Installation
### If you have not yet installed git enter in your terminal:
```console
sudo apt install git
```
### If you have not used git on your machine before:
```console
git config --global user.name "FirstName LastName"
git config --global user.email "FirstName.LastName@student.uibk.ac.at"
```
### After that:
```console
git clone https://github.com/alex-ebi/astro_scripts_uibk.git
python3 -m pip install ./astro_scripts_uibk
```
### Build Documentation:
```console
cd astro_scripts_uibk
python3 -m pip install --upgrade -r requirements-dev.txt
make -C docs html
firefox docs/build/html/index.html 
```
*(or use your favourite browser)*


## Python Template
Using [Scientific Python Cookiecutter](https://nsls-ii.github.io/scientific-python-cookiecutter/philosophy.html) 
as template with minor changes for university GitLab.

## GitLab Settings
Sensible settings to ensure code quality:

### Protect Master Branch
Settings > Repository > Protected Branches

"Nobody should be able to push straight to master"

### Merge Checks
Settings > General > Merge Requests > Merge Checks

"Pipelines should be okay and all discussions resolved"

### Repository Status and Reports
Settings > CI / CD > General pipelines

Pipeline status and coverage report

Requires test coverage parsing (with `^TOTAL.+?(\d+\%)$`)


## Unused Python Developer Tools 
* Pylint: code analysis tool which helps to enforce coding standards
* Bandit: analyzes code to find common security issues
* Black: formats Python code without compromise
* isort: formats imports by sorting alphabetically and separating into sections
* pip-tools: command line tools for managing packages and requirements
* tqdm: easy to use progress bar

[Even more dev tools](https://reposhub.com/python/learning-tutorial/ml-tooling-best-of-python-dev.html)
