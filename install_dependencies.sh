#!/bin/sh

# (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
# INRIA Rhone-Alpes
# STDP model : script automatically installing dependencies

UNAME=$(uname -a)
printf "Checking your os...\n$UNAME"

if echo $UNAME | grep -qi 'ubuntu' || echo $UNAME | grep -qi 'debian'; then
    sudo apt-get update
    sudo apt-get install gcc gfortran
    sudo apt-get install python-numpy python-scipy python-matplotlib
fi
if echo $UNAME | grep -qi 'fedora' || echo $UNAME | grep -qE 'fc[0-9]+'; then
    sudo yum install gcc gcc-gfortran
    sudo yum install redhat-rpm-config
    sudo yum install numpy scipy python-matplotlib
fi
if echo $UNAME | grep -qi 'arch' || echo $UNAME | grep -qi 'manjaro'; then
    sudo pacman -S gcc gcc-fortran
    sudo pacman -S python2 python2-numpy python2-scipy python2-matplotlib
fi
if echo $UNAME | grep -qi 'darwin'; then
    if type brew; then
        brew tap homebrew/science
        brew tap homebrew/python
        brew update && brew upgrade
        brew install gcc
        brew install python
        brew install numpy scipy matplotlib
    else
        echo 'Your OS is OSX, but "brew" command is not found. You can try to install dependencies manually'
    fi
fi
