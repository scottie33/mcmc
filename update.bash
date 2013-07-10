#!/bin/bash

wget https://github.com/scottie33/mcmc/archive/master.zip

unzip master.zip

cp -fr mcmc-master/* ./

rm -fr mcmc-master
rm master.zip


