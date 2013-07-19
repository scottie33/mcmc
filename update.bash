#!/bin/bash

wget https://github.com/scottie33/mcmc/archive/master.zip

if [ -f master.zip ]; then
	unzip master.zip
elif [ -f master ]; then
	unzip master
fi

cp -fr mcmc-master/* ./

rm -fr mcmc-master
rm master.zip


