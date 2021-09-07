#!/bin/bash

# Script for constructing a local homelette image with a valid modeller license
# Philipp Junk, 2021
# https://github.com/PhilippJunk/homelette

# set template image
export template_image="philippjunk/homelette_template:latest"

# define usage string
usage () {
	cat 1>&2 <<EOF
Usage: ${0} <modeller_key>

Creates a local image of homelette with all 
dependencies installed. 

Please be aware that this image contains your MODELLER
license key and should therefore not be shared.

A MODELLER license key can be acquired here:
https://salilab.org/modeller/registration.html

This script assumes that docker is installed on 
your system.
EOF
}

# check if modeller_key was given
if [ -z "$1" ] ; then usage ; exit 1 ; fi
modeller_key="${1}"


# check for template image, download if not found
if docker image inspect "${template_image}" > /dev/null 2>&1; then
	echo "Template image found"
else
	echo "Template image not found, downloading image"
	docker pull "${template_image}"
fi

# construct local image
# TODO figure out how to get MODELLER_VERSION from homelette_template
MODELLER_VERSION="10.1"
docker build --no-cache -t homelette:latest -<<EOF
FROM ${template_image}
# insert modeller license key
RUN cat /usr/lib/modeller${MODELLER_VERSION}/modlib/modeller/config.py | \
	sed "s/xxx/${modeller_key}/" > tmp_file && \
	mv tmp_file /usr/lib/modeller${MODELLER_VERSION}/modlib/modeller/config.py && \
	echo "Test if license key is accepted:" && \
	python3 -c "import modeller" 2>&1 /dev/null && \
	echo "License key successful integrated, local image build."
EOF



# print usage if building container was not successful
if [ $? != 0 ] ; then 
	echo "Invalid modeller license key: ${modeller_key}. Build unsuccessful." ; 
	exit 1 ; 
fi

# print success and warning about sharing image
cat <<EOF
homelette docker image successfully built.

!! Please be aware that this image should only be used locally. !!
It contains an active MODELLER license key and should therefore not 
be shared on DockerHub or other repositries!

To use the local image, either access the docker container yourself, or
use the docker/homelette.sh script that creates different access points
to the docker container. 

Type 
>>> ./homelette.sh -h
to learn more. 
EOF

exit 0
