#!/bin/bash

set -u

############################
# Run all tutorials as tests
############################

# set directory the script was run from as working dir
CURRENT_DIR=$( pwd )
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "${SCRIPT_DIR}"

# helper functions
function setup {
	mkdir -p tests
	cp -r data tests/
	cp Tutorial*.ipynb tests/
	cd tests
}

function cleanup {
	cd ..
	rm -r tests
	cd "${CURRENT_DIR}"
}

# setup
setup

# perform tests
for TUTORIAL_NOTEBOOK in $( ls Tutorial*.ipynb )
do
	jupyter nbconvert --to notebook --execute "${TUTORIAL_NOTEBOOK}"
	# if test fails, cleanup and exit
	result="$?"
	if [ "$result" -ne 0 ]; then echo "Tests unsuccessful." ; cleanup ; exit "$result" ;  fi
done

# cleanup
echo "Tests run sucessfully."
cleanup
exit 0
