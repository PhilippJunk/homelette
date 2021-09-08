#!/bin/bash

# Script for executing homelette container
# Philipp Junk, 2021

# usage 
usage () { 
	cat 1>&2 <<EOF
Usage: ${0} -m <mode> [-h] [-w <working_dir>] [-t <template_dir>] [-a <alignment_dir>] [-p <port>] [-s <script>]

Access to the local homelette:latest docker image with different modes
-m <mode>
	Mode
	One of "tutorial", "jupyterlab", "interactive" or "script". If no
	valid value is given, the scripts exits with an error. The default
	mode is "jupyterlab".
	tutorial: Opens the tutorials in an interactive jupyter lab session.
	jupyterlab: Opens an interactive jupyter lab session.
	interactive: Opens an interactive python interpreter.
	script: Executes a script given by [-s <script>] in the Python
	interpreter. -s is required in this mode!.

	In the modes "jupyterlab", "interactive" and "script", it is 
	possible configure the containers working directory with the flags
	-w, -t and -a.

	In the modes "tutorial" and "jupyterlab", it is possible to
	configure the port which Jupyter Lab uses with the -p flag.

-h
	Help
	Prints this message to screen.

-w <working_dir>
	Working Directory
	Connects the given working directory with the working directory
	in the homelette container (/home/). Flag requires an argument. 
	Flag works with modes "jupyterlab", "interactive" and "script"
	and is optional.

-t <template_dir>
	Template Directory
	Connects the given template directory with a template directory
	in the working directory of the homelette container
	(/home/templates/). Flag requires an argument. Flag works with
	modes "jupyterlab", "interactive" and "script" and is optional.

-a <alignment_dir>
	Alignment Directory
	Connects the given alignment directory with an alignment
	directory in the working directory of the homelette container
	(/home/alignments/). Flag requires an argument. Flag works with
	modes "jupyterlab", "interactive" and "script" and is optional.

-p <port>
	Port
	Changes the port that Jupyter Lab is connecting to on the host
	system. Standard port is 8888. Flag requires an argument. Flag
	works with modes "tutorial" and "jupyterlab" and is optional.

-s <script>
	Script
	Specifies the python script to be executes in "script" mode.
	Flag requires an argument. Flag only works with mode "script"
	and is required.
EOF
}

# check if image is present
IMAGE="homelette:latest"
docker image inspect "$IMAGE" >/dev/null 2>&1 || \
	echo "Could not find image '$IMAGE'. Please first create the local image with 'construct_local_image.sh'."

# argument parsing with getopts
# adapted from https://stackoverflow.com/a/16496491/7912251
MODE="jupyterlab"
WORKING_DIR=
TEMPLATE_DIR=
ALIGNMENT_DIR=
PORT=8888
SCRIPT=

while getopts ":m:w:t:a:p:s:h" o; do
	case "${o}" in
		h)
			usage; exit 0;
			;;
		m)
			MODE="${OPTARG}"
			;;
		w)
			WORKING_DIR=$( readlink -f "${OPTARG}" )
			;;
		t)
			TEMPLATE_DIR=$( readlink -f "${OPTARG}" )
			;;
		a)
			ALIGNMENT_DIR=$( readlink -f "${OPTARG}" )
			;;
		p)
			PORT="${OPTARG}"
			;;
		s)
			SCRIPT=$( readlink -f "${OPTARG}" )
			;;
		*)
			usage; exit 1;
			;;
	esac
done
shift "$((OPTIND-1))"

# wrapper for checking for open port
wrapper_port () {
	# check if port is open, if not exit
	if lsof -Pi :${PORT} -sTCP:LISTEN -t >/dev/null ; then
		echo "PORT $PORT occupied. Please choose another port (-p <PORT>)."; exit 1
	fi
}

# wrapper for capturing jupyter lab token and assembling URL
wrapper_jlab () {
	token=
	while [ -z "$token" ];
	do
		sleep 0.1
		# get name of last created container
		container_name=$( docker ps --filter ancestor="$IMAGE" --format "{{.Names}}" | head -1 )
		# get token from logs
		token=$( docker logs "$container_name" 2>&1 | grep token | head -1 | sed 's/^.*lab?token=//' )
	done

	# write output and open URL
	echo -e "\nOpen jupyter lab in your browser:\nhttp://127.0.0.1:${PORT}/lab?token=${token}\n\nTo stop this container from running, shut down jupyter lab or press cntrl-C in this terminal."
	xdg-open "http://127.0.0.1:${PORT}/lab?token=${token}"
}

# transform directories into volume information for docker
volumes=
if [ ! -z "${WORKING_DIR}" ] ; then volumes="${volumes} -v ${WORKING_DIR}:/home/" ; fi
if [ ! -z "${TEMPLATE_DIR}" ] ; then volumes="${volumes} -v ${TEMPLATE_DIR}:/home/templates" ; fi
if [ ! -z "${ALIGNMENT_DIR}" ] ; then volumes="${volumes} -v ${ALIGNMENT_DIR}:/home/alignments" ; fi

# implement different behaviour depending on mode
if [ "${MODE}" = "tutorials" ] ; then
	# check if port is open
	wrapper_port

	# start container
	arguments="--rm -p ${PORT}:8888 $IMAGE"
	docker run $arguments /bin/bash -c 'cd $TUTORIALS && jupyter lab --port=8888 --no-browser --ip=0.0.0.0 --allow-root' > /dev/null 2>&1 &

	# check exit status of container
	exit_code=$?
	if [ $exit_code -ne 0 ]; then echo "Docker exited wity non-zero exit status. Exit status : $exit_code"; exit $exit_code; fi

	# retrieve jupyter lab token
	wrapper_jlab

	# wait for process to finish
	wait

elif [ "${MODE}" = "jupyterlab" ] ; then
	# check if port is open
	wrapper_port

	# assemble arguments
	if [ -z "$volumes" ]; then
		arguments="--rm -p ${PORT}:8888 $IMAGE"
	else
		arguments="--rm -p ${PORT}:8888 $volumes $IMAGE"
	fi

	# start container
	docker run $arguments /bin/bash -c 'jupyter lab --port=8888 --no-browser --ip=0.0.0.0 --allow-root' > /dev/null 2>&1 &

	# check exit status of container
	exit_code=$?
	if [ $exit_code -ne 0 ]; then echo "Docker exited wity non-zero exit status. Exit status : $exit_code"; exit $exit_code; fi

	# retrieve jupyter lab token
	wrapper_jlab

	# wait for process to finish
	wait

elif [ "${MODE}" = "interactive" ] ; then
	# assemble arguments
	if [ -z "$volumes" ]; then
		arguments="-it --rm $IMAGE"
	else
		arguments="-it --rm $volumes $IMAGE"
	fi
	# start container
	docker run $arguments python3

elif [ "${MODE}" = "script" ] ; then
	# check if script is present
	if [ -z "$SCRIPT" ] || [ ! -f "$SCRIPT" ]; then
		echo -e "No valid script given for script mode.\n\n"
		usage ; exit 1
	fi
	# assemble arguments
	if [ -z "$volumes" ]; then
		arguments="-v ${SCRIPT}:/home/script.py --rm $IMAGE"
	else
		arguments="$volumes -v ${SCRIPT}:/home/script.py --rm $IMAGE"
	fi
	# start container
	docker run $arguments python3 script.py

else
	echo -e "Invalid mode: ${MODE}\n\n"
	usage ; exit 1
fi
