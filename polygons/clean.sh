#!/bin/bash

# Convenience script to clean up the output of rpg. Likely only works on Linux.

if [ "$1" == "" ]
then
	echo "Error: Expected a file to clean."
	exit 1
fi

sed -i "1 d" "$1"
sed -i "$ d" "$1"
