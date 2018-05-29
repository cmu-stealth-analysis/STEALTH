#!/bin/bash

# Silly little script that adds a common prefix to all lines in a file; useful for creating input file lists that start with the standard "root://cmseos.fnal.gov//store/..." prefix

sed -e "s|^|${1}|" ${2} > ${2}.tmp
mv ${2}.tmp ${2}
