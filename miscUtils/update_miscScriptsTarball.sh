#!/bin/bash

# Remove backup files
for backupFileName in `find miscScripts -type f -name "*~"`; do
    rm -v ${backupFileName}
done

# Remove tarball if it already exists
rm -f miscScripts.tar.gz

# Create backup
tar -cvzf miscScripts.tar.gz miscScripts/Makefile miscScripts/src/*.cc miscScripts/include/*.h
