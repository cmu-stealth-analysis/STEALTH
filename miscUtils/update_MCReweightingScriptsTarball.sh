#!/bin/bash

# Remove backup files
for backupFileName in `find MCReweightingScripts -type f -name "*~"`; do
    rm -v ${backupFileName}
done

# Remove tarball if it already exists
rm -f MCReweightingScripts.tar.gz

# Create backup
tar -cvzf MCReweightingScripts.tar.gz MCReweightingScripts/Makefile MCReweightingScripts/src/*.cc MCReweightingScripts/include/*.h
