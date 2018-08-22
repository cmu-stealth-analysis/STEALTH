#!/bin/bash

# Remove backup files
for file in `find eventSelection -type f -name "*~"`; do
    echo "Removing ${file}" && rm ${file}
done

# Remove tarball if it already exists
rm -f eventSelection.tar.gz

# Create backup
tar -cvzf eventSelection.tar.gz eventSelection/Makefile eventSelection/src/*.cpp eventSelection/include/*.h
