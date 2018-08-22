#!/bin/bash

tar -xvzf eventSelection.tar.gz && cd eventSelection && mkdir bin && mkdir obj && make clean && make && cd ..
