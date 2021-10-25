#!/bin/bash

tar -xvzf MCReweightingScripts.tar.gz && cd MCReweightingScripts && mkdir bin && mkdir obj && make clean && make && cd .. && rm -f MCReweightingScripts.tar.gz
