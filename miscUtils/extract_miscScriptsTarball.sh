#!/bin/bash

tar -xvzf miscScripts.tar.gz && cd miscScripts && mkdir bin && mkdir obj && make clean && make && cd .. && rm -f miscScripts.tar.gz
