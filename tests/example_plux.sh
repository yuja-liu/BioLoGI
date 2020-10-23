#!/bin/bash

if [ -z ${BIOLOGI_HOME+x} ]; then
    export BIOLOGI_HOME=~/Documents/graduation_project/YABMSS
fi

$BIOLOGI_HOME/biologic -c $BIOLOGI_HOME/data/BBa_T9002/BBa_T9002_config.json
