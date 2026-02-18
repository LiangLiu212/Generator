#!/bin/bash

cp setup/setup.tar.gz ../
cd ../
tar -xzf setup.tar.gz
source setup_env.sh
rm setup.tar.gz
cd Generator
