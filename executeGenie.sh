#!/bin/bash

echo "Running on All Data..."

./genie_analysis C12 2261 0 2000 &
./genie_analysis C12 2261 1 2000 &
./genie_analysis C12 4461 0 2000 &
./genie_analysis C12 4461 1 2000 &
./genie_analysis 56Fe 2261 0 2000 &
./genie_analysis 56Fe 2261 1 2000 &
./genie_analysis 56Fe 4461 0 2000 &
./genie_analysis 56Fe 4461 1 2000 

echo "Process Complete!"
