#!/bin/bash
sed -i '' -e 's/patch0/inlet/g' ./constant/polyMesh/boundary
sed -i '' -e 's/zone0/inlet/g' ./constant/polyMesh/boundary
rm -rf 0
cp -r 0.org 0
