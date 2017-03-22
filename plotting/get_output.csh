#!/bin/bash
grep 'MOM Day' $1 | awk -F'[, \t]*' '{print $3}' > day
grep 'MOM Day' $1 | awk -F'[, \t]*' '{print $6}' > energy
grep 'MOM Day' $1 | awk -F'[, \t]*' '{print $8}' > cfl
grep 'MOM Day' $1 | awk -F'[, \t]*' '{print $10}' > mass
grep 'MOM Day' $1 | awk -F'[, \t]*' '{print $14}' > temp
grep 'MOM Day' $1 | awk -F'[, \t]*' '{print $12}' > salt
