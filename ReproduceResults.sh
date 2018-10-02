#!/bin/bash
#Runs the executables for RGBLP and CPLEX to timing files (in ./timing)
for r in {1..5}; do
    for x in {4..16}; do
        for y in {4..16}; do
            BATCH=$((2**$x))
            SIZE=$((2**$y))
            echo "batch" $BATCH " size "$SIZE
            FILENAME="benchmarks/"$SIZE
            echo $FILENAME
            ./x64/Release/CPLEX.exe $FILENAME $BATCH
            ./x64/Release/LP.exe $FILENAME $BATCH
            echo ""
        done
    done
done