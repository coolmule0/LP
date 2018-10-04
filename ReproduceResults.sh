#!/bin/bash
#Runs the executables for RGBLP and CPLEX to timing files (in ./timing)
for r in {1..3}; do
    for x in {4..14}; do
        for y in {4..16}; do
            BATCH=$((2**$x))
            SIZE=$((2**$y))
            echo "-------"
            echo "batch" $BATCH " size "$SIZE
            FILENAME="benchmarks/"$SIZE
            
            #echo "CPLEX"
            #./x64/Release/CPLEX.exe $FILENAME $BATCH > /dev/null
            
            #echo "LP" 
            #./x64/Release/LP.exe $FILENAME $BATCH > /dev/null
            
            echo "GLPK"
            ./x64/Release/GLPK.exe $FILENAME 1 $BATCH 1 1 > /dev/null
            
            echo "mGLPK"
            ./x64/Release/mGLPK.exe $FILENAME 1 $BATCH 1 1 > /dev/null

            #echo "GR"
            #./x64/Release/GurungRay.exe $FILENAME 1 $BATCH 1 1 > /dev/null
            
            echo ""
        done
    done
done