#!/bin/bash
#Runs the executables for RGBLP and CPLEX to timing files (in ./timing)

#Location of executables
EXELOC = "./x64/Release/"
#Exetention (e.g. .exe if windows)
EXETENSION = "";

for r in {1..30}; do
    GRFAIL = 20;
    for x in {4..14}; do
        for y in {4..16}; do
            BATCH=$((2**$x))
            SIZE=$((2**$y))
            echo "-------"
            echo "batch" $BATCH " size "$SIZE
            FILENAME="benchmarks/"$SIZE
            
            echo "CLP"
            EXELOC+"CLP"+EXETENSION $FILENAME $BATCH >/dev/null

            echo "CPLEX"
            EXELOC+"CPLEX"+EXETENSION $FILENAME $BATCH > /dev/null
            
            echo "LP" 
            EXELOC+"LP"+EXETENSION $FILENAME $BATCH > /dev/null
            
            echo "mGLPK"
            EXELOC+"GLPK"+EXETENSION $FILENAME 1 $BATCH 1 1 > /dev/null
            
            #GurungRay model tends to fail at certain problem sizes. Ensure that if it fails, dont bother running it for larger problems as it will still fail.
            if(($GR < $y))
	    	echo "GR"
            	if ! ./x64/Release/GurungRay $FILENAME 1 $BATCH 1 1 > /dev/null ; then
                    (($GR = $y))
                fi
            fi
            echo ""
        done
    done
done
