#!/bin/bash
#Runs the executables for RGBLP and CPLEX to timing files (in ./timing)
for r in {1..3}; do
    GRFAIL = 20;
    for x in {4..14}; do
        for y in {4..16}; do
            BATCH=$((2**$x))
            SIZE=$((2**$y))
            echo "-------"
            echo "batch" $BATCH " size "$SIZE
            FILENAME="benchmarks/"$SIZE
            
	    echo "CLP"
	    ./x64/Release/CLP $FILENAME $BATCH >/dev/null

            #echo "CPLEX"
            #./x64/Release/CPLEX $FILENAME $BATCH > /dev/null
            
            echo "LP" 
            ./x64/Release/LP $FILENAME $BATCH > /dev/null
            
            echo "mGLPK"
            ./x64/Release/GLPK $FILENAME 1 $BATCH 1 1 > /dev/null

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



for r in {1..300}; do
    #amount GR fails at
    GRFAIL=20
     for x in {4..14}; do
         for y in {4..16}; do
             BATCH=$((2**$x))
             echo "batch" $BATCH " size "$SIZE
             FILENAME="benchmarks/"$SIZE

           echo "CLP"
            ./x64/Release/CLP $FILENAME $BATCH > /dev/null

             #echo "CPLEX"
             #./x64/Release/CPLEX.exe $FILENAME $BATCH > /dev/null


            echo "LP"
            ./x64/Release/LP $FILENAME $BATCH > /dev/null

             echo "mGLPK"
            ./x64/Release/GLPK $FILENAME 1 $BATCH 1 1 > /dev/null

            if(($GRFAIL > $y)); then
              echo "GR"
              if ! ./x64/Release/GurungRay $FILENAME 1 $BATCH 1 1 > /dev/null; then
                ((GRFAIL = $y))
              fi
            fi

             echo ""
         done
     done
-done
