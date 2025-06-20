#!/bin/bash

# This bash script accompanies the TestPaperCommandLineVertexSimulationTest. This script will execute the test TestPaperCommandLineVertexSimulation.
# This script is set up to take in 10 specific varaibles for Lambda and Gamma from a CSV file.
# The script expects the output csv to only contain the values for Lambda and Gamma and the headers above each value :
# Note: this script expects that values are given in a comman delimited format
# Note this bash script should read all numbers in the bash script regardless of length
# So one could begin with a single set of parameters (low fidelity) then supply 10 in the following runs (high fidelity)

# -----   Begining of CSV ------
#       Lambda     Gamma
#          1         10
#          2         9
#          3         8
#          4         7
#          5         6
#          6         5
#          7         4
#          8         3
#          9         2
#          10        1
# -------   End of CSV --------

# Here we will create a loop that goes through each row and takes the given values for Lambda and Gamma


while IFS="," read -r rec_column1 rec_column2
do
     echo "Lambda: $rec_column1"
     echo "Gamma: $rec_column2"
     ## Running simulation with read in parameters
     ##~/bui/projects/BayesianTissueProject/test/TestPaperCommandLineVertexSimulation -opt1 $rec_column1 -opt2 $rec_column2 &
     echo ""
done < <(tail -n +2 ~/Chaste/projects/BayesianTissueProject/ExampleCommandLineCSV.csv)
# Echos that the simulations have all been ran, this does not mean they did not error.
echo All done