# g++ func.cpp && ./a.out 1 1 50 50 0 0 > gnu.gp && tail -n 1 gnu.gp && gnuplot gnu.gp

g++ func.cpp \
&& \
echo -e -n "$\mu = $ \n begin{center} \n begin{tabular}{|c|c|c|c|c|} \n \hline \n $ tau backslash h$ & 0.1 & 0.01 & 0.001 & 0.0001" > RESIDUALS.txt
echo \\\\ >> RESIDUALS.txt
echo "\hline" >> RESIDUALS.txt
for N in 10 100 1000 10000
do
        #for i in {1..4}
        for i in 5
        do
                for M in 10 100 1000 10000
                do
                        if [[ $M == 10 ]]
                        then
                        echo "scale=5; 1 / $N" | bc                  >> RESIDUALS.txt
                        #result=$(echo "scale=5; 1 / $N" | bc); echo "$result" >> RESIDUALS.txt
                        fi
                        
                        echo -n "&  "                                   >> RESIDUALS.txt
                        echo -n "$(echo $i | ./a.out 1 1 $M $N 0 7)  "  >> RESIDUALS.txt
                        #if [[ $M -ne 10000 ]]
                        #then
                        #        echo -n "&  "                           >> RESIDUALS.txt
                        #else
                        #        echo ""                                 >> RESIDUALS.txt
                        #fi
                done
                echo \\\\                                                >> RESIDUALS.txt
        done
        echo "\hline"                                                    >> RESIDUALS.txt
done
echo -e "\n end{tabular} \n end{center}"                                 >> RESIDUALS.txt
echo "Скрипт закончен"

cat RESIDUALS.txt
