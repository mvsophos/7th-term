# g++ func.cpp && ./a.out 1 1 50 50 0 0 > gnu.gp && tail -n 1 gnu.gp && gnuplot gnu.gp

g++ func.cpp \
&& \
echo "Normas"                                                           >  RESIDUALS.txt
echo ""                                                                 >> RESIDUALS.txt
for N in 10 100 1000 10000
do
        for i in {1..4}
        do
                for M in 10 100 1000 10000
                do
                        echo -n "$(echo $i | ./a.out 1 1 $M $N 1 7)  "  >> RESIDUALS.txt
                        if [[ $M -ne 10000 ]]
                        then
                                echo -n "&  "                           >> RESIDUALS.txt
                        else
                                echo ""                                 >> RESIDUALS.txt
                        fi
                done
        done
        echo ""                                                         >> RESIDUALS.txt
done
echo "Скрипт закончен"
