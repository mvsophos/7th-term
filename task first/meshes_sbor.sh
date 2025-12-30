# Подаем на вход параметры, которые необходимы, чтобы определить какой файл надо закачать. Далее выкачиваем файл в программу и считаем нормы, с помощью функций, которые уже были написаны в func.cpp. Имя файла можно задать как аргумент самой функции
# В файле всегда ровно 100 чисел

init_M=$1
mult=$2
rezhim=0

filename="common_latex_tablic_for_better_M.txt"
echo "" > $filename

g++ calc_normas.cpp -o calc_normas.out \
&& \
if [ 0 -eq 1 ]; then
#for C in 1 10 100 ; do
for C in 10 100; do
        for mu in 0.1 0.01 0.001 ; do
                echo -E "\centering    \begin{tabular}{ |c|c| } \hline    & $ C = $C, \; \mu = $mu $\\\\ \hline" >> $filename
                for level in 1 2 3 -1 ; do
                        for rezhim in 0                         # эта штука, на режиме Cr
                        do
                                if [ $level == -1 ]; then
                                        echo -E " $ v - u $  & " >> $filename
                                else
                                        echo -E " $ v - v^$level $  & " >> $filename
                                fi
                                # эти строки задают параметры, для передачи файла как аргумента
                                if  [ $level == -1 ]; then      M=$init_M;           N=$(($M * $mult))
                                elif [ $level == 0 ]; then      M=$init_M;           N=$(($M * $mult))
                                elif [ $level == 1 ]; then      M=$(($init_M * 3));  N=$(($M * $mult))
                                elif [ $level == 2 ]; then      M=$(($init_M * 9));  N=$(($M * $mult))
                                else                            M=$(($init_M * 27)); N=$(($M * $mult))
                                fi
                                
                                if [ $rezhim == 0 ]; then
                                        input_file="file---$C-$mu---$M-$N---($level).txt"
                                        input_base_file="file---$C-$mu---$init_M-$(($init_M * $mult))---(0).txt"
                                else
                                        input_file="file------$mu---$M-$N---($level).txt"
                                        input_base_file="file------$mu---$init_M-$(($init_M * $mult))---(0).txt"
                                fi
                                (echo $M | ./calc_normas.out  $input_file  $input_base_file)   >>  $filename
                                echo -E "\\\\ \hline" >> $filename
                        done
                done
                echo -E "\end{tabular}" >> $filename
                echo "" >> $filename
        done
        echo "" >> $filename
        echo "" >> $filename
done
echo "Скрипт для Cr выполнен"
fi

if [ 0 == 0 ]; then
for C in 1 ; do
        for mu in 0.1 0.01 0.001 ; do
                echo -E "\centering    \begin{tabular}{ |c|c| } \hline    & $ p=\rho^\gamma, \; \mu = $mu $\\\\ \hline" >> $filename
                for level in 1 2 3 -1 ; do
                        for rezhim in 1
                        do
                                if [ $level == -1 ]; then
                                        echo -E " $ v - u $  & " >> $filename
                                else
                                        echo -E " $ v - v^$level $  & " >> $filename
                                fi
                                # эти строки задают параметры, для передачи файла как аргумента
                                if  [ $level == -1 ]; then      M=$init_M;           N=$(($M * $mult))
                                elif [ $level == 0 ]; then      M=$init_M;           N=$(($M * $mult))
                                elif [ $level == 1 ]; then      M=$(($init_M * 3));  N=$(($M * $mult))
                                elif [ $level == 2 ]; then      M=$(($init_M * 9));  N=$(($M * $mult))
                                else                            M=$(($init_M * 27)); N=$(($M * $mult))
                                fi
                                
                                if [ $rezhim == 0 ]; then
                                        input_file="file---$C-$mu---$M-$N---($level).txt"
                                        input_base_file="file---$C-$mu---$init_M-$init_M---(0).txt"
                                else
                                        input_file="file------$mu---$M-$N---($level).txt"
                                        input_base_file="file------$mu---$init_M-$init_M---(0).txt"
                                fi
                                (echo $M | ./calc_normas.out  $input_file  $input_base_file)   >>  $filename
                                echo -E "\\\\ \hline" >> $filename
                        done
                done
                echo -E "\end{tabular}" >> $filename
                echo "" >> $filename
        done
        echo "" >> $filename
        echo "" >> $filename
done
fi

echo "Скрипт закончен полностью"
