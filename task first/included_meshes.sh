# g++ func.cpp && ./a.out 1 1 50 50 0 0 > gnu.gp && tail -n 1 gnu.gp && gnuplot gnu.gp

# вот как это работает: я делаю файлы, каждый со своим названием в зависимости от M, N, C, mu, режима и номера задачи и уровень вложенной сетки. В них скриптом записываются все что на вложенных сетках получается, то есть в итоге там 5 файлов (0, 1, 2, 3 и -1). Далее питон-скрипт (или еще один баш скрипт) пробежится по всем файлам и уже из них составит сами матрицы, чтобы записать их в латехе. M и N возьму одинаковыми и для начала равными 100.
# необходимо написать еще одну программу, которая будет читать файлы, основываясь на введенных аргументах и будет выводить числа равные норме (точнее сразу в виде таблицы латеха).

init_M=$1
rezhim=0
mult=$2
num_of_task=0

echo "$((init_M*mult))"

g++ func.cpp \
&& \
#for C in 1 10 100
for C in 100
do
        for mu in 0.1 0.01 0.001
        do
                for level in -1 0 1 2 3 
                do
                        for rezhim in 1                         # 0 но вообще эта штука, на режиме Cr
                        do
                                if  [ $level == -1 ]; then      M=$init_M;           N=$((M*mult))
                                elif [ $level == 0 ]; then      M=$init_M;           N=$((M*mult))
                                elif [ $level == 1 ]; then      M=$(($init_M * 3));  N=$((M*mult))
                                elif [ $level == 2 ]; then      M=$(($init_M * 9));  N=$((M*mult))
                                else                            M=$(($init_M * 27)); N=$((M*mult))
                                fi
                                
                                if [ $rezhim == 0 ]; then
                                        filename="file---$C-$mu---$M-$N---($level).txt"
                                else
                                        filename="file------$mu---$M-$N---($level).txt"
                                fi
                                echo ""                                         >   $filename
                                ./a.out  $C  $mu  $M  $N  $rezhim  $level       >>  $filename
                        done
                done
                echo "Выполнено для параметров $C $mu"
        done
done

./meshes_sbor.sh $1 $2

echo "Скрипт закончен"
