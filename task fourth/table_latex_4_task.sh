# вот как вызывать для параметров C=10 и mu = 0.001
# table_latex_4_task.sh 10 0.001

C=$1
mu=$2
filename="TABLE_rho=$2.txt"

echo -E "\begin{table}[htbp] \centering \begin{tabular}{|c * {7}{|c}|} \hline \multicolumn{8}{|c|}{$ C=$C, \; \mu = $mu$} \\\\ \hline \diagbox{$\rho$}{u} & 1 & 2 & 3 & 4 & 5 & 6 & 7 \\\\ \hline"                                                    >  $filename
for rho in {1..7}
do
        for u in {1..7}
        do
                if [[ $u == 1 ]]
                then
                        echo -n "$rho "                                                         >> $filename
                fi
                echo -n "& $(echo -1 | ./a.out 1000 100000 1      $u $rho      $C $mu 10) "      >> $filename
        done
        echo -E "\\\\ \hline"                                                                   >> $filename
done
echo -E "\end{tabular}"                                                                         >> $filename
echo -E "\end{table}"                                                                           >> $filename

echo "для параметров C=$C и mu=$mu скрипт окончен"

