#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {

    if (argc != 6) {
        cout << "FUCK OFF" << endl;
        return -1;
    }

    int M;
    if (sscanf(argv[5], "%d", &M) != 1) cout << "HUYNA A NE M" << endl;

    ifstream file_1(argv[1]);
    
    if (!file_1) {
        cerr << "Ошибка открытия файла!" << endl;
        return 1;
    }

    vector<double> numbers_of_obsch((istream_iterator<double>(file_1)), istream_iterator<double>());
    
    for (int yu : {1, 4, 13}) {


        // Чтение всех чисел из файла в вектор одной строкой
        int id_file = -1;
        if (yu == 1)        id_file = 2;
        else if (yu == 4)   id_file = 3;
        else if (yu == 13)  id_file = 4;
        else                cout << "Hernya s yu" << endl;
        ifstream file(argv[id_file]);
        vector<double> numbers((istream_iterator<double>(file)), istream_iterator<double>());

        int M = numbers.size ();
        int shag = 0;

        switch (yu) {
        case 1:
            shag = 3;
            break;
        case 4:
            shag = 9;
            break;
        case 13:
            shag = 27;
            break;
        default:
            cout << "Xuynu vvel" << endl;
            break;
        }

        double norma_C = 0;
        double norma_L2 = 0;
        double h = 1 / (int)pow(3, id_file - 1) / M;
        for (int i = yu; i < M; i += shag) {
            norma_C = fmax(fabs(numbers_of_obsch[i] - numbers[i]) , norma_C);
            norma_L2 += h * fabs(numbers_of_obsch[i] - numbers[i]);
        }

        printf("%le  %le\n", norma_C, norma_L2);
        //cout << norma_C << "  " << norma_L2 << endl;
    }
    
    return 0;
}