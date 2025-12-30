#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

double C = 1, mu = 0.1;
double con = 1;

double rho (double t, double x)
{
//return exp (t) * (cos (M_PI * x * con) + 1.5);
	return 0;
}

double L2_norm (std::vector<double> v, double h, int st, int M)
{
double scal = 0.;
//int size = static_cast<int> (v.size ());
for (int i = st; i < M - 1; i++) 		scal += v[i] * v[i];
for (int i = 0; i < st; i++) 			scal += (1./2. * v[i] * v[i]);
scal += (v[M - 1] * v[M - 1] / 2.);
return sqrt (h * scal);
}

double W2_1h_norm (std::vector<double> v, double h, int st, int M)
{
double first = L2_norm (v, h, st, M);
double second = 0.;
//int size = static_cast<int> (v.size ());
for (int i = 1; i < M; i++)
	{
	second += (v[i] - v[i - 1]) * (v[i] - v[i - 1]);
	}
second /= h;
return sqrt (first * first + second);
}

void residuals(std::vector<double> H, double time, double h, int M, double &r1, double &r2, double &r3) {
	r1 = fabs(H[0]);
	double buf;
	for (int i = 1; i < M; i++) {
		buf = H[i] - rho (time, (i + 0.5) * h);
		r1 = fmax(r1, fabs(buf));
		r2 += buf * buf;
	}
	r2 = sqrt(r2 * h);
	r3 = W2_1h_norm (H, h, 0, M);
}

void substitute(std::vector<double> &V, std::vector<double> &V_level, int M) {
	for (int i = 0; i < M; i++) V_level[i] -= V[i];
}



int main(int argc, char *argv[]) {
	// задаем начальные данные из введенных значений
	int M, N, rezhim, mode;	// M = шагов по пространству, N = шагов по времени, mode = 0 или 1
	double tau;			// h - шаг по пространству, tau - шфг по времени
	double x = 1, time = 1;			// x - длина отрезка (считаем что он равен 1), time - временной отрезок
	int level;				// уровень вложенности. 0 = обычный режим, и это 0, 1 = вложенность и это с 1, 2 начать с 4, а 3 начать с 13. -1 это будет как-бы точное решение

	std::cin >> M;
	//std::cin >> N;

	//printf("%d %d\n", M, N);

	std::ifstream file(argv[1]);
	if (!file) {
		printf("Файл 1 не открылся\n");
		return 1;
	}
	std::vector<double> arr(M);
	for (int i = 0; i < M && file >> arr[i]; i++);
	file.close();

	//printf("%s\n", argv[1]);

	std::ifstream file__2(argv[2]);
	if (!file__2) {
		printf("Файл 2 не открылся\n");
		return 1;
	}
	std::vector<double> arr__2(M);
	for (int i = 0; i < M && file__2 >> arr__2[i]; i++);
	file__2.close();

	double h = x / M;
	double r1, r2, r3;

	substitute(arr, arr__2, M);

	residuals(arr__2, time, h, M, r1, r2, r3);

	printf(" \\makecell[tl] { %.3le \\\\ %.3le \\\\ %.3le } ", r1, r2, r3);

	return 0;
}