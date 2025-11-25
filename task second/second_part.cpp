#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
// для записи в файл
#include <iostream>
#include <fstream>
#include <unistd.h>

#define eps 1.e-15

enum class mode_my_mode {
	rho_v_stepeny, C_rho
};
enum class number_task {
	first, second
};

double gamma_ = 1.4;
double C = 10;
double mu = 0.01;

double zero(double, double, int) {
	return 0;
}

bool is_equal (double x, double y)
{
return fabs (x - y) < 1e-14;
}

namespace {
double u(double time, double x, number_task mode) {
	if (mode == number_task::first) return 0;
	else {
		if (x < 4.5 || x > 5.5) return 0;
		else return 1;
	}
}

double r(double time, double x, number_task mode) {
	if (mode == number_task::first) {
		if (x < 4.5 || x > 5.5) return 1;
		else return 2;
	}
	else return 1;
}

double f_1(double, double) {
	return 0;
}

double f_2(double, double) {
	return 0;
}
}



// функция задающая сами значения массива, которые используются для задания матрицы
void set_array_H (double *H, double time, double h, int M, number_task mode) {
	for (int i = 0; i < M; i++) { H[i] = r(time, (i + 0.5) * h, mode); }
}

void set_array_V (double *V, double time, double h, int M, number_task mode) {
	for (int i = 0; i <= M; i++) { V[i] = u(time, i * h, mode); }
}

void print_array(double *array, int n) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i]);
	printf("\n");
}
void print_array__V(double *array, double time, int n, number_task mode) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i] - u(time, i / (n - 1), mode));
	printf("\n");
}
void print_array__H(double *array, double time, int n, number_task mode) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i] - r(time, (i + 0.5) / n, mode));
	printf("\n");
}

void print_H_like_gnu(double *z, double time_x, double h, int M, number_task mode) {
	for (int i = 0; i < M; i++) {
		if ((i % (M / 100) == 0) || (i == M - 1)) printf("%lf  %lf  %lf \n", time_x, i * h, z[i]);
	}
	printf("\n");
}

// функция правильная, проверена работает для 0 в левом верхнем углу матрицы
int solver (double *p_half, double *phalf, double *phalf1, double *pprev, int N){ // pprev - это правая часть, и вместе с тем туда пишется решение
	double buf;
	bool flag = false;
	if ((fabs(phalf[0]) < eps) && (fabs(phalf1[0]) < eps)) printf("\nСистема нерешаема!\n\n");
	else {
		if (fabs(phalf[0]) < eps) {
			flag = true;
			//buf = phalf[0];  phalf[0] = p_half[1]; p_half[1] = buf;
			phalf[0] = p_half[1]; p_half[1] = 0;
			buf = phalf1[0]; phalf1[0] = phalf[1]; phalf[1]  = buf;
			
			buf = pprev[0];  pprev[0] = pprev[1];  pprev[1] = buf;

			buf = phalf1[1];		// запомнили элемент третий элемент в первой строке
			phalf1[1] = 0;
		}
		
		if (fabs(phalf[0]) < eps) {
			printf("Начальный массив имееет фигню\n");
			return -1;
		}
		else {
			phalf1[0] /= phalf[0];
			pprev[0] /= phalf[0];
			if (flag) buf /= phalf[0];
			phalf[0] = 1;
		}

		for (int i = 1; i <= N - 1; i++){
			// p_half[i] -= p_half[i] * phalf[i - 1];
			pprev[i]  -= p_half[i] * pprev[i - 1];
			phalf[i]  -= p_half[i] * phalf1[i - 1];
			p_half[i] = 0;
			// сейчас элемент p_half[i] = 0
			if (fabs(phalf[i]) < eps) {
				printf("Нееееет %d\n", i);
				return -2;
			}
			else {
				phalf1[i] /= phalf[i];
				pprev[i]  /= phalf[i];
				phalf[i]  =  1;
			}
		}
		for (int i = N - 1; i >= 1; i--){
			// pprev[i] = f[i]; // запоминаем плотность на прошлом шаге
			pprev[i - 1] -= pprev[i] * phalf1[i - 1];
			phalf1[i - 1] = 0;
		}

		if (flag == true) {
			pprev[0] -= pprev[2] * buf;
		}
	}
	return 0;
}



// функция задающая элементы матрицы под, на и над главной диагональю, n - это размер матрицы
int set_matrix_H (double *H_under, double *H_main, double *H_above, double *V, double *right_part, double *H_prev, double tau, double h, int timestep, int M, number_task mode) {
	H_above[0] = 0.5 * (V[1] - fabs(V[1])) / h * tau;
	H_main[0]  = 1 + 0.5 * (V[1] + fabs(V[1])) / h * tau;
	H_under[0] = 0;
	right_part[0] = H_prev[0] + f_1((timestep) * tau, 0.5 * h) * tau;

	for (int i = 1; i <= M - 2; i++) {
		H_above[i] = 0.5 * 1 / h * (V[i + 1] - fabs(V[i + 1])) * tau;
		H_main[i]  = 1 +  0.5 * 1 / h * (V[i + 1] + fabs(V[i + 1]) - V[i] + fabs(V[i])) * tau;
		H_under[i] = (-0.5) * 1 / h * (V[i] + fabs(V[i])) * tau;

		right_part[i] = H_prev[i] + f_1((timestep) * tau, (i + 0.5) * h) * tau;
	}

	H_above[M - 1] = 0;
	H_main[M - 1]  = 1 + 0.5 * (-V[M - 1] + fabs(V[M - 1])) / h * tau;
	H_under[M - 1] = (-0.5) * (V[M - 1] + fabs(V[M - 1])) / h * tau;
	right_part[M - 1] = H_prev[M - 1] + f_1((timestep) * tau, (M - 0.5) * h) * tau;

	// далее решаем прогонкой, а в правой части получаем решение
	int result = solver(H_under, H_main, H_above, right_part, M);
	for (int iter = 0; iter < M; iter++)
		if (right_part[iter] < eps) right_part[iter] = 0; //printf("Пипяу, есть отрицательная плотность, этого не может быть\n %d это шаг во времени, а значение %le\n", timestep, right_part[iter]);
	switch (result) {
	case  0:
		break;
	case -1:
	case -2:
		printf("Все плохо\n");
		return -1;
	}
	//right_part[0] = H_prev[0];
	//right_part[M - 1] = H_prev[M - 1];
	//print_array(right_part, M);
	// далее пишем 0 в края
	return 0;
}

int set_matrix_V (double *V_under, double *V_main, double *V_above, double *H, double *right_part, double *H_prev, double *V_prev, double tau, double h, int M, int time_step, number_task mode, mode_my_mode rezhym) {
	V_above[0] = 0;
	V_main[0] = 1;
	V_under[0] = 0;
	//right_part[0] = V_prev[0];
	right_part[0] = 0;
	for (int i = 1; i < M; i++) {
		// если сумма равна 0, то надо задать коэффициент 1, а числа слева и справа равны 0
		if (fabs(H[i] + H[i - 1]) < eps) {
			V_above[i] = 0;
			V_main[i] = 1;
			V_under[i] = 0;
			right_part[i] = 0;
		}
		else {
			V_above[i] =  0.5 * (H[i] + H[i - 1]) * (1 / h) * 0.5 * (V_prev[i] - fabs(V_prev[i]))  -  mu * 1 / (h * h);
			V_main[i]  =  0.5 * (H[i] + H[i - 1]) * (1 / tau  +  (1 / h) * fabs(V_prev[i]))        +  2 * mu * 1 / (h * h);
			V_under[i] = -0.5 * (H[i] + H[i - 1]) * (1 / h) * 0.5 * (V_prev[i] + fabs(V_prev[i]))  -  mu * 1 / (h * h);

			if (rezhym == mode_my_mode::rho_v_stepeny) {
				if (H[i] < 0) {
					printf("Плохое говно в %d\n", i);
				}
				else if (H[i - 1] < 0) {
					printf("Плохое говно в %d\n", i - 1);
				}
				right_part[i] = 0.5 * (H[i] + H[i - 1]) * V_prev[i] / tau  +  0.5 * 1 * (H[i] + H[i - 1]) * f_2((time_step) * tau, i * h) - (gamma_ / (gamma_ - 1)) * (1 / h) * 0.5 * (H[i] + H[i - 1]) * (pow(H[i], gamma_ - 1) - pow(H[i - 1], gamma_ - 1));
			}
			else if (rezhym == mode_my_mode::C_rho) {
				right_part[i] = 0.5 * (H[i] + H[i - 1]) * V_prev[i] / tau  +  0.5 * 1 * (H[i] + H[i - 1]) * f_2((time_step) * tau, i * h) - C * (1 / h) * (H[i] - H[i - 1]);
			}
		}
	}
	V_above[M] = 0;
	V_main[M] = 1;
	V_under[M] = 0;
	//right_part[M] = V_prev[M];
	right_part[0] = 0;

	// тут решаем
	int result = solver(V_under, V_main, V_above, right_part, M + 1);
	switch (result) {
	case  0:
		break;
	case -1:
	case -2:
		printf("Все плохо\n");
		return -1;
	}
	//right_part[0] = 0; right_part[M] = 0;
	return 0;
}

double L2_norm (double *v, double h, int st, int M)
{
double scal = 0.;
//int size = static_cast<int> (v.size ());
for (int i = st; i < M - 1; i++) 		scal += v[i] * v[i];
for (int i = 0; i < st; i++) 			scal += (1./2. * v[i] * v[i]);
scal += (v[M - 1] * v[M - 1] / 2.);
return sqrt (h * scal);
}
double W2_1h_norm (double *v, double h, int st, int M)
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


void residuals(double (*rho_func)(double, double, number_task), double *H, double time, double h, int M, number_task mode, double &r1, double &r2, double &r3, double &r4) {
	r1 = fabs(H[0] - rho_func(time, 0.5 * h, mode));
	double buf;
	for (int i = 1; i < M; i++) {
		buf = H[i] - rho_func(time, (i + 0.5) * h, mode);
		r1 = fmax(r1, fabs(buf));
		r2 += buf * buf;

	}
	r2 = sqrt(r2 * h);
	r3 = L2_norm (H, h, 0, M);
	r4 = W2_1h_norm (H, h, 0, M);
}

int main(int argc, char *argv[]) {
	// задаем начальные данные из введенных значений
	int M, N, rezhim, mode_of_task;	// M = шагов по пространству, N = шагов по времени, mode = 0 или 1
	double h, tau;			        // h - шаг по пространству, tau - шфг по времени
	double x, time;			        // x - длина отрезка (считаем что он равен 1), time - временной отрезок
	double length_on_space = 10.0, length_on_time = 1.0;
	number_task mode;
	mode_my_mode rezh;

	printf("# Usage: ./a.out  x  time  M  N  <rezhim = 0 or 1>  <mode_of_task = 1 or 2>\n");

	if (!(argc == 7 
		&& sscanf(argv[1], "%lf", &x) == 1
		&& sscanf(argv[2], "%lf", &time) == 1
		&& sscanf(argv[3], "%d",  &M) == 1
		&& sscanf(argv[4], "%d",  &N) == 1
		&& sscanf(argv[5], "%d",  &rezhim) == 1
		&& sscanf(argv[6], "%d",  &mode_of_task) == 1
		&& x > eps && time > eps && N >= 1 && M >= 1 
		&& fabs(rezhim - 0.5) < 0.6)){ return -1; }

	//x = length_on_space;
	//time = length_on_time;
	//rezhim = 1;

	int kakaya_norma = 0;
	std::cin >> kakaya_norma;
	//kakaya_norma = 0;
	
	h = x / M;
	tau = time / N;

	switch (rezhim) {
	case 0:
		rezh = mode_my_mode::C_rho;
		break;
	case 1:
		rezh = mode_my_mode::rho_v_stepeny;
		break;
	default:
		break;
	}

	switch (mode_of_task) {
	case 1:
		mode = number_task::first;
		break;
	case 2:
		mode = number_task::second;
		break;
	default:
		printf("Укуси себя за палец\n");
		break;
	}

	double res1 = 0, res2 = 0, res3 = 0, res4 = 0, timing = 0;
	double *H_prev = new double[M];
	double *V_prev = new double[M + 1];
	double *H = new double[M];
	double *V = new double[M + 1];
	double *H_under = new double[M];
	double *H_main  = new double[M];
	double *H_above = new double[M];
	double *V_under = new double[M + 1];
	double *V_main  = new double[M + 1];
	double *V_above = new double[M + 1];

	timing = clock();

	// задаем начальные массивы, они будут служить "предыдущими" массивами
	set_array_H(H_prev, 0, h, M, mode);
	set_array_V(V_prev, 0, h, M + 1, mode);

	if (kakaya_norma < 0) {
		double r1 = 1, buf;	int j = 0;

		set_matrix_H(H_under, H_main, H_above, V, H, H_prev, tau, h, 0, M, mode);

		r1 = fabs(H[0] - H_prev[0]);
		for (int i = 1; i < M; i++) { buf = H[i] - H_prev[i];		r1 = fmax(r1, fabs(buf)); }

		set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, 0, mode, rezh);
		memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
		j += 1;

		do {			// хз почему надо сделать как минимум два шага, так как после первого она нулевая, соответственно надо по-нормальному
			set_matrix_H(H_under, H_main, H_above, V, H, H_prev, tau, h, j, M, mode);
			
			r1 = fabs(H[4] - H_prev[4]);
			for (int i = 5; i < M - 4; i++) { buf = H[i] - H_prev[i];		r1 = fmax(r1, fabs(buf)); }

			set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, j, mode, rezh);
			memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
			j += 1;
		} while (r1 > 0.001);
		printf("%d %le\n", j, r1);
	}
	else if (kakaya_norma == 1 || kakaya_norma == 2 || kakaya_norma == 3 || kakaya_norma == 4 || kakaya_norma == 5) {
		for (int i = 0; i < N; i++) {
			//set_array_H(H, i * tau, h, M, mode);
			int result = set_matrix_H(H_under, H_main, H_above, V_prev, H, H_prev, tau, h, i, M, mode);
			if (result != 0) {
				printf("Плохое решение плотности на шаге %d\n", i);
				return -1;
			}
			//set_array_V(V, i * tau, h, M + 1, mode);
			result = set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, i, mode, rezh);
			if (result != 0) {
				printf("Плохое решение скорости на шаге %d\n", i);
				return -2;
			}
			memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
		}
		residuals(r, H, N * tau, h, M, mode, res1, res2, res3, res4);
		//residuals(u, V, N * tau, h, M + 1, mode, res1, res2, res3, res4);
		switch (kakaya_norma) {
		case 1:
			printf("%le", res1);
			break;
		case 2:
			printf("%le", res2);
			break;
		case 3:
			printf("%le", res3);
			break;
		case 4:
			printf("%le", res4);
			break;
		case 5:
			printf("%0.3le,  %0.3le", res1, res2);
			break;
		}
	}
	else if (kakaya_norma == 0) {				// РИСОВАНИЕ КАРТИНКИ
		printf("# h = %lf и tau = %lf\n", h, tau);
		printf("# Введенное время и пространство ни на что не влияют\n");

		printf("set terminal png size 1000, 600 \n" \
				"set encoding utf8 \n" \
				"set output 'gnuplot_heatmap_%d.png' \n" \
				"set palette viridis \n" \
				"set colorbox vertical default \n" \
				"set rmargin at screen 0.9 \n" \
				"set yrange [0:10] \n" \
				"unset border \n" \
				"$POV << EOD\n", mode_of_task);

		for (int i = 0; i < N; i++) {
			set_matrix_H(H_under, H_main, H_above, V, H, H_prev, tau, h, i, M, mode);
			if (i % (N / 400) == 0) print_H_like_gnu(H, (i + 1) * tau, h, M, mode);
			set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, i, mode, rezh);
			memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
		}

		printf("EOD\n");

		if (rezh == mode_my_mode::C_rho) printf("set bmargin 8 \n" \
												"set label \"C = %d,  \u00B5 = %g\" at screen 0.5, screen 0.05 center font \"Open Sans, 20\" \n", (int)C, mu);
		else printf("set bmargin 8 \n" \
					"set label \"\u00B5 = %g\" at screen 0.5, screen 0.05 center font \"Open Sans, 20\" \n", mu);

		printf("set xlabel \"ВРЕМЯ\" textcolor rgb \"red\"\nset ylabel \"ПРОСТРАНСТВО\" \n" \
				"plot $POV using 1:2:3 with image \n");

		printf("# время = %lf\n", (clock() - timing) / CLOCKS_PER_SEC);
	}
	else {											// РИСОВАНИЕ ГРАФИКА ПОВЕРХНОСТИ
		printf("$POV << EOD\n");

		for (int i = 0; i < N; i++) {
			set_matrix_H(H_under, H_main, H_above, V, H, H_prev, tau, h, i, M, mode);
			if (i % (N / 400) == 0) print_H_like_gnu(H, (i + 1) * tau, h, M, mode);
			set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, i, mode, rezh);
			memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
		}

		printf("EOD\n");
		printf("set xlabel \"ВРЕМЯ\" textcolor rgb \"red\"\nset ylabel \"ПРОСТРАНСТВО\"\n");
		printf("set title \"График поверхности\"\n");
		printf("splot $POV using 1:2:3 with lines\n");
		printf("pause -1\n");
	}
	return 0;
}