#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#define eps 1e-15

enum class mode_my_mode {
	rho_v_stepeny, C_rho
};

double gamma_ = 1.4;
double C = 1;
double mu = 0.001;

double u(double time, double x) {
	//return cos(2 * M_PI * time) * sin(4 * M_PI * x);
	//return cos(M_PI * time) * sin(M_PI * x);		// было вот это для f
	return 1;
}

double r(double time, double x) {
	//return exp(time) * (cos(3 * M_PI * x) + 1.5);
	//return exp(time) * (cos(M_PI * x) + 1.5);		// было вот это для f
	//return x + 1 - time;
	return 1;
}

double f_0(double time, double x) {
	return r(time, x) * (1 + M_PI * cos(M_PI * time) * cos(M_PI * x)) - M_PI * u(time, x) * exp(time) * sin(M_PI * x);
}

double f(double time, double x) {
	/* return (-gamma_ * M_PI * exp(time) * sin(M_PI * x) + 
			mu * M_PI * M_PI * cos(M_PI * time) * sin(M_PI * x) +
			M_PI * exp(time) * (cos(M_PI * x) + 1.5) * 
			cos(M_PI * time) * cos(M_PI * x) * cos(M_PI * time) * sin(M_PI * x) -
			M_PI * exp(time) * sin(M_PI * time) * sin(M_PI * x) * (cos(M_PI * x) + 1.5)) / 
		(exp(time) * (cos(M_PI * x) + 1.5)); */
	//return C / (x + 1 - time);
	return 0;
}

/* inline double f_exp(double t, double x) {

	const double pi = M_PI;
	const double sin_px = sin(pi * x);
	const double cos_px = cos(pi * x);
	const double sin_pt = sin(pi * t);
	const double cos_pt = cos(pi * t);
	const double exp_t = exp(t);
	const double rho = r(t, x); // = exp_t * (cos_px + 1.5)


	// 1. Конвективный член: (r u_t + r u u_x) / r
	double term1 = pi * sin_px * (-sin_pt + cos_pt * cos_pt * cos_px);

	// 2. Градиент давления: -γ ∂p/∂x / r
	double term2 = 0.0;
	if (gamma_ > 1.0 + 1e-12) {  // Защита от деления на ноль
		const double exp_gamma_minus_1_t = exp((gamma_ - 1.0) * t);
		const double r_pow_gamma_minus_2 = pow(cos_px + 1.5, gamma_ - 2.0);
		term2 = (C * gamma_ * gamma_ * pi / (gamma_ - 1.0)) * exp_gamma_minus_1_t * r_pow_gamma_minus_2 * sin_px;
	}

	// 3. Вязкостный член: μ u_xx / r
	double term3 = mu * pi * pi * cos_pt * sin_px / rho;

	// Сумма всех членов
	return term1 + term2 + term3;
} */

/* double f_linear(double t, double x) {
	double pi = M_PI;
	double sin_px = sin(pi * x);
	double cos_px = cos(pi * x);
	double sin_pt = sin(pi * t);
	double cos_pt = cos(pi * t);
	double exp_t = exp(t);
	double rho = r(t, x);  // = exp_t * (cos_px + 1.5)

	// 1. Конвективный член (r u_t + r u u_x) / r
	double term1 = pi * sin_px * (-sin_pt + cos_pt * cos_pt * cos_px);

	// 2. Градиент давления (только если gamma_ > 1)
	double term2 = 0.0;
	if (gamma_ > 1.0 + eps) {
		double exp_gamma_minus_1_t = exp((gamma_ - 1.0) * t);
		double r_pow_gamma_minus_2 = pow(cos_px + 1.5, gamma_ - 2.0);
		term2 = (C * gamma_ * gamma_ * pi / (gamma_ - 1.0)) * exp_gamma_minus_1_t * r_pow_gamma_minus_2 * sin_px;
	}

	// 3. Вязкостный член (μ u_xx)
	double term3 = mu * pi * pi * cos_pt * sin_px / rho;

	return term1 + term2 + term3;
} */

// функция задающая сами значения массива, которые используются для задания матрицы
void set_array_H (double *H, double (*r)(double, double), double time, double h, int M) {
	for (int i = 0; i < M; i++) {
		H[i] = r(time, (i + 0.5) * h);
	}
}

void set_array_V (double *V, double (*u)(double, double), double time, double h, int M) {
	for (int i = 0; i <= M; i++) {
		V[i] = u(time, i * h);
	}
}

void print_array(double *array, int n) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i]);
	printf("\n");
}
void print_array__V(double *array, double time, int n) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i] - u(time, i / (n - 1)));
	printf("\n");
}
void print_array__H(double *array, double time, int n) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i] - r(time, (i + 0.5) / n));
	printf("\n");
}

// функция правильная, проверена работает для 0 в левом верхнем углу матрицы
void solver (double *p_half, double *phalf, double *phalf1, double *pprev, int N){ // pprev - это правая часть, и вместе с тем туда пишется решение
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
		
		phalf1[0] /= phalf[0];
		pprev[0] /= phalf[0];
		if (flag) buf /= phalf[0];
		phalf[0] = 1;

		for (int i = 1; i <= N - 1; i++){
			// p_half[i] -= p_half[i] * phalf[i - 1];
			pprev[i]  -= p_half[i] * pprev[i - 1];
			phalf[i]  -= p_half[i] * phalf1[i - 1];
			p_half[i] = 0;
			// сейчас элемент p_half[i] = 0
			phalf1[i] /= phalf[i];
			pprev[i]  /= phalf[i];
			phalf[i]  =  1;
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
}



// функция задающая элементы матрицы под, на и над главной диагональю, n - это размер матрицы
int set_matrix_H (double *H_under, double *H_main, double *H_above, double *V, double *right_part, double *H_prev, double tau, double h, int M) {
	for (int i = 0; i < M; i++) {
		H_above[i] = (i == M - 1) ? 0 : 0.5 * tau / h * (V[i + 1] - fabs(V[i + 1]));
		H_main[i]  = 1 + 0.5 * tau / h * (V[i + 1] + fabs(V[i + 1]) - V[i] + fabs(V[i]));
		H_under[i] = (i == 0) ? 0 : (-0.5) * tau / h * (V[i] + fabs(V[i]));

		right_part[i] = H_prev[i];
	}

	// далее решаем прогонкой, а в правой части получаем решение
	solver(H_under, H_main, H_above, right_part, M);
	//print_array(right_part, M);
	// далее пишем 0 в края
	return 0;
}

int set_matrix_V (double *V_under, double *V_main, double *V_above, double *H, double *right_part, double *H_prev, double *V_prev, double (*f)(double, double), double tau, double h, int M, int time_step, mode_my_mode mode) {
	for (int i = 1; i < M; i++) {
		// если сумма равна 0, то надо задать коэффициент 1, а числа слева и справа равны 0
		if (fabs(H[i] + H[i - 1]) < eps) {
			V_above[i] = 0;
			V_main[i] = 1;
			V_under[i] = 0;
			right_part[i] = 0;
		}
		else {
			V_above[i] =  0.5 * (H[i] + H[i - 1]) * (tau / h) * 0.5 * (V_prev[i] - fabs(V_prev[i]))  -  mu * tau / (h * h);
			V_main[i]  =  0.5 * (H[i] + H[i - 1]) * (1  +  (tau / h) * fabs(V_prev[i]))              +  2 * mu * tau / (h * h);
			V_under[i] = -0.5 * (H[i] + H[i - 1]) * (tau / h) * 0.5 * (V_prev[i] + fabs(V_prev[i]))  -  mu * tau / (h * h);

			right_part[i] = 0.5 * (H[i] + H[i - 1]) * V_prev[i]  +  0.5 * tau * (H[i] + H[i - 1]) * f((time_step + 1) * tau, i * h);			////////////////////////////////// тут надо правильно установить функцию f
			
			if (mode == mode_my_mode::rho_v_stepeny) {
				right_part[i] -= (gamma_ / (gamma_ - 1)) * (tau / h) * 0.5 * (H[i] + H[i - 1]) * (pow(H[i], gamma_ - 1) - pow(H[i - 1], gamma_ - 1));
			}
			else if (mode == mode_my_mode::C_rho) {
				right_part[i] -= C * (tau / h) * (H[i] - H[i - 1]);
			}
		}
	}
	// тут решаем
	solver(V_under + 1, V_main + 1, V_above + 1, right_part + 1, M - 1);
	right_part[0] = 0; right_part[M] = 0;
	return 0;
}

void residuals(double (*rho)(double, double), double *H, double time, double h, int M, double &r1) {
	r1 = fabs(H[0] - rho(time, 0.5 * h));
	for (int i = 1; i < M; i++) {
		r1 = fmax(r1, fabs(H[i] - rho(time, (i + 0.5) * h)));
	}
}

int main(int argc, char *argv[]) {
	// задаем начальные данные из введенных значений
	int M, N, mode;			// M = шагов по пространству, N = шагов по времени, mode = 0 или 1
	double h, tau;			// h - шаг по пространству, tau - шфг по времени
	double x, time;			// x - длина отрезка (считаем что он равен 1), time - временной отрезок
	mode_my_mode rezh;

	if (!(argc == 6 
		&& sscanf(argv[1], "%lf", &x) == 1
		&& sscanf(argv[2], "%lf", &time) == 1
		&& sscanf(argv[3], "%d",  &M) == 1
		&& sscanf(argv[4], "%d",  &N) == 1
		&& sscanf(argv[5], "%d",  &mode) == 1
		&& x > eps && time > eps && N >= 1 && M >= 1 
		&& fabs(mode - 0.5) < 0.6))
		{
			printf("Usage: ./a.out  x  time  M  N  mode\n");
			return -1;
		}


	
	h = x / M;
	tau = time / N;
	printf("h = %lf и tau = %lf\n", h, tau);

	switch (mode) {
	case 0:
		rezh = mode_my_mode::C_rho;
		break;
	case 1:
		rezh = mode_my_mode::rho_v_stepeny;
		break;
	default:
		break;
	}

	double res = 0, timing = 0;
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
	set_array_H(H_prev, r, 0, h, M);
	set_array_V(V_prev, u, 0, h, M + 1);



	// в цикле функциями от set_matrix решаем и запоминаем решение, оно понадобится дальше
	for (int i = 0; i < N; i++) {
		set_matrix_H(H_under, H_main, H_above, V_prev, H, H_prev, tau, h, M);
		/* print_array__H(H, i * tau, M);
		printf("\n"); */
		set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, f, tau, h, M, i, rezh);
		memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
	}

	print_array__H(H, time, M);

	residuals(r, H, time, h, M, res);

	printf("разность = %lf, время = %lf\n", res, (clock() - timing) / CLOCKS_PER_SEC);
	return 0;
}