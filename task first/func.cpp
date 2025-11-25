#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
// для записи в файл
#include <iostream>
#include <fstream>

#define eps 1.e-15

namespace {
enum class mode_my_mode {
	rho_v_stepeny, C_rho
};

double gamma_ = 1.4;
double C = 1;
double mu = 0.1;

double zero(double, double, int) {
	return 0;
}

double con = 1;

double u_head (double t, double x)
{
  return cos (2. * M_PI * t) * sin (M_PI * x * x * con);
}
double rho (double t, double x)
{
  return exp (t) * (cos (M_PI * x * con) + 1.5);
}
double ddx_u (double t, double x)
{
  return M_PI * x * cos (2. * M_PI * t) * cos (M_PI * x * x * con) * 2 * con;
}
double ddx_rho (double t, double x)
{
  return -M_PI * exp (t) * sin (M_PI * x * con) * con;
}
double ddx_rho_u (double t, double x)
{
  return ddx_u (t, x) * rho (t, x) + ddx_rho (t, x) * u_head (t, x);
}
double ddt_u (double t, double x)
{
  return -2. * M_PI * sin (2. * M_PI * t) * sin (M_PI * x * x * con);
}
double d_2dx_2_u (double t, double x)
{
  /* double constant = M_PI * cos (2. * M_PI * t) * 2 * con;
  double first = cos (M_PI * x * x * con);
  double second = M_PI * x * x * sin (M_PI * x * x * con) * 2 * con; */
  return cos (2 * M_PI * t) * (-sin (M_PI * x * x * con) * 2 * x * 2 * x * con * con * M_PI * M_PI       +       cos (M_PI * x * x * con) * 2 * M_PI * con);
  //return constant * first - constant * second;
}



double f0 (double t, double x)
{
  return rho (t, x) + ddx_rho_u (t, x);
}

double f (double t, double x, double mu)
{
 return (rho (t, x) * ddt_u (t, x) + rho (t, x) * u_head (t, x) * ddx_u (t, x) +
     gamma_ * std::pow (rho (t, x), gamma_ - 1) * ddx_rho (t, x) -
     -mu * d_2dx_2_u (t, x)) / rho (t, x);
}
double f_linear (double t, double x, double mu) {
	return (rho(t, x) * ddt_u(t, x)  +  rho(t, x) * u_head(t, x) * ddx_u(t, x)  +  C * ddx_rho(t, x) - mu * d_2dx_2_u(t, x)) / rho(t, x);
}
bool is_equal (double x, double y)
{
  return fabs (x - y) < 1e-14;
}

double u(double time, double x, int mode) {
	switch (mode) {
	case 0:
		return cos(M_PI * time) * sin(M_PI * x);
	case 1:
		return 0;
	case 2:
		return 1 - 2 * time;
	case 3:
		return 1;
	case 4:
		return 0.5 - x;
	case 5:
		return 0.5 - x;
	case 6:
		return 0.5;
	case 7:
		return u_head(time, x);
	case 8:
		return exp(time) * cos(2 * M_PI * x);

	default:
		return 0;
		break;
	}
	//return cos(2 * M_PI * time) * sin(4 * M_PI * x);
}

double r(double time, double x, int mode) {
	switch (mode) {
	case 0:
		return exp(time) * (cos(M_PI * x) + 1.5);
	case 1:
		return 1 + x;
	case 2:
		return 2 - time;
	case 3:
		return 1 + time * x;
	case 4:
		return 1;
	case 5:
		return 1;
	case 6:
		return (0.2 <= x && x <= 0.4) ? 1 : 0;
	case 7:
		return rho(time, x);
	case 8:
		return exp(2 * time) * (cos(3 * M_PI * x) + 1.5);

	default:
		return 0;
		break;
	}
	//return exp(time) * (cos(3 * M_PI * x) + 1.5);
}

double f_1(double time, double x, int mode) {
	switch (mode) {
	case 0:
			return r(time, x, mode) * (1 + M_PI * cos(M_PI * time) * cos(M_PI * x)) - M_PI * u(time, x, mode) * exp(time) * sin(M_PI * x);
	case 1:
		return 0;
	case 2:
		return -1;
	case 3:
		return x + time;
	case 4:
		return -1;
	case 5:
		return -1;
	case 6:
		return 0;
	case 7:
		return f0(time, x);
	case 8:
		return 2 * exp(2 * time) * (cos(3 * M_PI * x) + 1.5) - exp(time) * (2 * M_PI * sin(2 * M_PI * x) * exp(2 * time) * (cos(3 * M_PI * x) + 1.5)) - exp(3 * time) * 3 * M_PI * cos(2 * M_PI * x) * sin(3 * M_PI * x);

	default:
		return 0;
		break;
	}
}

double f_2(double time, double x, int mode, mode_my_mode rezhym) {
	switch (mode) {
	case 0:
		if (rezhym == mode_my_mode::C_rho) return (-gamma_ * M_PI * exp(time) * sin(M_PI * x) + 
			mu * M_PI * M_PI * cos(M_PI * time) * sin(M_PI * x) +
			M_PI * exp(time) * (cos(M_PI * x) + 1.5) * 
			cos(M_PI * time) * cos(M_PI * x) * cos(M_PI * time) * sin(M_PI * x) -
			M_PI * exp(time) * sin(M_PI * time) * sin(M_PI * x) * (cos(M_PI * x) + 1.5)) / 
			(exp(time) * (cos(M_PI * x) + 1.5));
		else ;
	case 1:
		if (rezhym == mode_my_mode::C_rho) return C / (1 + x);
		else return gamma_ / pow(1 + x, gamma_ - 1);
	case 2:
		return (4 * time - 5) / (2 - time);
	case 3:
		if (rezhym == mode_my_mode::C_rho) return (x + time + C * time) / (time * x + 1);
		else return (x + time + gamma_ * time * pow(1 + time * x, gamma_ - 1)) / (time * x + 1);
	case 4:
		return (-1 + 2 * x - 2 * mu) / 1;
	case 5:				// этот кейс такой же, но считается по другой формуле (по второй формуле)
		return (x - 0.5 - 2 * mu) / 1;
	case 6:
		return 0;
	case 7:
		if (rezhym == mode_my_mode::rho_v_stepeny) return f(time, x, mu);
		else return f_linear(time, x, mu);
	case 8:
		if (rezhym == mode_my_mode::C_rho) return ( 	(cos(3 * M_PI * x) + 1.5) * exp(time) * cos(2 * M_PI * x) * exp(2 * time) 
						- 2 * M_PI * exp(4 * time) * cos(2 * M_PI * x) * (cos(3 * M_PI * x) + 1.5) * sin(2 * M_PI * x) 
						+ C * 3 * M_PI * (-sin(3 * M_PI * x)) * exp(2 * time) 
						+ mu * 4 * M_PI * M_PI * cos(2 * M_PI * x) * exp(time) ) 
						/ (exp(2 * time) * (cos(3 * M_PI * x) + 1.5));
		else return ( 	(cos(3 * M_PI * x) + 1.5) * exp(time) * cos(2 * M_PI * x) * exp(2 * time) 
						- 2 * M_PI * exp(4 * time) * cos(2 * M_PI * x) * (cos(3 * M_PI * x) + 1.5) * sin(2 * M_PI * x) 
						+ exp(2 * time) * gamma_ * pow((cos(3 * M_PI * x) + 1.5), gamma_ - 1) * 3 * M_PI * (-sin(3 * M_PI * x)) 
						+ mu * 4 * M_PI * M_PI * cos(2 * M_PI * x) * exp(time) ) 
						/ (exp(2 * time) * (cos(3 * M_PI * x) + 1.5));

	default:
		return 0;
		break;
	}
}
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
void set_array_H (double *H, double time, double h, int M, int mode) {
	for (int i = 0; i < M; i++) {
		H[i] = r(time, (i + 0.5) * h, mode);
	}
}

void set_array_V (double *V, double time, double h, int M, int mode) {
	for (int i = 0; i <= M; i++) {
		V[i] = u(time, i * h, mode);
	}
}

void print_array(double *array, int n) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i]);
	printf("\n");
}
void print_array__V(double *array, double time, int n, int mode) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i] - u(time, i / (n - 1), mode));
	printf("\n");
}
void print_array__H(double *array, double time, int n, int mode) {
	for (int i = 0; i < n; i++) printf(" %3.7le ", array[i] - r(time, (i + 0.5) / n, mode));
	printf("\n");
}

void print_res_H_like_gnu(double *z, double (*f)(double, double, int), double time_x, double h, int M, int mode) {
	for (int i = 0; i < M; i++) {
		//printf("%lf  %lf  %lf \n", time_x, i + h / 2, fabs(z[i] - f(time_x, i + h / 2, mode)));
		printf("%lf  %lf  %lf \n", time_x, (i + 0.5) * h, z[i] - f(time_x, (i + 0.5) * h, mode));
	}
	printf("\n");
}

void print_res_V_like_gnu(double *z, double (*f)(double, double, int), double time_x, double h, int M, int mode) {
	for (int i = 0; i < M; i++) {
		//printf("%lf  %lf  %lf \n", time_x, i + h / 2, fabs(z[i] - f(time_x, i + h / 2, mode)));
		printf("%lf  %lf  %lf \n", time_x, i * h, z[i] - f(time_x, i * h, mode));
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
int set_matrix_H (double *H_under, double *H_main, double *H_above, double *V, double *right_part, double *H_prev, double tau, double h, int timestep, int M, int mode) {
	H_above[0] = 0.5 * (V[1] - fabs(V[1])) / h * tau;
	H_main[0]  = 1 + 0.5 * (V[1] + fabs(V[1])) / h * tau;
	H_under[0] = 0;
	right_part[0] = H_prev[0] + f_1((timestep) * tau, 0.5 * h, mode) * tau;

	for (int i = 1; i <= M - 2; i++) {
		H_above[i] = 0.5 * 1 / h * (V[i + 1] - fabs(V[i + 1])) * tau;
		H_main[i]  = 1 +  0.5 * 1 / h * (V[i + 1] + fabs(V[i + 1]) - V[i] + fabs(V[i])) * tau;
		H_under[i] = (-0.5) * 1 / h * (V[i] + fabs(V[i])) * tau;

		right_part[i] = H_prev[i] + f_1((timestep) * tau, (i + 0.5) * h, mode) * tau;
	}

	H_above[M - 1] = 0;
	H_main[M - 1]  = 1 + 0.5 * (-V[M - 1] + fabs(V[M - 1])) / h * tau;
	H_under[M - 1] = (-0.5) * (V[M - 1] + fabs(V[M - 1])) / h * tau;
	right_part[M - 1] = H_prev[M - 1] + f_1((timestep) * tau, (M - 0.5) * h, mode) * tau;

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

int set_matrix_V (double *V_under, double *V_main, double *V_above, double *H, double *right_part, double *H_prev, double *V_prev, double tau, double h, int M, int time_step, int mode,mode_my_mode rezhym) {
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
				right_part[i] = 0.5 * (H[i] + H[i - 1]) * V_prev[i] / tau  +  0.5 * 1 * (H[i] + H[i - 1]) * f_2((time_step) * tau, i * h, mode, mode_my_mode::rho_v_stepeny) - (gamma_ / (gamma_ - 1)) * (1 / h) * 0.5 * (H[i] + H[i - 1]) * (pow(H[i], gamma_ - 1) - pow(H[i - 1], gamma_ - 1));
			}
			else if (rezhym == mode_my_mode::C_rho) {
				right_part[i] = 0.5 * (H[i] + H[i - 1]) * V_prev[i] / tau  +  0.5 * 1 * (H[i] + H[i - 1]) * f_2((time_step) * tau, i * h, mode, mode_my_mode::C_rho) - C * (1 / h) * (H[i] - H[i - 1]);
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


void residuals(double (*rho_func)(double, double, int), double *H, double time, double h, int M, int mode, double &r1, double &r2, double &r3, double &r4) {
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
	int M, N, rezhim, mode;	// M = шагов по пространству, N = шагов по времени, mode = 0 или 1
	double h, tau;			// h - шаг по пространству, tau - шфг по времени
	double x, time;			// x - длина отрезка (считаем что он равен 1), time - временной отрезок
	double length_on_space = 1.0, length_on_time = 1.0;
	mode_my_mode rezh;

	if (!(argc == 7 
		&& sscanf(argv[1], "%lf", &x) == 1
		&& sscanf(argv[2], "%lf", &time) == 1
		&& sscanf(argv[3], "%d",  &M) == 1
		&& sscanf(argv[4], "%d",  &N) == 1
		&& sscanf(argv[5], "%d",  &rezhim) == 1
		&& sscanf(argv[6], "%d",  &mode) == 1
		&& x > eps && time > eps && N >= 1 && M >= 1 
		&& fabs(rezhim - 0.5) < 0.6))
		{
			printf("# Usage: ./a.out  x  time  M  N  rezhim  mode\n");
			return -1;
		}

	//x = length_on_space;
	//time = length_on_time;
	//rezhim = 1;

	int kakaya_norma = 0;
	std::cin >> kakaya_norma;
	
	//char stroim_graphic_of = 'V';		// эта переменная для того чтобы понять правильно ли задаются преобразования для V
	char stroim_graphic_of = 'H';
	
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

	
	if (kakaya_norma == 1 || kakaya_norma == 2 || kakaya_norma == 3 || kakaya_norma == 4 || kakaya_norma == 5) {
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
	else
		switch (stroim_graphic_of) {
		case 'V':
		{
			printf("# h = %lf и tau = %lf\n", h, tau);
			printf("# Введенное время и пространство ни на что не влияют\n");
			printf("\n\n$POV << EOD\n");
			//print_res_V_like_gnu(V_prev, u, 0, h, M, mode);

			// в цикле функциями от set_matrix решаем и запоминаем решение, оно понадобится дальше
			for (int i = 1; i < N; i++) {
				set_array_H(H, i * tau, h, M, mode);
				//set_matrix_H(H_under, H_main, H_above, V_prev, H, H_prev, tau, h, i, M, mode);

				//set_array_V(V, i * tau, h, M + 1, mode);
				set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, i, mode, rezh);
				print_res_V_like_gnu(V, u, (i + 1) * tau, h, M, mode);
				memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
			}
			residuals(u, V, N * tau, h, M, mode, res1, res2, res3, res4);

			printf("EOD\n");
			printf("set xlabel \"ВРЕМЯ\" textcolor rgb \"red\"\nset ylabel \"ПРОСТРАНСТВО\"\n");
			printf("set title \"График поверхности\"\n");
			printf("splot $POV using 1:2:3 with lines\n");
			printf("pause -1\n");
			//printf("# Итоговый ответ:  ");	print_array(H, M);
			printf("# разность = %lf, время = %lf\n", res1, (clock() - timing) / CLOCKS_PER_SEC);

			break;
		}
		case 'H':
		{
			printf("# h = %lf и tau = %lf\n", h, tau);
			printf("# Введенное время и пространство ни на что не влияют\n");
			printf("\n\n$POV << EOD\n");
			//print_res_H_like_gnu(H_prev, r, 0, h, M, mode);

			// в цикле функциями от set_matrix решаем и запоминаем решение, оно понадобится дальше
			for (int i = 0; i < N; i++) {
				//set_array_H(H, i * tau, h, M, mode);
				set_matrix_H(H_under, H_main, H_above, V_prev, H, H_prev, tau, h, i, M, mode);
				print_res_H_like_gnu(H, r, (i + 1) * tau, h, M, mode);

				set_array_V(V, i * tau, h, M + 1, mode);
				//set_matrix_V(V_under, V_main, V_above, H, V, H_prev, V_prev, tau, h, M, i, mode, rezh);
				memcpy(H_prev, H, M * sizeof(double)); memcpy(V_prev, V, (M + 1) * sizeof(double));
			}
			residuals(r, H, N * tau, h, M, mode, res1, res2, res3, res4);

			printf("EOD\n");
			printf("set xlabel \"ВРЕМЯ\" textcolor rgb \"red\"\nset ylabel \"ПРОСТРАНСТВО\"\n");
			printf("set title \"График поверхности\"\n");
			printf("splot $POV using 1:2:3 with lines\n");
			printf("pause -1\n");
			//printf("# Итоговый ответ:  ");	print_array(H, M);
			printf("# разность = %lf, время = %lf\n", res1, (clock() - timing) / CLOCKS_PER_SEC);

			break;
		}
		}
	return 0;
}