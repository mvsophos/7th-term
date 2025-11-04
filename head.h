#include <iostream>
#include <cmath>
/*
double u_head (double t, double x)
{
  return cos (2. * M_PI * t) * sin (M_PI * x * x / 100.);
}
double rho (double t, double x)
{
  return exp (t) * (cos (M_PI * x / 10.) + 1.5);
}
double ddx_u (double t, double x)
{
  return M_PI * x * cos (2. * M_PI * t) * cos (M_PI * x * x / 100.) / 50.;
}
double ddx_rho (double t, double x)
{
  return -M_PI * exp (t) * sin (M_PI  * x / 10.) / 10.;
}
double ddx_rho_u (double t, double x)
{
  return ddx_u (t, x) * rho (t, x) + ddx_rho (t, x) * u_head (t, x);
}
double ddt_u (double t, double x)
{
  return -2. * M_PI * sin (2. * M_PI * t) * sin (M_PI * x * x / 100.);
}
double d_2dx_2_u (double t, double x)
{
  double constant = M_PI * cos (2. * M_PI * t) / 50.;
  double first = cos (M_PI * x * x / 100.);
  double second = M_PI * x * x * sin (M_PI * x * x / 100.) / 50.;

  return constant * first - constant * second;
}
*/

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

double f (double t, double x, double mu, double gamma)
{
 return (rho (t, x) * ddt_u (t, x) + rho (t, x) * u_head (t, x) * ddx_u (t, x) +
     gamma * std::pow (rho (t, x), gamma - 1) * ddx_rho (t, x) -
     -mu * d_2dx_2_u (t, x)) / rho (t, x);
}

bool is_equal (double x, double y)
{
  return fabs (x - y) < 1e-14;
}