#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) :
_df(data_file)
{
   std::cout << "Read functions of initial condition, source term and boundary conditions." << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;
}


// double Function::Initial_condition(const double x, const double y) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       return 15;
//    }
//    else if (this->_df->Get_scenario() == "advection_hom_neumann")
//    {
//       double x0(0.5), y0(0.5), a(100);
//       return 0.5*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)));
//    }
//    else { return Exact_solution(x,y,0); }
// }

double Function::Initial_condition1(const double x, const double y) const
{
   if (x<_df->Get_xb()){return _df->Get_h0g();}
   else {return _df->Get_h0d();}
}

double Function::Initial_condition2(const double x, const double y) const
{
   return 0;
}

double Function::Initial_condition3(const double x, const double y) const
{
   return 0;
}

// double Function::Source_term(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_hom_neumann")
//    {
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double x0(0.5), y0(0.5), a(5);
//       return 12.5*M_PI*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))*cos(0.25*M_PI*t)-this->_mu*(50.0*pow(a, 2)*pow(-2*x0+2*x, 2)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+50.0*pow(a, 2)*pow(-2*y0+2*y, 2)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-100.0*a*(-2*x0+2*x)*pow(x-xmax, 2)*(2*x-2*xmin)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-50.0*a*(-2*x0+2*x)*pow(x-xmin, 2)*(2*x-2*xmax)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-a*(-2*x0+2*x)*pow(x-xmin, 2)*(100.0*x-100.0*xmax)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-100.0*a*(-2*y0+2*y)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*(2*y-2*ymin)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-100.0*a*(-2*y0+2*y)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymin, 2)*(2*y-2*ymax)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-200.0*a*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(x-xmin, 2)*(2*y-2*ymax)*(2*y-2*ymin)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+2*(2*x-2*xmin)*(100.0*x-100.0*xmax)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2))));
//    }
//    else if (this->_df->Get_scenario() == "diffusion_all_BC")
//    {
//       double x0(0.5), y0(0.5), a(5);
//       return 0.125*M_PI*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)))*cos(0.25*M_PI*t)-this->_mu*(pow(a, 2)*pow(2*x-2*x0, 2)*(0.5*sin(0.25*M_PI*t) + 0.5)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2))) + pow(a, 2)*pow(2*y-2*y0, 2)*(0.5*sin(0.25*M_PI*t) + 0.5)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)))-4*a*(0.5*sin(0.25*M_PI*t) + 0.5)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2))));
//    }
//    else if (this->_df->Get_scenario() == "advection_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_advection_all_BC")
//    {
//       double x0(0.9), y0(0.1), a(50);
//       return-this->_mu*(pow(a, 2)*pow(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0), 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + pow(a, 2)*pow(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0), 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(288*pow(t, 2)*pow(y, 2) + (20*t*pow(y, 2) + 1)*(40*t*pow(y, 2) + 2))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(3200*pow(t, 2)*pow(x, 2)*pow(y, 2) + 80*t*x*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x + 1))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)))) - this->_mu*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + 20*a*x*pow(y, 2)*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - 12*a*x*y*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(40*x*pow(y, 2)*(20*t*x*pow(y, 2) + x - x0) - 24*x*y*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + 12*x*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - 20*pow(y, 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//    }
//    else { std::cout << "ST - This choice does not exist!" << std::endl; exit(0);}
// }

double Function::Source_term1(const double x, const double y, const double t) const
{
   return 0;
}

double Function::Source_term2(const double x, const double y, const double t) const
{
   return 0;
}

double Function::Source_term3(const double x, const double y, const double t) const
{
   return 0;
}

// double Function::Neumann_Function(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       double P(_df->Get_P_kettle()/0.2), R(0.015), S(2*M_PI*R), lambda(0.598);
//       return P/(S*lambda);
//    }
//    else if (this->_df->Get_scenario() == "diffusion_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_all_BC")
//    {
//       double x0(0.5), y0(0.5), a(5);
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double diffx = -a*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//       double diffy = -a*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//       double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
//       double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
//       return diffx*nx+diffy*ny;
//    }
//    else if (this->_df->Get_scenario() == "advection_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_advection_all_BC")
//    {
//       double x0(0.9), y0(0.1), a(50);
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double diffx = -a*(-12*t*y*(-6*t*x*y + y - y0) + (20*t*pow(y, 2) + 2)*(10*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
//       double diffy = -a*(40*t*x*y*(10*t*x*pow(y, 2) + x - x0) + (-12*t*x + 2)*(-6*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
//       double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
//       double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
//       return diffx*nx+diffy*ny;
//    }
//    else { std::cout << "NBC - This choice does not exist!" << std::endl; exit(0);}
// }

// double Function::Neumann_Function2(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       double P(_df->Get_P_kettle()/0.2), R(0.015), S(2*M_PI*R), lambda(0.598);
//       return P/(S*lambda);
//    }
//    else if (this->_df->Get_scenario() == "diffusion_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_all_BC")
//    {
//       double x0(0.5), y0(0.5), a(5);
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double diffx = -a*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//       double diffy = -a*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//       double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
//       double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
//       return diffx*nx+diffy*ny;
//    }
//    else if (this->_df->Get_scenario() == "advection_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_advection_all_BC")
//    {
//       double x0(0.9), y0(0.1), a(50);
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double diffx = -a*(-12*t*y*(-6*t*x*y + y - y0) + (20*t*pow(y, 2) + 2)*(10*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
//       double diffy = -a*(40*t*x*y*(10*t*x*pow(y, 2) + x - x0) + (-12*t*x + 2)*(-6*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
//       double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
//       double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
//       return diffx*nx+diffy*ny;
//    }
//    else { std::cout << "NBC - This choice does not exist!" << std::endl; exit(0);}
// }

// double Function::Neumann_Function3(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       double P(_df->Get_P_kettle()/0.2), R(0.015), S(2*M_PI*R), lambda(0.598);
//       return P/(S*lambda);
//    }
//    else if (this->_df->Get_scenario() == "diffusion_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_all_BC")
//    {
//       double x0(0.5), y0(0.5), a(5);
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double diffx = -a*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//       double diffy = -a*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
//       double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
//       double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
//       return diffx*nx+diffy*ny;
//    }
//    else if (this->_df->Get_scenario() == "advection_hom_neumann")
//    {
//       return 0;
//    }
//    else if (this->_df->Get_scenario() == "diffusion_advection_all_BC")
//    {
//       double x0(0.9), y0(0.1), a(50);
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double diffx = -a*(-12*t*y*(-6*t*x*y + y - y0) + (20*t*pow(y, 2) + 2)*(10*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
//       double diffy = -a*(40*t*x*y*(10*t*x*pow(y, 2) + x - x0) + (-12*t*x + 2)*(-6*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
//       double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
//       double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
//       return diffx*nx+diffy*ny;
//    }
//    else { std::cout << "NBC - This choice does not exist!" << std::endl; exit(0);}
// }

// double Function::Dirichlet_Function(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       return 15;
//    }
//    else
//    {
//       return Exact_solution(x,y,t);
//    }
// }

// double Function::Dirichlet_Function1(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "none")
//    {
//       // TODO for real case
//       return 15;
//    }
//    else
//    {
//       return Exact_solution(x,y,t);
//    }
// }

// double Function::Exact_solution(const double x, const double y, const double t) const
// {
//    if (this->_df->Get_scenario() == "diffusion_hom_neumann")
//    {
//       double xmin(0), ymin(0), xmax(1), ymax(1);
//       double x0(0.5), y0(0.5), a(5);
//       return 50.0*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)));
//    }
//    else if (this->_df->Get_scenario() == "diffusion_all_BC")
//    {
//       double x0(0.5), y0(0.5), a(5);
//       return 0.5*(sin(0.25*M_PI*t) + 1)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)));
//    }
//    else if (this->_df->Get_scenario() == "advection_hom_neumann")
//    {
//       return Initial_condition(x-Velocity_x(x,y,t)*t, y-Velocity_y(x,y,t)*t);
//    }
//    else if (this->_df->Get_scenario() == "diffusion_advection_all_BC")
//    {
//       double x0(0.9), y0(0.1), a(50);
//       double Vx = -20*y*y*x;
//       double Vy = 12*x*y;
//       return exp(-a*((x-Vx*t-x0)*(x-Vx*t-x0)+(y-Vy*t-y0)*(y-Vy*t-y0)))*exp(-this->_mu*t);
//    }
//    else { std::cout << "ES - This choice does not exist!" << std::endl; exit(0);}
// }


#define _FUNCTION_CPP
#endif
