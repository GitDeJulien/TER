#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* adv) :
_fin_vol(adv),_df(data_file), _sol1(adv->Initial_condition1()), _sol2(adv->Initial_condition2()), 
_sol3(adv->Initial_condition3()), _t(_df->Get_t0())
{
}

EulerScheme::EulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
   std::cout << "Build time scheme class." << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;
}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{
}

// Euler Explicite
void EulerScheme::Advance()
{
   Eigen::VectorXd F1, F2, F3;
   _fin_vol->Build_flux_mat_and_rhs(_t, _sol1, _sol2, _sol3);
   // _df->Adapt_dt(_df->Get_cfl()*_fin_vol->Get_min_d_b());
   F1 = _fin_vol->Get_F1();
   F2 = _fin_vol->Get_F2();
   F3 = _fin_vol->Get_F3();
   _sol1+=_df->Get_dt()*(_fin_vol->Source_term1(_t)-F1);
   _sol2+=_df->Get_dt()*(_fin_vol->Source_term2(_t)-F2);
   _sol3+=_df->Get_dt()*(_fin_vol->Source_term3(_t)-F3);
}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{
   // TODO
}

#define _TIME_SCHEME_CPP
#endif
