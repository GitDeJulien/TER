#ifndef _TIME_SCHEME_H

#include "FiniteVolume.h"

class TimeScheme
{
protected:
   // Pointeur vers la classe FiniteVolume
   FiniteVolume* _fin_vol;
   // Pointeur vers la classe FiniteVolume
   DataFile* _df;
   // Vecteur initial et vecteur solution
   // Eigen::VectorXd _sol;
   // Eigen::VectorXd _h;
   // Eigen::VectorXd _hu;
   // Eigen::VectorXd _hv;
   Eigen::VectorXd _sol1;
   Eigen::VectorXd _sol2;
   Eigen::VectorXd _sol3;

   // Time
   double _t;

public:
   // Constructeur par défaut
   TimeScheme(DataFile* data_file, FiniteVolume* adv);
   // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
   // alors le compilateur en génère un implicitement.
   virtual ~TimeScheme();
   // Enregistre la solution un fichier
   void Save_solution(int n) {_fin_vol->Save_sol(_sol1, n, "solution");};
   // Enregistre la solution un fichier
   void Save_solution(Eigen::VectorXd s, int n, std::string st) {_fin_vol->Save_sol(s, n, st);};
   // Une étape du schéma en temps
   virtual void Advance() = 0;
   // Permet de récupérer _sol
   // const Eigen::VectorXd & Get_sol() const {return _sol;};
   // const Eigen::VectorXd & Get_sol1() const {return _h;};
   // const Eigen::VectorXd & Get_sol2() const {return _hu;};
   // const Eigen::VectorXd & Get_sol3() const {return _hv;};
   const Eigen::VectorXd & Get_sol1() const {return _sol1;};
   const Eigen::VectorXd & Get_sol2() const {return _sol2;};
   const Eigen::VectorXd & Get_sol3() const {return _sol3;};
};

class EulerScheme : public TimeScheme
{
public:
   // Constructeur
   EulerScheme(DataFile* data_file, FiniteVolume* lap);
   // Une étape du schéma en temps
   void Advance();
};
// VectorXd FiniteVolume::Exact_solution(const double t)
// {
// 	VectorXd exactsol(this->_msh->Get_triangles().size());

// 	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
// 	{
// 		exactsol(i) = this->_fct->Exact_solution(this->_msh->Get_triangles_center()(i,0),
// 		this->_msh->Get_triangles_center()(i,1), t);
// 	}
// 	return exactsol;
// }

class ImplicitEulerScheme : public TimeScheme
{
private:
   Eigen::SparseLU<Eigen::SparseMatrix<double> > _solver_direct;
public:
   // Constructeur
   ImplicitEulerScheme(DataFile* data_file, FiniteVolume* lap);
   // Une étape du schéma en temps
   void Advance();
};

#define _TIME_SCHEME_H
#endif
