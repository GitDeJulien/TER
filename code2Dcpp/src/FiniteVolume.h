#ifndef _FINITEVOLUME_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "libraries/Data/Function.h"
#include "libraries/Mesh/Mesh2D.h"

class FiniteVolume {
private:
	// Pointeur de la classe Function (accès à la condition initiale,
	// la solution exacte et à la vitesse)
	const Function* _fct;
	// Pointeur de la classe DataFile pour récupérer toutes les
	// valeurs de paramètres
	const DataFile* _df;
	// Pointeur de la classe Mesh pour récupérer toutes les
	// données concernant le maillage
	const Mesh2D* _msh;

	// Vecteurs flux
	Eigen::VectorXd _F1;
	Eigen::VectorXd _F2;
	Eigen::VectorXd _F3;

	// Membre de droite pour les conditions aux bords
	// Eigen::VectorXd _BC_RHS;
	// Eigen::VectorXd _BC_RHS1;
	// Eigen::VectorXd _BC_RHS2;
	// Eigen::VectorXd _BC_RHS3;

	// Norme de u
	double _norm;
	double _norm1;
	double _norm2;
	double _norm3;

	// Minimum du rapport delta_k sur be
	double _min_d_b;

public:
	// Constructeur
	FiniteVolume(Function* fct, DataFile* data_file, Mesh2D* mesh);

	// Construit la matrice des flux et le membre de droite
	void Build_flux_mat_and_rhs(const double& t, const Eigen::VectorXd& _sol1, 
	const Eigen::VectorXd& _sol2, const Eigen::VectorXd& _sol3);

	// Renvoie la matrice qui permet d'obtenir le flux en faisant : M * sol
	// const Eigen::SparseMatrix<double>& Get_flux_matrix() const {return _mat_flux;};
	const Eigen::VectorXd& Get_F1() const {return _F1;};
	const Eigen::VectorXd& Get_F2() const {return _F2;};
	const Eigen::VectorXd& Get_F3() const {return _F3;};

	// Renvoie le vecteur M*sol = rhs
	// const Eigen::VectorXd& Get_BC_RHS() const {return _BC_RHS;};
	// const Eigen::VectorXd& Get_BC_RHS1() const {return _BC_RHS1;};
	// const Eigen::VectorXd& Get_BC_RHS2() const {return _BC_RHS2;};
	// const Eigen::VectorXd& Get_BC_RHS3() const {return _BC_RHS3;};

	// Renvoie la norme de u
	const double Get_norm() const {return _norm;};

	// Renvoie le min de d_b
	const double Get_min_d_b() const {return _min_d_b;};

	// Condition Initiale au centre des triangles
	Eigen::VectorXd Initial_condition1();
	Eigen::VectorXd Initial_condition2();
	Eigen::VectorXd Initial_condition3();

	// Terme source au centre des triangles
	Eigen::VectorXd Source_term1(double t);
	Eigen::VectorXd Source_term2(double t);
	Eigen::VectorXd Source_term3(double t);

	// Solution exacte au centre des triangles
	// Eigen::VectorXd Exact_solution(double t);

	// Sauvegarde la solution
	void Save_sol(const Eigen::VectorXd& sol1, int n, std::string st);
};


#define _FINITEVOLUME_H
#endif
