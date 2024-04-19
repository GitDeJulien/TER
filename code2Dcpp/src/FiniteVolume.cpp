#ifndef _FINITEVOLUME_CPP

#include "FiniteVolume.h"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh2D* mesh) :
_fct(function), _df(data_file), _msh(mesh)
{
	std::cout << "Build finite volume class." << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
void FiniteVolume::Build_flux_mat_and_rhs(const double& t, const Eigen::VectorXd& sol1, 
const Eigen::VectorXd& sol2, const Eigen::VectorXd& sol3)
{
	this->_F1.resize(this->_msh->Get_triangles().size());
	this->_F2.resize(this->_msh->Get_triangles().size());
	this->_F3.resize(this->_msh->Get_triangles().size());
	this->_F1.setZero();
	this->_F2.setZero();
	this->_F3.setZero();
	double d_b;
	_min_d_b=1000000000000000000;
	
	for (unsigned int i = 0; i < this->_msh->Get_edges().size(); i++)
	{
		// Mesh2D
		int t1 = _msh->Get_edges()[i].Get_T1(), t2 = _msh->Get_edges()[i].Get_T2();
		if(t1==-1){t1 = _msh->Get_edges()[i].Get_T2(), t2 = _msh->Get_edges()[i].Get_T1();}

		double xe = _msh->Get_edges_center()(i,0), ye = _msh->Get_edges_center()(i,1);
		double nx = _msh->Get_edges_normal()(i,0), ny = _msh->Get_edges_normal()(i,1);
		double le = _msh->Get_edges_length()[i];
		// nx=-nx; ny=-ny;
		double T1x = _msh->Get_triangles_center()(t1,0), T1y = _msh->Get_triangles_center()(t1,1);
		// ne normale sortante du triangle T1
		// le sens de ne doit etre le meme que celui du vecteur X1->Xe
		if ((xe-T1x)*nx+(ye-T1y)*ny<0) {nx=-nx; ny=-ny; cout<<"n->-n"<<endl;}

		if (t2==-1)
		{
			double A1 = _msh->Get_triangles_area()[t1];
			string BC = _msh->Get_edges()[i].Get_BC();
 
			if (BC=="Neumann")
			{
				// double h = sol1[t1];
				double h1 = sol1[t1];
				double hu1 = sol2[t1];
				double hv1 = sol3[t1];

				double Fi2 = le*_df->Get_g()*0.5*h1*h1*nx;
				double Fi3 = le*_df->Get_g()*0.5*h1*h1*ny;
				double be = max(abs((hu1*nx+hv1*ny)/h1+sqrt(_df->Get_g()*h1)),
								abs((hu1*nx+hv1*ny)/h1-sqrt(_df->Get_g()*h1)));

				d_b = A1/(_msh->Get_triangles_length()[t1]*be);
				if (d_b<_min_d_b) {_min_d_b=d_b;}
				
				// On ajoute les contributions du flux dans les mailles
				this->_F2[t1] += Fi2/A1;
				this->_F3[t1] += Fi3/A1;
			}
			else {cout<<"Nom CL <<"<<BC<<">> non valide"<<endl;}
		}
		else
		{
			// On récupère les paramétres des mailles
			double A1 = _msh->Get_triangles_area()[t1], A2 = _msh->Get_triangles_area()[t2];

			// Moyenne
			double  h1 = sol1[t1], h2 = sol1[t2];
			double hu1 = sol2[t1], hu2 = sol2[t2];
			double hv1 = sol3[t1], hv2 = sol3[t2];
			double huv1 = hu1*hv1/h1, huv2 = hu2*hv2/h2;

			double Fi1 = le*((hu1+hu2)*nx+(hv1+hv2)*ny)*0.5;
			double Fi2 = le*((hu1*hu1/h1+_df->Get_g()*0.5*h1*h1+
							  hu2*hu2/h2+_df->Get_g()*0.5*h2*h2)*nx+
							  (huv1+huv2)*ny)*0.5;
			double Fi3 = le*((huv1+huv2)*nx+
							(hv1*hv1/h1+_df->Get_g()*0.5*h1*h1+
							 hv2*hv2/h2+_df->Get_g()*0.5*h2*h2)*ny)*0.5;

			// Correction
			double be = max(max(abs((hu1*nx+hv1*ny)/h1+sqrt(_df->Get_g()*h1)),
								abs((hu1*nx+hv1*ny)/h1-sqrt(_df->Get_g()*h1))),
							max(abs((hu2*nx+hv2*ny)/h2+sqrt(_df->Get_g()*h2)),
								abs((hu2*nx+hv2*ny)/h2-sqrt(_df->Get_g()*h2))));

			d_b = A1/(_msh->Get_triangles_length()[t1]*be);
			if (d_b<_min_d_b) {_min_d_b=d_b;}
			d_b = A2/(_msh->Get_triangles_length()[t2]*be);
			if (d_b<_min_d_b) {_min_d_b=d_b;}

			// cout<<"be="<<be<<endl;

			Fi1+=-0.5*be*(h2-h1);
			Fi2+=-0.5*be*(hu2-hu1);
			Fi1+=-0.5*be*(hv2-hv1);
			
			// On ajoute les contributions du flux dans les mailles
			this->_F1[t1] += Fi1/A1; this->_F1[t2] += -Fi1/A2;
			this->_F2[t1] += Fi2/A1; this->_F2[t2] += -Fi2/A2;
			this->_F3[t1] += Fi3/A1; this->_F3[t2] += -Fi3/A2;
		}

	}
}


// --- Déjà implémenté ---
// Construit la condition initiale au centre des triangles
VectorXd FiniteVolume::Initial_condition1()
{
	VectorXd sol0(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sol0(i) = this->_fct->Initial_condition1(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1));
	}

	return sol0;
}

VectorXd FiniteVolume::Initial_condition2()
{
	VectorXd sol0(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sol0(i) = this->_fct->Initial_condition2(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1));
	}

	return sol0;
}

VectorXd FiniteVolume::Initial_condition3()
{
	VectorXd sol0(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sol0(i) = this->_fct->Initial_condition3(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1));
	}

	return sol0;
}

// Terme source au centre des triangles
VectorXd FiniteVolume::Source_term1(double t)
{
	VectorXd sourceterm(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = this->_fct->Source_term1(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}

	return sourceterm;
}

VectorXd FiniteVolume::Source_term2(double t)
{
	VectorXd sourceterm(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = this->_fct->Source_term2(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}

	return sourceterm;
}

VectorXd FiniteVolume::Source_term3(double t)
{
	VectorXd sourceterm(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = this->_fct->Source_term3(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}

	return sourceterm;
}

// Sauvegarde la solution
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol1, int n, std::string st)
{
	double norm1 = 0;
	for (unsigned int i = 0; i < sol1.rows(); i++)
	{
		norm1 += sol1(i)*sol1(i)*this->_msh->Get_triangles_area()[i];
	}
	norm1 = sqrt(norm1);
	this->_norm1 = norm1;

	if (st == "solution")
	{
		cout << "Norme de u1 = " << norm1 << endl;
	}

	string name_file = this->_df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
	unsigned int nb_vert = this->_msh->Get_vertices().size();
	assert(((long unsigned int)sol1.size() == this->_msh->Get_triangles().size())
	&& "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,sol1[i]) << endl;
	}
	solution << endl;

	//solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS CFL float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,this->_df->Get_dt()*fabs(sol1[i])/this->_msh->Get_triangles_length()(i)) << endl;
	}
	solution << endl;

	solution.close();
}

#define _FINITEVOLUME_CPP
#endif
