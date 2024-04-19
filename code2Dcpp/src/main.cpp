#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

   if (argc < 2)
   {
      cout << "Please, enter the name of your data file." << endl;
      exit(0);
   }
   const string data_file_name = argv[1];

   // ----------------------- Fichier de données --------------------------------
   DataFile* data_file = new DataFile(data_file_name);
   // ---------------------------------------------------------------------------

   // ------------------Définition du nombre d'itérations------------------------
   int nb_iterations = int(ceil((data_file->Get_tfinal()-data_file->Get_t0())
   /data_file->Get_dt()));
   data_file->Adapt_dt(data_file->Get_tfinal() / nb_iterations);
   // ---------------------------------------------------------------------------

   // ---------------------------- Résolution  ----------------------------------
   Mesh2D* mesh = new Mesh2D(data_file->Get_BC_ref(),data_file->Get_BC_type());
   mesh->Read_mesh(data_file->Get_mesh_name());
   Function* fct = new Function(data_file);
   FiniteVolume* fin_vol = new FiniteVolume(fct, data_file, mesh);
   TimeScheme* time_scheme = NULL;

   if (data_file->Get_scheme() == "ExplicitEuler")
   time_scheme = new EulerScheme(data_file, fin_vol);
   else
   time_scheme = new ImplicitEulerScheme(data_file, fin_vol);

   cout << "-------------------------------------------------" << endl;
   cout << "Search h, hu, hv such that : " << endl;
   // cout << "dt u + div (v u) - div(mu grad(u)) = f" << endl;
   cout << "dt h + div F1 = 0" << endl;
   cout << "dt hu + div F2 = 0" << endl;
   cout << "dt hv + div F3 = 0" << endl;
   cout << "-------------------------------------------------" << endl;

   // Démarrage du chrono
   auto start = chrono::high_resolution_clock::now();


   // Boucle temporelle classique
   
   double t(data_file->Get_t0());
   ofstream results;
   const string results_name = data_file->Get_results()+"/results.txt";
   results.open(results_name, ios::out);

   VectorXd sol1 = time_scheme->Get_sol1();
   // Calcul de la temperature moyenne et température max
   int I = size(sol1);
   double S = 0;
   double A = 0;
   double Tmax = 0;
   for (int i=0; i<I; i++)
   {
      S+=sol1[i]*mesh->Get_triangles_area()[i];
      A+=mesh->Get_triangles_area()[i];
      if (sol1[i]>Tmax)
      {
         Tmax = sol1[i];
      }
   }
   cout << "t = " << t << endl;
   cout <<"Temperature moyenne : "<<S/A<<endl;
   results<<t<<" "<<S/A<<endl;
   cout<<"Temperature maximale : "<<Tmax<<endl;
   cout << "Save initial condition " << endl;
   time_scheme->Save_solution(0);
   cout << "Time Loop" << endl;
   for (int n = 1; n <= nb_iterations and Tmax<100; n++) // Boucle en temps
   {
      time_scheme->Advance();
      t+=data_file->Get_dt();
      
      time_scheme->Save_solution(n);
      sol1 = time_scheme->Get_sol1();
      // Calcul de la temperature moyenne
      S = 0;
      A = 0;
      Tmax = 0;
      
      for (int i=0; i<I; i++)
      {
         S+=sol1[i]*mesh->Get_triangles_area()[i];
         A+=mesh->Get_triangles_area()[i];
         if (sol1[i]>Tmax)
         {
            Tmax = sol1[i];
         }
         
      }
      cout << "t = " << t << endl;
      cout<<"Temperature moyenne : "<<S/A<<endl;
      cout<<"Temperature maximale : "<<Tmax<<endl;
      results<<t<<" "<<S/A<<endl;
   }
   
   

   // Fin du chrono
   auto finish = chrono::high_resolution_clock::now();
   double t_it = chrono::duration_cast<chrono::seconds>(finish-start).count();
   // Affichage du résultat
   cout << "Cela a pris "<< t_it << " seconds" << endl;

   delete time_scheme;
   delete fin_vol;
   delete data_file;
   delete fct;

   return 0;
}
