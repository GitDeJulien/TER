// Rafinement du maillage
h = 0.1;

// Toutes les valeurs sont en mètre

// Longueur du bassin
L = 5; 
//Largeur du bassin
l = 5e-1;



// Définition des points
Point(1) = {-L, l, 0.0, h};
Point(2) = {L, l, 0.0, h};
Point(3) = {L, -l, 0.0, h};
Point(4) = {-L, -l, 0.0, h};
Point(5) = {0.0, l, 0.0, h};
Point(6) = {0.0, -l, 0.0, h};

//Point obstacle
//Point(8) = {1, -2e-1, 0.0, h};
//Point(9) = {1, 2e-1, 0.0, h};
//Point(10) = {2, 2e-1, 0.0, h};
//Point(11) = {2, -2e-1, 0.0, h};

// Définition des segments délimitant le domaine
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {3, 2};
Line(4) = {4, 1};
Line(5) = {6, 5};
Line(6) = {3, 6};
Line(7) = {6, 4};

//Obstacle
//Line(8) = {8, 9};
//Line(9) = {9, 10};
//Line(10) = {11, 10};
//Line(11) = {11, 8};

// Définition du bord du domaine
Curve Loop(1) = {4, 1, -5, 7};
Curve Loop(2) = {5, 2, -3, 6};

//Condition de bord
Physical Curve(1) = {1, 2, 6, 7}; //Mur
Physical Curve(2) = {3}; //Entrer
Physical Curve(3) = {4}; //Sortie

//Condition de surface
Physical Surface(10) = {1, 2, 3};


// Bord de l'obstacle
//Curve Loop(3) = {8, 9, -10, 11};

// Définition de la surface du domaine
Plane Surface(1) = {1};
Plane Surface(2) = {2};
//Plane Surface(3) = {3};


//Physical Curve(3) = {5};

// Définition bord et surface obtacle
//Physicla Curve(4) = {8, 9, 10, 11};





