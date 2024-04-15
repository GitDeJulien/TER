// Rafinement du maillage
h = 0.01;

// Toutes les valeurs sont en mètre

// Longueur du bassin
L = 5;
//Largeur du bassin
l = 3.7e-2;


// Définition des points
Point(1) = {-L, l, 0.0, h};
Point(2) = {L, l, 0.0, h};
Point(3) = {L, -l, 0.0, h};
Point(4) = {-L, -l, 0.0, h};
Point(5) = {0.0, l, 0.0, h};
Point(6) = {0.0, -l, 0.0, h};

//Point obstacle
Point(8) = {1-0.0125, -l, 0.0, h};
Point(9) = {1-0.0125, l, 0.0, h};
Point(10) = {1+0.0125, l, 0.0, h};
Point(11) = {1+0.0125, -l, 0.0, h};

// Définition des segments délimitant le domaine
Line(1) = {4, 1};
Line(2) = {1, 5};
Line(3) = {5, 6};
Line(4) = {6, 4};
Line(5) = {5, 9};
Line(6) = {9, 8};
Line(7) = {8, 6};
Line(8) = {9, 10};
Line(9) = {10, 2};
Line(10) = {3, 11};
Line(11) = {10, 11};
Line(12) = {11, 8};
Line(13) = {2, 3};


// Définition du bord du domaine
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-3, 5, 6, 7};
Curve Loop(3) = {-6, 8, 11, 12};
Curve Loop(4) = {-11, 9, 13, 10};


//Condition de bord
Physical Curve(1) = {2, 5, 8, 9, 10, 12, 7, 4};
Physical Curve(2) = {1, 13};

//Condition de surface
Physical Surface(10) = {1, 2};

// Définition de la surface du domaine
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
