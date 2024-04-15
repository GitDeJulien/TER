// Rafinement du maillage
h = 0.01;

// Toutes les valeurs sont en mètre

// Longueur du bassin
L = 5;
//Largeur du bassin
l = 3.7e-2;

// Définition des points
Point(1) = {-L, l, 0.0, h};
Point(2) = {1-0.1, l, 0.0, h};
Point(3) = {1+0.132, l-0.023, 0.0, h};
Point(4) = {1+0.2, l, 0.0, h};
Point(5) = {L, l, 0.0, h};
Point(6) = {L, -l, 0.0, h};
Point(9) = {1-0.1, -l, 0.0, h};
Point(8) = {1+0.132, -l+0.023, 0.0, h};
Point(7) = {1+0.2, -l, 0.0, h};
Point(10) = {-L, -l, 0.0, h};
Point(11) = {0.0, l, 0.0, h};
Point(12) = {0.0, -l, 0.0, h};

// Définition des segments délimitant le domaine
Line(1) = {1, 11};
Line(2) = {11, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 12};
Line(11) = {12, 10};
Line(12) = {10, 1};
Line(13) = {12, 11};

// Définition du bord du domaine
Curve Loop(1) = {1, -13, 11, 12};
Curve Loop(2) = {2, 3, 4, 5, 6, 7, 8, 9, 10,13};

//Condition de bord
Physical Curve(1) = {1, 2, 3, 4, 5, 7, 8, 9, 10, 11};
Physical Curve(2) = {6, 12};

//Condition de surface
Physical Surface(10) = {1, 2};

// Définition de la surface du domaine
Plane Surface(1) = {1};
Plane Surface(2) = {2};
