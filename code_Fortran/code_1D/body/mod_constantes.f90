module mod_constantes

    implicit none

    integer, parameter :: pr = 8
    real(pr), parameter :: PI = 4.*ATAN(1._pr)
    real(pr), parameter :: g = 9.80665
    real(pr), parameter :: h_debut = 2.
    real(pr), parameter :: h_fin = 1.

end module