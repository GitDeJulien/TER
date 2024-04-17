module mod_analyse_exp

    use mod_constantes
    use mod_initialize

    implicit none


contains 

subroutine evolution_capteur(file_in, file_out, imax,iter,dt, date, num)
        
    !-----Variables d'entrée-----!
    integer, intent(in)             :: imax, iter
    character(len=40),intent(in)    :: file_in, file_out, date, num
    real(PR),intent(in)             :: dt

    !-----Variables de sortie-----!
    character(len=40)            :: name_file_out

    !-----Variables locales-----!
    integer,dimension(:),allocatable  :: Pos_maillage
    real(pr), dimension(:), allocatable :: Pos_capteurs
    real,dimension(:),allocatable     :: H
    real(pr)                          :: a, xim1, t!, x
    integer                           :: k, km1, Nmax, j, l, iostatus, i, nb_capteurs, ligne

    open(unit = 2, file = file_in, action = "read")
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*) 
        read(2,*) 
        read(2,*)
        read(2,*) 
        read(2,*) nb_capteurs
    
        allocate(Pos_capteurs(0:nb_capteurs-1))

        print*, " nb capteurs = ",nb_capteurs

        ! print*, "A la position :"
        if ( nb_capteurs /= 0 ) then
            do i = 0, nb_capteurs-1
                read(2,*) Pos_capteurs(i) 
                ! print*, "x = ",Pos_capteurs(i) 
            end do
        end if
    close(2)


    allocate(Pos_maillage(0:nb_capteurs-1),H(0:nb_capteurs-1))

    open(10,file=file_out,action="read")

    Nmax=10000
    a=-40000             !Par sécurité
    xim1=0
    k=0

    !---------- Recherche du point du maillage le plus proche pour chaque capteurs ------!

    do i = 0,nb_capteurs-1
    !--- On considère que Pos_capteurs est classé de manière croissante ---!
        print*, "position du capteur x = ", Pos_capteurs(i)
        do while (Pos_capteurs(i)>a ) 
            xim1=a                          ! On stocke le point précédent !
            read(10,*) a                    ! On lit le point suivant !
            k = k + 1                       ! k compte la position du  point du maillage !
        end do
        Pos_maillage(i) = k
        ! print*, "Pour le capteur numéro ", i, ", le point de maillagee le plus proche est le numéro ", k
        !print*, "en x = ", a
        
    end do

    close(10)

    ! ---- Choix arbitraire de prendre le point juste d'après ----!

    ! ---- Ecriture de l'évolution des différents capteurs en fonction  

    name_file_out = "out/Hauteurs_num_"//trim(date)//"_exp_"//trim(num)//".txt"
    open(11,file=name_file_out,action="write")

    open(10,file=file_out,action="read")
    print*, file_out
    t = 0._PR
    ligne=0
    do l = 1,iter
        km1 = 0
        do i = 0, nb_capteurs -1
            k = Pos_maillage(i)

            do j = 1, k-1-km1
                read(10,*,iostat=iostatus)
                ! print*, j
            end do
            km1 = k
            ligne = ligne + j
            read(10,*,iostat=iostatus) a, H(i)
            print*, ligne, a, H(i)
            
        end do

        ! print*, "iter=", l, "h", H(0)

        ! write(11, '(F0.3,F12.6,F12.6,F12.6)') t, (H(j), j=0, nb_capteurs-1)
        write(11, '(F0.3,F12.6,F12.6,F12.6)') t, H(0), H(1), H(2)

        do j = k+1, imax + 4
            read(10,*,iostat=iostatus)
            if ( iostatus/=0 ) then
                print*, "Fin du fichier"
                exit
            end if
        end do
        ligne = ligne + imax+3-k+1
        t=t+dt
        
    end do
    
    close(11)
    close(10)

    print*, "Le nom du fichier des capteurs numériques est ",name_file_out
end subroutine evolution_capteur

subroutine convertion_tensions(date, num)
    implicit none
    !-----Variables d'entrée-----!
    character(len=40),intent(in)    :: date, num

    !-----Variables locales-----!
    real(pr), dimension(:), allocatable :: alpha, beta
    real(pr)                            :: t, h1, h2, h3, tmesure
    character(len=50)                   :: etalonnage, tensions, name_file_out
    integer                            :: iostatus, nb_capteurs, i

    etalonnage = "entry/exp_"//trim(date)//"/param/etalonnage.dat"
    tensions="entry/exp_"//trim(date)//"/mesures_capteurs/TER_exp"//trim(num)//".txt"
    name_file_out = "out/Hauteurs_EXP_"//trim(date)//"_exp_"//trim(num)//".txt"
    open(unit = 15, file = etalonnage , action = "read")
    open(unit = 16, file = tensions, action = "read", iostat=iostatus)
    open(unit = 17, file = name_file_out, action = "write") 
    
    nb_capteurs=3
    allocate(alpha(0:nb_capteurs-1), beta(0:nb_capteurs-1))
    do i =0, nb_capteurs-1
        read(15, *) alpha(i), beta(i)
    end do
    t=0
    do while (iostatus==0)
        read(16, *, iostat=iostatus) tmesure, h1, h2, h3
        write(17, *, iostat=iostatus) t, (h1-beta(0))/alpha(0), (h2-beta(1))/alpha(1), (h3-beta(2))/alpha(2)
        t=t+0.0002
    end do

    close(15)
    close(16)
    close(17)
    print*, "Le nom du fichier des capteurs expérimentaux est ",name_file_out
end subroutine convertion_tensions

end module 