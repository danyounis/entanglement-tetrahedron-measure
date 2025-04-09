! tetra_mini.f08
!
! Author: D. Younis
!         University of Rochester
!         Department of Physics
!
! Written: 12/26/2022
! Revised: 2/7/2023

program tetra_mini
    use prec
    use math
    use optimize
    use omp_lib

    implicit none

    integer, parameter :: xdim=3, adim=18
    complex(num), parameter :: Sigma_y(4,4) = &
        reshape(real((/0,0,0,-1,0,0,1,0,0,1,0,0,-1,0,0,0/), num), (/4,4/))
    real(num) :: feval, da, ftol
    character(len=24) :: deck_name
    character(len=20) :: Xeval_name, feval_name
    character(len=10) :: uid
    integer :: iu, nmins, itmax
    real(num) :: V(xdim) !={x,y,z} global

    if (num/=8) then
        print '(A)', 'ERROR: This code must be compiled using double-precision floats.'
        print '(A)', '[parameter num = kind(1.d0) = 8]'
        call exit(1)
    end if

    ! utilize all available threads
    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)
    ! parse input deck
    call get_command_argument(1,deck_name)
    call get_command_argument(2,Xeval_name)
    call read_deck()
    ! initialize (pseudo) Random Number Generator
    call init_RNG()

    ! Xeval_name = 'Xeval-xxxxxxxxxx.dat' for x in {0..9}
    uid = Xeval_name(7:16); feval_name = 'feval-'//uid//'.dat';

    ! cache {x,y,z} evaluation point
    open(newunit=iu, file=Xeval_name, form='unformatted')
    read(iu) V
    close(iu)

    ! calculate sigma_12
    feval = ObjectiveE1()

    ! output value
    open(newunit=iu, file=feval_name, form='unformatted')
    write(iu) feval
    close(iu)

    ! terminate
    call exit()

contains

! __________________________________________________________________________________________________
!
! OBJECTIVE FUNCTIONS
! __________________________________________________________________________________________________
!

! Interior (infimum) objective function (1).
! Input angles a(1:11)={11θ}, a(12:18)={7ϕ}.
real(num) function ObjectiveI1(a)
    implicit none
    real(num), intent(in) :: a(:)
    real(num) :: Cs(10), T(4), Denom, Delta
    Delta = 0.0_num

    ! calculate concurrences
    Cs = Concurrence(psi(a(1:11),a(12:18)))

    ! calculate tau's
    T(1) = Cs(1) - (Cs(5) + Cs(6) + Cs(7))
    T(2) = Cs(2) - (Cs(5) + Cs(8) + Cs(9))
    T(3) = Cs(3) - (Cs(6) + Cs(8) + Cs(10))
    T(4) = Cs(4) - (Cs(7) + Cs(9) + Cs(10))

    ! denominator term defined in notes, Eq. (13).
    Denom = V(3)*sum(T) + sum(Cs(5:10))
    if (Denom /= 0.0_num) Delta = -Numer(1,2,3,4,V,Cs)/Denom

    ObjectiveI1 = Cs(5) + Delta
    return
end function ObjectiveI1

! Exterior (supremum) objective function, Max of Min[I1].
! Must set array of parameters V(1:3)={x,y,z}.
real(num) function ObjectiveE1()
    implicit none
    integer :: n, im
    type(OptimizeND_NelderMead) :: H
    real(num), allocatable, save :: My(:), Mp(:,:)

    ! interior-step thread master variables
    if (.not.allocated(My)) &
        allocate(My(nmins), Mp(nmins,adim), source=0.0_num)

    ! execute interior minimization task 'nmins' times in-parallel
    !$OMP PARALLEL DO SHARED(My,Mp) PRIVATE(H,im)
    do n=1,nmins
        ! initialize & execute
        call H%create(func=ObjectiveI1, ndim=adim, ftol=ftol, itmax=itmax, warn=.false.)
        call init_opt(H)
        call H%minimize()
        ! save the best simplex value/vertex in the thread master array
        im = minloc(H%y,1)
        My(n) = H%y(im)
        Mp(n,:) = H%p(im,:)
        call H%destroy()
    end do
    !$OMP END PARALLEL DO

    ObjectiveE1 = minval(My)

    return
end function ObjectiveE1

! __________________________________________________________________________________________________
!
! MAIN CALCULATION FUNCTIONS
! __________________________________________________________________________________________________
!

! Calculate squared one-to-other & one-to-one concurrences given a 4-qubit state (s).
function Concurrence(s) result(Cs)
    implicit none
    real(num) :: Cs(10)
    complex(num), intent(in) :: s(0:1,0:1,0:1,0:1)

    real(num) :: Ca(4), Cb(6), WR(4)
    complex(num) :: rho_a(0:1,0:1), rho_b(4,4), tmp(0:1,0:1,0:1,0:1)
    integer :: nn,n,j,k,l

    ! temporary variables, LAPACK eigenvalue subroutine
    external :: ZGEES, SELECT
    logical :: BWORK(4)
    integer :: INFO, SDIM
    double precision :: RWORK(4)
    complex*16 :: A(4,4), W(4), VS(1), WORK(12)

    ! One-to-other, Ca:=C_i^2 for i={1..4}.
    do nn=1,4
    select case (nn)
    case (1)
        ForAll(j=0:1,k=0:1) & ! C_1^2
            rho_a(j,k) = sum(s(j,:,:,:)*conjg(s(k,:,:,:)))
    case (2)
        ForAll(j=0:1,k=0:1) & ! C_2^2
            rho_a(j,k) = sum(s(:,j,:,:)*conjg(s(:,k,:,:)))
    case (3)
        ForAll(j=0:1,k=0:1) & ! C_3^2
            rho_a(j,k) = sum(s(:,:,j,:)*conjg(s(:,:,k,:)))
    case (4)
        ForAll(j=0:1,k=0:1) & ! C_4^2
            rho_a(j,k) = sum(s(:,:,:,j)*conjg(s(:,:,:,k)))
    end select
        Ca(nn) = real(4.*(rho_a(0,0)*rho_a(1,1) - rho_a(0,1)*rho_a(1,0)), num)
    end do

    ! One-to-one, Cb:=C_{ij}^2 for {i,j}={1..4}, i/=j.
    do nn=1,6
    select case (nn)
    case (1)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! C_{12}^2
            tmp(n,j,k,l) = sum(s(n,j,:,:)*conjg(s(k,l,:,:)))
    case (2)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! C_{13}^2
            tmp(n,j,k,l) = sum(s(n,:,j,:)*conjg(s(k,:,l,:)))
    case (3)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! C_{14}^2
            tmp(n,j,k,l) = sum(s(n,:,:,j)*conjg(s(k,:,:,l)))
    case (4)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! C_{23}^2
            tmp(n,j,k,l) = sum(s(:,n,j,:)*conjg(s(:,k,l,:)))
    case (5)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! C_{24}^2
            tmp(n,j,k,l) = sum(s(:,n,:,j)*conjg(s(:,k,:,l)))
    case (6)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! C_{34}^2
            tmp(n,j,k,l) = sum(s(:,:,n,j)*conjg(s(:,:,k,l)))
    end select
        rho_b = to4x4(tmp)
        A = matmul(rho_b, matmul(Sigma_y, matmul(conjg(rho_b), Sigma_y)))
        CALL ZGEES('N','N',UNUSED,4,A,4,SDIM,W,VS,1,WORK,12,RWORK,BWORK,INFO)
        ! if (INFO/=0) call exit(INFO)
        WR = sort4(real(W,num))
        where (WR<0.0_num) WR=0.0_num
        Cb(nn) = max(0.0_num, sqrt(WR(4))-sqrt(WR(3))-sqrt(WR(2))-sqrt(WR(1)))**2
    end do

    Cs = (/Ca,Cb/)
    return
end function Concurrence

! 4-qubit state parametrized by 11 (t)heta and 7 (p)hi angles.
function psi(t,p) result(s)
    implicit none
    complex(num) :: s(0:1,0:1,0:1,0:1)
    real(num), intent(in) :: t(11), p(7)
    s = (0.0_num, 0.0_num)
    s(0,0,0,0) = cos(t(1))
    s(0,1,0,0) = s(0,0,0,0)*tan(t(1))*cos(t(2))
    s(0,1,0,1) = s(0,1,0,0)*tan(t(2))*cos(t(3))
    s(0,1,1,0) = s(0,1,0,1)*tan(t(3))*cos(t(4))
    s(1,0,0,0) = s(0,1,1,0)*tan(t(4))*cos(t(5))
    s(1,0,0,1) = s(1,0,0,0)*tan(t(5))*cos(t(6))*exp(i*p(1))
    s(1,0,1,0) = s(1,0,0,1)*tan(t(6))*cos(t(7))*exp(i*(p(2)-p(1)))
    s(1,0,1,1) = s(1,0,1,0)*tan(t(7))*cos(t(8))*exp(i*(p(3)-p(2)))
    s(1,1,0,0) = s(1,0,1,1)*tan(t(8))*cos(t(9))*exp(i*(p(4)-p(3)))
    s(1,1,0,1) = s(1,1,0,0)*tan(t(9))*cos(t(10))*exp(i*(p(5)-p(4)))
    s(1,1,1,0) = s(1,1,0,1)*tan(t(10))*cos(t(11))*exp(i*(p(6)-p(5)))
    s(1,1,1,1) = s(1,1,1,0)*tan(t(11))*exp(i*(p(7)-p(6)))
    return
end function psi

! Numerator term defined in notes, Eq. (14).
! Input integers (i,j,k,l), array X={x,y,z}, and concurrence array Cs.
function Numer(i,j,k,l,X,Cs) result(r)
    implicit none
    real(num) :: r
    integer, intent(in) :: i,j,k,l
    real(num), intent(in) :: X(3), Cs(10)

    real(num) :: Ti, Tj, Tk, Tl
    real(num) :: Cij, Cik, Cil, Cjk, Cjl, Ckl

    Cij = getConc(Cs,i,j)
    Cik = getConc(Cs,i,k)
    Cil = getConc(Cs,i,l)
    Cjk = getConc(Cs,j,k)
    Cjl = getConc(Cs,j,l)
    Ckl = getConc(Cs,k,l)

    Ti = Cs(i) - (Cij + Cik + Cil)
    Tj = Cs(j) - (Cij + Cjk + Cjl)
    Tk = Cs(k) - (Cik + Cjk + Ckl)
    Tl = Cs(l) - (Cil + Cjl + Ckl)

    r = Cij*Ckl - 0.5*Cik*Cjl - 0.5*Cil*Cjk &
      - (X(3)/3.0_num)*(Ti**2 + Tj**2) + (X(3)/6.0_num)*(Tk**2 + Tl**2) &
      + X(1)*X(3)*Ti*Tj + X(3)*(X(1)+1.0_num)*Tk*Tl - 0.5*X(3)*(X(1)+1.0_num)*(Ti+Tj)*(Tk+Tl) &
      + X(2)*Cij*(Ti+Tj) + (X(2)+0.5)*Ckl*(Tk+Tl) &
      - 0.5*(X(2)+1.0_num)*(Ti*(Cik+Cil) + Tj*(Cjk+Cjl)) &
      - 0.5*X(2)*(Tk*(Cik+Cjk) + Tl*(Cil+Cjl)) &
      + X(2)*(Ti*(Cjk+Cjl) + Tj*(Cik+Cil)) &
      + (X(2)+0.5)*(Tk*(Cil+Cjl) + Tl*(Cik+Cjk)) &
      - (2*X(2)+1.0_num)*Ckl*(Ti+Tj) - (2*X(2)+0.5)*Cij*(Tk+Tl)

    return
end function Numer

! Fetch the appropriate one-to-one concurrence value from the Cs array.
function getConc(Cs,a,b) result(r)
    implicit none
    real(num) :: r
    real(num), intent(in) :: Cs(10)
    integer, intent(in) :: a,b
    character(len=2) :: ab
    write(ab,'(I1,I1)') a,b
    select case (ab)
    case ('12','21')
        r = Cs(5)
    case ('13','31')
        r = Cs(6)
    case ('14','41')
        r = Cs(7)
    case ('23','32')
        r = Cs(8)
    case ('24','42')
        r = Cs(9)
    case ('34','43')
        r = Cs(10)
    end select
    return
end function getConc

! Order 4-qubit state components in 4-by-4 matrix form.
!  |0000>, |0001>, |0010>, |0011>
!  |0100>, |0101>, |0110>, |0111>
!  |1000>, |1001>, |1010>, |1011>
!  |1100>, |1101>, |1110>, |1111>
function to4x4(A) result(B)
    implicit none
    complex(num) :: B(4,4)
    complex(num), intent(in) :: A(0:1,0:1,0:1,0:1)
    B(:,1) = (/A(0,0,0,0), A(0,1,0,0), A(1,0,0,0), A(1,1,0,0)/)
    B(:,2) = (/A(0,0,0,1), A(0,1,0,1), A(1,0,0,1), A(1,1,0,1)/)
    B(:,3) = (/A(0,0,1,0), A(0,1,1,0), A(1,0,1,0), A(1,1,1,0)/)
    B(:,4) = (/A(0,0,1,1), A(0,1,1,1), A(1,0,1,1), A(1,1,1,1)/)
    return
end function to4x4

! Sort a 4D array in ascending order,
! x->y such that y(4)>y(3)>y(2)>y(1).
function sort4(x) result(y)
    implicit none
    real(num) :: y(4)
    real(num), intent(in) :: x(4)
    real(num) :: l(2), m(2), h(2)
    if (x(1)<x(2)) then
        l(1)=x(1); h(1)=x(2);
    else
        l(1)=x(2); h(1)=x(1);
    end if
    if (x(3)<x(4)) then
        l(2)=x(3); h(2)=x(4);
    else
        l(2)=x(4); h(2)=x(3);
    end if
    if (l(1)<l(2)) then
        y(1)=l(1); m(1)=l(2);
    else
        y(1)=l(2); m(1)=l(1);
    end if
    if (h(1)>h(2)) then
        y(4)=h(1); m(2)=h(2);
    else
        y(4)=h(2); m(2)=h(1);
    end if
    if (m(1)<m(2)) then
        y(2:3) = (/m(1),m(2)/)
    else
        y(2:3) = (/m(2),m(1)/)
    end if
    return
end function sort4

! Dummy function for LAPACK's eigenvalue subroutine.
LOGICAL FUNCTION UNUSED(X,Y)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X,Y
    UNUSED = .TRUE.
    RETURN
END FUNCTION UNUSED

! __________________________________________________________________________________________________
!
! SUPPORTING ROUTINES
! __________________________________________________________________________________________________
!

! Initialize an interior minimization task by generating a random simplex.
subroutine init_opt(this)
    implicit none
    class(OptimizeND_NelderMead), intent(inout) :: this
    real(num) :: rand
    integer :: n, sgn

    this%iter = 0

    ! 1st vertex, random guess in [0,pi)
    call random_number(this%p(1,:))
    this%p(1,:) = pi*this%p(1,:)
    this%y(1) = this%func(this%p(1,:))

    ! other vertices obtained by taking
    ! unit da steps in each direction
    do n=2,this%ndim+1
        this%p(n,:) = this%p(1,:)
        sgn = 1
        call random_number(rand)
        if (rand < 0.5) sgn = -1
        this%p(n,n-1) = this%p(n,n-1) + sgn*da
        this%y(n) = this%func(this%p(n,:))
    end do

    return
end subroutine init_opt

! __________________________________________________________________________________________________
!
! DATA I/O ROUTINES
! __________________________________________________________________________________________________
!

! Read input deck parameters.
subroutine read_deck
    implicit none
    character(len=128) :: label
    real(num) :: itmax_t
    integer :: iu
    open(newunit=iu, file=trim(adjustl(deck_name)), status='old', action='read')
    read(iu,*)!--------------------------------------------------
    read(iu,*)!PARAMETERS FOR INTERIOR MINIMIZATION (OVER 11θ/7φ)
    read(iu,*)!--------------------------------------------------
    read(iu,*) label, nmins
    read(iu,*) label, itmax_t; itmax = int(itmax_t);
    read(iu,*) label, ftol
    read(iu,*) label, da; da = pi*da;
    close(iu)
    return
end subroutine read_deck

end program tetra_mini
