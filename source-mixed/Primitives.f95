! 4-qubit state parametrized by 11 (t)heta and 7 (p)hi angles.
subroutine psi(t,p,s)
    implicit none
    complex, parameter :: i = (0,1)
    integer, parameter :: num = selected_real_kind(15,307)
    real(num), intent(in) :: t(11), p(7)
    complex(num), intent(out) :: s(0:1,0:1,0:1,0:1)
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
end subroutine psi

! Calculate squared one-to-other & one-to-one concurrences given a 4-qubit state (s).
subroutine Concurrence(s,Cs)
    implicit none
    integer, parameter :: num = selected_real_kind(15,307)
    complex, parameter :: i = (0,1)
    complex(num), parameter :: Sigma_y(4,4) = &
        reshape(real((/0,0,0,-1,0,0,1,0,0,1,0,0,-1,0,0,0/), num), (/4,4/))

    real(num), intent(out) :: Cs(10)
    complex(num), intent(in) :: s(0:1,0:1,0:1,0:1)

    real(num) :: Ca(4), Cb(6), WR(4)
    complex(num) :: rho_a(0:1,0:1), rho_b(4,4), tmp(0:1,0:1,0:1,0:1)
    integer :: nn,n,j,k,l

    ! temporary variables, LAPACK eigenvalue subroutine
    external :: ZGEES, SELECT
    logical :: BWORK(4), UNUSED
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
        call to4x4(tmp,rho_b)
        A = matmul(rho_b, matmul(Sigma_y, matmul(conjg(rho_b), Sigma_y)))
        CALL ZGEES('N','N',UNUSED,4,A,4,SDIM,W,VS,1,WORK,12,RWORK,BWORK,INFO)
        if (INFO/=0) then
            print '(A)', 'ERROR: In subroutine ZGEES.'
            print '(A,I3)', 'CODE ', INFO
            call exit(1)
        end if
        call sort4(real(W,num),WR)
        where (WR<0.0_num) WR=0.0_num
        Cb(nn) = max(0.0_num, sqrt(WR(4))-sqrt(WR(3))-sqrt(WR(2))-sqrt(WR(1)))**2
    end do

    Cs = (/Ca,Cb/)
    return
end subroutine Concurrence

! Calculate squared one-to-other & two-to-two concurrences given a 4-qubit state (s).
subroutine Concurrence_set2(s,Cs)
    implicit none
    integer, parameter :: num = selected_real_kind(15,307)
    complex, parameter :: i = (0,1)
    complex(num), parameter :: Sigma_y(4,4) = &
        reshape(real((/0,0,0,-1,0,0,1,0,0,1,0,0,-1,0,0,0/), num), (/4,4/))

    real(num), intent(out) :: Cs(7)
    complex(num), intent(in) :: s(0:1,0:1,0:1,0:1)

    real(num) :: Ca(4), Cb(3)
    complex(num) :: rho_a(0:1,0:1), rho_b(4,4), rho_b2(4,4), tmp(0:1,0:1,0:1,0:1)
    integer :: nn,n,j,k,l

    ! One-to-other, Ca:=C_{i(jkl)}^2 for i={1..4} and {j,k,l}/=i.
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

    ! Two-to-two, Cb:=C_{(ij)(kl)}^2 for {i,j}={1..4}, i/=j and {k,l}/={i,j}.
    do nn=1,3
    select case (nn)
    case (1)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! rho_12
            tmp(n,j,k,l) = sum(s(n,j,:,:)*conjg(s(k,l,:,:)))
    case (2)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! rho_13
            tmp(n,j,k,l) = sum(s(n,:,j,:)*conjg(s(k,:,l,:)))
    case (3)
        ForAll(n=0:1,j=0:1,k=0:1,l=0:1) & ! rho_14
            tmp(n,j,k,l) = sum(s(n,:,:,j)*conjg(s(k,:,:,l)))
    end select
        call to4x4(tmp,rho_b)
        rho_b2 = matmul(rho_b,rho_b)
        Cb(nn) = real(2.*(1.-rho_b2(1,1)-rho_b2(2,2)-rho_b2(3,3)-rho_b2(4,4)), num)
    end do

    Cs = (/Ca,Cb/)
    return
end subroutine Concurrence_set2

! Numerator term defined in notes, Eq. (14).
! Input integers (i,j,k,l), array X={x,y,z}, and concurrence array Cs.
subroutine Numer(i,j,k,l,X,Cs,r)
    implicit none
    integer, parameter :: num = selected_real_kind(15,307)
    real(num), intent(out) :: r
    integer, intent(in) :: i,j,k,l
    real(num), intent(in) :: X(3), Cs(10)

    real(num) :: Ti, Tj, Tk, Tl
    real(num) :: Cij, Cik, Cil, Cjk, Cjl, Ckl

    call getConc(Cs,i,j,Cij)
    call getConc(Cs,i,k,Cik)
    call getConc(Cs,i,l,Cil)
    call getConc(Cs,j,k,Cjk)
    call getConc(Cs,j,l,Cjl)
    call getConc(Cs,k,l,Ckl)

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
end subroutine Numer

! Interior (infimum) objective function (1).
! Input angles a(1:11)={11θ}, a(12:18)={7ϕ}, and array X={x,y,z}.
subroutine ObjectiveI1(a,X,r)
    implicit none
    real(8), intent(out) :: r
    real(8), intent(in) :: a(:), X(:)
    real(8) :: Cs(10), T(4), Denom, Delta, tmp
    complex(8) :: s(0:1,0:1,0:1,0:1)

    ! calculate concurrences
    call psi(a(1:11),a(12:18),s)
    call Concurrence(s,Cs)

    ! calculate tau's
    T(1) = Cs(1) - (Cs(5) + Cs(6) + Cs(7))
    T(2) = Cs(2) - (Cs(5) + Cs(8) + Cs(9))
    T(3) = Cs(3) - (Cs(6) + Cs(8) + Cs(10))
    T(4) = Cs(4) - (Cs(7) + Cs(9) + Cs(10))

    ! denominator term defined in notes, Eq. (13).
    Denom = X(3)*sum(T) + sum(Cs(5:10))
    if (Denom /= 0.0_8) then
        call Numer(1,2,3,4,X,Cs,tmp)
        Delta = -tmp/Denom
    else
        Delta = 0.0_8
    end if

    r = Cs(5) + Delta
    return
end subroutine ObjectiveI1

! Return vector of sigma_ij's for (i,j) in {1..4} with i<j.
! Input angles a(1:11)={11θ}, a(12:18)={7ϕ}, and array X={x,y,z}.
subroutine sigmas(a,X,r)
    implicit none
    real(8), intent(out) :: r(6)
    real(8), intent(in) :: a(:), X(:)
    real(8) :: Cs(10), T(4), Denom, Delta(6), tmp
    complex(8) :: s(0:1,0:1,0:1,0:1)
    integer :: n

    ! calculate concurrences
    call psi(a(1:11),a(12:18),s)
    call Concurrence(s,Cs)

    ! calculate tau's
    T(1) = Cs(1) - (Cs(5) + Cs(6) + Cs(7))
    T(2) = Cs(2) - (Cs(5) + Cs(8) + Cs(9))
    T(3) = Cs(3) - (Cs(6) + Cs(8) + Cs(10))
    T(4) = Cs(4) - (Cs(7) + Cs(9) + Cs(10))

    ! denominator term defined in notes, Eq. (13).
    Denom = X(3)*sum(T) + sum(Cs(5:10))
    if (Denom /= 0.0_8) then
        call Numer(1,2,3,4,X,Cs,tmp); Delta(1) = -tmp/Denom;
        call Numer(1,3,2,4,X,Cs,tmp); Delta(2) = -tmp/Denom;
        call Numer(1,4,2,3,X,Cs,tmp); Delta(3) = -tmp/Denom;
        call Numer(2,3,1,4,X,Cs,tmp); Delta(4) = -tmp/Denom;
        call Numer(2,4,1,3,X,Cs,tmp); Delta(5) = -tmp/Denom;
        call Numer(3,4,1,2,X,Cs,tmp); Delta(6) = -tmp/Denom;
    else
        Delta = 0.0_8
    end if

    ! r is ordered as: (σ12, σ13, σ14, σ23, σ24, σ34)
    ForAll(n=1:6) r(n) = Cs(n+4) + Delta(n)

    return
end subroutine sigmas

! Fetch the appropriate one-to-one concurrence value from the Cs array.
subroutine getConc(Cs,a,b,r)
    implicit none
    integer, parameter :: num = selected_real_kind(15,307)
    real(num), intent(out) :: r
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
end subroutine getConc

! Sort a 4D array in ascending order,
! x->y such that y(4)>y(3)>y(2)>y(1).
subroutine sort4(x,y)
    implicit none
    integer, parameter :: num = selected_real_kind(15,307)
    real(num), intent(out) :: y(4)
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
end subroutine sort4

! Order 4-qubit state components in 4-by-4 matrix form.
!  |0000>, |0001>, |0010>, |0011>
!  |0100>, |0101>, |0110>, |0111>
!  |1000>, |1001>, |1010>, |1011>
!  |1100>, |1101>, |1110>, |1111>
subroutine to4x4(A,B)
    implicit none
    integer, parameter :: num = selected_real_kind(15,307)
    complex(num), intent(out) :: B(4,4)
    complex(num), intent(in) :: A(0:1,0:1,0:1,0:1)
    B(:,1) = (/A(0,0,0,0), A(0,1,0,0), A(1,0,0,0), A(1,1,0,0)/)
    B(:,2) = (/A(0,0,0,1), A(0,1,0,1), A(1,0,0,1), A(1,1,0,1)/)
    B(:,3) = (/A(0,0,1,0), A(0,1,1,0), A(1,0,1,0), A(1,1,1,0)/)
    B(:,4) = (/A(0,0,1,1), A(0,1,1,1), A(1,0,1,1), A(1,1,1,1)/)
    return
end subroutine to4x4
