program compute_bandstructure

    use fem
    use class_linkedlist
    use potential


    implicit none
    
    include 'mpif.h'

    type :: matrix

        integer,dimension(:),allocatable :: col
        complex(8),dimension(:),allocatable :: A,S
        complex(8),dimension(:,:),allocatable :: AS

    end type



    integer :: dummy,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,Ne,Nn,Nlocal,Nquadrature,Nn_prime,Nxyz,Ne_sub,Ne_sub_remainder
    integer :: g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,neigh_size
    double precision :: dummy1,dummy2,dummy3,dummy4
    integer :: i,k,l,ii,ll,m,n,jj,kk,Nk,Nk1,Nk2,Nk3,Nk4,l1,l2,l3,l4,l5,l6,l7,mm,k_prime
    logical :: mtest
    character(len=1) :: dummyc
    character(len=2),dimension(:),allocatable :: element

    double precision :: r00,tmp00,cut


    double precision :: A_temp, V_temp, c_potential, x_real, y_real, x_imag, y_imag, lc, pi, pot,ik,k2
    complex(8) :: x,y,j
    complex(8),dimension(:),allocatable :: xy,ab,ab00

    double precision,dimension(:,:),allocatable :: Vhxc_g,nq_g,Vi_g,Vext_g
    double precision :: Vhxc_g_tmp

    double precision,dimension(:,:),allocatable :: kpoints,kpoints_sweep
    
    integer, dimension(:), allocatable :: index,list,color,color_prime,color1,color2,color3
    
    double precision, dimension(:), allocatable :: gp, basis_m,basis_n, location,location1,kx,ky,kz,linspace

    double precision, dimension(:,:), allocatable :: del_m,del_n,potential999
    type(matrix), dimension(:), allocatable :: row,row0,row00
    type(LinkedList) ,dimension(:),allocatable:: neigh,neigh_temp,neigh1,neigh2,neigh3,neigh4,neigh5,neigh6,neigh7

    integer, dimension(:,:), allocatable :: ele_prime
    double precision, dimension(:,:), allocatable :: point_prime,point_vi


    double precision :: latticeconstant,d,temp
    integer :: Ml

    integer :: index1,index2
    integer,dimension(:),allocatable :: Zc

    character(len=10) :: file_id
    character(len=50) :: file_name

    double precision :: latticeconstant1,latticeconstant2,lcx,lcz
    
    
     double precision, dimension(:,:), allocatable :: Ek
     integer, dimension(:), allocatable :: IA, JA, IB, JB, IA00
     complex(8), dimension(:), allocatable :: A,B,Amatrix_temp,Bmatrix_temp
    
    
    ! input parameters for FEAST

     character(len=1) :: UPLO='F' 

     integer, dimension(64) :: fpm
     integer :: M0,M0_1,M0_2,M0_plus ! search subspace dimension
     double precision :: Emin, Emax, Emin_1,Emax_1,Emin_2,Emax_2,Emin_plus,Emax_plus ! search interval
    ! output variables for FEAST
     double precision, dimension(:), allocatable :: E, res
     double precision, dimension(:,:), allocatable :: psi
     double precision :: epsout
     integer :: loop, info, M00, M00_plus
     
     
     !!!! MPI parameters 
     
    !  integer :: code, nb_procs, rank,NEW_COMM_WORLD
     integer :: code,rank,nb_procs,lrank,lnb_procs,color_mpi,key,NEW_COMM_WORLD
     complex(8) :: sum1,sum2
     character(len=3) :: cnL3
     integer :: nL3
     double precision :: start_time,finish_time

     integer, dimension(MPI_STATUS_SIZE) :: mpi_recv_status
     integer :: M000
     integer, parameter :: tag=100


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_INIT(code)

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code )
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code )


  start_time  = MPI_Wtime()




    




  if (rank==0) print *, "++++++++++++ load files ++++++++++++"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Mesh definition
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    open(10,file='Si_primitivecell_crown_new.1_p3_new.node',status='old')
    read(10,*) Nn 
    allocate(point(1:Nn,1:3))
    allocate(color(1:Nn))
    do i=1,Nn
    read(10,*) dummy,point(i,1),point(i,2),point(i,3)!,color(i)!,color1(i),color2(i),color3(i)
    enddo
    close(10)

    point=point/0.529177249d0 ! convert to atomic units

    


    open(10,file='Si_primitivecell_crown_new.1_p3_new.ele',status='old')
    read(10,*) Ne
    allocate(ele(1:Ne,1:20))
    do i=1,Ne
    read(10,*) dummy,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20
    ele(i,1)=g1+1
    ele(i,2)=g2+1
    ele(i,3)=g3+1
    ele(i,4)=g4+1
    ele(i,5)=g5+1
    ele(i,6)=g6+1
    ele(i,7)=g7+1
    ele(i,8)=g8+1
    ele(i,9)=g9+1
    ele(i,10)=g10+1
    ele(i,11)=g11+1
    ele(i,12)=g12+1
    ele(i,13)=g13+1
    ele(i,14)=g14+1
    ele(i,15)=g15+1
    ele(i,16)=g16+1
    ele(i,17)=g17+1
    ele(i,18)=g18+1
    ele(i,19)=g19+1
    ele(i,20)=g20+1
    end do

    close(10)

    

    open(10,file='Si_primitivecell_crown_new.1_p3_new_boundary.node',status='old')
    read(10,*) Nn_prime
    
    allocate(point_prime(1:Nn_prime,1:3))
    allocate(color_prime(1:Nn_prime))
    
    do i=1,Nn_prime
    read(10,*) dummy,point_prime(i,1),point_prime(i,2),point_prime(i,3)!,color_prime(i)
    enddo
    close(10)

    point_prime=point_prime/0.529177249d0 ! convert to atomic units


    open(10,file='Si_primitivecell_crown_new.1_p3_new_boundary.ele',status='old')
    read(10,*) Ne
    allocate(ele_prime(1:Ne,1:20))

    do i=1,Ne
    read(10,*) dummy,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20
    ele_prime(i,1)=g1+1
    ele_prime(i,2)=g2+1
    ele_prime(i,3)=g3+1
    ele_prime(i,4)=g4+1
    ele_prime(i,5)=g5+1
    ele_prime(i,6)=g6+1
    ele_prime(i,7)=g7+1
    ele_prime(i,8)=g8+1
    ele_prime(i,9)=g9+1
    ele_prime(i,10)=g10+1
    ele_prime(i,11)=g11+1
    ele_prime(i,12)=g12+1
    ele_prime(i,13)=g13+1
    ele_prime(i,14)=g14+1
    ele_prime(i,15)=g15+1
    ele_prime(i,16)=g16+1
    ele_prime(i,17)=g17+1
    ele_prime(i,18)=g18+1
    ele_prime(i,19)=g19+1
    ele_prime(i,20)=g20+1
    end do

    close(10)



    









    open(10,file='Si_primitivecell_crown_new.1_p3_gp_nessie.v',status='old')  !!!! Si
    allocate(Vhxc_g(1:Ne,1:11))
    do i=1,Ne
        do k=1,11
        read(10,*) Vhxc_g(i,k)!,Vi_g(i,k)
        enddo
    enddo
    close(10)



    open(10,file='Si_primitivecell_crown_new.1_p3_gp_nessie.vext',status='old')  !!!! Si
    allocate(Vext_g(1:Ne,1:11))
    do i=1,Ne
        do k=1,11
        read(10,*) Vext_g(i,k)!Vhxc_g_tmp
        enddo
    enddo
    close(10)

    Vext_g=Vext_g/27.2113961d0


    Vhxc_g=Vhxc_g/27.2113961d0




    open(10,file='Si26H42.xyz',status='old')  !!!! Si
    read(10,*) Nxyz
    read(10,*) dummyc
    allocate(point_vi(1:Nxyz,1:3))
    allocate(Zc(1:Nxyz))
    allocate(element(1:Nxyz))
    do i=1,Nxyz
        read(10,*) element(i),point_vi(i,1),point_vi(i,2),point_vi(i,3)!,Zc(i)
        ! read(10,*) point_vi(i,1),point_vi(i,2),point_vi(i,3)
    enddo
    close(10)

    point_vi=point_vi/0.529177249d0





    Nquadrature=11
    allocate(location(1:3))


    ! l=0

    ! cut=3.0d0!1.95d0

    ! do i=1,Ne

    !     call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    !     call gaussian_integral

    !     do k=1,Nquadrature

    !         location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

    !         tmp00=0.0d0

    !         do m=1,Nxyz

    !             r00=sqrt((location(1)-point_vi(m,1))**2+(location(2)-point_vi(m,2))**2+(location(3)-point_vi(m,3))**2)

    !             if (r00<cut) then
    !                 tmp00=tmp00+(1.0d0-(r00/cut)**8)**3
    !             else
    !                 tmp00=tmp00+0.0d0
    !             end if

    !         enddo

    !         if (tmp00>1.0d0) then
    !             tmp00=1.0d0
    !             l=l+1
    !         end if

    !         Vext_g(i,k)=tmp00*Vext_g(i,k)

    !     enddo
    ! enddo

    ! ! if (rank==0) print *,l

    ! Vhxc_g=Vhxc_g+Vext_g





    




    


latticeconstant= 5.4437024d0 ! 5.59d0 for NaCl ! 5.4437d0 for Si / 3.5671d0 for C
latticeconstant=latticeconstant/0.529177249d0 ! convert to atomic units

Ml=1        ! How many neighboring cell that need to be accountted
d=4.0d0     ! Atom radius by angstrom
d=d/0.529177249d0

pi=4.0d0*atan(1.0d0)

Nlocal=20
Nquadrature=11



j=cmplx(0.0d0,1.0d0,8)









! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!! Construct FE matrix without MPI parallelization !!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print *,"+++++++++++++"

! allocate(del_m(1:Nlocal,1:3))
! allocate(del_n(1:Nlocal,1:3))
! allocate(basis_m(1:Nlocal))
! allocate(basis_n(1:Nlocal))
! allocate(location(1:3))
! allocate(gp(1:3))


! allocate(neigh(Nn_prime))
! ! allocate(neigh(Nn))

! ll=0

! do i=1,Ne

!     ll=ll+1

!     if (ll==Ne/3) then
!       print *,ll
!     else if (ll==2*Ne/3) then
!       print *,LL
!     else if (ll==Ne) then
!       print *,ll
!     end if

!     call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))    
!     call mass_mat
!     call gaussian_integral

!     do m=1,Nlocal
        
!         do n=1,Nlocal
        
!             A_temp=0
!             V_temp=0
!             ik=0
!             k2=0
            
!             do k=1,Nquadrature
            
            
!                 !del_m = basis_p2_del(gpoint(k,:))
!                 !del_n = basis_p2_del(gpoint(k,:))
!                 !basis_m = basis_p2(gpoint(k,:))
!                 !basis_n = basis_p2(gpoint(k,:))
                
!                 del_m = basis_p3_del(gpoint(k,:))
!                 del_n = basis_p3_del(gpoint(k,:))
!                 basis_m = basis_p3(gpoint(k,:))
!                 basis_n = basis_p3(gpoint(k,:))
                
!                 A_temp = A_temp + dot_product(matmul(Jit,del_m(m,:)),matmul(Jit,del_n(n,:)))*gweight(k)

!                 ! location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
!                 ! c_potential = coulombpotential(location)
!                 ! c_potential = latticepotential(location,latticeconstant,Ml,d)
!                 ! V_temp = V_temp + c_potential*basis_m(m)*basis_n(n)*gweight(k)
!                 ! V_temp = V_temp + potential999(i,k)*basis_m(m)*basis_n(n)*gweight(k)


!                 ! temp = Vhxc_g(i,k)+Vi_g(i,k)!+Vx(i,k)+Vc(i,k)
                
!                 ! if (temp>0.0d0) then
!                 !     V_temp = V_temp + 0.0d0
!                 ! else
!                 !     V_temp = V_temp + temp*basis_m(m)*basis_n(n)*gweight(k)
!                 ! ! V_temp = V_temp + (Vhxc_g(i,k)+Vi_g(i,k)-7.0d0)*basis_m(m)*basis_n(n)*gweight(k)
!                 ! endif



!                 temp=Vhxc_g(i,k)
!                 location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
        
!                 do ii=1,Nxyz
! !temp=temp-float(Zc(ii))/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
!                     ! if (ii<=26) then !!! 5ppp
!                     ! if (ii<=18) then !!! 3ppp
!                     ! if (ii<=16) then !!! 3ppp
!                     if (element(ii)=='C') then
!  temp = temp-6.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     
! !                    else if (element(ii)=='O') then
! ! temp = temp-8.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
!                      else
!  temp = temp-1.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     
!                      end if
!                 enddo


!                 V_temp = V_temp + temp*basis_m(m)*basis_n(n)*gweight(k)


!             enddo

!             x_real=0.5d0*A_temp*(Jdet)/6.0d0+V_temp*(Jdet)/6.0d0
!             y_real=Jdet*p3_matrix(m,n)

!             x_imag=0.0d0
!             y_imag=0.0d0
            
!             x=cmplx(x_real,x_imag,8)
!             y=cmplx(y_real,y_imag,8)


!             ! x=0.5d0*A_temp*(Jdet)/6.0d0+V_temp*(Jdet)/6.0d0
!             ! y=Jdet*p2_matrix(m,n)
            
!             mtest=.true.
!             if (neigh(ele_prime(i,m))%find(ele_prime(i,n)))  mtest=.false.
!             if (mtest) call neigh(ele_prime(i,m))%insert(ele_prime(i,n))
!             call neigh(ele_prime(i,m))%insertAB(ele_prime(i,n),x,y)
!             ! if (neigh(ele(i,m))%find(ele(i,n)))  mtest=.false.
!             ! if (mtest) call neigh(ele(i,m))%insert(ele(i,n))
!             ! call neigh(ele(i,m))%insertAB(ele(i,n),x,y)


!             ! A(ele(i,m),ele(i,n))=A(ele(i,m),ele(i,n))+0.5d0*A_temp*(Jdet)/6.0d0+V_temp*(Jdet)/6.0d0
!             ! S(ele(i,m),ele(i,n))=S(ele(i,m),ele(i,n))+Jdet*p2_matrix(m,n)
!         enddo
!     enddo
! enddo


! neigh_size = 0

! do i=1,Nn_prime
!   do ii=1,neigh(i)%size
!         neigh_size=neigh_size+1
!   enddo
! enddo


! allocate(row(1:Nn_prime))
! ! allocate(row(1:Nn))
! allocate(ab(1:2))

! ll=0
! do i=1,Nn_prime
! ! do i=1,Nn
!     allocate(row(i)%col(1:neigh(i)%size))
!     allocate(row(i)%A(1:neigh(i)%size))
!     allocate(row(i)%S(1:neigh(i)%size))
!     do ii=1,neigh(i)%size
!         row(i)%col(ii)=neigh(i)%get(ii)
!         ab=neigh(i)%getAB(ii)
!         row(i)%A(ii)=ab(1)
!         row(i)%S(ii)=ab(2)        
!         ! ll=ll+1
!     enddo
! enddo

! print *,"++++++ Constructed FE matrix without MPI parallelization ++++++"

! ! stop


! ! allocate(xy(1:2))
! ! ! open(12,file='C_primitivecell_crown.1_p20_AS.txt',action="write",status='replace')
! ! open(12,file='C_unitcell.1_p2_AS.txt',action="write",status='replace')
! ! do i=1,Nn_prime
! ! ! do i=1,Nn
! !     do l=1,neigh(i)%size
! !         xy=neigh(i)%getab(l)
! ! write(12,*) i-1,neigh(i)%get(l)-1,real(xy(1)),aimag(xy(1)),real(xy(2)),aimag(xy(2))
! ! end do
! ! end do
! ! close(12)

! ! stop


















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!   Construct FE matrix with MPI parallelization  !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (rank==0) print *,"++++++ Construct FE Hamiltonian matrix using MPI parallelization ++++++"

! print *,"+++++++++++++"

allocate(del_m(1:Nlocal,1:3))
allocate(del_n(1:Nlocal,1:3))
allocate(basis_m(1:Nlocal))
allocate(basis_n(1:Nlocal))
! allocate(location(1:3))
allocate(gp(1:3))


allocate(neigh(Nn_prime))


do i=1,Ne

    do m=1,Nlocal
        
        do n=1,Nlocal

            mtest=.true.
            if (neigh(ele_prime(i,m))%find(ele_prime(i,n)))  mtest=.false.
            if (mtest) call neigh(ele_prime(i,m))%insert(ele_prime(i,n))


        enddo
    enddo
enddo






Ne_sub = Ne/nb_procs
Ne_sub_remainder = mod(Ne,nb_procs)



do i=rank*Ne_sub+1,rank*Ne_sub+Ne_sub

    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))    
    call mass_mat
    call gaussian_integral

    do m=1,Nlocal
        
        do n=1,Nlocal

            A_temp=0
            V_temp=0
            ik=0
            k2=0
            
            do k=1,Nquadrature
            
            
                ! del_m = basis_p2_del(gpoint(k,:))
                ! del_n = basis_p2_del(gpoint(k,:))
                ! basis_m = basis_p2(gpoint(k,:))
                ! basis_n = basis_p2(gpoint(k,:))
                
                del_m = basis_p3_del(gpoint(k,:))
                del_n = basis_p3_del(gpoint(k,:))
                basis_m = basis_p3(gpoint(k,:))
                basis_n = basis_p3(gpoint(k,:))
                
                A_temp = A_temp + dot_product(matmul(Jit,del_m(m,:)),matmul(Jit,del_n(n,:)))*gweight(k)




                temp=Vhxc_g(i,k)
                location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
        
                do ii=1,Nxyz
!temp=temp-float(Zc(ii))/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                    ! if (ii<=26) then !!! 5ppp
                    ! if (ii<=18) then !!! 3ppp
                    ! if (ii<=16) then !!! 3ppp
                    if (element(ii)=='Si') then
 temp = temp-14.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     
!                    else if (element(ii)=='O') then
! temp = temp-8.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     else
 temp = temp-1.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     
                     end if
                enddo


                V_temp = V_temp + temp*basis_m(m)*basis_n(n)*gweight(k)


            enddo

            x_real=0.5d0*A_temp*(Jdet)/6.0d0+V_temp*(Jdet)/6.0d0
            y_real=Jdet*p3_matrix(m,n)

            x_imag=0.0d0
            y_imag=0.0d0
            
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)


            call neigh(ele_prime(i,m))%insertAB(ele_prime(i,n),x,y)

            ! row(m)%A(n)=row(m)%A(n)+x
            ! row(m)%S(n)=row(m)%S(n)+y

        enddo
    enddo
enddo


if (rank==nb_procs-1) then

    ! print *,rank*Ne_sub+Ne_sub+1,rank*Ne_sub+Ne_sub+Ne_sub_remainder

do i=rank*Ne_sub+Ne_sub+1,rank*Ne_sub+Ne_sub+Ne_sub_remainder

    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))    
    call mass_mat
    call gaussian_integral

    do m=1,Nlocal
        
        do n=1,Nlocal

            A_temp=0
            V_temp=0
            ik=0
            k2=0
            
            do k=1,Nquadrature
            
            
                ! del_m = basis_p2_del(gpoint(k,:))
                ! del_n = basis_p2_del(gpoint(k,:))
                ! basis_m = basis_p2(gpoint(k,:))
                ! basis_n = basis_p2(gpoint(k,:))
                
                del_m = basis_p3_del(gpoint(k,:))
                del_n = basis_p3_del(gpoint(k,:))
                basis_m = basis_p3(gpoint(k,:))
                basis_n = basis_p3(gpoint(k,:))
                
                A_temp = A_temp + dot_product(matmul(Jit,del_m(m,:)),matmul(Jit,del_n(n,:)))*gweight(k)

                ! location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                ! c_potential = coulombpotential(location)
                ! c_potential = latticepotential(location,latticeconstant,Ml,d)
                ! V_temp = V_temp + c_potential*basis_m(m)*basis_n(n)*gweight(k)
                ! V_temp = V_temp + potential999(i,k)*basis_m(m)*basis_n(n)*gweight(k)


                ! temp = Vhxc_g(i,k)+Vi_g(i,k)!+Vx(i,k)+Vc(i,k)
                
                ! if (temp>0.0d0) then
                !     V_temp = V_temp + 0.0d0
                ! else
                !     V_temp = V_temp + temp*basis_m(m)*basis_n(n)*gweight(k)
                ! ! V_temp = V_temp + (Vhxc_g(i,k)+Vi_g(i,k)-7.0d0)*basis_m(m)*basis_n(n)*gweight(k)
                ! endif



                temp=Vhxc_g(i,k)
                location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
        
                do ii=1,Nxyz
!temp=temp-float(Zc(ii))/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                    ! if (ii<=26) then !!! 5ppp
                    ! if (ii<=18) then !!! 3ppp
                    ! if (ii<=16) then !!! 3ppp
                    if (element(ii)=='Si') then
 temp = temp-14.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     
!                    else if (element(ii)=='O') then
! temp = temp-8.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     else
 temp = temp-1.0d0/sqrt((location(1)-point_vi(ii,1))**2+(location(2)-point_vi(ii,2))**2+(location(3)-point_vi(ii,3))**2)
                     
                     end if
                enddo


                V_temp = V_temp + temp*basis_m(m)*basis_n(n)*gweight(k)


            enddo

            x_real=0.5d0*A_temp*(Jdet)/6.0d0+V_temp*(Jdet)/6.0d0
            y_real=Jdet*p3_matrix(m,n)

            x_imag=0.0d0
            y_imag=0.0d0
            
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)


            call neigh(ele_prime(i,m))%insertAB(ele_prime(i,n),x,y)

            ! row(m)%A(n)=row(m)%A(n)+x
            ! row(m)%S(n)=row(m)%S(n)+y

        enddo
    enddo
enddo


end if



call MPI_BARRIER(MPI_COMM_WORLD ,code)

neigh_size = 0

do i=1,Nn_prime
  do ii=1,neigh(i)%size
        neigh_size=neigh_size+1
  enddo
enddo




allocate(row(1:Nn_prime))
! allocate(row(1:Nn))
allocate(ab(1:2))

ll=0
do i=1,Nn_prime
! do i=1,Nn
    allocate(row(i)%col(1:neigh(i)%size))
    allocate(row(i)%A(1:neigh(i)%size))
    allocate(row(i)%S(1:neigh(i)%size))
    allocate(row(i)%AS(1:neigh(i)%size,1:2))
    do ii=1,neigh(i)%size
        row(i)%col(ii)=neigh(i)%get(ii)
        ! ab=neigh(i)%getAB(ii)
        ! row(i)%A(ii)=ab(1)
        ! row(i)%S(ii)=ab(2)        
        ! ll=ll+1
        ! row(i)%A(ii)=cmplx(0.0d0,0.0d0,8)
        ! row(i)%S(ii)=cmplx(0.0d0,0.0d0,8) 
    enddo
enddo





allocate(IA(1:Nn_prime+1))
allocate(JA(1:neigh_size))
allocate(A(1:neigh_size))
allocate(B(1:neigh_size))
allocate(Amatrix_temp(1:neigh_size))
allocate(Bmatrix_temp(1:neigh_size))


do i=1,Nn_prime
    do ii=1,neigh(i)%size
call MPI_ALLREDUCE(neigh(i)%getAB(ii),row(i)%AS(ii,1:2),2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)
    enddo
enddo








































if (rank==0) print *,"++++++ Compute the first Brillouin zone high-symmetry k points ++++++"


allocate(kpoints(1:6,1:3))


kpoints(1,:)=(/0.0d0,0.0d0,0.0d0/)    !!! Gamma 
! kpoints(2,:)=(/2.0d0/3.0d0,0.0d0,0.0d0/)     !!! X
kpoints(2,:)=(/1.0d0,0.0d0,0.0d0/)     !!! X
kpoints(3,:)=(/0.5d0,0.5d0,0.5d0/)    !!! L
! kpoints(3,:)=(/1.0d0/2.0d0,0.0d0,1.0d0/(2.0d0*sqrt(3.0d0))/)    !!! L
kpoints(4,:)=(/0.5d0,1.0d0,0.0d0/)    !!! W
!kpoints(5,:)=(/0.25d0,1.0d0,0.25d0/)  !!! U
!kpoints(5,:)=(/1.0d0/3.0d0,sqrt(3.0d0)/3.0d0,0.0d0/)  !!! U
kpoints(6,:)=(/0.75d0,0.75d0,0.0d0/)  !!! K

! kpoints(6,:)=(/0.5d0,0.0d0,0.5d0/)  !!! K



lc=latticeconstant



Nk1 = 7
Nk2 = 7
Nk3 = 3
Nk4 = 5
Nk = Nk1+Nk2+Nk3+Nk4-2

allocate(kx(1:Nk))
allocate(ky(1:Nk))
allocate(kz(1:Nk))


!!!! L -> Gamma !!!!
do i=1,Nk1
    kx(i)=kpoints(3,1)+(kpoints(1,1)-kpoints(3,1))/(dfloat(Nk1-1))*dfloat((i-1))
    ky(i)=kpoints(3,2)+(kpoints(1,2)-kpoints(3,2))/(dfloat(Nk1-1))*dfloat((i-1))
    kz(i)=kpoints(3,3)+(kpoints(1,3)-kpoints(3,3))/(dfloat(Nk1-1))*dfloat((i-1))
enddo
!!!! Gamma -> X !!!!
do i=1,Nk2
    kx(Nk1-1+i)=kpoints(1,1)+(kpoints(2,1)-kpoints(1,1))/(dfloat(Nk2-1))*dfloat((i-1))
    ky(Nk1-1+i)=kpoints(1,2)+(kpoints(2,2)-kpoints(1,2))/(dfloat(Nk2-1))*dfloat((i-1))
    kz(Nk1-1+i)=kpoints(1,3)+(kpoints(2,3)-kpoints(1,3))/(dfloat(Nk2-1))*dfloat((i-1))
enddo
!!!! X -> W !!!!
do i=1,Nk3
    kx(Nk1+Nk2-1+i)=kpoints(4,1)+(kpoints(6,1)-kpoints(4,1))/(dfloat(Nk3-1))*dfloat((i-1))
    ky(Nk1+Nk2-1+i)=kpoints(4,2)+(kpoints(6,2)-kpoints(4,2))/(dfloat(Nk3-1))*dfloat((i-1))
    kz(Nk1+Nk2-1+i)=kpoints(4,3)+(kpoints(6,3)-kpoints(4,3))/(dfloat(Nk3-1))*dfloat((i-1))
enddo
!!!! W -> K !!!!
do i=1,Nk4
    kx(Nk1+Nk2+Nk3-2+i)=kpoints(6,1)+(kpoints(1,1)-kpoints(6,1))/(dfloat(Nk4-1))*dfloat((i-1))
    ky(Nk1+Nk2+Nk3-2+i)=kpoints(6,2)+(kpoints(1,2)-kpoints(6,2))/(dfloat(Nk4-1))*dfloat((i-1))
    kz(Nk1+Nk2+Nk3-2+i)=kpoints(6,3)+(kpoints(1,3)-kpoints(6,3))/(dfloat(Nk4-1))*dfloat((i-1))
enddo




kx=kx*2.0d0*pi/lc!*2.0d0
ky=ky*2.0d0*pi/lc
kz=kz*2.0d0*pi/lc!*2.0d0



allocate(linspace(1:11))

linspace(1:11)=(/-0.5d0,-0.4d0,-0.3d0,-0.2d0,-0.1d0,0.0d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0/)















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST with L2 !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

M0=60
Emin=-80.0d0!-15.0d0
Emax= 0.5d0!-0.1d0

M0_1=2
Emin_1=-105.0d0
Emax_1= 0.0d0

M0_2=2
Emin_2=-39.0d0
Emax_2=-35.0d0

M0_plus=20
Emin_plus=0.0d0
Emax_plus=0.5d0




!!! Allocate memory for eigenvalues. eigenvectors, residual 
allocate(E(M0), psi(Nn_prime,M0), res(M0))


! if (rank==0) allocate(Ek(1:M0+M0_plus,1:13))
allocate(Ek(1:M0+M0_plus,1:13))
!allocate(Ek(1:M0,1:143))



do k=1,13


if (rank==0) print *, "---------- cycle",k,"started ----------"


allocate(row0(1:Nn_prime))
do i=1,Nn_prime
! allocate(row0(1:Nn))
! do i=1,Nn
    allocate(row0(i)%col(1:size(row(i)%col)))
    allocate(row0(i)%A(1:size(row(i)%col)))
    allocate(row0(i)%S(1:size(row(i)%col)))
    allocate(row0(i)%AS(1:size(row(i)%col),1:2))
    do ii=1,size(row(i)%col)
        row0(i)%A(ii)=cmplx(0.0d0,0.0d0,8)
        row0(i)%S(ii)=cmplx(0.0d0,0.0d0,8)
    enddo
enddo



if (rank/=nb_procs-1) then


do i=rank*Ne_sub+1,rank*Ne_sub+Ne_sub

    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    call mass_mat
    call gaussian_integral
        
        

    do m=1,Nlocal
        do n=1,Nlocal
        
        
        
            ik=0.0d0
            k2=0.0d0
            do l=1,Nquadrature
                ! del_m = basis_p2_del(gpoint(l,:))
                ! del_n = basis_p2_del(gpoint(l,:))
                ! basis_m = basis_p2(gpoint(l,:))
                ! basis_n = basis_p2(gpoint(l,:))
                
                del_m = basis_p3_del(gpoint(l,:))
                del_n = basis_p3_del(gpoint(l,:))
                basis_m = basis_p3(gpoint(l,:))
                basis_n = basis_p3(gpoint(l,:))

                ik = ik + (dot_product((/kx(k),ky(k),kz(k)/),matmul(Jit,del_n(n,:)))*basis_m(m)&
                          -dot_product((/kx(k),ky(k),kz(k)/),matmul(Jit,del_m(m,:)))*basis_n(n))*gweight(l)

                k2 = k2 + (kx(k)**2+ky(k)**2+kz(k)**2)*basis_m(m)*basis_n(n)*gweight(l)
            
            enddo


            x_real=0.5d0*k2*(Jdet)/6.0d0
            y_real=0.0d0

            x_imag=0.5d0*ik*(Jdet)/6.0d0
            y_imag=0.0d0

            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)

       
            index1=neigh(ele_prime(i,m))%findindex(ele_prime(i,n))
            row0(ele_prime(i,m))%A(index1)=row0(ele_prime(i,m))%A(index1)+x
            row0(ele_prime(i,m))%S(index1)=row0(ele_prime(i,m))%S(index1)+y
            ! index1=neigh(ele(i,m))%findindex(ele(i,n))
            ! row0(ele(i,m))%A(index1)=row0(ele(i,m))%A(index1)+x
            ! row0(ele(i,m))%S(index1)=row0(ele(i,m))%S(index1)+y

            ! call neigh2(ele_prime(i,m))%insertAB(ele_prime(i,n),x,y)

        enddo
    enddo
enddo




else



do i=rank*Ne_sub+Ne_sub+1,rank*Ne_sub+Ne_sub+Ne_sub_remainder
    


    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    call mass_mat
    call gaussian_integral
        
        

    do m=1,Nlocal
        do n=1,Nlocal
        
        
        
            ik=0.0d0
            k2=0.0d0
            do l=1,Nquadrature
                ! del_m = basis_p2_del(gpoint(l,:))
                ! del_n = basis_p2_del(gpoint(l,:))
                ! basis_m = basis_p2(gpoint(l,:))
                ! basis_n = basis_p2(gpoint(l,:))
                
                del_m = basis_p3_del(gpoint(l,:))
                del_n = basis_p3_del(gpoint(l,:))
                basis_m = basis_p3(gpoint(l,:))
                basis_n = basis_p3(gpoint(l,:))

                ik = ik + (dot_product((/kx(k),ky(k),kz(k)/),matmul(Jit,del_n(n,:)))*basis_m(m)&
                          -dot_product((/kx(k),ky(k),kz(k)/),matmul(Jit,del_m(m,:)))*basis_n(n))*gweight(l)

                k2 = k2 + (kx(k)**2+ky(k)**2+kz(k)**2)*basis_m(m)*basis_n(n)*gweight(l)
            
            enddo


            x_real=0.5d0*k2*(Jdet)/6.0d0
            y_real=0.0d0

            x_imag=0.5d0*ik*(Jdet)/6.0d0
            y_imag=0.0d0

            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)

       
            index1=neigh(ele_prime(i,m))%findindex(ele_prime(i,n))
            row0(ele_prime(i,m))%A(index1)=row0(ele_prime(i,m))%A(index1)+x
            row0(ele_prime(i,m))%S(index1)=row0(ele_prime(i,m))%S(index1)+y

            ! index1=neigh(ele(i,m))%findindex(ele(i,n))
            ! row0(ele(i,m))%A(index1)=row0(ele(i,m))%A(index1)+x
            ! row0(ele(i,m))%S(index1)=row0(ele(i,m))%S(index1)+y

            ! call neigh2(ele_prime(i,m))%insertAB(ele_prime(i,n),x,y)

        enddo
    enddo
enddo

end if



 ll=0
 IA(1)=1
 do jj=1,Nn_prime
 ! do jj=1,Nn
 IA(jj+1)=IA(jj)+size(row(jj)%col)
     do kk=1,size(row(jj)%col)
       ll=ll+1
       JA(ll)=row(jj)%col(kk)
       Amatrix_temp(ll)=row0(jj)%A(kk)+row(jj)%AS(kk,1)/nb_procs
       Bmatrix_temp(ll)=row0(jj)%S(kk)+row(jj)%AS(kk,2)/nb_procs
     enddo
 enddo






call MPI_ALLREDUCE(Amatrix_temp,A,neigh_size,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(Bmatrix_temp,B,neigh_size,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)








 call feastinit(fpm)
 !fpm(1)=1
 M0_1=30
 call zfeast_hcsrgv(UPLO, Nn_prime, A, IA, JA, B, IA, JA, fpm, epsout, loop, Emin_1, Emax_1, M0_1, E, psi, M00, res, info)

!  if (rank==0) print *, 'FEAST-info---------',info
!  if (rank==0) print *, 'FEAST-eigenvalue--- ',E(1)!,E(2),E(3),E(4)!,E(5),E(6),E(7),E(8),E(9),E(10)

!  Ek(1,k)=E(1)
!  Ek(:,(k_prime-1)*13+k)=E(1:40)
!  Ek(:,k-6)=E(1:40)
!  M0_2=2
!  call zfeast_hcsrgv(UPLO,Nn_prime,A,IA,JA,B,IA,JA,fpm,epsout,loop,Emin_2,Emax_2,M0_2,E,psi,M00,res,info)
 
!  if (rank==0) print *, 'FEAST-info---------',info
!  if (rank==0) print *, 'FEAST-eigenvalue---',E(1)!,E(2),E(3),E(4)
 
!  Ek(2,k)=E(1)
 
!  M0_1=40
! !  call zfeast_hcsrgv(UPLO,Nn_prime,A,IA,JA,B,IA,JA,fpm,epsout,loop,Emin,Emax,M0,E,psi,M00,res,info)
!  call zfeast_hcsrgv(UPLO, Nn_prime, A, IA, JA, B, IA, JA, fpm, epsout, loop, Emin_1, Emax_1, M0_1, E, psi, M00, res, info)
 




 if (rank==0) then
 print *, 'FEAST-info---------',info
 print *, 'FEAST-eigenvalues--',E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),E(12),E(13),E(14),E(15)
 Ek(1:M0,k)=E(1:M0)
 end if
 
 !print *,M00
 
 !stop
 
 !call zfeast_hcsrgv(UPLO,Nn_prime,A,IA,JA,B,IA,JA,fpm,epsout,loop,Emin_plus,Emax_plus,M0_plus,E,psi,M00_plus,res,info)

 !print *, 'FEAST-info---------',info
 !print *, 'FEAST-eigenvalue--- ',E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10)

 !Ek(M00+1:M00_plus,k)=E(1:M00_plus)
 !Ek(1:M00_plus,k)=E(1:M00_plus)

 !print *,M00_plus

 !stop
 
 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!! pFEAST !!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call pfeastinit(fpm,MPI_COMM_WORLD,nL3)
!  call pzfeast_hcsrgv(UPLO,Nn_prime,A,IA,JA,B,IA,JA,fpm,epsout,loop,Emin,Emax,M0,E,psi,M00,res,info)
 
!  !call pdfeast_scsrgv(UPLO,this%dH%N,this%dH%sa,this%dH%isa,&
!  !              this%dH%jsa,this%dS%sa,this%dS%isa,this%dS%jsa,this%fpm,epsout,&
!  !              loop,this%slice(k)%Emin,this%slice(k)%Emax,this%slice(k)%M0,this%slice(k)%E,&
!  !              this%slice(k)%X,this%slice(k)%M,this%slice(k)%res,info)
 
 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   if (rank==0) print *,'FEAST OUTPUT INFO',info
!   if ((info/=0).and.(rank==0)) print *,'PF90sparse_pdfeast_scsrgv -- failed'
!   if ((info==0).and.(rank==0)) then
!      print *,'PF90sparse_pdfeast_scsrgv -- success'
!      print *,'*************************************************'
!      print *,'************** REPORT ***************************'
!      print *,'*************************************************'
!      print *,'Eigenvalues/Residuals (inside interval)'
!      do i=1,M
!         print *,i,E(i),res(i)
!      enddo
!   endif
 
 
 
 



! enddo

!call MPI_FINALIZE(code)

!stop




















! 999 continue

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!! pFEAST with L123 !!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!! Definition of the two intervals 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!  if (rank<=nb_procs/2-1) then
!     color_mpi=1 ! first interval
!     else
!     color_mpi=2 ! second interval
!     endif 
    
!     !!!!!!!!!!!!!!!!! create new_mpi_comm_world
!     key=0
!     call MPI_COMM_SPLIT(MPI_COMM_WORLD,color_mpi,key,NEW_COMM_WORLD,code)
!     call MPI_COMM_RANK(NEW_COMM_WORLD,lrank,code) ! local rank
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!     !!! search interval [Emin,Emax] including M eigenpairs
!     if (color_mpi==1) then !! 1st interval
!      Emin=-15.0d0
!      Emax=-7.0d0 
!      M0=4
!     elseif(color_mpi==2) then !! 2nd interval
!      Emin=-7.0d0
!      Emax=-0.1d0 
!      M0=40
!     endif

!     ! print *,rank,lrank,color_mpi,key

!     ! call MPI_BARRIER(NEW_COMM_WORLD ,code)
!     ! call MPI_BARRIER(MPI_COMM_WORLD ,code)

!     ! stop





!     !!!!!!!!!!!!!! RUN INTERVALS in PARALLEL


!     ! !!!!!!!!!!!!! ALLOCATE VARIABLE 
!     ! allocate(E(1:M0))     ! Eigenvalue
!     ! allocate(psi(1:Nn_prime,1:M0)) ! Eigenvectors
!     ! allocate(res(1:M0))   ! Residual 

!     !!!!!!!!!!!!!  FEAST
!     call pfeastinit(fpm,NEW_COMM_WORLD,nL3)

!     fpm(1)=-color_mpi !! print info the files feast1.log for contour1 and feast2.log for contour 2 
    
!     call pzfeast_hcsrgv(UPLO,Nn_prime,A,IA,JA,B,IA,JA,fpm,epsout,loop,Emin,Emax,M0,E,psi,M00,res,info)

!     ! print *, 'FEAST-info---------',info
!     ! print *, 'FEAST-eigenvalue---',E(1),E(2),E(3),E(4),E(5),E(6),E(7),E(8),E(9),E(10),E(11),E(12),E(13)

!     ! call MPI_BARRIER(MPI_COMM_WORLD ,code)

    

!     if ((info==0).and.(lrank==0)) then
!         print *,'Eigenvalues/Residuals - inside the interval ',color_mpi
!         do i=1,M00
!            print *,i,E(i),res(i)
!         enddo

!         ! if (color_mpi==1) then
!         if (rank==0) then
!             Ek(1:M00,k)=E(1:M00)
!             call MPI_SEND(M00,1,MPI_INTEGER,4,tag,MPI_COMM_WORLD,code)
!             call MPI_SEND(Ek,size(Ek),MPI_DOUBLE,4,tag,MPI_COMM_WORLD,code)
!             call MPI_RECV(Ek,size(Ek),MPI_DOUBLE,4,tag,MPI_COMM_WORLD,mpi_recv_status,code)
!             ! print *,Ek(1:10,1)
!             ! call MPI_BCAST(Ek,size(Ek),MPI_DOUBLE,rank,MPI_COMM_WORLD,code)
!         endif
!         ! endif

!         ! if (color_mpi==2) then
!         if (rank==4) then
!             call MPI_RECV(M000,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,mpi_recv_status,code)
!             call MPI_RECV(Ek,size(Ek),MPI_DOUBLE,0,tag,MPI_COMM_WORLD,mpi_recv_status,code)
!             ! print *,Ek(1:2,1)
!             Ek((M000+1):(M000+M0),k)=E(1:M0)
!             ! print *,Ek(1:10,1)
!             call MPI_SEND(Ek,size(Ek),MPI_DOUBLE,0,tag,MPI_COMM_WORLD,code)
!             ! call MPI_BCAST(Ek,size(Ek),MPI_DOUBLE,rank,MPI_COMM_WORLD,code)
!         endif
!         ! endif
!         ! if ((rank==nb_procs/2-1).and.(color_mpi==1)) Ek(1:M00,k)=E(1:M00)

!         ! if ((rank==nb_procs/2+1).and.(color_mpi==2)) then
!         !     ! Ek(3:M0,k)=E(1:(M0-2))
!         !     print *,E(1:4)
!         ! end if
!      endif

!     !  stop




    call MPI_BARRIER(MPI_COMM_WORLD ,code)

    
    do i=1,Nn_prime
    ! allocate(row0(1:Nn))
    ! do i=1,Nn
        deallocate(row0(i)%col)
        deallocate(row0(i)%A)
        deallocate(row0(i)%S)
        deallocate(row0(i)%AS)
        ! do ii=1,neigh(i)%size
        !     row0(i)%A(ii)=cmplx(/0.0d0,0.0d0,8/)
        !     row0(i)%S(ii)=cmplx(/0.0d0,0.0d0,8/)
        ! enddo
    enddo
    deallocate(row0)

    
    ! if (rank==0) print *, "-----------------------------------------------"

enddo


if (rank==0) then


 !open(12,file='NaCl_primitive_new.1_p3_new_E_test.txt',action="write",status='replace')
 !open(12,file='ppp_huasuan_unitcell.1_p3_new_E_test.txt',action="write",status='replace')
!  open(12,file='graphene_new_primitive_itat.1_p3_E_test.txt',action="write",status='replace')
!  open(12,file='diamond_primitivecell_crown.1_p2_new_C26H42_E_test.txt',action="write",status='replace')
 !open(12,file='cnt_primitivecell.1_p3_new_E_test.txt',action="write",status='replace')
!  open(12,file='Si_primitivecell_crown_new2.1_p3_Ebands.txt',action="write",status='replace')
    open(12,file='Si_primitivecell_crown_new.1_p3_Si26H42_new1_E_test.txt',action="write",status='replace')
 do i=1,size(Ek(:,1))
 write(12,*) Ek(i,1),Ek(i,2),Ek(i,3),Ek(i,4),Ek(i,5),Ek(i,6),Ek(i,7),Ek(i,8),Ek(i,9),Ek(i,10),Ek(i,11),Ek(i,12),Ek(i,13)&
         ,Ek(i,14),Ek(i,15),Ek(i,16),Ek(i,17),Ek(i,18),Ek(i,19),Ek(i,20)!,Ek(i,21),Ek(i,22),Ek(i,23),Ek(i,24),Ek(i,25),Ek(i,26)!&
         !,Ek(i,27),Ek(i,28),Ek(i,29),Ek(i,30),Ek(i,31),Ek(i,32),Ek(i,33),Ek(i,34),Ek(i,35),Ek(i,36)
 end do
 close(12)

 print *, "bandstructure calculation finished"
end if

call MPI_BARRIER(MPI_COMM_WORLD ,code)

finish_time  = MPI_Wtime()

if (rank==0) print *, finish_time-start_time, "seconds"

call MPI_FINALIZE(code)







end program compute_bandstructure
