module potential
implicit none


    double precision :: pi0
    double precision, dimension(:,:), allocatable :: Vx,Vc

    

    contains


    subroutine exchangepotential(n)

        double precision,dimension(:,:),allocatable :: n
        double precision :: pi


        pi=4.0d0*atan(1.0d0)

        Vx(:,:) = -(3.0d0/pi)**(1.0d0/3.0d0)*n(:,:)**(1.0d0/3.0d0)
    
    end subroutine exchangepotential



    subroutine correlationpotential(n)

        double precision,dimension(:,:),allocatable :: n
        double precision :: rs,eps,pi,Ec
        integer :: i,k

        pi=4.0d0*atan(1.0d0)

        do i=1,size(n(:,1))
            do k=1,11

            eps = 0.0d0
            if (n(i,k)==0.0d0) eps=1E-10

            rs = (3.0d0/4.0d0/pi/(n(i,k)+eps))**(1.0d0/3.0d0)

            if (rs<1.0d0) then

    Ec = n(i,k)*(0.0311d0*log(rs)-0.0480d0+0.002d0*rs*log(rs)-0.0116d0*rs)
    Vc(i,k) = 0.0311d0*log(rs)-0.0480d0-1.0d0/3.0d0*0.0311d0+2.0d0/3.0d0*0.0002d0*rs*log(rs)&
                                                            +1.0d0/3.0d0*(-2.0d0*0.0116d0-0.002d0)*rs

            else if (rs>=1.0d0) then

    Ec = n(i,k)*(-0.1423d0)/(1.0d0+1.0529d0*sqrt(rs)+0.3334d0*rs)
    Vc(i,k) = Ec/(n(i,k)+eps)*(1.0d0+7.0d0/6.0d0*1.0529d0*sqrt(rs)+4.0d0/3.0d0*0.3334d0*rs)/&
                                                            (1.0d0+1.0529d0*sqrt(rs)+0.3334d0*rs)

            end if

        enddo
            
        enddo


    end subroutine correlationpotential

    function coulombpotential(location1) result(coulombpotential1)

        ! type(potential) :: this
        double precision, dimension(1:3) :: location1
        double precision :: coulombpotential1

        coulombpotential1 = -1/sqrt(location1(1)**2+location1(2)**2+location1(3)**2)

    end function coulombpotential


    function coulombpotential01(location11,location2,d0) result(coulombpotential011)

        ! type(potential) :: this
        double precision, dimension(1:3) :: location11,location2
        double precision :: coulombpotential011,d0,dd

        dd=sqrt((location11(1)-location2(1))**2+(location11(2)-location2(2))**2+(location11(3)-location2(3))**2)

        coulombpotential011 = -48.0d0*exp(-1.0d0/4.0d0*dd**2) ! 14 for Si // 6 for C


    end function coulombpotential01


    function coulombpotential02(location11,location2,d0) result(coulombpotential021)

        ! type(potential) :: this
        double precision, dimension(1:3) :: location11,location2
        double precision :: coulombpotential021,d0,dd,c

        dd=sqrt((location11(1)-location2(1))**2+(location11(2)-location2(2))**2+(location11(3)-location2(3))**2)

        c=(1.42d0/1.5446d0)**2

        coulombpotential021 = 0.3898d0*exp(-0.12126d0*dd**2*c)+1.091859d0*exp(-1.9148*dd**2*c)-3.63491d0*exp(-0.60078*dd**2*c)


    end function coulombpotential02












    function coulombpotential0(location11,location2,d0) result(coulombpotential01)

        ! type(potential) :: this
        double precision, dimension(1:3) :: location11,location2
        double precision :: coulombpotential01,d0,dd

        dd=sqrt((location11(1)-location2(1))**2+(location11(2)-location2(2))**2+(location11(3)-location2(3))**2)
        if (dd>d0) then
            coulombpotential01=0.0d0
        else if (dd<=d0) then
            coulombpotential01 = -6/dd ! 14 for Si // 6 for C
            ! coulombpotential01 = -1/dd
        end if

    end function coulombpotential0

    function latticepotential(location111,aa,Ml0,d) result(latticepotential1)

        double precision,dimension(1:3) :: location111
        double precision :: latticepotential1,aa,d
        double precision,dimension(:,:),allocatable :: cxyz
        integer :: ii,jj,kk,ll,Nxyz,Ml0
        character(len = 8) :: dummyc


        open(10,file='C_lattice_xyz_unitcell.txt',status='old')
        read(10,*) Nxyz        
        allocate(cxyz(1:Nxyz,1:3))
        do ii=1,Nxyz
        read(10,*) dummyc,cxyz(ii,1),cxyz(ii,2),cxyz(ii,3)
        enddo
        close(10)

        cxyz=cxyz/0.529177249d0

        latticepotential1=0.0d0
        do ll=1,Nxyz
        do ii=-Ml0,Ml0
            do jj=-Ml0,Ml0
                do kk=-Ml0,Ml0
                    latticepotential1=latticepotential1&
                    +coulombpotential02(location111,(/cxyz(ll,1)+aa*ii,cxyz(ll,2)+aa*jj,cxyz(ll,3)+aa*kk/),d)
                                ! +coulombpotential0(location111,(/cxyz(ll,1)+aa*ii,cxyz(ll,2)+aa*jj,cxyz(ll,3)+aa*kk/),d)!&
                                !+coulombpotential00(location111,(/cxyz(ll,1)+aa*ii,cxyz(ll,2)+aa*jj,cxyz(ll,3)+aa*kk/),d)
                enddo
            enddo
        enddo
        enddo

    end function latticepotential


    function coulombpotential00(location11,location2,d0) result(coulombpotential001)

        ! type(potential) :: this
        double precision, dimension(1:3) :: location11,location2
        double precision :: coulombpotential001,d0,dd,pi

        pi=4.0d0*atan(1.0d0)

        dd=sqrt((location11(1)-location2(1))**2+(location11(2)-location2(2))**2+(location11(3)-location2(3))**2)
        if (dd>d0) then
            coulombpotential001=0.0d0
        else if (dd<=d0) then
        coulombpotential001=(1.0d0-(2.0d0*6.0d0**2.0d0*dd**2+2.0d0*6.0d0*dd+1.0d0)*exp(-2.0d0*6.0d0*dd))*1.0d0/dd*2.0d0&
                                    -(1.0d0+2.0d0*6.0d0*dd)*exp(-2.0d0*6.0d0*dd)*2.0d0
        end if

    end function coulombpotential00




    function latticepotential00(location111,aa,Ml0,d,lc) result(latticepotential1)

        double precision,dimension(1:3) :: location111
        double precision :: latticepotential1,aa,d,lc
        double precision,dimension(:,:),allocatable :: cxyz
        integer :: ii,jj,kk,ll,Nxyz,Ml0,nn
        character(len = 8) :: dummyc


        open(10,file='C_lattice_xyz_primitivecell.txt',status='old')
        read(10,*) Nxyz       
        allocate(cxyz(1:(Nxyz*63),1:3))
        do ii=1,Nxyz
        read(10,*) dummyc,cxyz(ii,1),cxyz(ii,2),cxyz(ii,3)
        enddo
        close(10)
        
        cxyz=cxyz/0.529177249d0

        nn=3
        do ii=-1,1,2
            do jj=-1,1,2
                cxyz(nn,1)=cxyz(1,1)+lc*ii
                cxyz(nn,2)=cxyz(1,2)+lc*jj
                cxyz(nn,3)=cxyz(1,3)
                nn=nn+1
                cxyz(nn,1)=cxyz(2,1)+lc*ii
                cxyz(nn,2)=cxyz(2,2)+lc*jj
                cxyz(nn,3)=cxyz(2,3)
                nn=nn+1
            enddo
        enddo

        do ii=-1,1,2
            do jj=-1,1,2
                cxyz(nn,1)=cxyz(1,1)+lc*ii
                cxyz(nn,3)=cxyz(1,3)+lc*jj
                cxyz(nn,2)=cxyz(1,2)
                nn=nn+1
                cxyz(nn,1)=cxyz(2,1)+lc*ii
                cxyz(nn,3)=cxyz(2,3)+lc*jj
                cxyz(nn,2)=cxyz(2,2)
                nn=nn+1
            enddo
        enddo

        do ii=-1,1,2
            do jj=-1,1,2
                cxyz(nn,2)=cxyz(1,2)+lc*ii
                cxyz(nn,3)=cxyz(1,3)+lc*jj
                cxyz(nn,1)=cxyz(1,1)
                nn=nn+1
                cxyz(nn,2)=cxyz(2,2)+lc*ii
                cxyz(nn,3)=cxyz(2,3)+lc*jj
                cxyz(nn,1)=cxyz(2,1)
                nn=nn+1
            enddo
        enddo


        loop1 : do ii=-2,2,2
            loop2 : do jj=-2,2,2
                    if ((ii==0).and.(jj==0)) cycle
                    do kk=-2,2,2
                    cxyz(nn,1)=cxyz(1,1)+lc*ii
                    cxyz(nn,2)=cxyz(1,2)+lc*jj
                    cxyz(nn,3)=cxyz(1,3)+lc*kk
                    nn=nn+1
                    cxyz(nn,1)=cxyz(2,1)+lc*ii
                    cxyz(nn,2)=cxyz(2,2)+lc*jj
                    cxyz(nn,3)=cxyz(2,3)+lc*kk
                    nn=nn+1
                    enddo
            enddo loop2
        enddo loop1


        loop3 : do ii=-2,2,1
                    if (ii==0) cycle
            loop4 : do jj=-2,2,1
                    if (abs(ii)==abs(jj)) cycle
                    if (jj==0) cycle
                    do kk=-1,1,2
                    cxyz(nn,1)=cxyz(1,1)+lc*ii
                    cxyz(nn,2)=cxyz(1,2)+lc*jj
                    cxyz(nn,3)=cxyz(1,3)+lc*kk
                    nn=nn+1
                    cxyz(nn,1)=cxyz(2,1)+lc*ii
                    cxyz(nn,2)=cxyz(2,2)+lc*jj
                    cxyz(nn,3)=cxyz(2,3)+lc*kk
                    nn=nn+1
                    enddo
            enddo loop4
        enddo loop3


        do ii=-1,1,2
            do jj=-1,1,2
                do kk=-2,2,4
                    cxyz(nn,1)=cxyz(1,1)+lc*ii
                    cxyz(nn,2)=cxyz(1,2)+lc*jj
                    cxyz(nn,3)=cxyz(1,3)+lc*kk
                    nn=nn+1
                    cxyz(nn,1)=cxyz(2,1)+lc*ii
                    cxyz(nn,2)=cxyz(2,2)+lc*jj
                    cxyz(nn,3)=cxyz(2,3)+lc*kk
                    nn=nn+1
                enddo
            enddo
        enddo


        cxyz(nn,1)=cxyz(1,1)
        cxyz(nn,2)=cxyz(1,2)
        cxyz(nn,3)=cxyz(1,3)+lc*(-2.0d0)
        nn=nn+1
        cxyz(nn,1)=cxyz(2,1)
        cxyz(nn,2)=cxyz(2,2)
        cxyz(nn,3)=cxyz(2,3)+lc*(-2.0d0)
        nn=nn+1
        cxyz(nn,1)=cxyz(1,1)
        cxyz(nn,2)=cxyz(1,2)
        cxyz(nn,3)=cxyz(1,3)+lc*(2.0d0)
        nn=nn+1
        cxyz(nn,1)=cxyz(2,1)
        cxyz(nn,2)=cxyz(2,2)
        cxyz(nn,3)=cxyz(2,3)+lc*(2.0d0)
        nn=nn+1

        

        latticepotential1=0.0d0
        do ll=1,Nxyz*63
        latticepotential1=latticepotential1&
                    +coulombpotential0(location111,(/cxyz(ll,1),cxyz(ll,2),cxyz(ll,3)/),d)!&
                    !+coulombpotential00(location111,(/cxyz(ll,1)+aa*ii,cxyz(ll,2)+aa*jj,cxyz(ll,3)+aa*kk/),d)
        enddo

    end function latticepotential00





end module potential
