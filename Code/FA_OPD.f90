       PROGRAM FA_OPD 
!!!===============================================================================
!!!    The FA_OPD.f90 is used to calculate the optical properties of fractal
!!!    aggregtes with small (compared to incident wavelegth) spheres, which is
!!!    achieved by interplation of a large database.
!!!
!!!    INPUT:  
!!!     (1). Shaple parameters:   DF: Fractal Dimension  
!!!                                   (KF is fixed to be 1.2 in current version)
!!!                               D_MON: Monomer Diameter
!!!     (2). Size/Size Distribution N_PSD:  1. Single-size
!!!                                         2. Size Distribution
!!!          Particle Size N_Unit:  1. Monomer Number
!!!                                 2. Equivalent Volume Diameter
!!!
!!!          N_PSD=1: Read one variable
!!!                   N_Unit=1: NN: Monomer number
!!!                   N_Unit=2: D_EFF: Equivalent volume diameter
!!!                             
!!!          N_PSD=2: Read two variables
!!!                   N_Unit=1: NN: Monomer number
!!!                             SIGMA: Standard Deviation - Lognormal Distribution
!!!                   N_Unit=2: CMD: Equivalent Volume Diameter
!!!                             SIGMA: Standard Deviation - Lognormal Distribution
!!! 
!!!     (3). Incident Wavelegth:  WL: Wavlength        
!!!     (4). Refractive Index:    RE: Real Part
!!!                               IM: Imaginary Part
!!!
!!!    INPUT DATABASE:
!!!     DF, RE, IM, XM, N, AREA, QEXT, QSCA, QABS, SSA, G, P11(0:1:180)
!!!   
!!!
!!!    OUTPUT:
!!!      (1). Extinction, Absorption Cross Section 
!!!      (2). Mass Extinction, Absorption Cross Section
!!!      (3). Single Scattering Albedo
!!!      (4). Asymmetry Factor
!!!      (5). Phase Function P11
!!!           with scattering angles from 0 to 180 degree in step of 1 degree  
!!!
!!!===============================================================================
!!!    Last Updated: 14 Mar, 2019 By Ms Shiwen Teng
!!!    Current Range of parameters:
!!!       Fractal Dimension:       [1.80:2.80]
!!!       Monomer Size Parameter:  [0.05:0.50]
!!!       Monomer Number:          [   1:3000]
!!!       Real Part of RI:         [1.05:2.00]
!!!       Imaginary Part of RI:    [0.00:1.00]
!!!    The database range will be kept extending for different applications
!!!===============================================================================
!!!    This Fractal_Aggregate_Optical_Property_Database is developed by:
!!!    Dr. Chao Liu
!!!    School of Atmospheric Physics, NUIST, Nanjing, China
!!!    E-mail: chao_liu@nuist.edu; liuchao1225@gmail.com
!!!===============================================================================
!!!    Reference:
!!!    C. Liu, X. Xu, Y. Yin, M. Schnaiter, Y. L. Yung, 2018: 
!!!    Black carbon aggregate: A database for optical properties,
!!!    submitted to J. Quant Spectrosc Radiat Transfer. 
!!!===============================================================================


       DOUBLE PRECISION  DF0(3),XM0(10),RE0(6),IM0(10),NN0(135),XA0(135)
       DOUBLE PRECISION  SCA(186,3,10,6,10,135),D_WEIGHT(135),AREA(10,135)
       DOUBLE PRECISION  NORM,ASY,DANG,SANG(181),CANG(181),TEMP(40)
       DOUBLE PRECISION  BAREA,BVOL,BEXT,BSCA,BABS,BSSA,BASY,BP(181),MEC,MAC
       DOUBLE PRECISION  WV,RE,IM,XM,XA,DF,YQ(186,5,135)
       DOUBLE PRECISION  ISPLINE
       DOUBLE PRECISION  D_MON,D_UP,D_DOWN,CMD,SIGMA,D_EFF
       DOUBLE PRECISION  AREA_EFF,VOL_EFF,WEIGHT
       INTEGER N_DF,N_RE,N_IM,N_XM,N_N,NN
       INTEGER I_DF,I_RE,I_IM,I_XM,I_N 
       INTEGER N_PSD,N_Unit

!!!    READ INPUT FILE
       OPEN(UNIT=1,FILE='input.txt',STATUS='OLD')
       READ(1,*) !SKIP FIRST LINE
       READ(1,*) DF
       READ(1,*) D_MON
       READ(1,*) N_PSD
       READ(1,*) N_Unit
       IF ((N_PSD .EQ. 1) .AND. (N_Unit .EQ. 1)) THEN
         READ(1,*) NN
         D_EFF = DBLE(NN)**(1.D0/3.D0) * D_MON
       ELSEIF ((N_PSD .EQ. 1) .AND. (N_Unit .EQ. 2)) THEN
         READ(1,*) D_EFF
       ELSEIF ((N_PSD .EQ. 2) .AND. (N_Unit .EQ. 1)) THEN
         READ(1,*) NN, SIGMA
         CMD = DBLE(NN)**(1.D0/3.D0) * D_MON
       ELSEIF ((N_PSD .EQ. 2) .AND. (N_Unit .EQ. 2)) THEN
         READ(1,*) CMD, SIGMA
       ENDIF
       READ(1,*) WV
       READ(1,*) RE, IM  
       CLOSE(1)

       XM = 2.0D0 * DACOS(0.0D0) * D_MON / WV
       D_UP   = DBLE(3000)**(1.D0/3.D0) * D_MON
       D_DOWN = DBLE(1)**(1.D0/3.D0) * D_MON

!!!    CHECK INPUT
       IF ((DF .LT. 1.8) .OR. (DF .GT. 2.8)) THEN
         WRITE(*,*) 'Error! Input out of range!'
         WRITE(*,*) 'DF must be larger than 1.8 and smaller than 2.8!'
         STOP
       ENDIF
       IF ((XM .LT. 0.05) .OR. (XM .GT. 0.5)) THEN
         WRITE(*,*) 'Error! Input out of range!'
         WRITE(*,*) 'Pi*D_MON/WV must be larger than 0.05 and smaller than 0.5!'
         STOP
       ENDIF
       IF ((RE .LT. 1.05) .OR. (RE .GT. 2.0)) THEN
         WRITE(*,*) 'Error! Input out of range!'
         WRITE(*,*)  'Real Part of Refractive Index must be larger than 1.05 and smaller than 2.0!'
         STOP
       ENDIF
       IF ((IM .LT. 0.0) .OR. (IM .GT. 1.0)) THEN
         WRITE(*,*) 'Error! Input out of range!'
         WRITE(*,*) 'Imaginary Part of Refractive Index must be larger than 0 and smaller than 1.0!'
         STOP
       ENDIF
       IF (N_PSD .EQ. 1) THEN
         IF ((D_EFF .LT. D_DOWN) .OR. (D_EFF .GT. D_UP)) THEN
           WRITE(*,*) 'Error! Particle size out of range!'
           STOP
         ENDIF
       ELSEIF (N_PSD .EQ. 2) THEN
         IF ((CMD .LT. D_DOWN) .OR. (CMD .GT. D_UP)) THEN
           WRITE(*,*) 'Error! Particle size out of range!'
           STOP
         ENDIF
       ENDIF 


!!!    READ THE DATABASE
       I_DF = 3
       I_IM = 10
       I_RE = 6
       I_XM = 10
       I_N  = 135
       DENSITY = 1.8D-12   !! in unit of gram per micron^3: g/um^3

       DF0 = (/1.8D0, 2.3D0, 2.8D0/)
       XM0 = (/0.05D0, 0.10D0, 0.15D0, 0.20D0, 0.25D0, &
               0.30D0, 0.35D0, 0.40D0, 0.45D0, 0.50D0/)
       RE0 = (/1.05D0, 1.2D0, 1.4D0, 1.6D0, 1.8D0, 2.0D0/)
       IM0 = (/0.00D0, 0.0001D0, 0.001D0, 0.01D0, 0.1D0, &
               0.20D0,    0.4D0,   0.6D0,  0.8D0, 1.0D0/)

       DANG = DACOS(0.D0)/90.D0
       DO I=1,181
          SANG(I) = DSIN(DBLE(I-1)*DANG)
          CANG(I) = DCOS(DBLE(I-1)*DANG)
       ENDDO

       OPEN(UNIT=1,FILE='FA_1p8.DAT',STATUS='OLD')
       OPEN(UNIT=2,FILE='FA_2p3.DAT',STATUS='OLD')
       OPEN(UNIT=3,FILE='FA_2p8.DAT',STATUS='OLD')

       DO N_DF=1,I_DF
       DO N_RE=1,I_RE
       DO N_IM=1,I_IM
          DO N_XM=1,I_XM
          DO N_N=1,I_N
             READ(N_DF,*) TEMP(1:4),NN0(N_N),AREA(N_XM,N_N),SCA(1:186,N_DF,N_XM,N_RE,N_IM,N_N)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO

       CLOSE(1)
       CLOSE(2)
       CLOSE(3)

       D_WEIGHT(1) = DBLE(NN0(2))**(1.D0/3.D0)*D_MON-DBLE(NN0(1))**(1.D0/3.D0)*D_MON
       DO N_N = 2, I_N-1
         D_WEIGHT(N_N)=DBLE(NN0(N_N+1))**(1.D0/3.D0)*D_MON-DBLE(NN0(N_N-1))**(1.D0/3.D0)*D_MON 
       ENDDO
       D_WEIGHT(I_N)=DBLE(NN0(I_N))**(1.D0/3.D0)*D_MON-DBLE(NN0(I_N-1))**(1.D0/3.D0)*D_MON 


!!!    BEGINNING CALCULATIONS
       OPEN(UNIT=1,FILE='output.txt',STATUS='UNKNOWN')

       DO N_N = 1,I_N
         DO N_DF = 1,I_DF
           DO N_XM = 1,I_XM
             DO N_IM = 1,I_IM

!!!             INTERPOLATE FOR RE
                DO I=1,186
                   YQ(I,1,N_IM) = ISPLINE(RE,6,RE0,SCA(I,N_DF,N_XM,1:I_RE,N_IM,N_N))
                ENDDO
             ENDDO

!!!          INTERPOLATE FOR IM
             DO I=1,186
                YQ(I,2,N_XM) = ISPLINE(IM,10,IM0,YQ(I,1,1:I_IM))
             ENDDO
           ENDDO

!!!        INTERPOLATE FOR XM
           DO I=1,186
              YQ(I,3,N_DF) = ISPLINE(XM,10,XM0,YQ(I,2,1:I_XM))
           ENDDO
         ENDDO

!!!      INTERPOLATE FOR DF 
         DO I=1,186
            YQ(I,4,N_N) = ISPLINE(DF,3,DF0,YQ(I,3,1:I_DF))
         ENDDO
       ENDDO

       BNORM=0.0D0
       BAREA=0.0D0
       BVOL=0.0D0
       BEXT=0.0D0
       BSCA=0.0D0
       BABS=0.0D0
       BSSA=0.0D0
       BASY=0.0D0
       BP(1:181)=0.0D0


       IF (N_PSD .EQ. 1) THEN
         XA = 2.0D0 * DACOS(0.0D0) * D_EFF / WV
         DO N_N = 1, I_N
            XA0(N_N) = XM * DBLE(NN0(N_N))**(1.D0/3.D0)
         ENDDO
         AREA_EFF  = DACOS(-1.D0)*(D_EFF*D_EFF*0.25D0)
         VOL_EFF  = DACOS(-1.D0)*(D_EFF*D_EFF*D_EFF)/6.D0

!!!      INTERPOLATE FOR XA
         DO I=1,186
            YQ(I,5,1) = ISPLINE(XA,135,XA0,YQ(I,4,1:I_N))
         ENDDO

         BEXT = YQ(1,5,1) * AREA_EFF
         BSCA = YQ(2,5,1) * AREA_EFF
         BABS = YQ(3,5,1) * AREA_EFF

         MEC  = BEXT / VOL_EFF / DENSITY * 1.0D-12
         MAC  = BABS / VOL_EFF / DENSITY * 1.0D-12
         BSSA = BSCA / BEXT
      
         NORM = 0.0D0
         ASY  = 0.0D0
         DO I = 6,186
            NORM = NORM + YQ(I,5,1)*SANG(I-5)*0.5D0*DANG
             ASY = ASY  + YQ(I,5,1)*CANG(I-5)*SANG(I-5)*0.5D0*DANG
         ENDDO
         BASY = ASY/NORM

         BP(1:181) = YQ(6:186,5,1)/NORM


       ELSEIF (N_PSD .EQ. 2) THEN
         CMD   = DLOG(CMD)
         SIGMA = DLOG(SIGMA)
       
!!!      BULK PROPERTIES FOR GIVEN DISTRIBUTION
         DO N_N = 1, I_N
            D_EFF = DBLE(NN0(N_N))**(1.D0/3.D0) * D_MON
            AREA_EFF = DACOS(-1.D0)*(D_EFF*D_EFF*0.25D0) 
            VOL_EFF =  DACOS(-1.D0)*(D_EFF*D_EFF*D_EFF)/6.D0

            NORM = 0.0D0
            ASY  = 0.0D0
            DO I = 6,186
               NORM = NORM + YQ(I,4,N_N)*SANG(I-5)*0.5D0*DANG
                ASY = ASY  + YQ(I,4,N_N)*CANG(I-5)*SANG(I-5)*0.5D0*DANG
            ENDDO
            YQ(5,4,N_N) = ASY/NORM
            YQ(6:186,4,N_N) = YQ(6:186,4,N_N)/NORM

            CALL LOGNORMAL(CMD,SIGMA,D_EFF,WEIGHT)
            WEIGHT = WEIGHT * D_WEIGHT(N_N)
          
            BNORM = BNORM + WEIGHT 
            BAREA = BAREA + WEIGHT * AREA_EFF
            BVOL  = BVOL  + WEIGHT * VOL_EFF
            BEXT  = BEXT  + WEIGHT * AREA_EFF * YQ(1,4,N_N)
            BSCA  = BSCA  + WEIGHT * AREA_EFF * YQ(2,4,N_N)
            BABS  = BABS  + WEIGHT * AREA_EFF * YQ(3,4,N_N)
            BASY  = BASY  + WEIGHT * AREA_EFF * YQ(2,4,N_N) * YQ(5,4,N_N)
            BP(1:181) = BP(1:181)  + WEIGHT * AREA_EFF * YQ(2,4,N_N) * YQ(6:186,4,N_N)
          ENDDO

          MEC   = BEXT / BVOL / DENSITY * 1.0D-12
          MAC   = BABS / BVOL / DENSITY * 1.0D-12
          BSSA  = BSCA / BEXT
          BASY  = BASY / BSCA
          BP(1:181) = BP(1:181) / BSCA

        ENDIF

!!!     OUTPUT
        WRITE(1,'(A)') 'Scattering Properties Based on Fractal_Aggregate_Optical_Property_Database'
        WRITE(1,'(A)') 'Developed by Dr. Chao Liu, School of Atmospheric Physics, NUIST, 2019'             
        WRITE(1,'(A)') '' 
        WRITE(1,'(A)') 'Input:'  
        WRITE(1,'(A,F5.1)') 'Fractal Dimension:', DF
        WRITE(1,'(A,F6.2,A)') 'Monomer Diameter:', D_MON, ' um'
        IF (N_PSD .EQ. 1) WRITE(1,'(A,F6.2)') 'Equivalent Volume Diameter:',D_EFF
        IF (N_PSD .EQ. 2) WRITE(1,'(A,2F6.2)') 'CMD and SIGMA for Lognormal Size Distribution:', DEXP(CMD), DEXP(SIGMA)
        WRITE(1,'(A,F6.2,A)') 'Wavelength:',WV,' um'
        WRITE(1,'(A,F6.2,A,F6.2,A)') 'Refractive Index:',RE,'+',IM,'i'
        WRITE(1,'(A)') ''
        WRITE(1,'(A)') 'Output:'
        WRITE(1,'(A,E20.6,A)') 'Extinction Cross Section:',BEXT,' um^2'
        WRITE(1,'(A,E15.6,A)') 'Mass Extinction Cross Section:',MEC,' m^2/g'
        WRITE(1,'(A,E20.6,A)') 'Absorption Cross Section:',BABS,' um^2'
        WRITE(1,'(A,E15.6,A)') 'Mass Absorption Cross Section:',MAC,' m^2/g'
        WRITE(1,'(A,E20.6)') 'Single Scattering Albedo:',BSSA
        WRITE(1,'(A,E28.6)') 'Asymmetry Factor:', BASY
        WRITE(1,'(A)') 'Normalized Phase Function:'
        DO I = 1,181
           WRITE(1,'(I5,E20.8)') (I-1), BP(I)
        ENDDO
        CLOSE(1)


        END PROGRAM FA_OPD



        Function ispline(u, n, x, y)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
        implicit none
        double precision ispline
        integer n
        double precision  u, x(n), y(n), b(n), c(n), d(n)
        integer i, j, k
        double precision dx

        call spline(x,y,b,c,d,n)

! if u is ouside the x() interval take a boundary value (left or right)
        if(u <= x(1)) then
          ispline = y(1)
          return
        end if
        if(u >= x(n)) then
          ispline = y(n)
          return
        end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(u < x(k)) then
            j=k
          else
            i=k
          end if
        end do
!*
!  evaluate spline interpolation
!*
        dx = u - x(i)
        ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
        end function ispline



        subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
        implicit none
        integer n
        double precision x(n), y(n), b(n), c(n), d(n)
        integer i, j, gap
        double precision h

        gap = n-1
! check input
        if ( n < 2 ) return
        if ( n < 3 ) then
           b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
           c(1) = 0.
           d(1) = 0.
           b(2) = b(1)
           c(2) = 0.
           d(2) = 0.
           return
        end if
!
! step 1: preparation
!
        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, gap
           d(i) = x(i+1) - x(i)
           b(i) = 2.0*(d(i-1) + d(i))
           c(i+1) = (y(i+1) - y(i))/d(i)
           c(i) = c(i+1) - c(i)
        end do
!
! step 2: end conditions 
!
        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = 0.0
        c(n) = 0.0
        if(n /= 3) then
          c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
          c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
          c(1) = c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
        end if
!
! step 3: forward elimination 
!
        do i = 2, n
           h = d(i-1)/b(i-1)
           b(i) = b(i) - h*d(i-1)
           c(i) = c(i) - h*c(i-1)
        end do
!
! step 4: back substitution
!
        c(n) = c(n)/b(n)
        do j = 1, gap
           i = n-j
           c(i) = (c(i) - d(i)*c(i+1))/b(i)
        end do
!
! step 5: compute spline coefficients
!
        b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
        do i = 1, gap
           b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
           d(i) = (c(i+1) - c(i))/d(i)
           c(i) = 3.*c(i)
        end do
        c(n) = 3.0*c(n)
        d(n) = d(n-1)
        end subroutine spline

        SUBROUTINE LOGNORMAL(A,B,X,Y)
        DOUBLE PRECISION  A,B,X,Y
!     A: mu      LOG(D)
!     B: sigma   LOG(SIGMA)
!     Y = dN/dD  
!     the LOG in the equations can be changed as LOG10, where all variables should be changed
!     It include A, B and LOG(X) in the equation.
        DOUBLE PRECISION PI
        PI = DACOS(-1.0D0)
        Y  = DEXP(-0.50D0*(DLOG(X)-A)**2.0D0/B/B)/X/B/DSQRT(2.0D0*PI)
        RETURN
        END SUBROUTINE LOGNORMAL
