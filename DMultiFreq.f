      program DMultiFreq

! As this code seems to be all working now, I'll call it 1.0
! Making this all work as expected requires a couple other pieces
! Bash: NewMultiShell, FortranRunner, DataCollector
! 	FortranRunner2 also exists, with 40 phi points instead of 20
! Text: P, Q, MultiInStart
!	From there, NewMultiShell will end up making MultiIn0
! Python: PQ_Plotter
!	PQ_Plotter works for both MultiFrequency and Gating
c****************************************************************
	  IMPLICIT REAL*8 (A-H,O-Z)
	  IMPLICIT INTEGER*4 (I-N)
          REAL*8 :: mdrive1,mdrive2
	  dimension par(20)

c-------------------------------------------------MPI stuff START
      REAL*8, ALLOCATABLE :: SAMPLE(:),SAMPLEB1(:),SAMPLEB2(:)
      CHARACTER*8, ALLOCATABLE :: SAMPNAME(:) 

      include 'mpif.h'
      INTEGER,  ALLOCATABLE :: SEED(:),OLDSEED(:),NSEED(:)
      double precision starttime,endtime
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid,ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr)
!-------------------------------------------------MPI stuff END
! Notes for the curious, 7/28/19
! This code is in full working condition but is a bit ugly
! Many things aren't labeled well, comments aren't great
! A much better example is DGating or DLattice
! Eventually this will be rewritten to follow their style
!****************************************************************

      input=33
	  open(unit=nf_in,file='input.dat',status='old')
c	num_atom
      read(nf_in,*) num_atom
c   gammaP
      read(nf_in,*) gamma_in
c   ntstep
      read(nf_in,*) ntstep
c   phi
      read(nf_in,*) phi
c   U_0
      read(nf_in,*) U_in
c   dt0
      read(nf_in,*) dt0
c   idum
      read(nf_in,*) idum
c   jdum
      read(nf_in,*) jdum
c   F_0
      read(nf_in,*) F0_in

      read(nf_in,*) pqratio
	  
	  
      close(nf_in)
c-------------------------------------------------
c     INITIALLIZE RANDOM NUMBER GENERATOR:
c      CALL SYSTEM_CLOCK(COUNT=idum)
      CALL RANDOM_SEED(SIZE=n)
      ALLOCATE(SEED(n),NSEED(n))
c     Set maximum seeds according to the intel manual:
      seedmax1=2147483562
      seedmax2=2147483398
      nseed(1)=nint(real(seedmax1)/2.0/numprocs)
      nseed(2)=nint(real(seedmax2)/2.0/numprocs)
      DO i=1,n
         SEED(i)=idum+nseed(i)*myid
      enddo
      CALL RANDOM_SEED(PUT=SEED)
      DEALLOCATE(SEED)
c-------------------------------------------------
C     SAMPLE VARIABLES: (BLOCK DATA AVERAGE)

c     Setting NBLOCK to number of processors
      nblock=numprocs

      NSAMPLE=2
      ALLOCATE(SAMPLE(NSAMPLE),SAMPLEB1(NSAMPLE))
      ALLOCATE(SAMPLEB2(NSAMPLE),SAMPNAME(NSAMPLE))
      SAMPNAME(1)='<v>    :'
      SAMPNAME(2)='Diff   :'

      DO I=1,NSAMPLE
         SAMPLEB1(I)=0.D0
         SAMPLEB2(I)=0.D0
      ENDDO

c****************************************************************
!      nf=51 222  format(a9,20(f20.13))
!      open(unit=nf,file='p_list.dat',status='unknown')
      nstate = 1
	  state=dble(nstate)
      oneninth = (1.d0/9.d0)
      deltaP = 0.d0
	  
      dt=dt0/gamma_in

      tsample=dt*ntstep
      
      F0=F0_in/2.d0
      gammaP = gamma_in/2.d0

      U_overER = U_in
      U_0 = U_in / 2.d0
      


      w_1 = dsqrt(U_overER) * .74 ! As per figure 4.9
	  pi= 3.141592653589793d0
	  !period=2.d0*pi/w_d
      w_2 = (pqratio) * w_1
  
      A = 1.d0
      B = .3d0

      
!     Write parameters to output file
      if(myid.eq.0)then
         nf_ou=10
         open(unit=nf_ou,file='output.dat',status='unknown')
         write(nf_ou,*) 'dt = ',dt
         write(nf_ou,*)'U0/Er=',U_overER
         write(nf_ou,*)'tsample=',tsample 
         write(nf_ou,*)'period=',period 	
         write(nf_ou,*) 'phi = ',phi
         write(nf_ou,*) 'F_0 = ',F0_in
         write(nf_ou,*) 'nblock=',nblock
	 write(nf_ou,*) 'pqratio=',pqratio
         close(nf_ou)
      endif


!      CALL CPU_TIME ( time_begin)
      starttime=MPI_WTIME()
	  
c****************************************************************
c     Transition time for adiabatic turn on

      ramptime = 10.d0

!		 rampup=(1.d0 - exp(-t/ramptime))

c****************************************************************


      xmean=0.d0
      xmean2=0.d0
      np=0
      vmean=0.d0

c     Start the party

      do i=1,num_atom

! Set the initial conditions, momentum and position random
! Mostly just to make sure not every single atom is starting the exact same way
        z0=gasdev(idum)
        p0=gasdev(idum)
        t0=0.d0
! Now actually say these are the variables of the atom
        p = p0
        z = z0
	t = t0
	
        ! This is the loop for running through the entire time sequence
        do j=1,ntstep

! Use a variable in place os dcos(z*2) every time
         cosz2=dcos(2.d0*z)
! This is diffusion from the spontaneous emission
         Dpp = (gammaP/90.d0)*(35.d0 + state*7.d0*cosz2)
! Actual kick of momentum from emission
         kick = ((2.d0*Dpp*dt)**(0.5d0))*gasdev(idum)
! Split the multi-frequency drive into multiple parts because it's really long
	 sinw1 = dsin(w_1 * t)
	 sin2w1= dsin(2.d0 * w_1 * t)
	 cosw1 = dcos(w_1 * t)
	 cos2w1= dcos(2.d0 * w_1 * t)
	 sinw2 = dsin((w_2 * t) + phi)
	 cosw2 = dcos((w_2 * t) + phi)
         mdrive1 = F0 * w_1 * sinw2 *((A*cosw1)+(2.d0 * B * cos2w1))
	 mdrive2 = F0 * w_2 * cosw2 * ((A * sinw1) + (B * sin2w1))
! Combine the two parts back together
	 mdrive = mdrive1 + mdrive2
! Force is a negative, equation pulled directly from Brown thesis
! deltaP is actual change in the momentum of the atom
         deltaP = (U_0*sin(2.d0*z)*dt)*state + kick - mdrive*dt


	! We need to ask if the atom changes wells
	 call RANDOM_NUMBER(rand)

         jumprate = (oneninth)*dt*gammaP*(1.d0 + state*cosz2)
	! Now compare random number to probability of changing
         if (rand.lt.jumprate)then !you change wells
!		switch the state
	    nstate = - nstate
		! Make sure that state is a double for use in the above equations
	    state=dble(nstate)
!           print*, 'switch'
	 endif

	! Now update the variables of the atom
         p = p + deltaP
         z = z + (p*dt)
	 t = t + dt 
	
	
	
        end do

c       Drift and diffusion calculation
! np is only keeping index on things
        np=np+1
        xd = z - z0
        timed = t - t0
        vmean=vmean+xd/timed
        xmean=xmean+xd
        xmean2=xmean2+xd*xd    
	  
      end do
	  
c     -------------------------------------------------------
c     SAMPLE VARIABLES:

      vmean=vmean/dble(np) 
      xmean=xmean/dble(np)
      diff=(xmean2/dble(np)-xmean*xmean)/2.d0/timed
      SAMPLE(1)=vmean
      SAMPLE(2)=diff	  

C     STORE SAMPLE VARIABLES:
c     Gather data from all processes
      ndim=nsample
      call MPI_REDUCE(sample,sampleb1,ndim,MPI_DOUBLE_PRECISION,
     &       MPI_SUM,0,MPI_COMM_WORLD,ierr)
      do k=1,nsample
         sample(k)=sample(k)*sample(k)
      enddo
      call MPI_REDUCE(sample,sampleb2,ndim,MPI_DOUBLE_PRECISION,
     &       MPI_SUM,0,MPI_COMM_WORLD,ierr)
      
      if(myid.eq.0)then

       open(unit=nf_ou,file='output.dat',status='old',position='append')

c     Drift and effective diffussion

         write(nf_ou,*)'Block data: '

C     SAMPLE VARIABLES:
         DO K=1,NSAMPLE
            SAMPLEB1(K)=SAMPLEB1(K)/DBLE(NBLOCK)
            SAMPLEB2(K)=(SAMPLEB2(K)-DBLE(NBLOCK)*SAMPLEB1(K)
     &           *SAMPLEB1(K))/DBLE(NBLOCK-1)
            write(nf_ou,222)SAMPNAME(K),SAMPLEB1(K)
            write(nf_ou,222)'    (dis)',DSQRT(SAMPLEB2(K))
         ENDDO

c      Peclet: 
       pi= 3.141592653589793d0
       pe=sampleb1(1)*pi/sampleb1(2)
       write(nf_ou,222)'Peclet : ',pe

c      CALL CPU_TIME ( time_end)
         endtime=MPI_WTIME()
         write(nf_ou,*)
     &        'total time spent in root ID:'
     &        ,endtime-starttime,' seconds ('
     &        ,(endtime-starttime)/60.,' minutes)'
c      write(nf_ou,*)
c      write(nf_ou,*)'total CPU time spent:'
c     &     ,time_end-time_begin,' seconds'



         close(nf_ou)
      endif

 222  format(a9,20(f20.13))
      call MPI_FINALIZE(ierr)
      end



      
****************************************************************
* Algorithms to generate random numbers(Ref. W.H.Press et Al.  *
* "NUMERICAL RECIPES".Cambridge U.P. 1986)                     *
* Function GASDEV generates a Gaussian Distribution            *
* of zero mean and variance unity.(Box-Mueller formula)        *
****************************************************************
      FUNCTION GASDEV(IDUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      REAL*8 GASDEV,X(2)
      SAVE ISET,GSET
      DATA ISET/0/
      IF(ISET.EQ.0) THEN
 1       CALL RANDOM_NUMBER(X)
         V1=2.0d0*X(1)-1.d0
         V2=2.0d0*X(2)-1.d0
         R=V1**2+V2**2
      IF(R.GE.1.0d0) GO TO 1
         FAC=DSQRT(-2.0d0*DLOG(R)/R)
         GSET=V1*FAC
         GASDEV=V2*FAC
         ISET=1
      ELSE
         GASDEV=GSET
         ISET=0
      ENDIF
      RETURN
      END
      
	  
      
         
	  
		
