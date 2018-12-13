c23456789012345678901234567890123456789012345678901234567890123456789012
c **********************************************************************
c (C) Copyright 2009-2019 by Ralph C Aldredge, PhD, PE                 
c  All Rights Reserved.                           		       
c Licensed under the GNU Affero General Public License v3.0
c **********************************************************************

c SLAPCONS

c*****************************************
c NEED TO UPDATE THE departure subroutine!
c*****************************************

c This program computes the evolution of advected propagating surfaces.
c It uses a fully three-dimensional semi-Lagrangian scheme based on 
c transient interpolation modelling (TIM) with the upwind-biased 
c twenty-point grid-point stencil of Leonard, MacVean and Lock (1993), 
c and is 2nd order in time and 3rd order in space; as described in 
c R.C.  Aldredge, “Semi‐Lagrangian advection–propagation (SLAP) scheme
c   for three-dimensional interface tracking,” Journal of Computational 
c   Physics, vol. 229, pp. 4683-4702 (2010)

c Correction still needed: dt used to compute cx_ave should not be the same
c every time step.

        program slap

        implicit double precision (a-h,o-z)
        character*160 fileplt
        character*80 str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
        character*80 str11,str12,str13,str14,str15,str16,str17,str18
        character*80 str19,str20,str21,str22,str23
	dimension matcon1(256,16),matcon2(12,2)

        allocatable :: igmin(:,:),igmax(:,:)
        allocatable :: x(:),y(:),z(:)
        allocatable :: uprmvar(:)
        allocatable :: mask(:,:,:),indx(:,:)
c	allocatable :: u(:,:,:)
c	allocatable :: v(:,:,:)
c    	allocatable :: w(:,:,:)
        allocatable :: cx(:,:,:)
     	allocatable :: cy(:,:,:)
     	allocatable :: cz(:,:,:)
        allocatable :: cx_old(:,:,:)
     	allocatable :: cy_old(:,:,:)
     	allocatable :: cz_old(:,:,:)
        allocatable :: cx_old2(:,:,:)
     	allocatable :: cy_old2(:,:,:)
     	allocatable :: cz_old2(:,:,:)
        allocatable :: cx_ave(:,:,:)
     	allocatable :: cy_ave(:,:,:)
     	allocatable :: cz_ave(:,:,:)
        allocatable :: cx_tmp(:,:,:)
     	allocatable :: cy_tmp(:,:,:)
     	allocatable :: cz_tmp(:,:,:)
        allocatable :: dep(:,:),dephalf(:,:)
	allocatable :: cfn(:,:,:)
	allocatable :: g(:,:,:)
	allocatable :: gradg(:,:,:)
	allocatable :: gradgold(:,:,:)
 	allocatable :: gradgold2(:,:,:)
	allocatable :: gradg_ave(:,:,:)
	allocatable :: gtransf(:,:,:)
	allocatable :: xtransf(:)

c23456789012345678901234567890123456789012345678901234567890123456789012

        common /rparam/small,pi
        common /iparam/iquad
        common /isodat/matcon1,matcon2
        common /surfave/grdave0,grdave1
	common /tubedat/twa,twb,gjump,nomask
	common /global/mstep,fileplt
	common /reindat/irk,igrd,cfl_re
	common /veldat/umax,sl,umod,tke,kr,nmodes,iso,ifreq,cfl
	common /grid/icush,imin,imax,izone
	common /kechk/uprime_rel
	
	open (unit=4,file='slapdat.txt',status='old')
 	read(4,*) cfl	     ! CFL number
 	read(4,*) sl	     ! normal propagation speed
        read(4,*) npltkp     ! number of data output events
        read(4,*) ireplt     ! if 1, replotting to gdat.txt occurs
                             !   otherwise, to gdat1.txt, gdat2.txt, ...
        read(4,*) limit      ! 1 => Shape-Preservation Algorithm (SPA) used
			     ! SHOULD BE ZERO (otherwise, get spurious results
			     ! for large-amplitude velocity fluctuations)
        read(4,*) itc	     ! number of iterations for determining departure point
        read(4,*) kr         ! ratio (kx/ky=kx/kz) of transverse to
                             !   longitudinal length scales
        read(4,*) iso        ! 1 => isotropic advection field when kr = 1
                             ! otherwise (when kr.ne.1) iso=1 => equal advection 
                             ! amplitudes but unequal long and trans length scales.
                             ! Note that results for (i.eq.1) and (i.ne.1) are the 
                             ! same when k=2
        read(4,*) iquad      ! 0 => computation on the full-period domain 
                             !   (y,z = [-0.5,0.5]); otherwise on a quadrant of the 
			     !   full domain (y,z = [-0.25,0.25])
			     ! Note that iquad is forced to 0 if nmodes > 1
        read(4,*) pltstrt    ! beginning of time2pi range for plotting 
                             ! (requires ireplt=0)
        read(4,*) pltstop    ! end of time2pi range for plotting
                             ! (requires ireplt=0)
        read(4,*) irein	     ! number of time steps between reinitialization application 
			     ! (should be >= 2; get spurious results at low-res otherwise))
        read(4,*) nordg	     ! order of spatial differentiation in wenogradts
			     ! subroutine (where gradient magnitude is calculated
			     ! for reinitialization subroutines)
        read(4,*) nreinit    ! number of reinitialization iterations
			     ! (nreinit=0) => don't reinitialize
        read(4,*) iscf	     ! 1 => reinitialize with 'sub-cell fix'
        read(4,*) irk	     ! 0 (1) => 1st order (RK) reinitialization
        read(4,*) igrd	     ! 0 (~0) => ENO (WENO) gradients
        read(4,*) wind	     ! mean velocity field in x direction (units of Sl0)
        read(4,*) ifreq	     ! 0 => coordinate transformation to eliminate time
			     ! dependence of velocity field
        read(4,*) cfl_re     ! cfl for reinitialization
        read(4,*) nt0	     ! no. of grid points per vel-field length scale [-0.5,0.5]
			     ! along each coord axis (must be even); if iquad = 1
			     ! then comp on only 1+nt0/2 grid points with same resolution
        read(4,*) itimelev   ! number of time levels used for time averaging
			     ! (should be 2 or 3)
        read(4,*) nmodes     ! number of wavenumber modes in velocity field
	read(4,*) norda	     ! order of spatial differentiation for gradient in advection 
			     ! scheme: 4 or 5 => WENO; otherwise 4rd-order upwind (faster)
			     ! value of 0 recommended
        read(4,*) lrunge     ! 0 (1) => vel-ave collocation (runge-kutta integration) to get 
			     ! characteristic departure point; the code is slower and no more 
			     ! accurate with lrunge=1 than with lrunge=0 (disappointing); static 
			     ! departure-point calculation (only at t=0) not yet implemented for 
			     ! lrunge=1

	close (4)

c23456789012345678901234567890123456789012345678901234567890123456789012

c The even modes present when nmodes > 1 (i.e., k=2,4,...) break symmetry
c and thus require computation on the full domain (iquad=0)
        if (nmodes.gt.1) then
	   iquad=0
	elseif (iquad.ne.0) then
	   iquad=1
	endif
        if ((itimelev.ne.2).and.(itimelev.ne.3)) itimelev=2
c       if ((itimelev.eq.2).and.(irein.eq.1)) irein=2

	icush=3; itwa=3; itwb=6
c	itwb=nt0/2-icush

 	dx=1.0d0/dfloat(nt0); dy=dx; dz=dx
 	twa=dfloat(itwa)*dx
 	twb=dfloat(itwb)*dx
 	gjump=2.0d0*twb

        if ((norda.eq.4).or.(norda.eq.5)) limit=0
	lflag=0

	nx=2*(icush+itwb); ny=nt0/(1+iquad)+iquad; nz=ny
	nxp=nx+1; nyp=ny+1; nzp=nz+1
        nptot=nx*ny*nz

        allocate (igmin(ny,nz),igmax(ny,nz))
        allocate (x(nx+1),y(ny+1),z(nz+1))
        allocate (uprmvar(50))
	allocate (g(-2:nx+3,-2:ny+3,-2:nz+3))

c23456789012345678901234567890123456789012345678901234567890123456789012

c Initialization
c-------------------------------------------------------------------------
        call isoinit(matcon1,matcon2)

        small=1.0d-16
        pi=dacos(-1.0d0)
        ireplt0=ireplt
        sl0=sl
	wind0=wind

        if ((kr.ne.1).and.(kr.ne.2)) then
           write (6,*) 'ERROR: kr must be either 1 or 2; kr = ',kr
           stop
        endif

c23456789012345678901234567890123456789012345678901234567890123456789012

	time=0.0d0
        ixave=0
        igauss=0
        maxstp=1000000
	ip=nx/5+1
        yc=0.15d0
        zc=0.15d0
        ibcv=0     ! ibcv=0 => periodic; constant slope otherwise
                   ! ibcv=4 => periodic in x, constant slope in y & z
        ibcg=2     ! ibcg=(0,2) => periodic, jump-periodic in x 
                   !    periodic in y & z
                   ! ibcg=3 => jump-periodic in x, constant slope in y & z
                   ! ibcg.ne.(0,2,3) => constant slope in x, y & z

c23456789012345678901234567890123456789012345678901234567890123456789012

	do 100 i=1,nxp
  	   x(i)=0.5d0*dfloat(2*(i-1)-nx)*dx
	   if (dabs(x(i)).lt.small) x(i)=0.0d0
  100	continue
        do 150 j=1,nyp
c 	   y(j)=dfloat(j-1-ny/2)*dy
c 	   y(j)=-0.5d0+dfloat(j-1)*dy
  	   y(j)=-0.5d0/dfloat(1+iquad)+dfloat(j-1)*dy
	   if (dabs(y(j)).lt.small) y(j)=0.0d0
  150	continue
        do 175 k=1,nzp
c 	   z(k)=dfloat(k-1-nz/2)*dz
c 	   z(k)=-0.5d0+dfloat(k-1)*dz
  	   z(k)=-0.5d0/dfloat(1+iquad)+dfloat(k-1)*dz
	   if (dabs(z(k)).lt.small) z(k)=0.0d0
  175	continue

        do 200 k=1,nz
        do 200 j=1,ny
        do 200 i=1,nx
c Linear (variation with x only)
           if (dabs(x(i)).le.twb) then
  	      g(i,j,k)=x(i)
 	   else
  	      g(i,j,k)=dsign(twb,x(i))
 	   endif
	   if (dabs(g(i,j,k)).lt.small) g(i,j,k)=0.0d0
  200	continue
  	call bound(nx,ny,nz,g,ibcg)

 	
c23456789012345678901234567890123456789012345678901234567890123456789012
 	
        nustrt=5 ; nustop=5

        do 220 i=1,nustop
           uprmvar(i)=2.0d0*dfloat(i)
c          uprmvar(i)=1+0.5*dfloat(i)
c          uprmvar(i)=0.1*dfloat(i)
c          uprmvar(i)=0.01*dfloat(i)
c          uprmvar(i)=0.08*dfloat(i)
c          uprmvar(i)=10.0*dfloat(i)
c          uprmvar(i)=5.0*dfloat(i)
  220	continue

	cfl0 = cfl
        do 2000 mustep=nustrt,nustop
	
c	if (uprmvar(mustep).le.0.08) then
c	   cfl = 0.05
c	else
 	   cfl = cfl0
c	endif

        iplot=0; itransf=0
	if (lflag.eq.1) deallocate (cx_ave,cy_ave,cz_ave)

	nx=2*(icush+itwb); nxp=nx+1; nptot=nx*ny*nz; 
	ishft=(nt0-nx)/2; lflag=0
        deallocate (x,g)
	allocate (x(nxp))
        allocate (g(-2:nx+3,-2:ny+3,-2:nz+3))
 	if (mustep.eq.nustrt) then
	   allocate (gradg(-2:nx+3,-2:ny+3,-2:nz+3)); gradg=0.0d0
           allocate (gradgold(-2:nx+3,-2:ny+3,-2:nz+3)); gradgold=0.0d0
           allocate (mask(nx,ny,nz),indx(nptot,4)); mask=0; indx=0
 	   allocate (cfn(-2:nx+3,-2:ny+3,-2:nz+3)); cfn=0.0d0
 	endif
	do 225 i=1,nxp
  	   x(i)=0.5d0*dfloat(2*(i-1)-nx)*dx
	   if (dabs(x(i)).lt.small) x(i)=0.0d0
  225	continue

        do 230 k=1,nz
        do 230 j=1,ny
        do 230 i=1,nx
           if (dabs(x(i)).le.twb) then
  	      g(i,j,k)=x(i)
 	   else
  	      g(i,j,k)=dsign(twb,x(i))
 	   endif
	   if (dabs(g(i,j,k)).lt.small) g(i,j,k)=0.0d0
  230	continue
  	call bound(nx,ny,nz,g,ibcg)
  	call tubes(nx,ny,nz,nptot,ntc,mask,indx,cfn,isotube1,
     &	   nisotubes,g,ibcg,ibcv,0,mstep)
	
	timepre1=0.0d0; timepre=0.0d0; time=0.0d0
        uprime=uprmvar(mustep)

c23456789012345678901234567890123456789012345678901234567890123456789012

 	adiv=0.0d0; adiv_dt=0.0d0; adivomega=0.0d0
	fm=2.0d0*pi

	goto 233

        do 232 n3=1,nmodes
c          rn3=dfloat(n3)
           rn3=dfloat(n3)*fm
        do 232 n2=1,nmodes
c          rn2=dfloat(n2)
           rn2=dfloat(n2)*fm
        do 232 n1=1,nmodes
c          rn1=dfloat(n1)
           rn1=dfloat(n1)*fm
	   rnmag=dsqrt(rn1**2.0d0+rn2**2.0d0+rn3**2.0d0)
	   ratio=rnmag/rn1

c	   rnmax=dmax1(rn1,rn2,rn3)
cc	   rnfac=(rnmax/dsqrt(0.5*(rnmag**2.0d0-rnmax**2.0d0)))**8.0d0
c	   rnfac=(rnmax/dsqrt(0.5*(rnmag**2.0d0-rnmax**2.0d0)))**15.0d0
c          ai=rnmag**(-5.0d0/6.0d0)/rnfac
cc	   ai = rnmag**(-5.0d0/6.0d0)

c          epsfac = fm/8.0d0
c          a12 = (dabs(rn1-rn2) + epsfac)/dsqrt(rn1**2.0d0 + rn2**2.0d0)
c          a13 = (dabs(rn1-rn3) + epsfac)/dsqrt(rn1**2.0d0 + rn3**2.0d0)
c          a23 = (dabs(rn2-rn3) + epsfac)/dsqrt(rn2**2.0d0 + rn3**2.0d0)
c          wnfac = a12*a13*a23*dsqrt(8.0d0)
   	   ai = rnmag**(-5.0d0/6.0d0)*wnfac

   	   ai = rnmag**(-5.0d0/6.0d0)
 	   adiv=adiv+((ratio-3.0d0/ratio)*ai/rnmag)**2.0d0
	   dum=(1.0d0-3.0d0*(rn1/rnmag)**2.0d0)/rn1
 	   adiv_dt=adiv_dt+dabs(dum*ai)

 	   adivomega = adivomega + ((rn2/rn3-rn3/rn2)*ai)**2.0d0
c	   adivomega = adivomega + (((rn2/rn3-rn3/rn2)*ai)**2.0d0)
c    &	      /(rn2*rn3)
c	   adivomega = adivomega + (((rn2/rn3-rn3/rn2)*ai/rnmag)**2.0d0)
c	   adivomega = adivomega + (((rn2/rn3-rn3/rn2)*ai)**2.0d0)
c    &	      /(rn2**2+rn3**2)

  232	continue

  233	continue
        do 234 n1=1,nmodes
           adiv = adiv + 1
  234	continue

c       sigma=dsqrt(adiv/8.0d0)
c       sigma=dsqrt(adivomega/8.0d0)

	if (nmodes.gt.1) then
c          sigma=dsqrt(adivomega/8.0d0)/dfloat(nmodes) ! for enstrophy 
c          sigma=dsqrt(adivomega)/(2.0d0*fm) ! for enstrophy 
c          sigma=dsqrt(adivomega/8.0d0)/(4.0d0*pi) ! for enstrophy 
c          sigma=dsqrt(adiv/8.0d0)/dfloat(nmodes)**3.0d0 ! for intensity correlation
           sigma=dsqrt(adiv/8.0d0) ! for intensity correlation
	else
c          sigma=0.5d0 ! for both correlations
           sigma=dsqrt(1.0d0/8.0d0) ! for both correlations
	endif

        sigma=dsqrt(adiv/8.0d0)

	amplt=uprime/sigma
        write(6,*) 'amplt: ',amplt,' uprime: ',uprime,' sigma: ',
     &     sigma,' adiv: ',adiv

c	adiv_dt=adiv_dt/sigma

        sl=sl0/amplt
        wind=wind0/amplt
c       sl=sl0/uprime
c       wind=wind0/uprime

        fac2pi=5.0d0
        nplot=npltkp

c Defining the plotting and stopping times

        if (uprime.ge.10.0d0) then
           fac2pi=0.21d0
        else if (uprime.ge.4.0d0) then
           fac2pi=0.31d0
        else if (uprime.gt.3.0d0) then
           fac2pi=0.5d0
        elseif (uprime.gt.2.0d0) then
c          fac2pi=0.75d0
           fac2pi=1.0d0
        elseif (uprime.ge.1.0d0) then
           fac2pi=1.5d0
        elseif (uprime.ge.0.5) then
c          fac2pi=5.0d0
c          fac2pi=10.0d0
           fac2pi=3.0d0
        elseif (uprime.ge.0.05) then
           fac2pi=10.0d0
        elseif (uprime.ge.0.04) then
           fac2pi=16.0d0
        elseif (uprime.ge.0.03) then
           fac2pi=20.0d0
        elseif (uprime.ge.0.02) then
           fac2pi=40.0d0
        elseif (uprime.ge.0.01) then
           fac2pi=80.0d0
        else
           fac2pi=5
        endif

c       fac2pi=3.0d0*fac2pi

c23456789012345678901234567890123456789012345678901234567890123456789012

cc      stime=dfloat(2/kr)*pi*fac2pi
cc      stime=2.0d0*pi*fac2pi*amplt

        stime=fac2pi*amplt
c       stime=fac2pi*amplt*sigma

cc      delplt=2.0d0*pi*amplt/dfloat(kr*nplot)
        delplt=amplt/dfloat(kr*nplot)
c       delplt=uprime/dfloat(kr*nplot)
        nextplot=-2
        write (6,*) 'stime=',stime,'delplt=',delplt

c Creating the output-file name
        if (uprime.lt.10.0d0) then
           write (str1,'(f6.4)') uprime
        else
           write (str1,'(f6.3)') uprime
        endif
        if (fac2pi.lt.10.0d0) then
           write (str3,'(f4.2)') fac2pi
        else
           write (str3,'(f4.1)') fac2pi
        endif
        write (str2,'(I10)') kr
        write (str4,'(f4.2)') cfl
        write (str5,'(I10)') nt0
        write (str6,'(I10)') ny
        write (str7,'(I10)') limit
        write (str8,'(I10)') itc
        write (str9,'(I10)') irein
        write (str10,'(I10)') nordg
        write (str11,'(I10)') nreinit
        write (str12,'(I10)') iscf
        write (str13,'(I10)') iso
        write (str14,'(I10)') irk
        if (igrd.eq.0) then
	   str15='e'
	else
	   str15='w'
	endif	
        write (str16,'(f4.2)') wind0
        write (str17,'(I10)') ifreq
        write (str18,'(f4.2)') cfl_re
        write (str19,'(I10)') iquad
        write (str20,'(I10)') nmodes
        write (str21,'(I10)') itimelev
        if ((norda.eq.4).or.(norda.eq.5)) then
           write (str22,'(I10)') norda
           str22=adjustl(str22)
	   str22='w'//trim(str22)
	else
	   str22='u4'
	endif	
        write (str23,'(I10)') lrunge
        str1=adjustl(str1)
        str2=adjustl(str2)
        str3=adjustl(str3)
        str4=adjustl(str4)
        str5=adjustl(str5)
        str6=adjustl(str6)
        str7=adjustl(str7)
        str8=adjustl(str8)
        str9=adjustl(str9)
        str10=adjustl(str10)
        str11=adjustl(str11)
        str12=adjustl(str12)
        str13=adjustl(str13)
        str14=adjustl(str14)
        str15=adjustl(str15)
        str16=adjustl(str16)
        str17=adjustl(str17)
        str18=adjustl(str18)
        str19=adjustl(str19)
        str20=adjustl(str20)
        str21=adjustl(str21)
        str23=adjustl(str23)

c23456789012345678901234567890123456789012345678901234567890123456789012

        fileplt='u'//trim(str1)//'k'//trim(str2)//
     &     't'//trim(str3)//'c'//trim(str4)//'n'//trim(str5)//
     &     'x'//trim(str6)//'l'//trim(str7)//'it'//trim(str8)//
     &     trim(str15)//trim(str10)//'tlv'//trim(str21)//
     &	   'rk'//trim(str14)//'r'//trim(str11)//'rf'//trim(str9)//
     &     'rc'//trim(str18)//'sc'//trim(str12)//'w'//trim(str16)//
     &     'f'//trim(str17)//'q'//trim(str19)//'m'//trim(str20)//
     &	   trim(str22)//'lr'//trim(str23)//'.txt'

        write (6,*) fileplt
c       call plot(nx,ny,nz,x,y,z,g,0,ireplt,fileplt)
	
c	call timestamp
        call isosurfts(nt0,nx,ny,nz,x,y,z,nptot,ntc,indx,mask,
     &     isosurf1,nisosurf,g,area,localmax,0,cx_ave,cy_ave,
     &	   cz_ave,igmin,igmax)
c	write (6,*) 'area = ',area
c	call timestamp
c	stop

        open(unit=7,file=fileplt,status='unknown',access='append')
 	write (7,9300) 0,0,size(x),0,0,0,0,x(1),x(nxp),0.0d0,0.0d0,
     &	   1.0d0,area,0
        close(7)
c	write (6,*) 'plot # = ',0,'time/2*pi =',0,'area =',area
 	write (6,*) 'plot # = ',0,'time =',0,'area =',area
        nextplt=-1

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 1000 mstep=1,maxstp
c       do 1000 mstep=1,1
	
	if ((mod(mstep,irein).ne.0).or.(nreinit.eq.0)) then
	   irein_now=0
	else
	   irein_now=1
	endif

	igmin=nxp; igmax=0
        do 235 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
	   if (indx(nt,4).eq.1) then
	      igmin(j,k)=min(igmin(j,k),i)
	      igmax(j,k)=max(igmax(j,k),i)
	   endif
  235	continue

	if (itransf.ne.1) goto 300
	   
c Determination of nx...
c	igmin=nxp; igmax=0
c       do 235 nt=1,ntc
c          i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c	   if (indx(nt,4).eq.1) then
c	      igmin(j,k)=min(igmin(j,k),i)
c	      igmax(j,k)=max(igmax(j,k),i)
c	   endif
c 235	continue

  	imin = minval(igmin); imax = maxval(igmax); izone=imax-imin
 	nxold=nx; nx=2*icush+imax-imin; nxp=nx+1; nptot=nx*ny*nz

c Transfering data from old-time-step variables to new-time-step 
c variables (with new nx)...

        allocate(xtransf(1:nxp))
	call transf1(x,xtransf,nxold,nx,ishft)
        deallocate (x); allocate (x(nxp))
	x=xtransf; deallocate (xtransf)

c23456789012345678901234567890123456789012345678901234567890123456789012

        deallocate (mask,indx,cfn)
        allocate (mask(nx,ny,nz),indx(nptot,4)); mask=0; indx=0
 	allocate (cfn(-2:nx+3,-2:ny+3,-2:nz+3)); cfn=0.0d0
	if (itimelev.eq.3) then
     	   allocate (gradgold2(-2:nx+3,-2:ny+3,-2:nz+3))
	   gradgold2=0.0d0
	endif
 	allocate (gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)); gradg_ave=0.0d0

c23456789012345678901234567890123456789012345678901234567890123456789012
c
c Uncomment the following lines if a formula that depends on the maximumum
c velocity amplitude is used.
c       allocate (u(-2:nx+3,-2:ny+3,-2:nz+3)); u=0
c       allocate (v(-2:nx+3,-2:ny+3,-2:nz+3)); v=0
c       allocate (w(-2:nx+3,-2:ny+3,-2:nz+3)); w=0
c       if (ifreq.eq.0) then
c          trel=0.0d0
c	   call vfield(nx,ny,nz,dx,dy,dz,x,u,v,w,trel,wind-sl)
c	else
c          trel=time/amplt
c    	   call vfield(nx,ny,nz,dx,dy,dy,x,u,v,w,trel,wind)
c 	endif
c Move this line to a location after which (u,v,w) are no longer needed
c      	deallocate (u,v,w)
c
c23456789012345678901234567890123456789012345678901234567890123456789012

c	call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
c    &     gradg_ave,cx_ave,cy_ave,cz_ave)
c	localmx1=localmax
c	if (localmx1.gt.0) then
c	   write (6,*) 'BEFORE GTRANSF: ',localmx1,' local maxima'
c	endif

  	allocate(gtransf(-2:nx+3,-2:ny+3,-2:nz+3))
 	call transf2(g,gtransf,nxold,nx,ny,nz,1)
       	deallocate (g); allocate (g(-2:nx+3,-2:ny+3,-2:nz+3))
  	g=gtransf
cc	deallocate (gtransf)

c	call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,grad,
c    &     gradg_ave,cx_ave,cy_ave,cz_ave)
c	localmx2=localmax
c	if (localmx2.ne.localmx1) then
c	   write (6,*) 'GTRANSF: local max change from ',
c    &        localmx1,' to ',localmx2
c	endif

  	call tubes(nx,ny,nz,nptot,ntc,mask,indx,cfn,isotube1,
     &	   nisotubes,g,ibcg,ibcv,0,mstep)

c 	allocate(gtransf(-2:nx+3,-2:ny+3,-2:nz+3))
 	call transf2(gradgold,gtransf,nxold,nx,ny,nz,0)
       	deallocate (gradgold);
	allocate (gradgold(-2:nx+3,-2:ny+3,-2:nz+3))
  	gradgold=gtransf
c	deallocate (gtransf)
 	
c 	allocate(gtransf(-2:nx+3,-2:ny+3,-2:nz+3))
 	call transf2(gradg,gtransf,nxold,nx,ny,nz,0)
       	deallocate (gradg);
	allocate (gradg(-2:nx+3,-2:ny+3,-2:nz+3))
  	gradg=gtransf; deallocate (gtransf)

	itransf=0
  300	continue

	if (itimelev.eq.3) gradgold2=gradgold
	gradgold=gradg

	if (igrd.eq.0) then
 	   call enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &	   nord,g,g,gradg,ibcv,'evol',0)
	else 
 	   call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &        nord,g,g,gradg,ibcv,'evol',10)
c	   call gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,g,gradg,ibcv,'evol',10)
	endif

	if (mstep.eq.1) then
	   if (itimelev.eq.3) gradgold2=gradg
	   gradgold=gradg

	endif

c23456789012345678901234567890123456789012345678901234567890123456789012


c time derivatives will need to be redirived if dt varies between time steps
c       dt=dmin1(dx,dy,dz)*cfl/umax
c       dt=cfl*dx*amplt/ceiling(sl0+amplt*dsqrt(3.0d0))
cc      dt=cfl*dx/ceiling(wind+dsqrt(3.0d0)) 
c       dt=cfl*dx/(wind+adiv*dsqrt(3.0d0)) 

	if ((ifreq.eq.0).and.(mstep.eq.1)) then
           dt=cfl*dx/(adiv_dt*dsqrt(3.0d0))
	   if (nmodes.eq.1) dt = cfl*dx/dsqrt(3.0d0)

c	The number 13.7 is an estimate of the maximum velocity in the 
c	excitation field (umax), which would need to be determined before
c	calculating dt with the formula below. With the appropriate
c	value for umax replacing 13.7, the formula below will give a larger 
c	CFL number than the one above (better).
c          dt=cfl*dx/(13.7*dsqrt(3.0d0))

	elseif (mstep.eq.1) then
 	   dt=cfl*dx/(wind+adiv_dt*dsqrt(3.0d0))
	   if (nmodes.eq.1) dt = cfl*dx/(wind+dsqrt(3.0d0))
	endif

c	write (6,*) mstep, time+dt, dt, adiv_dt, cfl*dx/adiv_dt, 
c    &	   cfl*dx/(wind+adiv_dt*dsqrt(3.0d0)),dx

c	if (mstep.eq.1) then
c	   dt0=dt
c	elseif (ifreq.eq.0) then
c	  cx_ave=dt*cx_ave/dt0 
c	  cy_ave=dt*cy_ave/dt0 
c	  cz_ave=dt*cz_ave/dt0 
c	endif

	timepre2=timepre1; timepre1=timepre; timepre=time

	if (lrunge.ne.1) then

	if ((ifreq.eq.0).and.(lflag.eq.0)) then
           allocate (cx_ave(-2:nt0+3,-2:ny+3,-2:nz+3)); cx_ave=0.0d0
     	   allocate (cy_ave(-2:nt0+3,-2:ny+3,-2:nz+3)); cy_ave=0.0d0
     	   allocate (cz_ave(-2:nt0+3,-2:ny+3,-2:nz+3)); cz_ave=0.0d0
	   trel=0.0d0; uconst=wind-sl
 	   call cfieldtot(nt0,ny,nz,dx,dy,dz,dt,trel,uconst,
     &	      cx_ave,cy_ave,cz_ave,0.0d0,1) ! Evaluation at x(i,j,k)

	   write (6,*) 'number of modes = ',nmodes,
     &	      '   uprime_rel = ',uprime_rel

	elseif (ifreq.ne.0) then
           allocate (cx_ave(-2:nx+3,-2:ny+3,-2:nz+3)); cx_ave=0.0d0
     	   allocate (cy_ave(-2:nx+3,-2:ny+3,-2:nz+3)); cy_ave=0.0d0
     	   allocate (cz_ave(-2:nx+3,-2:ny+3,-2:nz+3)); cz_ave=0.0d0
           allocate (cx(-2:nx+3,-2:ny+3,-2:nz+3)); cx=0.0d0
     	   allocate (cy(-2:nx+3,-2:ny+3,-2:nz+3)); cy=0.0d0
     	   allocate (cz(-2:nx+3,-2:ny+3,-2:nz+3)); cy=0.0d0
           allocate (cx_old(-2:nx+3,-2:ny+3,-2:nz+3)); cx_old=0.0d0
     	   allocate (cy_old(-2:nx+3,-2:ny+3,-2:nz+3)); cy_old=0.0d0
     	   allocate (cz_old(-2:nx+3,-2:ny+3,-2:nz+3)); cz_old=0.0d0
           allocate (cx_old2(-2:nx+3,-2:ny+3,-2:nz+3)); cx_old2=0.0d0
     	   allocate (cy_old2(-2:nx+3,-2:ny+3,-2:nz+3)); cy_old2=0.0d0
     	   allocate (cz_old2(-2:nx+3,-2:ny+3,-2:nz+3)); cz_old2=0.0d0

	   uconst=wind

           trel=timepre/amplt
c          trel=timepre/uprime
 	   call cfield(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel,uconst,
     &	      cx,cy,cz,0.0d0,nptot,ntc,indx,1)

           trel=timepre1/amplt
c          trel=timepre1/uprime
 	   call cfield(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel,uconst,
     &	      cx_old,cy_old,cz_old,0.0d0,nptot,ntc,indx,1)

           trel=timepre2/amplt
c          trel=timepre2/uprime
 	   call cfield(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel,uconst,
     &	      cx_old2,cy_old2,cz_old2,0.0d0,nptot,ntc,indx,1)

           cx_ave=(15.0d0*cx-10.0d0*cx_old+3.0d0*cx_old2)/8.0d0
           cy_ave=(15.0d0*cy-10.0d0*cy_old+3.0d0*cy_old2)/8.0d0
           cz_ave=(15.0d0*cz-10.0d0*cz_old+3.0d0*cz_old2)/8.0d0
           deallocate (cx_old,cy_old,cz_old,cx_old2,cy_old2,cz_old2)
           deallocate (cx,cy,cz)
	endif

	elseif (mstep.eq.1) then
           allocate (cx_ave(-2:nt0+3,-2:ny+3,-2:nz+3)); cx_ave=0.0d0
     	   allocate (cy_ave(-2:nt0+3,-2:ny+3,-2:nz+3)); cy_ave=0.0d0
     	   allocate (cz_ave(-2:nt0+3,-2:ny+3,-2:nz+3)); cz_ave=0.0d0
  	endif ! condtion that lrunge.ne.1

	if (itimelev.eq.3) then

c 3rd-order collocation of c at t + dt/2, 2nd-order approx to the time average;
c used to initiate the iteration below to find the fully 3rd-order average 
c (in time and space)--appropriate if reinitialization was not performed at
c the end of the previous time step (i.e., for irein >= 2).

           gradg_ave=(15.0d0*gradg-10.0d0*gradgold+3.0d0*gradgold2)
     &	      /8.0d0

c 2nd-order collocation of c at t + dt/2, 2nd-order approx to the time average
c Temperton & Staniforth-Q. J. R. Meteorol. SOC, vol. 113, pp. 1025-1039 (1987)
c Equations (39) & (40) therein. Only the most-recent and current times are 
c appropriate for time averaging when reinitialization was performed on the
c solution obtained during the previous time step (i.e., when irein = 1).

	else
           gradg_ave=(3.0d0*gradg-gradgold)/2.0d0
	endif

c	gradg_ave=gradg

c23456789012345678901234567890123456789012345678901234567890123456789012

cc This method (directly below) is used when the velocity field is not an
cc analytical function that can be directly evaluated at the departure point

cc Iteration a la Temperton & Staniforth-Q. J. R. Meteorol. SOC, 
cc vol. 113, pp. 1025-1039 (1987)-Equations (38), (40) & (41) therein

cc      allocate (cx_old2(-2:nx+3,-2:ny+3,-2:nz+3))
cc   	allocate (cy_old2(-2:nx+3,-2:ny+3,-2:nz+3))
cc   	allocate (cz_old2(-2:nx+3,-2:ny+3,-2:nz+3))

cc      do 590 itno=1,itc

c Initialize the on-grid time-integrated advection velocity (stored in 
c the dummy variables cx_old2, cy_old2 & cz_old2)

cc      if (itno.eq.1) then
cc         cx_old2=cx_ave; cy_old2=cy_ave; cz_old2=cz_ave
cc      endif

c Calculation of the interpolated [from x(i,j,k), y(i,j,k), z(i,j,k)] 
c departure-point time-integrated advection velocity (stored in cx_ave, 
c cy_ave & cz_ave)-[Aldredge, 61-124]

cc	if (itc.eq.1) then

cc	   call interp(nx,ny,nz,cx_old2,cy_old2,cz_old2,0.5d0,
cc   &        cx_old2,cx_ave,0)
cc	   call interp(nx,ny,nz,cx_old2,cy_old2,cz_old2,0.5d0,
cc   &        cy_old2,cy_ave,0)
cc	   call interp(nx,ny,nz,cx_old2,cy_old2,cz_old2,0.5d0,
cc   &        cz_old2,cz_ave,0)

cc	else

c23456789012345678901234567890123456789012345678901234567890123456789012

cc         allocate (cx_tmp(-2:nx+3,-2:ny+3,-2:nz+3))
cc   	   allocate (cy_tmp(-2:nx+3,-2:ny+3,-2:nz+3))
cc   	   allocate (cz_tmp(-2:nx+3,-2:ny+3,-2:nz+3))

cc	   call interp(nx,ny,nz,cx_ave,cy_ave,cz_ave,0.5d0,
cc   &        cx_old2,cx_tmp,0)
cc	   call interp(nx,ny,nz,cx_ave,cy_ave,cz_ave,0.5d0,
cc   &        cy_old2,cy_tmp,0)
cc	   call interp(nx,ny,nz,cx_ave,cy_ave,cz_ave,0.5d0,
cc   &        cz_old2,cz_tmp,0)
cc         cx_ave=cx_tmp; cy_ave=cy_tmp; cz_ave=cz_tmp

cc         deallocate (cx_tmp,cy_tmp,cz_tmp)

cc	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc This (below) does not seem to make any difference
cc 	call interp(nx,ny,nz,cx_ave,cy_ave,cz_ave,0.5d0,
cc    &     cx_old2,cx_ave,0)
cc 	call interp(nx,ny,nz,cx_ave,cy_ave,cz_ave,0.5d0,
cc    &     cy_old2,cy_ave,0)
cc 	call interp(nx,ny,nz,cx_ave,cy_ave,cz_ave,0.5d0,
cc    &     cz_old2,cz_ave,0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc590	continue

cc      deallocate (cx_old2,cy_old2,cz_old2)

c The methods (below) are more efficient but require an analytical 
c expression of the velocity field (specified in subroutine "cfield").

        allocate (dep(ntc,3)); dep=0.0d0
        allocate (dephalf(ntc,3)); dephalf=0.0d0

c23456789012345678901234567890123456789012345678901234567890123456789012

	if (lrunge.eq.1) then !  METHOD 1:

 	   call departure(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel,uconst,
     &	      dep,dephalf,nptot,ntc,indx)

	else ! METHOD 2:

c Iteration a la Temperton & Staniforth-Q. J. R. Meteorol. SOC, 
c vol. 113, pp. 1025-1039 (1987)-Equations (38), (40) & (41) therein

	   if (lflag.eq.0) then

c Evaluation of the vel field at x-dt*vel/2 (midway between x_ijk and 
c the departure point; needs to be done only once if ifreq==0 (no time 
c dependency of the velocity field)

	      if (ifreq.eq.0) then

c		 allocate (cx_tmp(-2:nx+3,-2:ny+3,-2:nz+3))
c		 allocate (cy_tmp(-2:nx+3,-2:ny+3,-2:nz+3))
c		 allocate (cz_tmp(-2:nx+3,-2:ny+3,-2:nz+3))
c		 cx_tmp=cx_ave; cy_tmp=cy_ave; cz_tmp=cz_ave
 
c (should give greater accuracy but appears to make the results worse)--10/21/15
                 do 590 itno=1,itc
 	            call cfieldtot(nt0,ny,nz,dx,dy,dz,dt,trel/2.0d0,
     &	               uconst,cx_ave,cy_ave,cz_ave,0.5d0,1)
cc   &	               uconst,cx_ave,cy_ave,cz_ave,1.0d0,0)
  590	         continue

c	         call cfieldtot(nt0,ny,nz,dx,dy,dz,dt,trel/2.0d0,
c    &	            uconst,cx_ave,cy_ave,cz_ave,1.0d0,0)
c		 cx_ave=0.5d0*(cx_ave+cx_tmp)
c		 cy_ave=0.5d0*(cy_ave+cy_tmp)
c		 cz_ave=0.5d0*(cz_ave+cz_tmp)
c       	 deallocate (cx_tmp,cy_tmp,cz_tmp)

	         lflag=1
	      else
                 do 600 itno=1,itc
 	   	    call cfield(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel/2.0d0,
     &	      	       uconst,cx_ave,cy_ave,cz_ave,0.5d0,nptot,ntc,indx,1)
  600	         continue
	      endif
	   endif

 	   do nt=1,ntc
 	      i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
	      if ((i+ishft).gt.0) then
		 i=mod((i+ishft-1),nt0)+1
	      else
 	         i=nt0+mod((i+ishft),nt0)
	      endif
	      dep(nt,1)=cx_ave(i,j,k); dephalf(nt,1)=cx_ave(i,j,k)/2.0d0
	      dep(nt,2)=cy_ave(i,j,k); dephalf(nt,2)=cy_ave(i,j,k)/2.0d0
	      dep(nt,3)=cz_ave(i,j,k); dephalf(nt,3)=cz_ave(i,j,k)/2.0d0

c	      if ((cx_tmp(i,j,k)*dep(nt,1)).lt.-small) then
c		 dep(nt,1)=0.0d0
c	      endif
c	      if ((cy_tmp(i,j,k)*dep(nt,2)).lt.-small) then
c		 dep(nt,2)=0.0d0
c	      endif
c	      if ((cz_tmp(i,j,k)*dep(nt,3)).lt.-small) then
c		 dep(nt,3)=0.0d0
c	      endif

c	      if (dabs(dep(nt,1)).lt.small) then
c		 dep(nt,1)=0.5d0*(cx_ave(i-1,j,k)+cx_ave(i+1,j,k))
c	      endif
c	      if (dabs(dep(nt,2)).lt.small) then
c		 dep(nt,2)=0.5d0*(cy_ave(i,j-1,k)+cy_ave(i,j+1,k))
c	      endif
c	      if (dabs(dep(nt,3)).lt.small) then
c		 dep(nt,3)=0.5d0*(cz_ave(i,j,k-1)+cz_ave(i,j,k+1))
c	      endif

c	      if (((cx_ave(i-1,j,k)*dep(nt,1)).lt.-small).or.
c    &		  ((cx_ave(i+1,j,k)*dep(nt,1)).lt.-small)) then
c		 dep(nt,1)=0.0d0
c	      endif
c	      if (((cy_ave(i,j-1,k)*dep(nt,2)).lt.-small).or.
c    &		  ((cy_ave(i,j+1,k)*dep(nt,2)).lt.-small)) then
c		 dep(nt,2)=0.0d0
c	      endif
c	      if (((cz_ave(i,j,k-1)*dep(nt,3)).lt.-small).or.
c    &		  ((cz_ave(i,j,k+1)*dep(nt,3)).lt.-small)) then
c		 dep(nt,3)=0.0d0
c	      endif

c	      if ((cx_ave(i-1,j,k)*cx_ave(i+1,j,k)).lt.-small) then
c		 dep(nt,1)=0.0d0
c	      endif
c	      if ((cy_ave(i,j-1,k)*cy_ave(i,j+1,k)).lt.-small) then
c		 dep(nt,2)=0.0d0
c	      endif
c	      if ((cz_ave(i,j,k-1)*cz_ave(i,j,k+1)).lt.-small) then
c		 dep(nt,3)=0.0d0
c	      endif
	   end do

	endif

c23456789012345678901234567890123456789012345678901234567890123456789012

c Evaluation of the grad at a location halfway between the 
c departure point and x(i,j,k):
c (Doesn't seem to improve the results)--10/21/15

c	if (itimelev.eq.3) then
c          gradgold2=gradg_ave
c 	   call interp(nx,ny,nz,nptot,ntc,indx,dephalf,
c    &        1.0d0,gradgold2,gradg_ave,0)
c	else
c          gradgold=gradg_ave
c 	   call interp(nx,ny,nz,nptot,ntc,indx,dephalf,
c    &        1.0d0,gradgold,gradg_ave,0)
c	endif

c Updating g: Interpolation from the on-grid value of g [at x(i,j,k), 
c y(i,j,k), z(i,j,k)] at t(n+1) to the get its value at the departure 
c point [xL(t(n)), yL(t(n)), zL(t(n))]

	time=time+dt
        time2pi=time/amplt
c       if ((ireplt0.eq.0).and.(time2pi.ge.pltstrt)
        if ((ireplt0.eq.0).and.(time.ge.pltstrt)
     &     .and.(time.le.pltstop)) then
c	   write (6,*) 'time: ',time,' dt: ',dt
           ireplt=0
        else
           ireplt=1
        endif

	call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &     gradg_ave,cx_ave,cy_ave,cz_ave)
	localmx1=localmax
	if (localmx1.gt.0) then
	   write (6,*) 'BEFORE INTERPG: ',localmx1,' local maxima'
	endif

	if ((norda.eq.4).or.(norda.eq.5)) then
 	   call interpgweno(nx,ny,nz,nptot,ntc,mask,indx,dep,
     &        1.0d0,g,ibcg,ibcv,limit,cfn,gradg_ave,dt,sl)
	else
 	   call interpg(nx,ny,nz,nptot,ntc,mask,indx,dep,x,ishft,
     &        1.0d0,g,ibcg,limit,cfn,gradg_ave,dt,sl,mstep,nt0,cx_ave,
     &	      cy_ave,cz_ave,nmodes)
c	   call interpg2(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
c    &        1.0d0,g,ibcg,limit,cfn,gradg_ave,dt,sl,mstep,nt0,cx_ave,
c    &	      cy_ave,cz_ave,nmodes)
c	   call interpg3(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
c    &        1.0d0,g,ibcg,limit,cfn,gradg_ave,dt,sl,mstep,nt0,cx_ave,
c    &	      cy_ave,cz_ave,nmodes)
c	   call interpg4(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
c    &        1.0d0,g,ibcg,limit,cfn,gradg_ave,dt,sl,mstep,nt0,cx_ave,
c    &	      cy_ave,cz_ave,nmodes)
c	   call interpg5(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
c    &        1.0d0,g,ibcg,limit,cfn,gradg_ave,dt,sl,mstep,nt0,cx_ave,
c    &	      cy_ave,cz_ave,nmodes)
c	   call interpg6(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
c    &        1.0d0,g,ibcg,limit,cfn,gradg_ave,dt,sl,mstep,nt0,cx_ave,
c    &	      cy_ave,cz_ave,nmodes)
	endif

	call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &     gradg_ave,cx_ave,cy_ave,cz_ave)
	localmx2=localmax
	if (localmx2.ne.localmx1) then
	   write (6,*) 'INTERPG6: local max change from ',
     &        localmx1,' to ',localmx2
	endif

 	call tubes(nx,ny,nz,nptot,ntc,mask,indx,cfn,isotube1,
     &	   nisotubes,g,ibcg,ibcv,1,mstep)
	itransf=1

        call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &     gradg_ave,cx_ave,cy_ave,cz_ave)
	localmx3=localmax
	if (localmx3.ne.localmx2) then
	   write (6,*) 'TUBES AFTER INTERPOLATION: local max change 
     &        from ',localmx2,' to ',localmx3
	endif

c23456789012345678901234567890123456789012345678901234567890123456789012

c Surface-area calculation, reinitialization and tube reconstruction

        call isosurfts(nt0,nx,ny,nz,x,y,z,nptot,ntc,indx,mask,
     &     isosurf1,nisosurf,g,area1,localmax,0,cx_ave,cy_ave,
     &	   cz_ave,igmin,igmax)

	if (irein_now.eq.1) then

	   call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &        gradg_ave,cx_ave,cy_ave,cz_ave)
	   localmx1=localmax

 	   call reinwenorkts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &	      iscf,nreinit,nordg,ibcg,ibcv,cfn,g)

	   call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &        gradg_ave,cx_ave,cy_ave,cz_ave)
	   localmx2=localmax
 	   if (localmx2.ne.localmx1) then
	      maskmin=minval(mask/gjump)
	      maskmax=maxval(mask/gjump)
              gmin=minval(g/gjump)
              gmax=maxval(g/gjump)
              gminL=minval(g(1,1:ny,1:nz)/gjump)
              gminR=minval(g(nx+1,1:ny,1:nz)/gjump)
              gmaxL=maxval(g(1,1:ny,1:nz)/gjump)
              gmaxR=maxval(g(nx+1,1:ny,1:nz)/gjump)
	      write (6,*) 'REINWENORKTS: local max change from ',
     &           localmx1,' to ',localmx2
	      write (6,21) maskmin,maskmax,gmin,gmax,gminL,gminR,
     &		 gmaxL,gmaxR
 21           format (2(1x,i3),6(1x,e18.12))
 	   endif

 	   call tubes(nx,ny,nz,nptot,ntc,mask,indx,cfn,isotube1,
     &	      nisotubes,g,ibcg,ibcv,1,mstep)
	   itransf=1

           call chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &        gradg_ave,cx_ave,cy_ave,cz_ave)
	   localmx3=localmax
	   if (localmx3.ne.localmx2) then
	      write (6,*) 'TUBES AFTER REINITIALIZATION: local max 
     &           change from ',localmx2,' to ',localmx3
	   endif

	endif

c       deallocate (dep,dephalf,gradg_ave)
	if (itimelev.eq.3) deallocate (gradgold2)
	if ((lrunge.ne.1).and.(lflag.eq.0)) then
	   deallocate (cx_ave,cy_ave,cz_ave)
	endif

c23456789012345678901234567890123456789012345678901234567890123456789012

c Printing the results 

c       if ((uprime.le.0.02).and.(mod(mstep,4).ne.0)) then
c          deallocate (dep,dephalf,gradg_ave)
c	   goto 700
c	endif

c	if ((mod(mstep,8).eq.0).or.(mstep.eq.300)) then
c	if (mstep.gt.0) then

        if (((mstep.gt.1).and.(dabs(mod(time,delplt)-dt).le.(0.5d0*dt)))
     &     .or.(mstep.ge.maxstp).or.(ireplt.eq.0)
     &     .or.((mod(mstep,2).eq.0).and.(delplt.lt.2.0d0*dt))
     &	   .or.(mstep.eq.1).or.(mstep.le.10)) then
cc   &	   .or.(mstep.eq.1)) then

 	   iplot=iplot+1

	   if (irein_now.eq.1) then
              call isosurfts(nt0,nx,ny,nz,x,y,z,nptot,ntc,indx,mask,
     &           isosurf1,nisosurf,g,area2,localmax,2,cx_ave,cy_ave,
     &		 cz_ave,igmin,igmax)
c	      write (6,*) "localmax: 2nd call ",localmax
              areachg=(area2-area1)/area1
	   else
	      area2=area1; areachg=0.0d0
	   endif

           open(unit=7,file=fileplt,status='unknown',
     &        access='append')

 	   write (7,9300) iplot,mstep,size(x),imin,
     &        imax,ishft,mod((imin+ishft),nt0+1),
     &	      x(1),x(nxp),time2pi,areachg,area1,area2,localmax
           close(7)
 	   write (6,9300) iplot,mstep,size(x),imin,
     &        imax,ishft,mod((imin+ishft),nt0+1),
     &	      x(1),x(nxp),time2pi,areachg,area1,area2,localmax

cc         if ((iplot.ge.58).and.(iplot.le.88)) then
c          if ((mstep.ge.180).and.(mstep.le.210)) then
           if (ireplt.eq.0) then
 	      call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &           nord,g,g,gradg,ibcv,'evol',10)
              call plot(nx,ny,nz,x,y,z,g,mstep,0,fileplt,'n')
              call plot(nt0,ny,nz,x,y,z,cx_ave,mstep,0,fileplt,'x')
              call plot(nt0,ny,nz,x,y,z,cy_ave,mstep,0,fileplt,'y')
              call plot(nt0,ny,nz,x,y,z,cz_ave,mstep,0,fileplt,'z')
           endif

        endif
        deallocate (dep,dephalf,gradg_ave)

 700    continue

c23456789012345678901234567890123456789012345678901234567890123456789012

c Error calculation
           if (igauss.eq.0) goto 900
           errsum=0.0d0
           errmax=0.0d0
           errmin=1.0d0/small
           do 800 k=1,nz
           do 800 j=1,ny
 	      y2=(y(j)+rad*dsin(2.0d0*pi*time))**2.0d0
 	      z2=(z(k)-rad*dcos(2.0d0*pi*time))**2.0d0
 	      temp=dexp(-0.5d0*(y2+z2)/sigma**2.0d0)
 	      err=dabs(g(1,j,k)-0.5d0*temp/(pi*sigma**2.0d0))
              errmax=dmax1(err,errmax)
              errmin=dmin1(err,errmin)
              errsum=errsum+err*dx**2.0d0
 800       continue
           open(unit=10,file='gerror.txt',status='unknown',
     &        access='append')
           if (iplot.eq.1) write (10,8500)
 	   write (10,8600) nx,dx,dt,cfl,errmin,errmax,errsum,cpu,time
 	   close(10)
           write (6,8500)
 	   write (6,8600) nx,dx,dt,cfl,errmin,errmax,errsum,cpu,time
 900       continue

 	   if ((dabs(time-stime).lt.(dt/2.0d0)).or.(mstep.ge.maxstp))
     &        goto 2000

 1000   continue        ! End of the timestep

 2000   continue        ! End of the loop through uprime values

c       call cpu_time(cpu)

 8500   format (2x,'nx',6x,'dx',9x,'dt',8x,'cfl',6x,'errmin',
     &     5x,'errmax',5x,'errsum',5x,'cputime',5x,'time')
 8600   format (1x,i3,8(1x,e10.4))
 9000   format (1x,i5,4(1x,e11.5),i5)
 9100   format (1x,5(1x,f20.15))
c9200   format (1x,2(i4,1x),2(e15.7,1x))
 9200   format (1x,2(i6,1x),5(e18.12,1x))
 9300   format (2(i6,1x),i4,1x,4(i5,1x),2x,6(e18.12,2x),i5)
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine bound(nx,ny,nz,g,ibc)
        implicit double precision (a-h,o-z)
        dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	common /tubedat/twa,twb,gjump,nomask
        common /iparam/iquad

c 1 => computation on a quadrant of the full-period domain; otherwise 
c   on the full-period domain

        if (iquad.eq.0) goto 400

c Periodic (ibc=0) or jump periodic (ibc=2) in x and periodic in y & z

	if ((ibc.eq.0).or.(ibc.eq.2)) then
	   bjump=gjump*dfloat(ibc/2)
	   do 50 k=1,nz
	   do 50 j=1,ny
              g(0,j,k)=-0.5d0*bjump
              g(-1,j,k)=-0.5d0*bjump
              g(-2,j,k)=-0.5d0*bjump
              g(nx+1,j,k)=0.5d0*bjump
              g(nx+2,j,k)=0.5d0*bjump
              g(nx+3,j,k)=0.5d0*bjump
  50	   continue

	   do 100 k=1,nz
	   do 100 i=-2,nx+3
              g(i,0,k)=g(i,2,k)
              g(i,-1,k)=g(i,3,k)
              g(i,-2,k)=g(i,4,k)
              g(i,ny+1,k)=g(i,ny-1,k)
              g(i,ny+2,k)=g(i,ny-2,k)
              g(i,ny+3,k)=g(i,ny-3,k)
  100	   continue

	   do 200 j=-2,ny+3
	   do 200 i=-2,nx+3
              g(i,j,0)=g(i,j,2)
              g(i,j,-1)=g(i,j,3)
              g(i,j,-2)=g(i,j,4)
              g(i,j,nz+1)=g(i,j,nz-1)
              g(i,j,nz+2)=g(i,j,nz-2)
              g(i,j,nz+3)=g(i,j,nz-3)
  200	   continue

        else
           write (6,*) 'ERROR in the bound subroutine'
           write (6,*) 'ibcv must be 0 and ibcg either 0 or 2 if 
     &        iquad = 1'
           write (6,*) 'ibc = ', ibc
           stop
        endif
        goto 900

  400	   continue

	if ((ibc.eq.0).or.(ibc.eq.2)) then
	   bjump=gjump*dfloat(ibc/2)
	   do 450 k=1,nz
	   do 450 j=1,ny
              g(0,j,k)=-0.5d0*bjump
              g(-1,j,k)=-0.5d0*bjump
              g(-2,j,k)=-0.5d0*bjump
              g(nx+1,j,k)=0.5d0*bjump
              g(nx+2,j,k)=0.5d0*bjump
              g(nx+3,j,k)=0.5d0*bjump
  450	   continue
	   do 500 k=1,nz
	   do 500 i=-2,nx+3
              g(i,0,k)=g(i,ny,k)
              g(i,-1,k)=g(i,ny-1,k)
              g(i,-2,k)=g(i,ny-2,k)
              g(i,ny+1,k)=g(i,1,k)
              g(i,ny+2,k)=g(i,2,k)
              g(i,ny+3,k)=g(i,3,k)
  500	   continue
	   do 550 j=-2,ny+3
	   do 550 i=-2,nx+3
              g(i,j,0)=g(i,j,nz)
              g(i,j,-1)=g(i,j,nz-1)
              g(i,j,-2)=g(i,j,nz-2)
              g(i,j,nz+1)=g(i,j,1)
              g(i,j,nz+2)=g(i,j,2)
              g(i,j,nz+3)=g(i,j,3)
  550	   continue

c Periodic (ibc=3) or jump periodic (ibc=4) in x and constant slope in y & z:

	elseif ((ibc.eq.3).or.(ibc.eq.4)) then
	   bjump=gjump*dfloat(ibc-3)
	   do 600 k=1,nz
	   do 600 j=1,ny
              g(0,j,k)=g(nx,j,k)-bjump
              g(-1,j,k)=g(nx-1,j,k)-bjump
              g(-2,j,k)=g(nx-2,j,k)-bjump
              g(nx+1,j,k)=g(1,j,k)+bjump
              g(nx+2,j,k)=g(2,j,k)+bjump
              g(nx+3,j,k)=g(3,j,k)+bjump
  600	   continue
	   do 650 k=1,nz
	   do 650 i=-2,nx+3
              g(i,0,k)=2.0d0*g(i,1,k)-g(i,2,k)
              g(i,-1,k)=2.0d0*g(i,0,k)-g(i,1,k)
              g(i,-2,k)=2.0d0*g(i,-1,k)-g(i,0,k)
              g(i,ny+1,k)=2.0d0*g(i,ny,k)-g(i,ny-1,k)
              g(i,ny+2,k)=2.0d0*g(i,ny+1,k)-g(i,ny,k)
              g(i,ny+3,k)=2.0d0*g(i,ny+2,k)-g(i,ny+1,k)
  650	   continue
	   do 700 j=-2,ny+3
	   do 700 i=-2,nx+3
              g(i,j,0)=2.0d0*g(i,j,1)-g(i,j,2)
              g(i,j,-1)=2.0d0*g(i,j,0)-g(i,j,1)
              g(i,j,-2)=2.0d0*g(i,j,-1)-g(i,j,0)
              g(i,j,nz+1)=2.0d0*g(i,j,nz)-g(i,j,nz-1)
              g(i,j,nz+2)=2.0d0*g(i,j,nz+1)-g(i,j,nz)
              g(i,j,nz+3)=2.0d0*g(i,j,nz+2)-g(i,j,nz+1)
  700	   continue

c Constant slope in x, y & z:

	else
	   do 750 k=1,nz
	   do 750 j=1,ny
              g(0,j,k)=2.0d0*g(1,j,k)-g(2,j,k)
              g(-1,j,k)=2.0d0*g(0,j,k)-g(1,j,k)
              g(-2,j,k)=2.0d0*g(-1,j,k)-g(0,j,k)
              g(nx+1,j,k)=2.0d0*g(nx,j,k)-g(nx-1,j,k)
              g(nx+2,j,k)=2.0d0*g(nx+1,j,k)-g(nx,j,k)
              g(nx+3,j,k)=2.0d0*g(nx+2,j,k)-g(nx+1,j,k)
  750	   continue
	   do 800 k=1,nz
	   do 800 i=-2,nx+3
              g(i,0,k)=2.0d0*g(i,1,k)-g(i,2,k)
              g(i,-1,k)=2.0d0*g(i,0,k)-g(i,1,k)
              g(i,-2,k)=2.0d0*g(i,-1,k)-g(i,0,k)
              g(i,ny+1,k)=2.0d0*g(i,ny,k)-g(i,ny-1,k)
              g(i,ny+2,k)=2.0d0*g(i,ny+1,k)-g(i,ny,k)
              g(i,ny+3,k)=2.0d0*g(i,ny+2,k)-g(i,ny+1,k)
  800	   continue
	   do 850 j=-2,ny+3
	   do 850 i=-2,nx+3
              g(i,j,0)=2.0d0*g(i,j,1)-g(i,j,2)
              g(i,j,-1)=2.0d0*g(i,j,0)-g(i,j,1)
              g(i,j,-2)=2.0d0*g(i,j,-1)-g(i,j,0)
              g(i,j,nz+1)=2.0d0*g(i,j,nz)-g(i,j,nz-1)
              g(i,j,nz+2)=2.0d0*g(i,j,nz+1)-g(i,j,nz)
              g(i,j,nz+3)=2.0d0*g(i,j,nz+2)-g(i,j,nz+1)
  850	   continue
	endif

  900	continue
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012

	subroutine plot(nx,ny,nz,x,y,z,g,iplot,ireplt,pltname,dir)
        implicit double precision (a-h,o-z)
        character*80 str
        character*160 pltname,filename
	character dir
        dimension x(nx+1),y(ny+1),z(nz+1)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
        common /rparam/small,pi

	nxp=nx+1
	nyp=ny+1
	nzp=nz+1

        if ((iplot.eq.0).or.(ireplt.eq.0)) then
           write (str,'(I10)') iplot
	   if (dir.eq.'n') then
	      filename='g'//trim(adjustl(str))//trim(adjustl(pltname))
	   else
	      filename='u'//dir//trim(adjustl(str))
     &		 //trim(adjustl(pltname))
	   endif
        else   
           filename='g'//trim(adjustl(pltname))
	   write (6,*) filename
        endif

        open(unit=9,file=filename,status='replace')
	if (dir.ne.'n') then
	   write (9,6500) (-0.5d0+dfloat(i-1)/dfloat(nx),i=1,nxp)
	else
 	   write (9,6500) (x(i),i=1,nxp)
	endif
 	write (9,6500) (y(j),j=1,nyp)
 	write (9,6500) (z(k),k=1,nzp)
        do 400 i=1,nxp
 	write (9,6300)
c	if (dir.ne.'n') then
c	   write (9,6301) i,-0.5d0+dfloat(i-1)/dfloat(nx),y(1),y(nyp),
c    &	      z(1),z(nzp)
c	else
c	   write (9,6301) i,x(i),y(1),y(nyp),z(1),z(nzp)
c	endif
        do 300 j=1,nyp
 	   write (9,6500) (g(i,j,k),k=1,nzp)
c	   if ((i.eq.10).and.(j.eq.1)) then
c	      write (6,*) 'x: ',x(i),' y: ',y(j)
c	      write (6,6500) (g(i,j,k),k=1,nzp)
c	   endif
  300      continue
  400   continue

c	write (9,*) ""
c	write (9,6302) g(25,15,15),g(25,16,15),g(25,17,15),g(25,18,15)
c    &	   ,g(25,19,15)
c	write (9,6302) g(25,15,16),g(25,16,16),g(25,17,16),g(25,18,16)
c    &	   ,g(25,19,16)
c	write (9,6302) g(25,15,17),g(25,16,17),g(25,17,17),g(25,18,17)
c    &	   ,g(25,19,17)
c	write (9,6302) g(25,15,18),g(25,16,18),g(25,17,18),g(25,18,18)
c    &	   ,g(25,19,18)
c	write (9,6302) g(25,15,19),g(25,16,19),g(25,17,19),g(25,18,19)
c    &	   ,g(25,19,19)
c	write (9,*) ""
c	write (9,6302) g(26,15,15),g(26,16,15),g(26,17,15),g(26,18,15)
c    &	   ,g(26,19,15)
c	write (9,6302) g(26,15,16),g(26,16,16),g(26,17,16),g(26,18,16)
c    &	   ,g(26,19,16)
c	write (9,6302) g(26,15,17),g(26,16,17),g(26,17,17),g(26,18,17)
c    &	   ,g(26,19,17)
c	write (9,6302) g(26,15,18),g(26,16,18),g(26,17,18),g(26,18,18)
c    &	   ,g(26,19,18)
c	write (9,6302) g(26,15,19),g(26,16,19),g(26,17,19),g(26,18,19)
c    &	   ,g(26,19,19)
c	write (9,*) ""
c	write (9,6302) g(28,15,15),g(28,16,15),g(28,17,15),g(28,18,15)
c    &	   ,g(28,19,15)
c	write (9,6302) g(28,15,16),g(28,16,16),g(28,17,16),g(28,18,16)
c    &	   ,g(28,19,16)
c	write (9,6302) g(28,15,17),g(28,16,17),g(28,17,17),g(28,18,17)
c    &	   ,g(28,19,17)
c	write (9,6302) g(28,15,18),g(28,16,18),g(28,17,18),g(28,18,18)
c    &	   ,g(28,19,18)
c	write (9,6302) g(28,15,19),g(28,16,19),g(28,17,19),g(28,18,19)
c    &	   ,g(28,19,19)

 	close(9)
        open(unit=9,file='tmpfile.txt',status='replace')
        write (9,*) filename
 	close(9)
	return
 6300   format (' ')
 6301   format (i5,5(1x,e18.12))
 6302   format (5(1x,e18.12))
 6500   format (513(1x,e18.12))
	end

c23456789012345678901234567890123456789012345678901234567890123456789012

	subroutine plotnrm(nx,ny,nz,x,y,z,gnrm,iplot,ireplt,pltname)
        implicit double precision (a-h,o-z)
        character*80 filename,str
        character*4 pltname
        dimension x(nx+1),y(ny+1),z(nz+1)
	dimension gnrm(nx,ny,nz)
        common /rparam/small,pi
	nxp=nx+1
	nyp=ny+1
	nzp=nz+1
        if ((iplot.eq.0).or.(ireplt.eq.0)) then
           write (str,'(I10)') iplot
           str=adjustl(str)
           filename=pltname//trim(str)//'.txt'
        else   
           filename=pltname//'.txt'
        endif
        open(unit=9,file=filename,status='replace')
 	write (9,6500) (x(i),i=1,nxp)
 	write (9,6500) (y(j),j=1,nyp)
 	write (9,6500) (z(k),k=1,nzp)

        do 400 i=1,nx
	   write (9,6300)
           do 300 j=1,ny
 	      write (9,6500) (gnrm(i,j,k),k=1,nz),
     &           100.0d0*gnrm(i,j,1)
  300      continue
           write (9,6500) (100.0d0*gnrm(i,1,k),k=1,nz)
  400   continue
	write (9,6300)
        do 500 j=1,ny
 	   write (9,6500) (gnrm(1,j,k),k=1,nz),
     &        100.0d0*gnrm(1,j,1)
  500   continue
        write (9,6500) (100.0d0*gnrm(1,1,k),k=1,nz)
        
 	close(9)
        open(unit=9,file='tmpfile.txt',status='replace')
        write (9,*) filename
 	close(9)
	return
 6300   format (' ')
 6500   format (401(1x,e18.12))
	end

c23456789012345678901234567890123456789012345678901234567890123456789012
	
 	subroutine wenonorm(nx,ny,nz,dx,dy,dz,nord,
     &     gsgn,g,gnrmx,gnrmy,gnrmz,cx_ave,cy_ave,cz_ave,difftype)
        implicit double precision (a-h,o-z)
        character*4 difftype
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gsgn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnrmx(nx,ny,nz),gnrmy(nx,ny,nz),gnrmz(nx,ny,nz)
        dimension cx_ave(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cy_ave(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cz_ave(-2:nx+3,-2:ny+3,-2:nz+3)
        common /rparam/small,pi
        common /surfave/grdave0,grdave1
        grdave=0.0d0
c       eps=1.0d-6
        eps=small
        do 100 k=1,nz
        do 100 j=1,ny
        do 100 i=1,nx
           cxsgn=dsign(1.0d0,cx_ave(i,j,k))
           cysgn=dsign(1.0d0,cy_ave(i,j,k))
           czsgn=dsign(1.0d0,cz_ave(i,j,k))
	   if (cx_ave(i,j,k).eq.0.0d0) cxsgn=1.0d0
	   if (cy_ave(i,j,k).eq.0.0d0) cysgn=1.0d0
	   if (cz_ave(i,j,k).eq.0.0d0) czsgn=1.0d0
           gc=gsgn(i,j,k)

	   if (nord.eq.3) goto 50

           gxc=8.0d0*(g(i+1,j,k)-g(i-1,j,k))-(g(i+2,j,k)-g(i-2,j,k))
           gxc=gxc/12.0d0
           am=g(i-3,j,k)-2.0d0*g(i-2,j,k)+g(i-1,j,k)
           bm=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)
           cm=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
           dm=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)
           s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
           s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
           s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &        +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
           gxm=gxc-gweno

           ap=g(i+1,j,k)-2.0d0*g(i+2,j,k)+g(i+3,j,k)
           bp=dm; dp=bm; cp=cm
c          bp=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)
c          cp=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
c          dp=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)
           s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
           s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
           s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &        +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
           gxp=gxc+gweno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gyc=8.0d0*(g(i,j+1,k)-g(i,j-1,k))-(g(i,j+2,k)-g(i,j-2,k))
           gyc=gyc/12.0d0
           am=g(i,j-3,k)-2.0d0*g(i,j-2,k)+g(i,j-1,k)
           bm=g(i,j-2,k)-2.0d0*g(i,j-1,k)+g(i,j,k)
           cm=g(i,j-1,k)-2.0d0*g(i,j,k)+g(i,j+1,k)
           dm=g(i,j,k)-2.0d0*g(i,j+1,k)+g(i,j+2,k)
           s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
           s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
           s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &        +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
           gym=gyc-gweno

           ap=g(i,j+1,k)-2.0d0*g(i,j+2,k)+g(i,j+3,k)
           bp=dm; dp=bm; cp=cm
           s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
           s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
           s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &        +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
           gyp=gyc+gweno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gzc=8.0d0*(g(i,j,k+1)-g(i,j,k-1))-(g(i,j,k+2)-g(i,j,k-2))
           gzc=gzc/12.0d0
           am=g(i,j,k-3)-2.0d0*g(i,j,k-2)+g(i,j,k-1)
           bm=g(i,j,k-2)-2.0d0*g(i,j,k-1)+g(i,j,k)
           cm=g(i,j,k-1)-2.0d0*g(i,j,k)+g(i,j,k+1)
           dm=g(i,j,k)-2.0d0*g(i,j,k+1)+g(i,j,k+2)
           s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
           s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
           s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &        +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
           gzm=gzc-gweno

           ap=g(i,j,k+1)-2.0d0*g(i,j,k+2)+g(i,j,k+3)
           bp=dm; dp=bm; cp=cm
           s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
           s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
           s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &        +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
           gzp=gzc+gweno

	   goto 75
 50        continue

           gxc=g(i+1,j,k)-g(i-1,j,k)
           gyc=g(i,j+1,k)-g(i,j-1,k)
           gzc=g(i,j,k+1)-g(i,j,k-1)
	   gxm=0.5d0*(gxc-(g(i+1,j,k)
     &        +3.0d0*(g(i-1,j,k)-g(i,j,k))-g(i-2,j,k))/3.0d0)
	   gym=0.5d0*(gyc-(g(i,j+1,k)
     &        +3.0d0*(g(i,j-1,k)-g(i,j,k))-g(i,j-2,k))/3.0d0)
	   gzm=0.5d0*(gzc-(g(i,j,k+1)
     &        +3.0d0*(g(i,j,k-1)-g(i,j,k))-g(i,j,k-2))/3.0d0)
	   gxp=0.5d0*(gxc+(g(i-1,j,k)
     &        +3.0d0*(g(i+1,j,k)-g(i,j,k))-g(i+2,j,k))/3.0d0)
	   gyp=0.5d0*(gyc+(g(i,j-1,k)
     &        +3.0d0*(g(i,j+1,k)-g(i,j,k))-g(i,j+2,k))/3.0d0)
	   gzp=0.5d0*(gzc+(g(i,j,k-1)
     &        +3.0d0*(g(i,j,k+1)-g(i,j,k))-g(i,j,k+2))/3.0d0)

 75        continue

           if (dabs(gxm).lt.small) gxm=0.0d0
           if (dabs(gym).lt.small) gym=0.0d0
           if (dabs(gzm).lt.small) gzm=0.0d0
           if (dabs(gxp).lt.small) gxp=0.0d0
           if (dabs(gyp).lt.small) gyp=0.0d0
           if (dabs(gzp).lt.small) gzp=0.0d0
           gxm=gxm/dx
           gym=gym/dy
           gzm=gzm/dz
           gxp=gxp/dx
           gyp=gyp/dy
           gzp=gzp/dz

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           if (difftype.eq.'evol') then
              sgn=-1.0d0
           elseif (difftype.eq.'rein') then
              if (dabs(gc).le.dx) then
                 sgn=gc/dx+dsin(pi*gc/dx)/pi
              else
                 sgn=dsign(1.0d0,gc)
              endif
           else
              write (6,*) "wenonorm: difftype must be 'evol' or 'rein'"
              stop
           endif

           if (sgn.lt.0.0d0) then ! Propagation down the gradient

              if ((gxm+gxp).gt.0.0d0) then
                 gx=dmax1(0.0d0,gxp)
              elseif ((gxm+gxp).lt.0.0d0) then
                 gx=dmin1(0.0d0,gxm)
              else
		 gx=cxsgn*dmax1(0.0d0,gxp) ! necess. for prop. at local min.
c                gx=0.0d0
              endif
              if ((gym+gyp).gt.0.0d0) then
                 gy=dmax1(0.0d0,gyp)
              elseif ((gym+gyp).lt.0.0d0) then
                 gy=dmin1(0.0d0,gym)
              else
		 gy=cysgn*dmax1(0.0d0,gyp) ! necess. for prop. at local min.
c                gy=0.0d0
              endif
              if ((gzm+gzp).gt.0.0d0) then
                 gz=dmax1(0.0d0,gzp)
              elseif ((gzm+gzp).lt.0.0d0) then
                 gz=dmin1(0.0d0,gzm)
              else
		 gz=czsgn*dmax1(0.0d0,gzp) ! necess. for prop. at local min.
c                gz=0.0d0
              endif

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           else ! Propagation up the gradient

              if ((gxm+gxp).gt.0.0d0) then
                 gx=dmax1(0.0d0,gxm)
              elseif ((gxm+gxp).lt.0.0d0) then
                 gx=dmin1(0.0d0,gxp)
              else
		 gx=cxsgn*dmax1(0.0d0,gxm) ! necess. for prop. at local max.
c                gx=0.0d0
              endif
              if ((gym+gyp).gt.0.0d0) then
                 gy=dmax1(0.0d0,gym)
              elseif ((gym+gyp).lt.0.0d0) then
                 gy=dmin1(0.0d0,gyp)
              else
		 gy=cysgn*dmax1(0.0d0,gym) ! necess. for prop. at local max.
c                gy=0.0d0
              endif
              if ((gzm+gzp).gt.0.0d0) then
                 gz=dmax1(0.0d0,gzm)
              elseif ((gzm+gzp).lt.0.0d0) then
                 gz=dmin1(0.0d0,gzp)
              else
		 gz=czsgn*dmax1(0.0d0,gzm) ! necess. for prop. at local max.
c                gz=0.0d0
              endif

           endif

 	   grad=dsqrt(gx**2.0d0+gy**2.0d0+gz**2.0d0)
 	   if (grad.lt.small) then
c             write (6,*) 'wenonorm: grad0 = 0 at ',i,j,k,
c    &           'gradients= ',gx,gy,gz
 	      gnrmx(i,j,k)=0.0d0
 	      gnrmy(i,j,k)=0.0d0
 	      gnrmz(i,j,k)=0.0d0
 	   else
 	      gnrmx(i,j,k)=gx/grad
 	      gnrmy(i,j,k)=gy/grad
 	      gnrmz(i,j,k)=gz/grad
	   endif
	   if (dabs(gnrmx(i,j,k)).lt.small) gnrmx(i,j,k)=0.0d0
	   if (dabs(gnrmy(i,j,k)).lt.small) gnrmy(i,j,k)=0.0d0
 	   if (dabs(gnrmz(i,j,k)).lt.small) gnrmz(i,j,k)=0.0d0
           grdave=grdave+grad
  100   continue
        grdave0=grdave/dfloat(nx*ny*nz)
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012
	
 	subroutine gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &	   nord,gsgn,g,grad,ibcv,difftype,iflg)
        implicit double precision (a-h,o-z)
        character*4 difftype
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gsgn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension grad(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension indx(nptot,4),mask(nx,ny,nz)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

        eps=small

	grad=0.0d0
        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           gc=gsgn(i,j,k)-gjump*dfloat(mask(i,j,k))
           gxm1=(g(i,j,k)-g(i-1,j,k))/dx
           gxp1=(g(i+1,j,k)-g(i,j,k))/dx
           gym1=(g(i,j,k)-g(i,j-1,k))/dy
           gyp1=(g(i,j+1,k)-g(i,j,k))/dy
           gzm1=(g(i,j,k)-g(i,j,k-1))/dz
           gzp1=(g(i,j,k+1)-g(i,j,k))/dz

           if (dabs(gxm1).lt.small) gxm1=0.0d0
           if (dabs(gym1).lt.small) gym1=0.0d0
           if (dabs(gzm1).lt.small) gzm1=0.0d0
           if (dabs(gxp1).lt.small) gxp1=0.0d0
           if (dabs(gyp1).lt.small) gyp1=0.0d0
           if (dabs(gzp1).lt.small) gzp1=0.0d0

	   gxc=0.5d0*(g(i+1,j,k)-g(i-1,j,k))/dx
	   gyc=0.5d0*(g(i,j+1,k)-g(i,j-1,k))/dy
	   gzc=0.5d0*(g(i,j,k+1)-g(i,j,k-1))/dz

 
c23456789012345678901234567890123456789012345678901234567890123456789012
	
           if (difftype.eq.'evol') then
              sgn=-1.0d0
           elseif (difftype.eq.'rein') then
              if (dabs(gc).le.dx) then
                 sgn=gc/dx+dsin(pi*gc/dx)/pi
              else
                 sgn=dsign(1.0d0,gc)
              endif
           else
              write (6,*) "wenogradts: difftype must be 'evol' or 
     &		 'rein'"
              stop
           endif

           if (sgn.lt.0.0d0) then ! Propagation down the gradient
	      gx2=dmax1(dmax1(gxp1,0.0d0)**2.0d0,
     &		 dmin1(gxm1,0.0d0)**2.0d0)
	      gy2=dmax1(dmax1(gyp1,0.0d0)**2.0d0,
     &		 dmin1(gym1,0.0d0)**2.0d0)
	      gz2=dmax1(dmax1(gzp1,0.0d0)**2.0d0,
     &		 dmin1(gzm1,0.0d0)**2.0d0)
              grad(i,j,k)=dsqrt(gx2+gy2+gz2)

           elseif (sgn.gt.0.0d0) then ! Propagation up the gradient
	      gx2=dmax1(dmin1(gxp1,0.0d0)**2.0d0,
     &		 dmax1(gxm1,0.0d0)**2.0d0)
	      gy2=dmax1(dmin1(gyp1,0.0d0)**2.0d0,
     &		 dmax1(gym1,0.0d0)**2.0d0)
	      gz2=dmax1(dmin1(gzp1,0.0d0)**2.0d0,
     &		 dmax1(gzm1,0.0d0)**2.0d0)
              grad(i,j,k)=dsqrt(gx2+gy2+gz2)

 	   else
              gx2=gxc**2.0d0
              gy2=gyc**2.0d0
              gz2=gzc**2.0d0
              grad(i,j,k)=dsqrt(gx2+gy2+gz2)
           endif

  100   continue
        call bound(nx,ny,nz,grad,ibcv)
        return
        end


c23456789012345678901234567890123456789012345678901234567890123456789012
	
 	subroutine wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &	   nord,gsgn,g,grad,ibcv,difftype,iflg)
        implicit double precision (a-h,o-z)
        character*4 difftype
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gsgn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension grad(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension indx(nptot,4),mask(nx,ny,nz)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

        eps=small

	grad=0.0d0
        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           gc=gsgn(i,j,k)-gjump*dfloat(mask(i,j,k))

           if (nord.eq.3) goto 50

           gxc=8.0d0*(g(i+1,j,k)-g(i-1,j,k))-(g(i+2,j,k)-g(i-2,j,k))
           gxc=gxc/12.0d0
           am=g(i-3,j,k)-2.0d0*g(i-2,j,k)+g(i-1,j,k)
           bm=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)
           cm=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
           dm=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)
           s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
           s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
           s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &        +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
           gxm=gxc-gweno

           ap=g(i+1,j,k)-2.0d0*g(i+2,j,k)+g(i+3,j,k)
           bp=dm; dp=bm; cp=cm
c          bp=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)
c          cp=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
c          dp=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)
           s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
           s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
           s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &        +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
           gxp=gxc+gweno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gyc=8.0d0*(g(i,j+1,k)-g(i,j-1,k))-(g(i,j+2,k)-g(i,j-2,k))
           gyc=gyc/12.0d0
           am=g(i,j-3,k)-2.0d0*g(i,j-2,k)+g(i,j-1,k)
           bm=g(i,j-2,k)-2.0d0*g(i,j-1,k)+g(i,j,k)
           cm=g(i,j-1,k)-2.0d0*g(i,j,k)+g(i,j+1,k)
           dm=g(i,j,k)-2.0d0*g(i,j+1,k)+g(i,j+2,k)
           s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
           s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
           s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &        +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
           gym=gyc-gweno

           ap=g(i,j+1,k)-2.0d0*g(i,j+2,k)+g(i,j+3,k)
           bp=dm; dp=bm; cp=cm
           s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
           s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
           s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &        +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
           gyp=gyc+gweno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gzc=8.0d0*(g(i,j,k+1)-g(i,j,k-1))-(g(i,j,k+2)-g(i,j,k-2))
           gzc=gzc/12.0d0
           am=g(i,j,k-3)-2.0d0*g(i,j,k-2)+g(i,j,k-1)
           bm=g(i,j,k-2)-2.0d0*g(i,j,k-1)+g(i,j,k)
           cm=g(i,j,k-1)-2.0d0*g(i,j,k)+g(i,j,k+1)
           dm=g(i,j,k)-2.0d0*g(i,j,k+1)+g(i,j,k+2)
           s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
           s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
           s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &        +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
           gzm=gzc-gweno

           ap=g(i,j,k+1)-2.0d0*g(i,j,k+2)+g(i,j,k+3)
           bp=dm; dp=bm; cp=cm
           s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
           s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
           s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
           a0=1.0d0/(eps+s0)**2.0d0
           a1=6.0d0/(eps+s1)**2.0d0
           a2=3.0d0/(eps+s2)**2.0d0
           w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
           gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &        +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
           gzp=gzc+gweno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
	   gxc=gxc/dx
	   gyc=gyc/dy
	   gzc=gzc/dz
           goto 75

 50        continue

           gxc=g(i+1,j,k)-g(i-1,j,k)
           gyc=g(i,j+1,k)-g(i,j-1,k)
           gzc=g(i,j,k+1)-g(i,j,k-1)
           gxm=0.5d0*(gxc-(g(i+1,j,k)
     &        +3.0d0*(g(i-1,j,k)-g(i,j,k))-g(i-2,j,k))/3.0d0)
           gym=0.5d0*(gyc-(g(i,j+1,k)
     &        +3.0d0*(g(i,j-1,k)-g(i,j,k))-g(i,j-2,k))/3.0d0)
           gzm=0.5d0*(gzc-(g(i,j,k+1)
     &        +3.0d0*(g(i,j,k-1)-g(i,j,k))-g(i,j,k-2))/3.0d0)
           gxp=0.5d0*(gxc+(g(i-1,j,k)
     &        +3.0d0*(g(i+1,j,k)-g(i,j,k))-g(i+2,j,k))/3.0d0)
           gyp=0.5d0*(gyc+(g(i,j-1,k)
     &        +3.0d0*(g(i,j+1,k)-g(i,j,k))-g(i,j+2,k))/3.0d0)
           gzp=0.5d0*(gzc+(g(i,j,k-1)
     &        +3.0d0*(g(i,j,k+1)-g(i,j,k))-g(i,j,k+2))/3.0d0)
	   gxc=0.5d0*gxc/dx
	   gyc=0.5d0*gyc/dy
	   gzc=0.5d0*gzc/dz

 75        continue

           if (dabs(gxm).lt.small) gxm=0.0d0
           if (dabs(gym).lt.small) gym=0.0d0
           if (dabs(gzm).lt.small) gzm=0.0d0
           if (dabs(gxp).lt.small) gxp=0.0d0
           if (dabs(gyp).lt.small) gyp=0.0d0
           if (dabs(gzp).lt.small) gzp=0.0d0
           gxm=gxm/dx
           gym=gym/dy
           gzm=gzm/dz
           gxp=gxp/dx
           gyp=gyp/dy
           gzp=gzp/dz

           gxp1=(g(i+1,j,k)-g(i,j,k))/dx
           gxm1=(g(i,j,k)-g(i-1,j,k))/dx
           gyp1=(g(i,j+1,k)-g(i,j,k))/dy
           gym1=(g(i,j,k)-g(i,j-1,k))/dy
           gzp1=(g(i,j,k+1)-g(i,j,k))/dz
           gzm1=(g(i,j,k)-g(i,j,k-1))/dz
 
c23456789012345678901234567890123456789012345678901234567890123456789012
	
           if (difftype.eq.'evol') then
              sgn=-1.0d0
           elseif (difftype.eq.'rein') then
              if (dabs(gc).le.dx) then
                 sgn=gc/dx+dsin(pi*gc/dx)/pi
c                sgn=-(gc/dx+dsin(pi*gc/dx)/pi)
              else
                 sgn=dsign(1.0d0,gc)
c                sgn=dsign(-1.0d0,gc)
              endif
           else
              write (6,*) "wenogradts: difftype must be 'evol' or 
     &		 'rein'"
              stop
           endif

           if (sgn.lt.0.0d0) then ! Propagation down the gradient
              gx2=dmax1(dmax1(gxp,0.0d0)**2.0d0,dmin1(gxm,0.0d0)**2.0d0)
 	      if (gx2.eq.0.0d0) gx2=dmax1(dmax1(gxp1,0.0d0)**2.0d0,
     &		 dmin1(gxm1,0.0d0)**2.0d0)
c	      gx2s=dmax1(dmax1(gxp1,0.0d0)**2.0d0
c    &		 ,dmin1(gxm1,0.0d0)**2.0d0)
c	      if (gx2s.eq.0.0d0) gx2=gx2s

              gy2=dmax1(dmax1(gyp,0.0d0)**2.0d0,dmin1(gym,0.0d0)**2.0d0)
 	      if (gy2.eq.0.0d0) gy2=dmax1(dmax1(gyp1,0.0d0)**2.0d0,
     &		 dmin1(gym1,0.0d0)**2.0d0)
c	      gy2s=dmax1(dmax1(gyp1,0.0d0)**2.0d0
c    &		 ,dmin1(gym1,0.0d0)**2.0d0)
c	      if (gy2s.eq.0.0d0) gy2=gy2s

              gz2=dmax1(dmax1(gzp,0.0d0)**2.0d0,dmin1(gzm,0.0d0)**2.0d0)
 	      if (gz2.eq.0.0d0) gz2=dmax1(dmax1(gzp1,0.0d0)**2.0d0,
     &		 dmin1(gzm1,0.0d0)**2.0d0)
c	      gz2s=dmax1(dmax1(gzp1,0.0d0)**2.0d0
c    &		 ,dmin1(gzm1,0.0d0)**2.0d0)
c	      if (gz2s.eq.0.0d0) gz2=gz2s

              grad(i,j,k)=dsqrt(gx2+gy2+gz2)

           elseif (sgn.gt.0.0d0) then ! Propagation up the gradient
              gx2=dmax1(dmin1(gxp,0.0d0)**2.0d0,dmax1(gxm,0.0d0)**2.0d0)
	      if (gx2.eq.0.0d0) gx2=dmax1(dmin1(gxp1,0.0d0)**2.0d0,
     &		 dmax1(gxm1,0.0d0)**2.0d0)

              gy2=dmax1(dmin1(gyp,0.0d0)**2.0d0,dmax1(gym,0.0d0)**2.0d0)
	      if (gy2.eq.0.0d0) gy2=dmax1(dmin1(gyp1,0.0d0)**2.0d0,
     &		 dmax1(gym1,0.0d0)**2.0d0)

              gz2=dmax1(dmin1(gzp,0.0d0)**2.0d0,dmax1(gzm,0.0d0)**2.0d0)
	      if (gz2.eq.0.0d0) gz2=dmax1(dmin1(gzp1,0.0d0)**2.0d0,
     &		 dmax1(gzm1,0.0d0)**2.0d0)

              grad(i,j,k)=dsqrt(gx2+gy2+gz2)

 	   else
              gx2=gxc**2.0d0
              gy2=gyc**2.0d0
              gz2=gzc**2.0d0
              grad(i,j,k)=dsqrt(gx2+gy2+gz2)
           endif

c23456789012345678901234567890123456789012345678901234567890123456789012
	
c	   if ((grad(i,j,k).lt.small)
c    &	      .and.(((sgn*gxm).gt.0.0d0).or.((sgn*gxp).lt.0.0d0) 
c    &	        .or.((sgn*gym).gt.0.0d0).or.((sgn*gyp).lt.0.0d0) 
c    &	        .or.((sgn*gzm).gt.0.0d0).or.((sgn*gzp).lt.0.0d0))) then
c             write (6,2000) i,j,k,mask(i,j,k),indx(nt,4),mstep,
c    &           iflg,nt,sgn,grad(i,j,k),gx2,gy2,gz2
c	      write (6,3000) am,bm,cm,dm,am-2.0d0*bm+cm,bm-2.0d0*cm+dm
cc	      write (6,3000) gxm1,gxp1,gym1,gyp1,gzm1,gzp1
c             grad(i,j,k)=0.0d0
c	   elseif (grad(i,j,k).gt.1.0e2) then
c             write (6,4000) i,j,k,mask(i,j,k),indx(nt,4),mstep,
c    &           iflg,nt,sgn,grad(i,j,k),gx2,gy2,gz2
c          endif

  100   continue
        call bound(nx,ny,nz,grad,ibcv)
 
 2000   format ('ENOGRAD LOW: ',7(1x,i3),1x,i7,5(2x,e18.12))
 3000   format ('ENOGRAD: ',8(2x,e18.12))
 4000   format ('ENOGRAD BIG: ',7(1x,i3),1x,i7,5(2x,e18.12))
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &	   nord,gsgn,g,grad,ibcv,difftype,iflg)
        implicit double precision (a-h,o-z)
        character*4 difftype
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gsgn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension grad(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension indx(nptot,4),mask(nx,ny,nz)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

        eps=small

	grad=0.0d0
        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           gc=gsgn(i,j,k)-gjump*dfloat(mask(i,j,k))

           if (nord.eq.3) goto 50

           gxc=8.0d0*(g(i+1,j,k)-g(i-1,j,k))-(g(i+2,j,k)-g(i-2,j,k))
           gxc=gxc/12.0d0
           am=g(i-3,j,k)-2.0d0*g(i-2,j,k)+g(i-1,j,k)
           bm=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)
           cm=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
           dm=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)
	   phi=-(bm-2.0d0*cm+dm)/12.0d0

	   if ((dabs(bm).lt.dabs(cm))
     &	      .and.(dabs(am-bm).lt.dabs(bm-cm))) then
              geno=phi+(am-2.0d0*bm+cm)/3.0d0
	   elseif ((dabs(bm).gt.dabs(cm))
     &	      .and.(dabs(bm-cm).gt.dabs(cm-dm))) then
              geno=-phi
	   else
              geno=phi
	   endif
           gxm=gxc-geno

           ap=g(i+1,j,k)-2.0d0*g(i+2,j,k)+g(i+3,j,k)
           bp=dm; dp=bm; cp=cm
c          bp=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)
c          cp=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
c          dp=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)

	   phi=-(bp-2.0d0*cp+dp)/12.0d0
	   if ((dabs(bp).lt.dabs(cp))
     &	      .and.(dabs(ap-bp).lt.dabs(bp-cp))) then
              geno=phi+(ap-2.0d0*bp+cp)/3.0d0
	   elseif ((dabs(bp).gt.dabs(cp))
     &	      .and.(dabs(bp-cp).gt.dabs(cp-dp))) then
              geno=-phi
	   else
              geno=phi
	   endif
           gxp=gxc+geno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gyc=8.0d0*(g(i,j+1,k)-g(i,j-1,k))-(g(i,j+2,k)-g(i,j-2,k))
           gyc=gyc/12.0d0
           am=g(i,j-3,k)-2.0d0*g(i,j-2,k)+g(i,j-1,k)
           bm=g(i,j-2,k)-2.0d0*g(i,j-1,k)+g(i,j,k)
           cm=g(i,j-1,k)-2.0d0*g(i,j,k)+g(i,j+1,k)
           dm=g(i,j,k)-2.0d0*g(i,j+1,k)+g(i,j+2,k)

	   phi=-(bm-2.0d0*cm+dm)/12.0d0
	   if ((dabs(bm).lt.dabs(cm))
     &	      .and.(dabs(am-bm).lt.dabs(bm-cm))) then
              geno=phi+(am-2.0d0*bm+cm)/3.0d0
	   elseif ((dabs(bm).gt.dabs(cm))
     &	      .and.(dabs(bm-cm).gt.dabs(cm-dm))) then
              geno=-phi
	   else
              geno=phi
	   endif
           gym=gyc-geno

           ap=g(i,j+1,k)-2.0d0*g(i,j+2,k)+g(i,j+3,k)
           bp=dm; dp=bm; cp=cm

	   phi=-(bp-2.0d0*cp+dp)/12.0d0
	   if ((dabs(bp).lt.dabs(cp))
     &	      .and.(dabs(ap-bp).lt.dabs(bp-cp))) then
              geno=phi+(ap-2.0d0*bp+cp)/3.0d0
	   elseif ((dabs(bp).gt.dabs(cp))
     &	      .and.(dabs(bp-cp).gt.dabs(cp-dp))) then
              geno=-phi
	   else
              geno=phi
	   endif
           gyp=gyc+geno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gzc=8.0d0*(g(i,j,k+1)-g(i,j,k-1))-(g(i,j,k+2)-g(i,j,k-2))
           gzc=gzc/12.0d0
           am=g(i,j,k-3)-2.0d0*g(i,j,k-2)+g(i,j,k-1)
           bm=g(i,j,k-2)-2.0d0*g(i,j,k-1)+g(i,j,k)
           cm=g(i,j,k-1)-2.0d0*g(i,j,k)+g(i,j,k+1)
           dm=g(i,j,k)-2.0d0*g(i,j,k+1)+g(i,j,k+2)

	   phi=-(bm-2.0d0*cm+dm)/12.0d0
	   if ((dabs(bm).lt.dabs(cm))
     &	      .and.(dabs(am-bm).lt.dabs(bm-cm))) then
              geno=phi+(am-2.0d0*bm+cm)/3.0d0
	   elseif ((dabs(bm).gt.dabs(cm))
     &	      .and.(dabs(bm-cm).gt.dabs(cm-dm))) then
              geno=-phi
	   else
              geno=phi
	   endif
           gzm=gzc-geno

           ap=g(i,j,k+1)-2.0d0*g(i,j,k+2)+g(i,j,k+3)
           bp=dm; dp=bm; cp=cm

	   phi=-(bp-2.0d0*cp+dp)/12.0d0
	   if ((dabs(bp).lt.dabs(cp))
     &	      .and.(dabs(ap-bp).lt.dabs(bp-cp))) then
              geno=phi+(ap-2.0d0*bp+cp)/3.0d0
	   elseif ((dabs(bp).gt.dabs(cp))
     &	      .and.(dabs(bp-cp).gt.dabs(cp-dp))) then
              geno=-phi
	   else
              geno=phi
	   endif
           gzp=gzc+geno

c23456789012345678901234567890123456789012345678901234567890123456789012
	
	   gxc=gxc/dx
	   gyc=gyc/dy
	   gzc=gzc/dz
           goto 75

 50        continue

           gxc=g(i+1,j,k)-g(i-1,j,k)
           gyc=g(i,j+1,k)-g(i,j-1,k)
           gzc=g(i,j,k+1)-g(i,j,k-1)
           gxm=0.5d0*(gxc-(g(i+1,j,k)
     &        +3.0d0*(g(i-1,j,k)-g(i,j,k))-g(i-2,j,k))/3.0d0)
           gym=0.5d0*(gyc-(g(i,j+1,k)
     &        +3.0d0*(g(i,j-1,k)-g(i,j,k))-g(i,j-2,k))/3.0d0)
           gzm=0.5d0*(gzc-(g(i,j,k+1)
     &        +3.0d0*(g(i,j,k-1)-g(i,j,k))-g(i,j,k-2))/3.0d0)
           gxp=0.5d0*(gxc+(g(i-1,j,k)
     &        +3.0d0*(g(i+1,j,k)-g(i,j,k))-g(i+2,j,k))/3.0d0)
           gyp=0.5d0*(gyc+(g(i,j-1,k)
     &        +3.0d0*(g(i,j+1,k)-g(i,j,k))-g(i,j+2,k))/3.0d0)
           gzp=0.5d0*(gzc+(g(i,j,k-1)
     &        +3.0d0*(g(i,j,k+1)-g(i,j,k))-g(i,j,k+2))/3.0d0)
	   gxc=0.5d0*gxc/dx
	   gyc=0.5d0*gyc/dy
	   gzc=0.5d0*gzc/dz

 75        continue

           if (dabs(gxm).lt.small) gxm=0.0d0
           if (dabs(gym).lt.small) gym=0.0d0
           if (dabs(gzm).lt.small) gzm=0.0d0
           if (dabs(gxp).lt.small) gxp=0.0d0
           if (dabs(gyp).lt.small) gyp=0.0d0
           if (dabs(gzp).lt.small) gzp=0.0d0
           gxm=gxm/dx
           gym=gym/dy
           gzm=gzm/dz
           gxp=gxp/dx
           gyp=gyp/dy
           gzp=gzp/dz

           gxp1=(g(i+1,j,k)-g(i,j,k))/dx
           gxm1=(g(i,j,k)-g(i-1,j,k))/dx
           gyp1=(g(i,j+1,k)-g(i,j,k))/dy
           gym1=(g(i,j,k)-g(i,j-1,k))/dy
           gzp1=(g(i,j,k+1)-g(i,j,k))/dz
           gzm1=(g(i,j,k)-g(i,j,k-1))/dz
 
c23456789012345678901234567890123456789012345678901234567890123456789012
	
           if (difftype.eq.'evol') then
              sgn=-1.0d0
           elseif (difftype.eq.'rein') then
              if (dabs(gc).le.dx) then
                 sgn=gc/dx+dsin(pi*gc/dx)/pi
              else
                 sgn=dsign(1.0d0,gc)
              endif
           else
              write (6,*) "enogradts: difftype must be 'evol' or 'rein'"
              stop
           endif

           if (sgn.lt.0.0d0) then ! Propagation down the gradient
              gx2=dmax1(dmax1(gxp,0.0d0)**2.0d0,dmin1(gxm,0.0d0)**2.0d0)
	      if (gx2.eq.0.0d0) gx2=dmax1(dmax1(gxp1,0.0d0)**2.0d0,
     &		 dmin1(gxm1,0.0d0)**2.0d0)

              gy2=dmax1(dmax1(gyp,0.0d0)**2.0d0,dmin1(gym,0.0d0)**2.0d0)
	      if (gy2.eq.0.0d0) gy2=dmax1(dmax1(gyp1,0.0d0)**2.0d0,
     &		 dmin1(gym1,0.0d0)**2.0d0)

              gz2=dmax1(dmax1(gzp,0.0d0)**2.0d0,dmin1(gzm,0.0d0)**2.0d0)
	      if (gz2.eq.0.0d0) gz2=dmax1(dmax1(gzp1,0.0d0)**2.0d0,
     &		 dmin1(gzm1,0.0d0)**2.0d0)

              grad(i,j,k)=dsqrt(gx2+gy2+gz2)

           elseif (sgn.gt.0.0d0) then ! Propagation up the gradient
              gx2=dmax1(dmin1(gxp,0.0d0)**2.0d0,dmax1(gxm,0.0d0)**2.0d0)
	      if (gx2.eq.0.0d0) gx2=dmax1(dmin1(gxp1,0.0d0)**2.0d0,
     &		 dmax1(gxm1,0.0d0)**2.0d0)

              gy2=dmax1(dmin1(gyp,0.0d0)**2.0d0,dmax1(gym,0.0d0)**2.0d0)
	      if (gy2.eq.0.0d0) gy2=dmax1(dmin1(gyp1,0.0d0)**2.0d0,
     &		 dmax1(gym1,0.0d0)**2.0d0)

              gz2=dmax1(dmin1(gzp,0.0d0)**2.0d0,dmax1(gzm,0.0d0)**2.0d0)
	      if (gz2.eq.0.0d0) gz2=dmax1(dmin1(gzp1,0.0d0)**2.0d0,
     &		 dmax1(gzm1,0.0d0)**2.0d0)

              grad(i,j,k)=dsqrt(gx2+gy2+gz2)

 	   else
              gx2=gxc**2.0d0
              gy2=gyc**2.0d0
              gz2=gzc**2.0d0
              grad(i,j,k)=dsqrt(gx2+gy2+gz2)
           endif

c23456789012345678901234567890123456789012345678901234567890123456789012
	
 	   if ((grad(i,j,k).lt.small)
     &	      .and.(((sgn*gxm).gt.0.0d0).or.((sgn*gxp).lt.0.0d0) 
     &	        .or.((sgn*gym).gt.0.0d0).or.((sgn*gyp).lt.0.0d0) 
     &	        .or.((sgn*gzm).gt.0.0d0).or.((sgn*gzp).lt.0.0d0))) then
              write (6,2000) i,j,k,mask(i,j,k),indx(nt,4),mstep,
     &           iflg,nt,sgn,grad(i,j,k),gx2,gy2,gz2
	      write (6,3000) am,bm,cm,dm,am-2.0d0*bm+cm,bm-2.0d0*cm+dm,
     &		 geno
              grad(i,j,k)=0.0d0
 	   elseif (grad(i,j,k).gt.1.0e2) then
              write (6,4000) i,j,k,mask(i,j,k),indx(nt,4),mstep,
     &           iflg,nt,sgn,grad(i,j,k),gx2,gy2,gz2
           endif

  100   continue
        call bound(nx,ny,nz,grad,ibcv)
 
 2000   format ('ENOGRAD LOW: ',7(1x,i3),1x,i7,5(2x,e18.12))
 3000   format ('ENOGRAD: ',8(2x,e18.12))
 4000   format ('ENOGRAD BIG: ',7(1x,i3),1x,i7,5(2x,e18.12))
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012
	
 	subroutine interp(nx,ny,nz,nptot,ntc,indx,dep,
     &     distfac,g,gnew,limit)
        implicit double precision (a-h,o-z)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension indx(nptot,4),dep(ntc,3)
        common /rparam/small,pi

 	gnew=0.0d0
 	do 600 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 

c       do 600 k=1,nz
c       do 600 j=1,ny
c       do 600 i=1,nx

           isgn=idnint(dsign(1.0d0,dep(nt,1)))
           jsgn=idnint(dsign(1.0d0,dep(nt,2)))
           ksgn=idnint(dsign(1.0d0,dep(nt,3)))

c          isgn=1;jsgn=1;ksgn=1

           iu=i-isgn
           ju=j-jsgn
           ku=k-ksgn
           iu2=i-2*isgn
           ju2=j-2*jsgn
           ku2=k-2*ksgn
           id=i+isgn
           jd=j+jsgn
           kd=k+ksgn
	   if (dep(nt,1).eq.0.0d0) then
	      iu=i
	      id=i
	      iu2=i
           endif
	   if (dep(nt,2).eq.0.0d0) then
	      ju=j
	      jd=j
	      ju2=j
           endif
	   if (dep(nt,3).eq.0.0d0) then
	      ku=k
	      kd=k
	      ku2=k
           endif
           cxa=dabs(dep(nt,1))*distfac
           cya=dabs(dep(nt,2))*distfac
           cza=dabs(dep(nt,3))*distfac

	   gc=g(i,j,k)
	   gxu=g(iu,j,k)
	   gyu=g(i,ju,k)
	   gzu=g(i,j,ku)
	   gxd=g(id,j,k)
	   gyd=g(i,jd,k)
	   gzd=g(i,j,kd)
	   gxu2=g(iu2,j,k)
	   gyu2=g(i,ju2,k)
	   gzu2=g(i,j,ku2)
	   gxuyu=g(iu,ju,k)
	   gxuyd=g(iu,jd,k)
	   gxuzu=g(iu,j,ku)
	   gxuzd=g(iu,j,kd)
	   gyuxd=g(id,ju,k)
	   gyuzu=g(i,ju,ku)
	   gyuzd=g(i,ju,kd)
	   gzuxd=g(id,j,ku)
	   gzuyd=g(i,jd,ku)
	   guuu=g(iu,ju,ku)

	   grdx=gxd-gxu
	   grdy=gyd-gyu
	   grdz=gzd-gzu
	   grd2x=gxd-2.0d0*gc+gxu
	   grd2y=gyd-2.0d0*gc+gyu
	   grd2z=gzd-2.0d0*gc+gzu
	   grd2xu=gc-2.0d0*gxu+gxu2
	   grd2yu=gc-2.0d0*gyu+gyu2
	   grd2zu=gc-2.0d0*gzu+gzu2
 	   gtemp=-grdx*cxa-grdy*cya-grdz*cza
     &	      +grd2x*cxa**2.0d0+grd2y*cya**2.0d0+grd2z*cza**2.0d0
     &	      +(grd2x-grd2xu)*cxa*(1.0d0-cxa**2.0d0)/3.0d0
     &	      +(grd2y-grd2yu)*cya*(1.0d0-cya**2.0d0)/3.0d0
     &	      +(grd2z-grd2zu)*cza*(1.0d0-cza**2.0d0)/3.0d0
     &	      +(grd2x+gyu-gyuxd+gyd-gxuyd)*cxa*cya
     &	      +(grd2x+gzu-gzuxd+gzd-gxuzd)*cxa*cza
     &	      +(grd2z+gyu-gyuzd+gyd-gzuyd)*cza*cya
     &	      -(grd2x-(gyuxd-2.0d0*gyu+gxuyu))*cya*cxa**2.0d0
     &	      -(grd2y-(gxuyd-2.0d0*gxu+gxuyu))*cxa*cya**2.0d0
     &	      -(grd2x-(gzuxd-2.0d0*gzu+gxuzu))*cza*cxa**2.0d0
     &	      -(grd2z-(gxuzd-2.0d0*gxu+gxuzu))*cxa*cza**2.0d0
     &	      -(grd2y-(gzuyd-2.0d0*gzu+gyuzu))*cza*cya**2.0d0
     &	      -(grd2z-(gyuzd-2.0d0*gyu+gyuzu))*cya*cza**2.0d0
     &	      -2.0d0*(gc-gxu+gxuyu-gyu-(gzu-gxuzu+guuu-gyuzu))
     &        *cxa*cya*cza
	   if (dabs(gtemp).lt.small) gtemp=0.0d0
	   gnew(i,j,k)=gc+0.5d0*gtemp
c          if (limit.eq.1) then
c             gamin=dmin1(gc,gxu,gyu,gzu,gxuyu,gxuzu,gyuzu,guuu)
c             gamax=dmax1(gc,gxu,gyu,gzu,gxuyu,gxuzu,gyuzu,guuu)
c	      gnew(i,j,k)=dmin1(gnew(i,j,k),gamax)
c	      gnew(i,j,k)=dmax1(gnew(i,j,k),gamin)
c          endif
  600	continue
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine gface(gnu,gnc,gnd,gt1unc,gt1dnc,gt1und,
     &     gt2unc,gt2dnc,gt2und,gct12u,cn,ct1,ct2,gf)
        implicit double precision (a-h,o-z)
        common /param4/small,bjump,hlfwdt
        glin=0.5d0*(gnd+gnc)
        grdn=gnd-gnc
        crvn=gnd-2.0d0*gnc+gnu
        grdt1=gnc-gt1unc
        grdt2=gnc-gt2unc
        crvt1=gt1dnc-2.0d0*gnc+gt1unc
        crvt2=gt2dnc-2.0d0*gnc+gt2unc
        gtst1=gnd-gnc+gt1unc-gt1und
        gtst2=gnd-gnc+gt2unc-gt2und
        gtst12=gnc-gt1unc+gct12u-gt2unc
        acn=dabs(cn)
        act1=dabs(ct1)
        act2=dabs(ct2)
        gf=glin-0.5d0*acn*grdn-(1.0d0-acn**2)*crvn/6.0d0
     &     -0.5d0*(act1*grdt1+act2*grdt2)+act1*act2*gtst12/3.0d0
     &     -0.25d0*act1*((1.0d0-act1)*crvt1+(1.0d0-acn)*gtst1)
     &     -0.25d0*act2*((1.0d0-act2)*crvt2+(1.0d0-acn)*gtst2)
        if (dabs(gf).lt.small) gf=0.0d0

c	gf=gnc

        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine glim(gc,gd,gf,gcm1,gcp1,ctm1,ctp1,
     &     gcm2,gcp2,ctm2,ctp2,gdmin,gdmax,gfmin2,gfmax2)
        implicit double precision (a-h,o-z)
        gupmin=gc
        gupmax=gc
        if (ctm1.gt.0.0d0) then
           gupmin=dmin1(gupmin,gcm1)
           gupmax=dmax1(gupmax,gcm1)
        endif
        if (ctp1.lt.0.0d0) then
           gupmin=dmin1(gupmin,gcp1)
           gupmax=dmax1(gupmax,gcp1)
        endif
        if (ctm2.gt.0.0d0) then
           gupmin=dmin1(gupmin,gcm2)
           gupmax=dmax1(gupmax,gcm2)
        endif
        if (ctp2.lt.0.0d0) then
           gupmin=dmin1(gupmin,gcp2)
           gupmax=dmax1(gupmax,gcp2)
        endif
        gfmin1=dmin1(gd,gupmin)
        gfmax1=dmax1(gd,gupmax)
        gdmin=dmin1(gfmin1,gdmin)
        gdmax=dmax1(gfmax1,gdmax)
        gf=dmin1(gf,gfmax1)
        gf=dmax1(gf,gfmin1)
        gfmin2=dmin1(gf,gupmin)
        gfmax2=dmax1(gf,gupmax)
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012
c a la Thuburn (1996)...

 	subroutine interpg2(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
        dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)

        dimension cxl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cxb(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cxo(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cyl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cyb(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cyo(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension czl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension czb(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension czo(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension gmin(0:nx+1,0:ny+1,0:nz+1)
        dimension gmax(0:nx+1,0:ny+1,0:nz+1)
        dimension gl(nx+1,ny+1,nz+1)
        dimension gb(nx+1,ny+1,nz+1)
        dimension go(nx+1,ny+1,nz+1)
        dimension gnxl(nx,ny,nz),gnyl(nx,ny,nz),gnzl(nx,ny,nz)
        dimension gnxb(nx,ny,nz),gnyb(nx,ny,nz),gnzb(nx,ny,nz)
        dimension gnxo(nx,ny,nz),gnyo(nx,ny,nz),gnzo(nx,ny,nz)
        dimension gulmin(0:nx+1,0:ny+1,0:nz+1)
        dimension gulmax(0:nx+1,0:ny+1,0:nz+1)
        dimension gubmin(0:nx+1,0:ny+1,0:nz+1)
        dimension gubmax(0:nx+1,0:ny+1,0:nz+1)
        dimension guomin(0:nx+1,0:ny+1,0:nz+1)
        dimension guomax(0:nx+1,0:ny+1,0:nz+1)
        dimension ctot(nx,ny,nz)

        common /rparam/small,pi
        common /tubedat/twa,twb,gjump,nomask

        cxl=0.0d0; cyl=0.0d0; czl=0.0d0
        cxo=0.0d0; cyo=0.0d0; czo=0.0d0
        cxb=0.0d0; cyb=0.0d0; czb=0.0d0
	nxp=nx+1; nyp=ny+1; nzp=nz+1
        dx=1.0d0/dfloat(nt0)

        do 350 k=1,nzp
        do 350 j=1,nyp
        do 350 i=1,nxp
           if ((i+ishft).gt.0) then
              it0=mod((i+ishft-1),nt0)+1
           else
              it0=nt0+mod((i+ishft),nt0)
           endif
           cxl(i,j,k)=0.5d0*(cx(it0-1,j,k)+cx(it0,j,k))
           cyl(i,j,k)=0.5d0*(cy(it0-1,j,k)+cy(it0,j,k))
           czl(i,j,k)=0.5d0*(cz(it0-1,j,k)+cz(it0,j,k))
           cxo(i,j,k)=0.5d0*(cx(it0,j-1,k)+cx(it0,j,k))
           cyo(i,j,k)=0.5d0*(cy(it0,j-1,k)+cy(it0,j,k))
           czo(i,j,k)=0.5d0*(cz(it0,j-1,k)+cz(it0,j,k))
           cxb(i,j,k)=0.5d0*(cx(it0,j,k-1)+cx(it0,j,k))
           cyb(i,j,k)=0.5d0*(cy(it0,j,k-1)+cy(it0,j,k))
           czb(i,j,k)=0.5d0*(cz(it0,j,k-1)+cz(it0,j,k))
           if (dabs(cxl(i,j,k)).lt.small) cxl(i,j,k)=0.0d0
           if (dabs(cyl(i,j,k)).lt.small) cyl(i,j,k)=0.0d0
           if (dabs(czl(i,j,k)).lt.small) czl(i,j,k)=0.0d0
           if (dabs(cxo(i,j,k)).lt.small) cxo(i,j,k)=0.0d0
           if (dabs(cyo(i,j,k)).lt.small) cyo(i,j,k)=0.0d0
           if (dabs(czo(i,j,k)).lt.small) czo(i,j,k)=0.0d0
           if (dabs(cxb(i,j,k)).lt.small) cxb(i,j,k)=0.0d0
           if (dabs(cyb(i,j,k)).lt.small) cyb(i,j,k)=0.0d0
           if (dabs(czb(i,j,k)).lt.small) czb(i,j,k)=0.0d0
  350   continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 400 k=1,nzp
        do 400 j=1,nyp
        do 400 i=1,nxp
           if (dabs(cxl(i,j,k)).gt.small) then
              sgncxl=dsign(1.0d0,cxl(i,j,k))
           else
              ctmp=cxl(i,j+1,k)+cxl(i,j-1,k)+cxl(i,j,k+1)
     &		 +cxl(i,j,k-1)
              sgncxl=dsign(1.0d0,ctmp)
           endif
           if (dabs(cyl(i,j,k)).gt.small) then
              sgncyl=dsign(1.0d0,cyl(i,j,k))
           else
              ctmp=cyl(i,j+1,k)+cyl(i,j-1,k)+cyl(i,j,k+1)
     &		 +cyl(i,j,k-1)
              sgncyl=dsign(1.0d0,ctmp)
           endif
           if (dabs(czl(i,j,k)).gt.small) then
              sgnczl=dsign(1.0d0,czl(i,j,k))
           else
              ctmp=czl(i,j+1,k)+czl(i,j-1,k)+czl(i,j,k+1)
     &		 +czl(i,j,k-1)
              sgnczl=dsign(1.0d0,ctmp)
           endif
           iu=i-dnint(0.5d0+1.5d0*sgncxl)
           ic=i-dnint(0.5d0+0.5d0*sgncxl)
           id=i-dnint(0.5d0-0.5d0*sgncxl)
           ju=j-dnint(sgncyl)
           jd=j+dnint(sgncyl)
           ku=k-dnint(sgnczl)
           kd=k+dnint(sgnczl)
           giu=g(iu,j,k)
           gic=g(ic,j,k)
           gid=g(id,j,k)
           gicju=g(ic,ju,k)
           gicjd=g(ic,jd,k)
           gidju=g(id,ju,k)
           gicku=g(ic,j,ku)
           gickd=g(ic,j,kd)
           gidku=g(id,j,ku)
           gcjuku=g(ic,ju,ku)
           call gface(giu,gic,gid,gicju,gicjd,gidju,gicku,gickd,gidku,
     &        gcjuku,cxl(i,j,k),cyl(i,j,k),czl(i,j,k),gl(i,j,k))
           call glim(gic,gid,gl(i,j,k),
     &        g(ic,j-1,k),g(ic,j+1,k),cyo(ic,j,k),cyo(ic,j+1,k),
     &        g(ic,j,k-1),g(ic,j,k+1),czb(ic,j,k),czb(ic,j,k+1),
     &        gmin(id,j,k),gmax(id,j,k),gulmin(i,j,k),gulmax(i,j,k))

c23456789012345678901234567890123456789012345678901234567890123456789012

           if (dabs(cxo(i,j,k)).gt.small) then
              sgncxo=dsign(1.0d0,cxo(i,j,k))
           else
              ctmp=cxo(i+1,j,k)+cxo(i-1,j,k)+cxo(i,j,k+1)
     &		 +cxo(i,j,k-1)
              sgncxo=dsign(1.0d0,ctmp)
           endif
           if (dabs(cyo(i,j,k)).gt.small) then
              sgncyo=dsign(1.0d0,cyo(i,j,k))
           else
              ctmp=cyo(i+1,j,k)+cyo(i-1,j,k)+cyo(i,j,k+1)
     &		 +cyo(i,j,k-1)
              sgncyo=dsign(1.0d0,ctmp)
           endif
           if (dabs(czo(i,j,k)).gt.small) then
              sgnczo=dsign(1.0d0,czo(i,j,k))
           else
              ctmp=czo(i+1,j,k)+czo(i-1,j,k)+czo(i,j,k+1)
     &		 +czo(i,j,k-1)
              sgnczo=dsign(1.0d0,ctmp)
           endif
           ju=j-dnint(0.5d0+1.5d0*sgncyo)
           jc=j-dnint(0.5d0+0.5d0*sgncyo)
           jd=j-dnint(0.5d0-0.5d0*sgncyo)
           ku=k-dnint(sgnczo)
           kd=k+dnint(sgnczo)
           iu=i-dnint(sgncxo)
           id=i+dnint(sgncxo)
           gju=g(i,ju,k)
           gjc=g(i,jc,k)
           gjd=g(i,jd,k)
           giujc=g(iu,jc,k)
           gidjc=g(id,jc,k)
           giujd=g(iu,jd,k)
           gjcku=g(i,jc,ku)
           gjckd=g(i,jc,kd)
           gjdku=g(i,jd,ku)
           gciuku=g(iu,jc,ku)
           call gface(gju,gjc,gjd,giujc,gidjc,giujd,gjcku,gjckd,gjdku,
     &        gciuku,cyo(i,j,k),cxo(i,j,k),czo(i,j,k),go(i,j,k))
           call glim(gjc,gjd,go(i,j,k),
     &        g(i-1,jc,k),g(i+1,jc,k),cxl(i,jc,k),cxl(i+1,jc,k),
     &        g(i,jc,k-1),g(i,jc,k+1),czb(i,jc,k),czb(i,jc,k+1),
     &        gmin(i,jd,k),gmax(i,jd,k),guomin(i,j,k),guomax(i,j,k))

c23456789012345678901234567890123456789012345678901234567890123456789012

           if (dabs(cxb(i,j,k)).gt.small) then
              sgncxb=dsign(1.0d0,cxb(i,j,k))
           else
              ctmp=cxb(i+1,j,k)+cxb(i-1,j,k)+cxb(i,j+1,k)
     &		 +cxb(i,j-1,k)
              sgncxb=dsign(1.0d0,ctmp)
           endif
           if (dabs(cyb(i,j,k)).gt.small) then
              sgncyb=dsign(1.0d0,cyb(i,j,k))
           else
              ctmp=cyb(i+1,j,k)+cyb(i-1,j,k)+cyb(i,j+1,k)
     &		 +cyb(i,j-1,k)
              sgncyb=dsign(1.0d0,ctmp)
           endif
           if (dabs(czb(i,j,k)).gt.small) then
              sgnczb=dsign(1.0d0,czb(i,j,k))
           else
              ctmp=czb(i+1,j,k)+czb(i-1,j,k)+czb(i,j+1,k)
     &		 +czb(i,j-1,k)
              sgnczb=dsign(1.0d0,ctmp)
           endif
           ku=k-dnint(0.5d0+1.5d0*sgnczb)
           kc=k-dnint(0.5d0+0.5d0*sgnczb)
           kd=k-dnint(0.5d0-0.5d0*sgnczb)
           iu=i-dnint(sgncxb)
           id=i+dnint(sgncxb)
           ju=j-dnint(sgncyb)
           jd=j+dnint(sgncyb)
           gku=g(i,j,ku)
           gkc=g(i,j,kc)
           gkd=g(i,j,kd)
           giukc=g(iu,j,kc)
           gidkc=g(id,j,kc)
           giukd=g(iu,j,kd)
           gjukc=g(ic,j,ku)
           gjdkc=g(ic,j,kd)
           gjukd=g(id,j,ku)
           gciuju=g(ic,ju,kc)
           call gface(gku,gkc,gkd,giukc,gidkc,giukd,gjukc,gjdkc,gjukd,
     &        gciuju,czb(i,j,k),cxb(i,j,k),cyb(i,j,k),gb(i,j,k))
           call glim(gkc,gkd,gb(i,j,k),
     &        g(i-1,j,kc),g(i+1,j,kc),cxl(i,j,kc),cxl(i+1,j,kc),
     &        g(i,j-1,kc),g(i,j+1,kc),cyo(i,j,kc),cyo(i,j+1,kc),
     &        gmin(i,j,kd),gmax(i,j,kd),gubmin(i,j,k),gubmax(i,j,k))
  400   continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 450 k=1,nz
        do 450 j=1,ny
        do 450 i=1,nx
           clxi=dabs(cxl(i,j,k))+cxl(i,j,k)
           coyi=dabs(cyo(i,j,k))+cyo(i,j,k)
           cbzi=dabs(czb(i,j,k))+czb(i,j,k)
           clxo=dabs(cxl(i,j,k))-cxl(i,j,k)
           coyo=dabs(cyo(i,j,k))-cyo(i,j,k)
           cbzo=dabs(czb(i,j,k))-czb(i,j,k)
           crxi=dabs(cxl(i+1,j,k))-cxl(i+1,j,k)
           ciyi=dabs(cyo(i,j+1,k))-cyo(i,j+1,k)
           ctzi=dabs(czb(i,j,k+1))-czb(i,j,k+1)
           crxo=dabs(cxl(i+1,j,k))+cxl(i+1,j,k)
           ciyo=dabs(cyo(i,j+1,k))+cyo(i,j+1,k)
           ctzo=dabs(czb(i,j,k+1))+czb(i,j,k+1)
           citot=0.5d0*(clxi+crxi+coyi+ciyi+cbzi+ctzi)
           cotot=0.5d0*(clxo+crxo+coyo+ciyo+cbzo+ctzo)
           ctot(i,j,k)=citot-cotot
           flowin=0.5d0*(clxi*gulmin(i,j,k)+crxi*gulmin(i+1,j,k)
     &        +coyi*guomin(i,j,k)+ciyi*guomin(i,j+1,k)
     &        +cbzi*gubmin(i,j,k)+ctzi*gubmin(i,j,k+1))
           goutmx=(g(i,j,k)+flowin-gmin(i,j,k)*(1.0d0+ctot(i,j,k)))
     &        /cotot
           flowin=0.5d0*(clxi*gulmax(i,j,k)+crxi*gulmax(i+1,j,k)
     &        +coyi*guomax(i,j,k)+ciyi*guomax(i,j+1,k)
     &        +cbzi*gubmax(i,j,k)+ctzi*gubmax(i,j,k+1))
           goutmn=(g(i,j,k)+flowin-gmax(i,j,k)*(1.0d0+ctot(i,j,k)))
     &        /cotot

c23456789012345678901234567890123456789012345678901234567890123456789012

           if (cxl(i,j,k).lt.0.0d0) then
              gl(i,j,k)=dmin1(gl(i,j,k),goutmx)
              gl(i,j,k)=dmax1(gl(i,j,k),goutmn)
           endif
           if (cxl(i+1,j,k).gt.0.0d0) then
              gl(i+1,j,k)=dmin1(gl(i+1,j,k),goutmx)
              gl(i+1,j,k)=dmax1(gl(i+1,j,k),goutmn)
           endif
           if (cyo(i,j,k).lt.0.0d0) then
              go(i,j,k)=dmin1(go(i,j,k),goutmx)
              go(i,j,k)=dmax1(go(i,j,k),goutmn)
           endif
           if (cyo(i,j+1,k).gt.0.0d0) then
              go(i,j+1,k)=dmin1(go(i,j+1,k),goutmx)
              go(i,j+1,k)=dmax1(go(i,j+1,k),goutmn)
           endif
           if (czb(i,j,k).lt.0.0d0) then
              gb(i,j,k)=dmin1(gb(i,j,k),goutmx)
              gb(i,j,k)=dmax1(gb(i,j,k),goutmn)
           endif
           if (czb(i,j,k+1).gt.0.0d0) then
              gb(i,j,k+1)=dmin1(gb(i,j,k+1),goutmx)
              gb(i,j,k+1)=dmax1(gb(i,j,k+1),goutmn)
           endif
  450   continue

c       do 500 k=1,nz
c       do 500 j=1,ny
c       do 500 i=1,nx
        do 500 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 500
           g(i,j,k)=(g(i,j,k)
     &        +cxl(i,j,k)*gl(i,j,k)-cxl(i+1,j,k)*gl(i+1,j,k)
     &        +cyo(i,j,k)*go(i,j,k)-cyo(i,j+1,k)*go(i,j+1,k)
     &        +czb(i,j,k)*gb(i,j,k)-czb(i,j,k+1)*gb(i,j,k+1))
     &        /(1.0d0+ctot(i,j,k))
c          g(i,j,k)=g(i,j,k)+dt*sl*gradg_ave(i,j,k)
 	   g(i,j,k)=g(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
           if (dabs(g(i,j,k)).lt.small) g(i,j,k)=0.0d0
  500   continue

c       do 600 k=1,nz
c       do 600 j=1,ny
c       do 600 i=1,nx
cc      do 600 nt=1,ntc
cc         i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
cc         if (indx(nt,4).eq.2) then
c             giso=gjump*dfloat(mask(i,j,k))
c             if ((dabs(g(i,j,k)-giso).gt.twb).or.
c    &           (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
c                g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
c             endif
cc         endif
c 600   continue

        call bound(nx,ny,nz,g,ibcg)
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine interpg3(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)

        dimension flux1(nx,ny,nz)
        dimension flux2(nx,ny,nz)
        dimension flux3(nx,ny,nz)

        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

	flux1=0.0d0; flux2=0.0d0; flux3=0.0d0
        do 50 k=-2,nz+3
        do 50 j=-2,ny+3
        do 50 i=-2,nx+3
	   gnew(i,j,k)=g(i,j,k)
  50    continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 100

	   intsgn=idnint(dsign(1.0d0,dep(nt,1)))
           ca=dabs(dep(nt,1))*distfac
	   index=i

           indxu=index-intsgn
           if (ca.lt.small) indxu=index
	   flux1(i,j,k)=ca*(g(i,j,k)-g(indxu,j,k))
 	   gnew(i,j,k)=g(i,j,k)-flux1(i,j,k)
  100   continue
c       do 150 nt=1,ntc
c          i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c   	   g(i,j,k)=gnew(i,j,k)
c 150   continue

        do 200 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 200

	   intsgn=idnint(dsign(1.0d0,dep(nt,2)))
           ca=dabs(dep(nt,2))*distfac
	   index=j

           indxu=index-intsgn
           if (ca.lt.small) indxu=index
	   flux2(i,j,k)=ca*(gnew(i,j,k)-gnew(i,indxu,k))
 	   gnew(i,j,k)=g(i,j,k)-flux1(i,j,k)-flux2(i,j,k)
c	   flux2(i,j,k)=ca*(g(i,j,k)-g(i,indxu,k))
c	   gnew(i,j,k)=g(i,j,k)+ca*(g(i,j,k)-g(i,indxu,k))
  200   continue
c       do 250 nt=1,ntc
c          i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c	   g(i,j,k)=gnew(i,j,k)
c 250   continue

        do 300 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 300

	   intsgn=idnint(dsign(1.0d0,dep(nt,3)))
           ca=dabs(dep(nt,3))*distfac
	   index=k

           indxu=index-intsgn
           if (ca.lt.small) indxu=index
	   flux3(i,j,k)=ca*(gnew(i,j,k)-gnew(i,j,indxu))
 	   g(i,j,k)=g(i,j,k)-flux1(i,j,k)-flux2(i,j,k)-flux3(i,j,k)
     &	      +dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
c	   flux3(i,j,k)=ca*(g(i,j,k)-g(i,j,indxu))
c	   gnew(i,j,k)=g(i,j,k)+ca*(g(i,j,k)-g(i,j,indxu))
c	   gnew(i,j,k)=gnew(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
c	   g(i,j,k)=g(i,j,k)+flux1(i,j,k)+flux2(i,j,k)+flux3(i,j,k)
c    &	      +dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
  300   continue
c       do 350 nt=1,ntc
c          i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c	   g(i,j,k)=gnew(i,j,k)
c 350   continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 400 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).eq.2) then
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
           endif
  400   continue

        call bound(nx,ny,nz,g,ibcg)
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine interpg4(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)

        dimension cxl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cxb(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cxo(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cyl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cyb(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cyo(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension czl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension czb(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension czo(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension gl(nx+1,ny+1,nz+1)
        dimension gb(nx+1,ny+1,nz+1)
        dimension go(nx+1,ny+1,nz+1)

        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

        cxl=0.0d0; cyl=0.0d0; czl=0.0d0
        cxo=0.0d0; cyo=0.0d0; czo=0.0d0
        cxb=0.0d0; cyb=0.0d0; czb=0.0d0
	gl=0.0d0; gb=0.0d0; go=0.0d0 
        nxp=nx+1; nyp=ny+1; nzp=nz+1

        do 25 k=1,nzp
        do 25 j=1,nyp
        do 25 i=1,nxp
           if ((i+ishft).gt.0) then
              it0=mod((i+ishft-1),nt0)+1
           else
              it0=nt0+mod((i+ishft),nt0)
           endif
           cxl(i,j,k)=0.5d0*(cx(it0-1,j,k)+cx(it0,j,k))
           cyl(i,j,k)=0.5d0*(cy(it0-1,j,k)+cy(it0,j,k))
           czl(i,j,k)=0.5d0*(cz(it0-1,j,k)+cz(it0,j,k))
           cxo(i,j,k)=0.5d0*(cx(it0,j-1,k)+cx(it0,j,k))
           cyo(i,j,k)=0.5d0*(cy(it0,j-1,k)+cy(it0,j,k))
           czo(i,j,k)=0.5d0*(cz(it0,j-1,k)+cz(it0,j,k))
           cxb(i,j,k)=0.5d0*(cx(it0,j,k-1)+cx(it0,j,k))
           cyb(i,j,k)=0.5d0*(cy(it0,j,k-1)+cy(it0,j,k))
           czb(i,j,k)=0.5d0*(cz(it0,j,k-1)+cz(it0,j,k))
           if (dabs(cxl(i,j,k)).lt.small) cxl(i,j,k)=0.0d0
           if (dabs(cyl(i,j,k)).lt.small) cyl(i,j,k)=0.0d0
           if (dabs(czl(i,j,k)).lt.small) czl(i,j,k)=0.0d0
           if (dabs(cxo(i,j,k)).lt.small) cxo(i,j,k)=0.0d0
           if (dabs(cyo(i,j,k)).lt.small) cyo(i,j,k)=0.0d0
           if (dabs(czo(i,j,k)).lt.small) czo(i,j,k)=0.0d0
           if (dabs(cxb(i,j,k)).lt.small) cxb(i,j,k)=0.0d0
           if (dabs(cyb(i,j,k)).lt.small) cyb(i,j,k)=0.0d0
           if (dabs(czb(i,j,k)).lt.small) czb(i,j,k)=0.0d0
  25    continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 50 k=-2,nz+3
        do 50 j=-2,ny+3
        do 50 i=-2,nx+3
	   gnew(i,j,k)=g(i,j,k)
  50    continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 100

	   qp=(0.5d0-cxl(i,j,k))/(1.0d0+cxl(i+1,j,k)-cxl(i,j,k))
	   qm=(0.5d0+cxl(i,j,k))/(1.0d0+cxl(i,j,k)-cxl(i-1,j,k))
	   if (cxl(i,j,k)>0.0d0) then
 	      flux1=cxl(i,j,k)*g(i-1,j,k)
c	      flux1=cx(i-1,j,k)*g(i-1,j,k)
	   else
 	      flux1=cxl(i,j,k)*g(i,j,k)
c	      flux1=cx(i,j,k)*g(i,j,k)
	   endif
 	   flux1=flux1-0.125d0*(g(i,j,k)-g(i-1,j,k))
	   if (cxl(i+1,j,k)>0.0d0) then
 	      flux2=cxl(i+1,j,k)*g(i,j,k)
c	      flux2=cx(i,j,k)*g(i,j,k)
	   else
 	      flux2=cxl(i+1,j,k)*g(i+1,j,k)
c	      flux2=cx(i+1,j,k)*g(i+1,j,k)
	   endif
 	   flux2=flux2-0.125d0*(g(i+1,j,k)-g(i,j,k))
	   gnew(i,j,k)=g(i,j,k)+flux1-flux2
  100   continue
        do 150 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
    	   g(i,j,k)=gnew(i,j,k)
  150   continue

        do 200 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 200

	   if (cyo(i,j,k)>0.0d0) then
 	      flux1=cyo(i,j,k)*g(i,j-1,k)
c	      flux1=cy(i,j-1,k)*g(i,j-1,k)
	   else
 	      flux1=cyo(i,j,k)*g(i,j,k)
c	      flux1=cy(i,j,k)*g(i,j,k)
	   endif
 	   flux1=flux1-0.125d0*(g(i,j,k)-g(i,j-1,k))
	   if (cyo(i,j+1,k)>0.0d0) then
 	      flux2=cyo(i,j+1,k)*g(i,j,k)
c	      flux2=cy(i,j,k)*g(i,j,k)
	   else
 	      flux2=cyo(i,j+1,k)*g(i,j+1,k)
c	      flux2=cy(i,j+1,k)*g(i,j+1,k)
	   endif
 	   flux2=flux2-0.125d0*(g(i,j+1,k)-g(i,j,k))
	   gnew(i,j,k)=g(i,j,k)+flux1-flux2
  200   continue
        do 250 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
    	   g(i,j,k)=gnew(i,j,k)
  250   continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 300 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).ne.2) goto 300

	   if (czb(i,j,k)>0.0d0) then
 	      flux1=czb(i,j,k)*g(i,j,k-1)
c	      flux1=cz(i,j,k-1)*g(i,j,k-1)
	   else
 	      flux1=czb(i,j,k)*g(i,j,k)
c	      flux1=cz(i,j,k)*g(i,j,k)
	   endif
 	   flux1=flux1-0.125d0*(g(i,j,k)-g(i,j,k-1))
	   if (czb(i,j,k+1)>0.0d0) then
 	      flux2=czb(i,j,k+1)*g(i,j,k)
c	      flux2=cz(i,j,k)*g(i,j,k)
	   else
 	      flux2=czb(i,j,k+1)*g(i,j,k+1)
c	      flux2=cz(i,j,k+1)*g(i,j,k+1)
	   endif
 	   flux2=flux2-0.125d0*(g(i,j,k+1)-g(i,j,k))
	   gnew(i,j,k)=g(i,j,k)+flux1-flux2
 	   gnew(i,j,k)=gnew(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
  300   continue
        do 350 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
    	   g(i,j,k)=gnew(i,j,k)
  350   continue

        do 400 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).eq.2) then
c	   g(i,j,k)=g(i,j,k)
c    &	      +0.125d0*(g(i+1,j,k)-2.0d0*g(i,j,k)-g(i-1,j,k)
c    &	      +g(i,j+1,k)-2.0d0*g(i,j,k)-g(i,j-1,k)
c    &	      +g(i,j,k+1)-2.0d0*g(i,j,k)-g(i,j,k-1))
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
           endif
  400   continue

        call bound(nx,ny,nz,g,ibcg)
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012
c a la FCT...

 	subroutine interpg5(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)

        dimension cxa(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cya(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cza(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension flux(-2:nx+3,-2:ny+3,-2:nz+3)

        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

	flux=0.0d0; cxa=0.0d0; cya=0.0d0; cza=0.0d0
        nxp=nx+1; nyp=ny+1; nzp=nz+1

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 25 k=0,nzp
        do 25 j=0,nyp
        do 25 i=0,nxp
           if ((i+ishft).ge.0) then
              it0=mod((i+ishft-1),nt0)+1
           else
              it0=nt0+mod((i+ishft),nt0)
           endif
	   cabs=dabs(cx(it0,j,k))/2.0d0
	   if (cx(it0,j,k).ge.0.0d0) then
	      cxa(i,j,k)=(1.0d0-cabs)*cx(it0,j,k)+cabs*cx(it0+1,j,k)
	   else
	      cxa(i,j,k)=(1.0d0-cabs)*cx(it0,j,k)+cabs*cx(it0-1,j,k)
	   endif

c	   cxa(i,j,k)=cx(it0,j,k)

	   cabs=dabs(cy(it0,j,k))/2.0d0
	   if (cy(it0,j,k).ge.0.0d0) then
	      cya(i,j,k)=(1.0d0-cabs)*cy(it0,j,k)+cabs*cy(it0,j+1,k)
	   else
	      cya(i,j,k)=(1.0d0-cabs)*cy(it0,j,k)+cabs*cy(it0,j-1,k)
	   endif

c	   cya(i,j,k)=cy(it0,j,k)

	   cabs=dabs(cz(it0,j,k))/2.0d0
	   if (cz(it0,j,k).ge.0.0d0) then
	      cza(i,j,k)=(1.0d0-cabs)*cz(it0,j,k)+cabs*cz(it0,j,k+1)
	   else
	      cza(i,j,k)=(1.0d0-cabs)*cz(it0,j,k)+cabs*cz(it0,j,k-1)
	   endif

c	   cza(i,j,k)=cz(it0,j,k)

           if (dabs(cxa(i,j,k)).lt.small) cxa(i,j,k)=0.0d0
           if (dabs(cya(i,j,k)).lt.small) cya(i,j,k)=0.0d0
           if (dabs(cza(i,j,k)).lt.small) cza(i,j,k)=0.0d0
  25    continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 50 k=-2,nz+3
        do 50 j=-2,ny+3
        do 50 i=-2,nx+3
	   gnew(i,j,k)=g(i,j,k)
  50    continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 100 nt=1,ntc
c          if (indx(nt,4).ne.2) goto 100
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)

	   gc=g(i,j,k); gm=g(i-1,j,k); gp=g(i+1,j,k)
	   cc=cxa(i,j,k); cm=cxa(i-1,j,k); cp=cxa(i+1,j,k)

	   qp=(0.5d0-cc)/(1.0d0+cp-cc)
	   qm=(0.5d0+cc)/(1.0d0+cc-cm)
	   gnew(i,j,k)=0.5d0*((gp-gc)*qp**2.0d0-(gc-gm)*qm**2.0d0)
     &	      +(qm+qp)*gc
  100   continue

        do 125 k=0,nz
        do 125 j=0,ny
        do 125 i=0,nx
	   diff0=gnew(i,j,k)-gnew(i-1,j,k)
	   diff1=gnew(i+1,j,k)-gnew(i,j,k)
	   diff2=gnew(i+2,j,k)-gnew(i+1,j,k)
	   sgn=dsign(1.0d0,diff1)
	   flxmin=dmin1(diff0*sgn,0.125d0*dabs(diff1),diff2*sgn)
	   flux(i,j,k)=sgn*dmax1(0.0d0,flxmin)
  125   continue

        do 150 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 150

 	   gnew(i,j,k)=gnew(i,j,k)+flux(i-1,j,k)-flux(i,j,k)
    	   g(i,j,k)=gnew(i,j,k)
  150   continue

        do 200 nt=1,ntc
c          if (indx(nt,4).ne.2) goto 200
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)

           gc=g(i,j,k); gm=g(i,j-1,k); gp=g(i,j+1,k)
           cc=cya(i,j,k); cm=cya(i,j-1,k); cp=cya(i,j+1,k)

           qp=(0.5d0-cc)/(1.0d0+cp-cc)
           qm=(0.5d0+cc)/(1.0d0+cc-cm)
           gnew(i,j,k)=0.5d0*((gp-gc)*qp**2.0d0-(gc-gm)*qm**2.0d0)
     &        +(qm+qp)*gc
  200   continue

        do 225 k=0,nz
        do 225 j=0,ny
        do 225 i=0,nx
           diff0=gnew(i,j,k)-gnew(i,j-1,k)
           diff1=gnew(i,j+1,k)-gnew(i,j,k)
           diff2=gnew(i,j+2,k)-gnew(i,j+1,k)
           sgn=dsign(1.0d0,diff1)
           flxmin=dmin1(diff0*sgn,0.125d0*dabs(diff1),diff2*sgn)
           flux(i,j,k)=sgn*dmax1(0.0d0,flxmin)
  225   continue

        do 250 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 250

           gnew(i,j,k)=gnew(i,j,k)+flux(i,j-1,k)-flux(i,j,k)
           g(i,j,k)=gnew(i,j,k)
  250   continue

        do 300 nt=1,ntc
c          if (indx(nt,4).ne.2) goto 300
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)

           gc=g(i,j,k); gm=g(i,j,k-1); gp=g(i,j,k+1)
           cc=cza(i,j,k); cm=cza(i,j,k-1); cp=cza(i,j,k+1)

           qp=(0.5d0-cc)/(1.0d0+cp-cc)
           qm=(0.5d0+cc)/(1.0d0+cc-cm)
           gnew(i,j,k)=0.5d0*((gp-gc)*qp**2.0d0-(gc-gm)*qm**2.0d0)
     &        +(qm+qp)*gc
  300   continue

        do 325 k=0,nz
        do 325 j=0,ny
        do 325 i=0,nx
           diff0=gnew(i,j,k)-gnew(i,j,k-1)
           diff1=gnew(i,j,k+1)-gnew(i,j,k)
           diff2=gnew(i,j,k+2)-gnew(i,j,k+1)
           sgn=dsign(1.0d0,diff1)
           flxmin=dmin1(diff0*sgn,0.125d0*dabs(diff1),diff2*sgn)
           flux(i,j,k)=sgn*dmax1(0.0d0,flxmin)
  325   continue

        do 350 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 350

           gnew(i,j,k)=gnew(i,j,k)+flux(i,j,k-1)-flux(i,j,k)
c          gnew(i,j,k)=gnew(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
           gnew(i,j,k)=gnew(i,j,k)+dt*sl*gradg_ave(i,j,k)
           g(i,j,k)=gnew(i,j,k)
  350   continue

        do 400 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).eq.2) then
cc	   g(i,j,k)=g(i,j,k)
cc   &	      +0.125d0*(g(i+1,j,k)-2.0d0*g(i,j,k)-g(i-1,j,k)
cc   &	      +g(i,j+1,k)-2.0d0*g(i,j,k)-g(i,j-1,k)
cc   &	      +g(i,j,k+1)-2.0d0*g(i,j,k)-g(i,j,k-1))
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
c          endif
  400   continue

        call bound(nx,ny,nz,g,ibcg)
	return
	end


c23456789012345678901234567890123456789012345678901234567890123456789012
c a la FCT...(donor cell for the low-order method)

 	subroutine interpg6(nx,ny,nz,nptot,ntc,mask,indx,dep,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)

        dimension cxa(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cya(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension cza(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension flux(-2:nx+3,-2:ny+3,-2:nz+3)

        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

	flux=0.0d0; cxa=0.0d0; cya=0.0d0; cza=0.0d0
        nxp=nx+1; nyp=ny+1; nzp=nz+1

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 25 k=0,nzp
        do 25 j=0,nyp
        do 25 i=0,nxp
           if ((i+ishft).ge.0) then
              it0=mod((i+ishft-1),nt0)+1
           else
              it0=nt0+mod((i+ishft),nt0)
           endif
 	   cxa(i,j,k)=cx(it0,j,k)
 	   cya(i,j,k)=cy(it0,j,k)
 	   cza(i,j,k)=cz(it0,j,k)
           if (dabs(cxa(i,j,k)).lt.small) cxa(i,j,k)=0.0d0
           if (dabs(cya(i,j,k)).lt.small) cya(i,j,k)=0.0d0
           if (dabs(cza(i,j,k)).lt.small) cza(i,j,k)=0.0d0
  25    continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 50 k=-2,nz+3
        do 50 j=-2,ny+3
        do 50 i=-2,nx+3
	   gnew(i,j,k)=g(i,j,k)
  50    continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 100 nt=1,ntc
c          if (indx(nt,4).ne.2) goto 100
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)

	   gc=g(i,j,k); gm=g(i-1,j,k); gp=g(i+1,j,k)
	   cc=cxa(i,j,k); cm=cxa(i-1,j,k); cp=cxa(i+1,j,k)

	   if (cc.ge.0.0d0) then
c	      gnew(i,j,k)=gc+cm*gm-cc*gc
	      gnew(i,j,k)=gc-cc*(gc-gm)
	   else
c	      gnew(i,j,k)=gc+cc*gc-cp*gp
	      gnew(i,j,k)=gc-cc*(gp-gc)
	   endif

	   eps=dabs(cc)
	   eta=0.5d0*eps*(1.0d0-eps)
c	   gnew(i,j,k)=gnew(i,j,k)+0.125d0*(gp-2.0d0*gc+gm)
c	   gnew(i,j,k)=gnew(i,j,k)+eta*(gp-2.0d0*gc+gm)
  100   continue
        do 110 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).eq.2) then
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(gnew(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(gnew(i,j,k)-giso)-twb).lt.small)) then
                 gnew(i,j,k)=dsign(twb,(gnew(i,j,k)-giso))+giso
c             endif
           endif
  110   continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 125 k=0,nz
        do 125 j=0,ny
        do 125 i=0,nx
	   diff0=gnew(i,j,k)-gnew(i-1,j,k)
	   diff1=gnew(i+1,j,k)-gnew(i,j,k)
	   diff2=gnew(i+2,j,k)-gnew(i+1,j,k)
	   cc=dabs(cxa(i,j,k))
	   cp=dabs(cxa(i+1,j,k))

	   eps=dmax1(cc,cp)
	   eta=0.5d0*eps*(1.0d0-eps)
	   sgn=dsign(1.0d0,diff1)
c	   antidiff=eta*dabs(diff1)+0.125d0*diff1
 	   antidiff=eta*dabs(diff1)
	   flxmin=dmin1(diff0*sgn,antidiff,diff2*sgn)
	   flux(i,j,k)=sgn*dmax1(0.0d0,flxmin)
  125   continue

        do 150 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 150

 	   gnew(i,j,k)=gnew(i,j,k)+flux(i-1,j,k)-flux(i,j,k)
    	   g(i,j,k)=gnew(i,j,k)
  150   continue
        call bound(nx,ny,nz,g,ibcg)
        call bound(nx,ny,nz,gnew,ibcg)

        do 200 nt=1,ntc
c          if (indx(nt,4).ne.2) goto 200
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)

           gc=g(i,j,k); gm=g(i,j-1,k); gp=g(i,j+1,k)
           cc=cya(i,j,k); cm=cya(i,j-1,k); cp=cya(i,j+1,k)

	   if (cc.ge.0.0d0) then
c	      gnew(i,j,k)=gc+cm*gm-cc*gc
	      gnew(i,j,k)=gc-cc*(gc-gm)
	   else
c	      gnew(i,j,k)=gc+cc*gc-cp*gp
	      gnew(i,j,k)=gc-cc*(gp-gc)
	   endif

	   eps=dabs(cc)
	   eta=0.5d0*eps*(1.0d0-eps)
c	   gnew(i,j,k)=gnew(i,j,k)+0.125d0*(gp-2.0d0*gc+gm)
c	   gnew(i,j,k)=gnew(i,j,k)+eta*(gp-2.0d0*gc+gm)
  200   continue
        do 210 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).eq.2) then
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(gnew(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(gnew(i,j,k)-giso)-twb).lt.small)) then
                 gnew(i,j,k)=dsign(twb,(gnew(i,j,k)-giso))+giso
              endif
c          endif
  210   continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 225 k=0,nz
        do 225 j=0,ny
        do 225 i=0,nx
           diff0=gnew(i,j,k)-gnew(i,j-1,k)
           diff1=gnew(i,j+1,k)-gnew(i,j,k)
           diff2=gnew(i,j+2,k)-gnew(i,j+1,k)
	   cc=dabs(cya(i,j,k))
	   cp=dabs(cya(i,j+1,k))

	   eps=dmax1(cc,cp)
	   eta=0.5d0*eps*(1.0d0-eps)
           sgn=dsign(1.0d0,diff1)
c	   antidiff=eta*dabs(diff1)+0.125d0*diff1
 	   antidiff=eta*dabs(diff1)
	   flxmin=dmin1(diff0*sgn,antidiff,diff2*sgn)
           flux(i,j,k)=sgn*dmax1(0.0d0,flxmin)
  225   continue

        do 250 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 250

           gnew(i,j,k)=gnew(i,j,k)+flux(i,j-1,k)-flux(i,j,k)
           g(i,j,k)=gnew(i,j,k)
  250   continue
        call bound(nx,ny,nz,g,ibcg)
        call bound(nx,ny,nz,gnew,ibcg)

        do 300 nt=1,ntc
c          if (indx(nt,4).ne.2) goto 300
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)

           gc=g(i,j,k); gm=g(i,j,k-1); gp=g(i,j,k+1)
           cc=cza(i,j,k); cm=cza(i,j,k-1); cp=cza(i,j,k+1)

	   if (cc.ge.0.0d0) then
c	      gnew(i,j,k)=gc+cm*gm-cc*gc
	      gnew(i,j,k)=gc-cc*(gc-gm)
	   else
c	      gnew(i,j,k)=gc+cc*gc-cp*gp
	      gnew(i,j,k)=gc-cc*(gp-gc)
	   endif

	   eps=dabs(cc)
	   eta=0.5d0*eps*(1.0d0-eps)
c	   gnew(i,j,k)=gnew(i,j,k)+0.125d0*(gp-2.0d0*gc+gm)
c	   gnew(i,j,k)=gnew(i,j,k)+eta*(gp-2.0d0*gc+gm)
  300   continue
        do 310 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).eq.2) then
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(gnew(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(gnew(i,j,k)-giso)-twb).lt.small)) then
                 gnew(i,j,k)=dsign(twb,(gnew(i,j,k)-giso))+giso
              endif
c          endif
  310   continue
        call bound(nx,ny,nz,gnew,ibcg)

        do 325 k=0,nz
        do 325 j=0,ny
        do 325 i=0,nx
           diff0=gnew(i,j,k)-gnew(i,j,k-1)
           diff1=gnew(i,j,k+1)-gnew(i,j,k)
           diff2=gnew(i,j,k+2)-gnew(i,j,k+1)
	   cc=dabs(cza(i,j,k))
	   cp=dabs(cza(i,j,k+1))

	   eps=dmax1(cc,cp)
	   eta=0.5d0*eps*(1.0d0-eps)
           sgn=dsign(1.0d0,diff1)
c	   antidiff=eta*dabs(diff1)+0.125d0*diff1
 	   antidiff=eta*dabs(diff1)
	   flxmin=dmin1(diff0*sgn,antidiff,diff2*sgn)
           flux(i,j,k)=sgn*dmax1(0.0d0,flxmin)
  325   continue

        do 350 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 350

           gnew(i,j,k)=gnew(i,j,k)+flux(i,j,k-1)-flux(i,j,k)
c          gnew(i,j,k)=gnew(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
           gnew(i,j,k)=gnew(i,j,k)+dt*sl*gradg_ave(i,j,k)
           g(i,j,k)=gnew(i,j,k)
  350   continue
        call bound(nx,ny,nz,g,ibcg)
        call bound(nx,ny,nz,gnew,ibcg)

        do 400 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).eq.2) then
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
c          endif
  400   continue

        call bound(nx,ny,nz,g,ibcg)
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine interpg_old(nx,ny,nz,nptot,ntc,mask,indx,dep,xt,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
        dimension xt(nx+1)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cx_tmp(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy_tmp(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz_tmp(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

	cx_tmp=0.0d0; cy_tmp=0.0d0; cz_tmp=0.0d0 
	dx=1.0d0/dfloat(nt0)
	fm=2.0d0*pi

        do 100 k=-2,nz+3
        do 100 j=-2,ny+3
        do 100 i=-2,nx+3
	   if ((i+ishft).gt.0) then
	      ic=mod((i+ishft-1),nt0)+1
	   else
 	      ic=nt0+mod((i+ishft),nt0)
	   endif
  	   xvel=-0.5d0+dfloat(ic-1)*dx+0.5d0*dx
  	   yvel=-0.5d0+dfloat(j-1)*dx+0.5d0*dx
  	   zvel=-0.5d0+dfloat(k-1)*dx+0.5d0*dx
c 	   xvel=-0.5d0+dfloat(ic-1)*dx
c 	   yvel=-0.5d0+dfloat(j-1)*dx
c 	   zvel=-0.5d0+dfloat(k-1)*dx

           adiv=0.0d0
           utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0

           if (nmodes.gt.1) then
              do 50 n3=1,nmodes
                 rn3=dfloat(n3); wn3=fm*rn3
              do 50 n2=1,nmodes
                 rn2=dfloat(n2); wn2=fm*rn2
              do 50 n1=1,nmodes
                 if ((n1.eq.n2).and.(n2.eq.n3)) goto 50
                 rn1=dfloat(n1); wn1=fm*rn1
                 rnmag=dsqrt(rn1**2.0d0+rn2**2.0d0+rn3**2.0d0)
                 ai=rnmag**(-5.0d0/6.0d0)
                 dum1=(1.0d0-3.0d0*(rn1/rnmag)**2.0d0)/rn1
                 dum2=(1.0d0-3.0d0*(rn2/rnmag)**2.0d0)/rn2
                 dum3=(1.0d0-3.0d0*(rn3/rnmag)**2.0d0)/rn3
                 utmp=utmp+ai*dum1
     &              *dcos(wn1*xvel)*dsin(wn2*yvel)*dsin(wn3*zvel)
                 vtmp=vtmp+ai*dum2
     &              *dsin(wn1*xvel)*dcos(wn2*yvel)*dsin(wn3*zvel)
                 wtmp=wtmp+ai*dum3
     &              *dsin(wn1*xvel)*dsin(wn2*yvel)*dcos(wn3*zvel)
   50         continue
           else
              utmp =  dcos(fm*xvel)*dsin(fm*yvel)*dsin(fm*zvel)
              vtmp = -dsin(fm*xvel)*dcos(fm*yvel)*dsin(fm*zvel)
              wtmp = -dsin(fm*xvel)*dsin(fm*yvel)*dcos(fm*zvel)
           endif

           if (dabs(utmp).lt.small) utmp=0.0d0
           if (dabs(vtmp).lt.small) vtmp=0.0d0
           if (dabs(wtmp).lt.small) wtmp=0.0d0
           cx_tmp(ic,j,k)=dt*utmp/dx
           cy_tmp(ic,j,k)=dt*vtmp/dx
           cz_tmp(ic,j,k)=dt*wtmp/dx

c 	   cx_tmp(ic,j,k)=dt*dsin(fm*yvel)/dx
c 	   cx_tmp(ic,j,k)=dt*dsin(fm*yvel)*dsin(fm*zvel)/dx
c 	   cx_tmp(ic,j,k)=dt*dcos(fm*xvel)*dsin(fm*yvel)/dx
c 	   cx_tmp(ic,j,k)=dt*dcos(fm*xvel)*dsin(fm*zvel)/dx
  100	continue
c	cy_tmp=0.0d0; cz_tmp=0.0d0

c23456789012345678901234567890123456789012345678901234567890123456789012

cc	if (ipltdep.eq.0) then
c	if ((ipltdep.eq.1).or.(ipltdep.eq.50)) then
c	   write (numstr,'(I10)') ipltdep
c	   filename='t'//trim(adjustl(numstr))//'indxoutsc.txt'
c	   open(unit=8,file='t0indxoutsc.txt',status='replace')
c	   open(unit=9,file=filename,status='replace')
c	   dx = 1.0d0/dfloat(nt0)
c	   dy=dx; dz=dx
c	endif

  150	continue
        do 600 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c	   if (indx(nt,4).ne.2) goto 600

	   if ((i+ishft).gt.0) then
	      ic=mod((i+ishft-1),nt0)+1
	   else
 	      ic=nt0+mod((i+ishft),nt0)
	   endif
  	   xvel=-0.5d0+dfloat(ic-1)*dx
  	   yvel=-0.5d0+dfloat(j-1)*dx
  	   zvel=-0.5d0+dfloat(k-1)*dx

           isgn=idnint(dsign(1.0d0,dep(nt,1)))
           jsgn=idnint(dsign(1.0d0,dep(nt,2)))
           ksgn=idnint(dsign(1.0d0,dep(nt,3)))
c          isgn=idnint(dsign(1.0d0,cx_tmp(ic,j,k)))
c          jsgn=idnint(dsign(1.0d0,cy_tmp(ic,j,k)))
c          ksgn=idnint(dsign(1.0d0,cz_tmp(ic,j,k)))

c          isgn=1;jsgn=1;ksgn=1

           icu=ic-isgn
           iu=i-isgn
           ju=j-jsgn
           ku=k-ksgn
           iu2=i-2*isgn
           ju2=j-2*jsgn
           ku2=k-2*ksgn
           id=i+isgn
           jd=j+jsgn
           kd=k+ksgn
 	   if (dabs(dep(nt,1)).lt.small) then
c	   if (dabs(cx_tmp(ic,j,k)).lt.small) then
	      icu=ic
	      iu=i
	      id=i
	      iu2=i
           endif
 	   if (dabs(dep(nt,2)).lt.small) then
c	   if (dabs(cy_tmp(ic,j,k)).lt.small) then
	      ju=j
	      jd=j
	      ju2=j
           endif
 	   if (dabs(dep(nt,3)).lt.small) then
c	   if (dabs(cz_tmp(ic,j,k)).lt.small) then
	      ku=k
	      kd=k
	      ku2=k
           endif

c          cxa=dabs(cx_tmp(ic,j,k))*distfac
c          cya=dabs(cy_tmp(ic,j,k))*distfac
c          cza=dabs(cz_tmp(ic,j,k))*distfac
cc         cxa=dabs(cx_tmp(ic,j,k)+cx_tmp(icu,j,k))/2.0d0
cc         cya=dabs(cy_tmp(ic,j,k)+cy_tmp(icu,j,k))/2.0d0
cc         cza=dabs(cz_tmp(ic,j,k)+cz_tmp(icu,j,k))/2.0d0

           cxa=dabs(dep(nt,1))*distfac
           cya=dabs(dep(nt,2))*distfac
           cza=dabs(dep(nt,3))*distfac
c          cxa=cfn(i,j,k)*dabs(dep(nt,1))*distfac
c          cya=cfn(i,j,k)*dabs(dep(nt,2))*distfac
c          cza=cfn(i,j,k)*dabs(dep(nt,3))*distfac

c23456789012345678901234567890123456789012345678901234567890123456789012

c          if ((ipltdep.eq.1).or.(ipltdep.eq.50)) then
c             write (numstr,'(I10)') ipltdep
c             filename='t'//trim(adjustl(numstr))//'veldat.txt'
c             open(unit=7,file=filename,status='unknown',
c    &		 access='append')
c	      write(7,7120)ic,i,j,k,xt(i),xvel,yvel,zvel,dep(nt,1),
c    &		 cx_tmp(ic,j,k),dep(nt,2),cy_tmp(ic,j,k),dep(nt,3),
c    &		 cz_tmp(ic,j,k)
c             close(7)
c	   endif
c7120      format (4(1x,i3),1x,10(1x,e13.7))

c	   jfix=nt0/2+5; kfix=nt0/2+5; dx=1.0d0/dfloat(nt0)
c	   if ((j.eq.jfix).and.(k.eq.kfix)) then
c 	      x1=0.5d0*dfloat(2*(i-1)-nx)*dx
c 	      y1=-0.5d0+dfloat(j-1)*dx
c 	      z1=-0.5d0+dfloat(k-1)*dx
c 	      xx=-0.5d0+dfloat(i-1)*dx
c 	      yy=-0.5d0+dfloat(j-1)*dx
c 	      zz=-0.5d0+dfloat(k-1)*dx
c             open(unit=7,file="tdat",status='unknown',access='append')
c	      write(7,7122)i,j,k,x1,y1,z1,cxa,cya,cza,g(i,j,k),cfn(i,j,k)
c             close(7)
c	   endif
c7122      format (3(1x,i3),8(2x,e18.12))

c23456789012345678901234567890123456789012345678901234567890123456789012

c	   if ((ipltdep.eq.1).or.(ipltdep.eq.50)) then
cc         if (ipltdep.eq.1) then
c	      x= 0.5d0*dfloat(2*(i-1)-nx)*dx 
c	      y=-0.5d0+dfloat(j-1)*dy
c	      z=-0.5d0+dfloat(k-1)*dz
c	      is=i-1-nx/2
c             write (8,598) is,j,k,x,y,z,g(i,j,k),dep(nt,1),
c    &	         dep(nt,2),dep(nt,3)
c          endif
c 598      format (2(i5,',',1x),i5,',',3x,e18.12,6(',',1x,e18.12)) 

	   gc=g(i,j,k)
	   gxu=g(iu,j,k)
	   gyu=g(i,ju,k)
	   gzu=g(i,j,ku)
	   gxd=g(id,j,k)
	   gyd=g(i,jd,k)
	   gzd=g(i,j,kd)
	   gxu2=g(iu2,j,k)
	   gyu2=g(i,ju2,k)
	   gzu2=g(i,j,ku2)
	   gxuyu=g(iu,ju,k)
	   gxuyd=g(iu,jd,k)
	   gxuzu=g(iu,j,ku)
	   gxuzd=g(iu,j,kd)
	   gyuxd=g(id,ju,k)
	   gyuzu=g(i,ju,ku)
	   gyuzd=g(i,ju,kd)
	   gzuxd=g(id,j,ku)
	   gzuyd=g(i,jd,ku)
	   guuu=g(iu,ju,ku)

	   grdx=gxd-gxu
	   grdy=gyd-gyu
	   grdz=gzd-gzu
	   grd2x=gxd-2.0d0*gc+gxu
	   grd2y=gyd-2.0d0*gc+gyu
	   grd2z=gzd-2.0d0*gc+gzu
	   grd2xu=gc-2.0d0*gxu+gxu2
	   grd2yu=gc-2.0d0*gyu+gyu2
	   grd2zu=gc-2.0d0*gzu+gzu2

c	   if ((gc.ge.gxu).and.(gc.ge.gxd)) cxa = 0.0d0
c	   if ((gc.ge.gyu).and.(gc.ge.gyd)) cya = 0.0d0
c	   if ((gc.ge.gzu).and.(gc.ge.gzd)) cza = 0.0d0

c	   cza=0.0d0
c	   cya=0.0d0 
c	   cxa=0.0d0
c	   gtemp=-(gc-gxu)*cxa-(gc-gyu)*cya-(gc-gzu)*cza

 	   goto 550

	   cxm=0.0d0;cym=0.0d0;czm=0.0d0;cxp=0.0d0;cyp=0.0d0;czp=0.0d0

 	   cxm=0.5d0*(cx(ic-1,j,k)+cx(ic,j,k))
 	   cym=0.5d0*(cy(ic,j-1,k)+cy(ic,j,k))
 	   czm=0.5d0*(cz(ic,j,k-1)+cz(ic,j,k))
 	   cxp=0.5d0*(cx(ic+1,j,k)+cx(ic,j,k))
 	   cyp=0.5d0*(cy(ic,j+1,k)+cy(ic,j,k))
 	   czp=0.5d0*(cz(ic,j,k+1)+cz(ic,j,k))

c	   cxm=cx_tmp(ic-1,j,k)
c	   cym=cy_tmp(ic,j-1,k)
c	   czm=cz_tmp(ic,j,k-1)
c	   cxp=cx_tmp(ic,j,k)
c	   cyp=cy_tmp(ic,j,k)
c	   czp=cz_tmp(ic,j,k)

c	   cxm=0.5d0*(cx_tmp(ic-1,j,k)+cx_tmp(ic,j,k))
c	   cym=0.5d0*(cy_tmp(ic,j-1,k)+cy_tmp(ic,j,k))
c	   czm=0.5d0*(cz_tmp(ic,j,k-1)+cz_tmp(ic,j,k))
c	   cxp=0.5d0*(cx_tmp(ic+1,j,k)+cx_tmp(ic,j,k))
c	   cyp=0.5d0*(cy_tmp(ic,j+1,k)+cy_tmp(ic,j,k))
c	   czp=0.5d0*(cz_tmp(ic,j,k+1)+cz_tmp(ic,j,k))

	   if (dabs(cxm).lt.small) cxm=0.0d0
	   if (dabs(cym).lt.small) cym=0.0d0
	   if (dabs(czm).lt.small) czm=0.0d0
	   if (dabs(cxp).lt.small) cxp=0.0d0
	   if (dabs(cyp).lt.small) cyp=0.0d0
	   if (dabs(czp).lt.small) czp=0.0d0
	   if (cxm.ge.0.0d0) then
	      gxm=g(i-1,j,k)
	   else
	      gxm=g(i,j,k)
	   endif
	   if (cym.ge.0.0d0) then
	      gym=g(i,j-1,k)
	   else
	      gym=g(i,j,k)
	   endif
	   if (czm.ge.0.0d0) then
	      gzm=g(i,j,k-1)
	   else
	      gzm=g(i,j,k)
	   endif
	   if (cxp.le.0.0d0) then
	      gxp=g(i+1,j,k)
	   else
	      gxp=g(i,j,k)
	   endif
	   if (cyp.le.0.0d0) then
	      gyp=g(i,j+1,k)
	   else
	      gyp=g(i,j,k)
	   endif
	   if (czp.le.0.0d0) then
	      gzp=g(i,j,k+1)
	   else
	      gzp=g(i,j,k)
	   endif
c	   if ((cxm*cxp).lt.-small) then
c	      gtempx=0.0d0
c	   else
 	      gtempx=cxm*gxm-cxp*gxp
c	   endif
c	   if ((cym*cyp).lt.-small) then
c	      gtempy=0.0d0
c	   else
 	      gtempy=cym*gym-cyp*gyp
c	   endif
c	   if ((czm*czp).lt.-small) then
c	      gtempz=0.0d0
c	   else
 	      gtempz=czm*gzm-czp*gzp
c	   endif

 	   gtemp=gtempx+gtempy+gtempz

  550	   continue

c23456789012345678901234567890123456789012345678901234567890123456789012

 	   gtemp=-grdx*cxa-grdy*cya-grdz*cza
     &	      +grd2x*cxa**2.0d0+grd2y*cya**2.0d0+grd2z*cza**2.0d0
     &	      +(grd2x-grd2xu)*cxa*(1.0d0-cxa**2.0d0)/3.0d0
     &	      +(grd2y-grd2yu)*cya*(1.0d0-cya**2.0d0)/3.0d0
     &	      +(grd2z-grd2zu)*cza*(1.0d0-cza**2.0d0)/3.0d0
     &	      +(grd2x+gyu-gyuxd+gyd-gxuyd)*cxa*cya
     &	      +(grd2x+gzu-gzuxd+gzd-gxuzd)*cxa*cza
     &	      +(grd2z+gyu-gyuzd+gyd-gzuyd)*cza*cya
     &	      -(grd2x-(gyuxd-2.0d0*gyu+gxuyu))*cya*cxa**2.0d0
     &	      -(grd2y-(gxuyd-2.0d0*gxu+gxuyu))*cxa*cya**2.0d0
     &	      -(grd2x-(gzuxd-2.0d0*gzu+gxuzu))*cza*cxa**2.0d0
     &	      -(grd2z-(gxuzd-2.0d0*gxu+gxuzu))*cxa*cza**2.0d0
     &	      -(grd2y-(gzuyd-2.0d0*gzu+gyuzu))*cza*cya**2.0d0
     &	      -(grd2z-(gyuzd-2.0d0*gyu+gyuzu))*cya*cza**2.0d0
     &	      -2.0d0*(gc-gxu+gxuyu-gyu-(gzu-gxuzu+guuu-gyuzu))
     &        *cxa*cya*cza
	   if (dabs(gtemp).lt.small) gtemp=0.0d0
 	   gnew(i,j,k)=gc+0.5d0*gtemp
c	   gnew(i,j,k)=gc+gtemp
c          if (limit.eq.1) then
c             gamin=dmin1(gc,gxu,gyu,gzu,gxuyu,gxuzu,gyuzu,guuu)
c             gamax=dmax1(gc,gxu,gyu,gzu,gxuyu,gxuzu,gyuzu,guuu)
c	      gnew(i,j,k)=dmin1(gnew(i,j,k),gamax)
  	      gnew(i,j,k)=dmax1(gnew(i,j,k),gamin)
     &		 +dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
cc   &		 +dt*sl*gradg_ave(i,j,k)
c	   else
cc	      gnew(i,j,k)=gnew(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
 	      gnew(i,j,k)=gnew(i,j,k)+dt*sl*gradg_ave(i,j,k)
c          endif
c	   if (dabs(gnew(i,j,k)-gc).gt.small) write (6,*) 'g changed'

c          if (ipltdep.eq.0) then
c          if (ipltdep.eq.1) then
	   if ((ipltdep.eq.1).or.(ipltdep.eq.50)) then
	      x= 0.5d0*dfloat(2*(i-1)-nx)*dx 
	      y=-0.5d0+dfloat(j-1)*dy
	      z=-0.5d0+dfloat(k-1)*dz
	      is=i-1-nx/2
              write (9,599) is,j,k,x,y,z,gnew(i,j,k),dep(nt,1),
     &	         dep(nt,2),dep(nt,3)
           endif
  599      format (2(i5,',',1x),i5,',',3x,e18.12,6(',',1x,e18.12)) 

  600	continue
	if ((ipltdep.eq.1).or.(ipltdep.eq.50)) then
c       if (ipltdep.eq.0) then
	   close(8)
	   close(9)
        endif

        do 700 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
cc	   if (indx(nt,4).eq.2) then
	      g(i,j,k)=gnew(i,j,k)
 	      giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &		 (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
cc	   else
cc	      goto 700
cc	   endif
c	   gc = g(i,j,k)
c	   if ((dabs(gc).gt.small)
c    &	      .and.(gc.gt.g(i+1,j,k)).and.(gc.gt.g(i-1,j,k))
c    &	      .and.(gc.gt.g(i,j+1,k)).and.(gc.gt.g(i,j-1,k))
c    &	      .and.(gc.gt.g(i,j,k+1)).and.(gc.gt.g(i,j,k-1)))
c    &     then
c	      write (6,699) i,j,k,g(i-1,j,k),g(i,j,k),
c    &	         g(i+1,j,k),xt(1),xt(i),dep(nt,1),dep(nt,2),dep(nt,3)
c	   endif
c699    format (3(1x,i3),8(1x,e18.12))

c	   gcn = g(i,j,k)
c	   if ((gcn.gt.g(i-1,j,k)).and.(gcn.gt.g(i+1,j,k))
c    & 	      .and.(gcn.gt.g(i,j-1,k)).and.(gcn.gt.g(i,j+1,k))
c    &        .and.(gcn.gt.g(i,j,k-1)).and.(gcn.gt.g(i,j,k+1))) then
c	      g(i,j,k)
c    &		  =dmax1(g(i-1,j,k),g(i+1,j,k),g(i,j-1,k),g(i,j+1,k),
c    &		  g(i,j,k-1),g(i,j,k+1))
c	   endif
  700	continue

c23456789012345678901234567890123456789012345678901234567890123456789012

c	do it = 1,2
c       do 800 nt=1,ntc
c          i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c	   if ((indx(nt,4).eq.2).or.(indx(nt,4).eq.1)) then
c             gamin=dmin1(g(i,j,k),g(i-1,j,k),g(i+1,j,k),g(i,j-1,k),
c    &	         g(i,j+1,k),g(i,j,k-1),g(i,j,k+1))
c             gamax=dmax1(g(i,j,k),g(i-1,j,k),g(i+1,j,k),g(i,j-1,k),
c    &	         g(i,j+1,k),g(i,j,k-1),g(i,j,k+1))
c	      g(i,j,k)=dmin1(g(i,j,k),gamax)
c 	      g(i,j,k)=dmax1(g(i,j,k),gamin)
c	   endif
c 800	continue
c	end do

        call bound(nx,ny,nz,g,ibcg)
c       call bound(nx,ny,nz,gnew,ibcg)

	return
	end


c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine interpg(nx,ny,nz,nptot,ntc,mask,indx,dep,xt,ishft,
     &     distfac,g,ibcg,limit,cfn,gradg_ave,dt,sl,ipltdep,nt0,cx,cy,
     &	   cz,nmodes)
        implicit double precision (a-h,o-z)
	character*160 filename,numstr
        dimension xt(nx+1)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cx_tmp(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy_tmp(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz_tmp(-2:nt0+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

        do 600 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c	   if (indx(nt,4).ne.2) goto 600

           isgn=idnint(dsign(1.0d0,dep(nt,1)))
           jsgn=idnint(dsign(1.0d0,dep(nt,2)))
           ksgn=idnint(dsign(1.0d0,dep(nt,3)))
c          isgn=idnint(dsign(1.0d0,cx_tmp(ic,j,k)))
c          jsgn=idnint(dsign(1.0d0,cy_tmp(ic,j,k)))
c          ksgn=idnint(dsign(1.0d0,cz_tmp(ic,j,k)))

c          isgn=1;jsgn=1;ksgn=1

c23456789012345678901234567890123456789012345678901234567890123456789012

           icu=ic-isgn
           iu=i-isgn
           ju=j-jsgn
           ku=k-ksgn
           iu2=i-2*isgn
           ju2=j-2*jsgn
           ku2=k-2*ksgn
           id=i+isgn
           jd=j+jsgn
           kd=k+ksgn
 	   if (dabs(dep(nt,1)).lt.small) then
c	   if (dabs(cx_tmp(ic,j,k)).lt.small) then
	      icu=ic
	      iu=i
	      id=i
	      iu2=i
           endif
 	   if (dabs(dep(nt,2)).lt.small) then
c	   if (dabs(cy_tmp(ic,j,k)).lt.small) then
	      ju=j
	      jd=j
	      ju2=j
           endif
 	   if (dabs(dep(nt,3)).lt.small) then
c	   if (dabs(cz_tmp(ic,j,k)).lt.small) then
	      ku=k
	      kd=k
	      ku2=k
           endif

           cxa=dabs(dep(nt,1))*distfac
           cya=dabs(dep(nt,2))*distfac
           cza=dabs(dep(nt,3))*distfac
c          cxa=cfn(i,j,k)*dabs(dep(nt,1))*distfac
c          cya=cfn(i,j,k)*dabs(dep(nt,2))*distfac
c          cza=cfn(i,j,k)*dabs(dep(nt,3))*distfac

c23456789012345678901234567890123456789012345678901234567890123456789012

	   gc=g(i,j,k)
	   gxu=g(iu,j,k)
	   gyu=g(i,ju,k)
	   gzu=g(i,j,ku)
	   gxd=g(id,j,k)
	   gyd=g(i,jd,k)
	   gzd=g(i,j,kd)
	   gxu2=g(iu2,j,k)
	   gyu2=g(i,ju2,k)
	   gzu2=g(i,j,ku2)
	   gxuyu=g(iu,ju,k)
	   gxuyd=g(iu,jd,k)
	   gxuzu=g(iu,j,ku)
	   gxuzd=g(iu,j,kd)
	   gyuxd=g(id,ju,k)
	   gyuzu=g(i,ju,ku)
	   gyuzd=g(i,ju,kd)
	   gzuxd=g(id,j,ku)
	   gzuyd=g(i,jd,ku)
	   guuu=g(iu,ju,ku)

	   grdx=gxd-gxu
	   grdy=gyd-gyu
	   grdz=gzd-gzu
	   grd2x=gxd-2.0d0*gc+gxu
	   grd2y=gyd-2.0d0*gc+gyu
	   grd2z=gzd-2.0d0*gc+gzu
	   grd2xu=gc-2.0d0*gxu+gxu2
	   grd2yu=gc-2.0d0*gyu+gyu2
	   grd2zu=gc-2.0d0*gzu+gzu2

c23456789012345678901234567890123456789012345678901234567890123456789012

 	   gtemp=-grdx*cxa-grdy*cya-grdz*cza
     &	      +grd2x*cxa**2.0d0+grd2y*cya**2.0d0+grd2z*cza**2.0d0
     &	      +(grd2x-grd2xu)*cxa*(1.0d0-cxa**2.0d0)/3.0d0
     &	      +(grd2y-grd2yu)*cya*(1.0d0-cya**2.0d0)/3.0d0
     &	      +(grd2z-grd2zu)*cza*(1.0d0-cza**2.0d0)/3.0d0
     &	      +(grd2x+gyu-gyuxd+gyd-gxuyd)*cxa*cya
     &	      +(grd2x+gzu-gzuxd+gzd-gxuzd)*cxa*cza
     &	      +(grd2z+gyu-gyuzd+gyd-gzuyd)*cza*cya
     &	      -(grd2x-(gyuxd-2.0d0*gyu+gxuyu))*cya*cxa**2.0d0
     &	      -(grd2y-(gxuyd-2.0d0*gxu+gxuyu))*cxa*cya**2.0d0
     &	      -(grd2x-(gzuxd-2.0d0*gzu+gxuzu))*cza*cxa**2.0d0
     &	      -(grd2z-(gxuzd-2.0d0*gxu+gxuzu))*cxa*cza**2.0d0
     &	      -(grd2y-(gzuyd-2.0d0*gzu+gyuzu))*cza*cya**2.0d0
     &	      -(grd2z-(gyuzd-2.0d0*gyu+gyuzu))*cya*cza**2.0d0
     &	      -2.0d0*(gc-gxu+gxuyu-gyu-(gzu-gxuzu+guuu-gyuzu))
     &        *cxa*cya*cza
	   if (dabs(gtemp).lt.small) gtemp=0.0d0
 	   gnew(i,j,k)=gc+0.5d0*gtemp
c	   gnew(i,j,k)=gc+0.5d0*cfn(i,j,k)*gtemp


	   if ((dabs(grdx).lt.small).and.(dabs(grdy).lt.small)
     &	      .and.(dabs(grdz).lt.small)) then
	      gnew(i,j,k)=gc + dt*sl*gradg_ave(i,j,k)
 	   else
cc	      gnew(i,j,k)=gnew(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
 	      gnew(i,j,k)=gnew(i,j,k)+dt*sl*gradg_ave(i,j,k)
           endif
c	   if (dabs(gnew(i,j,k)-gc).gt.small) write (6,*) 'g changed'

  600	continue

        do 700 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
cc	   if (indx(nt,4).eq.2) then
	      g(i,j,k)=gnew(i,j,k)
 	      giso=gjump*dfloat(mask(i,j,k))
	      if (dabs(g(i,j,k)-giso).lt.small) then
		 g(i,j,k)=giso
	      endif
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &		 (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
  700	continue

        call bound(nx,ny,nz,g,ibcg)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,g,ibcg)
	endif

	return
	end


c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,g,ibcg)
        implicit double precision (a-h,o-z)
        dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4)
        common /tubedat/twa,twb,gjump,nomask
        common /rparam/small,pi

        do 600 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          if (indx(nt,4).ne.2) goto 600
 	   giso=gjump*dfloat(mask(i,j,k))

	   gc=g(i,j,k)-giso
	   gxm=g(i-1,j,k)-giso
           gym=g(i,j-1,k)-giso
           gzm=g(i,j,k-1)-giso
           gxp=g(i+1,j,k)-giso
           gyp=g(i,j+1,k)-giso
           gzp=g(i,j,k+1)-giso
	   gmax = dmax1(gxm,gym,gzm,gxp,gyp,gzp)
	   gmin = dmin1(gxm,gym,gzm,gxp,gyp,gzp)

	   if ((dabs(gxm-twb)<small).and.(dabs(gxp-twb)<small)
     &	      .and.(dabs(gym-twb)<small).and.(dabs(gyp-twb)<small)
     &	      .and.(dabs(gzm-twb)<small).and.(dabs(gzp-twb)<small)) then
	      g(i,j,k)=giso+twb
	   elseif ((dabs(gxm+twb)<small).and.(dabs(gxp+twb)<small)
     &        .and.(dabs(gym+twb)<small).and.(dabs(gyp+twb)<small)
     &        .and.(dabs(gzm+twb)<small).and.(dabs(gzp+twb)<small)) then
              g(i,j,k)=giso-twb
	   endif

	   if (gc > gmax) then 
	      g(i,j,k)=gmax+giso
	   elseif (gc < gmin) then
	      g(i,j,k)=gmin+giso
	   endif
	if (giso.ne.0.0d0) write(6,*) "giso = ",giso
  600	continue
        call bound(nx,ny,nz,g,ibcg)
	return
	end


c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine interpgweno(nx,ny,nz,nptot,ntc,mask,indx,dep,
     &     distfac,g,ibcg,ibcv,limit,cfn,gradg_ave,dt,sl)
        implicit double precision (a-h,o-z)
 	dimension gradDotCfl(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gnew(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gtmp(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gradg_ave(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension mask(nx,ny,nz),indx(nptot,4),dep(ntc,3)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

	gtmp=0.0d0

 	call wenoGradCfl(nx,ny,nz,nptot,ntc,indx,dep,
     &	   g,gradDotCfl,ibcv,distfac)
        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
 	   if (indx(nt,4).ne.2) goto 100
	   gtmp(i,j,k)=gradDotCfl(i,j,k)
c	   g(i,j,k)=g(i,j,k)-gtmp(i,j,k)*cfn(i,j,k)
	   g(i,j,k)=g(i,j,k)-gtmp(i,j,k)
  100	continue

c	goto 450

 	call wenoGradCfl(nx,ny,nz,nptot,ntc,indx,dep,
     &	   gtmp,gradDotCfl,ibcv,distfac)
        do 200 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
	   gtmp(i,j,k)=gradDotCfl(i,j,k)
c	   g(i,j,k)=g(i,j,k)+0.5d0*gtmp(i,j,k)*cfn(i,j,k)
	   g(i,j,k)=g(i,j,k)+0.5d0*gtmp(i,j,k)
  200	continue

 	call wenoGradCfl(nx,ny,nz,nptot,ntc,indx,dep,
     &	   gtmp,gradDotCfl,ibcv,distfac)
        do 300 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
	   gtmp(i,j,k)=gradDotCfl(i,j,k)
c	   g(i,j,k)=g(i,j,k)-gtmp(i,j,k)*cfn(i,j,k)/6.0d0
	   g(i,j,k)=g(i,j,k)-gtmp(i,j,k)/6.0d0
  300	continue

c23456789012345678901234567890123456789012345678901234567890123456789012

	if (norda.eq.5) then
 	   call wenoGradCfl(nx,ny,nz,nptot,ntc,indx,dep,
     &	      gtmp,gradDotCfl,ibcv,distfac)
           do 400 nt=1,ntc
              i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
	      gtmp(i,j,k)=gradDotCfl(i,j,k)
c	      g(i,j,k)=g(i,j,k)+gtmp(i,j,k)*cfn(i,j,k)/24.0d0
	      g(i,j,k)=g(i,j,k)+gtmp(i,j,k)/24.0d0
  400	   continue
	endif

  450	continue

        do 500 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
 	   if (indx(nt,4).ne.2) goto 500
	   gtmp(i,j,k)=gradDotCfl(i,j,k)
	   g(i,j,k)=g(i,j,k)+dt*sl*cfn(i,j,k)*gradg_ave(i,j,k)
	   if (dabs(g(i,j,k)).lt.small) g(i,j,k)=0.0d0
  500	continue

        do 700 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
           if (indx(nt,4).eq.2) then
              giso=gjump*dfloat(mask(i,j,k))
              if ((dabs(g(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(g(i,j,k)-giso)-twb).lt.1e-10)) then
                 g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
              endif
           else
              goto 700
           endif
  700   continue

        call bound(nx,ny,nz,g,ibcg)
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012
	
 	subroutine wenoGradCfl(nx,ny,nz,nptot,ntc,indx,dep,
     &	   g,gradDotCfl,ibcv,distfac)
        implicit double precision (a-h,o-z)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
 	dimension gradDotCfl(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension indx(nptot,4),dep(ntc,3)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask

        eps=small; gradDotCfl=0.0d0

        do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gxc=8.0d0*(g(i+1,j,k)-g(i-1,j,k))-(g(i+2,j,k)-g(i-2,j,k))
           gxc=gxc/12.0d0
           am=g(i-3,j,k)-2.0d0*g(i-2,j,k)+g(i-1,j,k)
           bm=g(i-2,j,k)-2.0d0*g(i-1,j,k)+g(i,j,k)
           cm=g(i-1,j,k)-2.0d0*g(i,j,k)+g(i+1,j,k)
           dm=g(i,j,k)-2.0d0*g(i+1,j,k)+g(i+2,j,k)

	   if (dep(nt,1).ge.0.0d0) then

              s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
              s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
              s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
              a0=1.0d0/(eps+s0)**2.0d0
              a1=6.0d0/(eps+s1)**2.0d0
              a2=3.0d0/(eps+s2)**2.0d0
              w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
              gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &           +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
              gxm=gxc-gweno
	      gradDotCfl(i,j,k)=distfac*dep(nt,1)*gxm

	   else

              ap=g(i+1,j,k)-2.0d0*g(i+2,j,k)+g(i+3,j,k)
              bp=dm; dp=bm; cp=cm
              s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
              s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
              s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
              a0=1.0d0/(eps+s0)**2.0d0
              a1=6.0d0/(eps+s1)**2.0d0
              a2=3.0d0/(eps+s2)**2.0d0
              w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
              gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &           +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
              gxp=gxc+gweno
	      gradDotCfl(i,j,k)=distfac*dep(nt,1)*gxp

	   endif

c	goto 100

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gyc=8.0d0*(g(i,j+1,k)-g(i,j-1,k))-(g(i,j+2,k)-g(i,j-2,k))
           gyc=gyc/12.0d0
           am=g(i,j-3,k)-2.0d0*g(i,j-2,k)+g(i,j-1,k)
           bm=g(i,j-2,k)-2.0d0*g(i,j-1,k)+g(i,j,k)
           cm=g(i,j-1,k)-2.0d0*g(i,j,k)+g(i,j+1,k)
           dm=g(i,j,k)-2.0d0*g(i,j+1,k)+g(i,j+2,k)

	   if (dep(nt,2).ge.0.0d0) then

              s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
              s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
              s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
              a0=1.0d0/(eps+s0)**2.0d0
              a1=6.0d0/(eps+s1)**2.0d0
              a2=3.0d0/(eps+s2)**2.0d0
              w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
              gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &           +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
              gym=gyc-gweno
	      gradDotCfl(i,j,k)=gradDotCfl(i,j,k)
     &		 +distfac*dep(nt,2)*gym

	   else

              ap=g(i,j+1,k)-2.0d0*g(i,j+2,k)+g(i,j+3,k)
              bp=dm; dp=bm; cp=cm
              s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
              s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
              s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
              a0=1.0d0/(eps+s0)**2.0d0
              a1=6.0d0/(eps+s1)**2.0d0
              a2=3.0d0/(eps+s2)**2.0d0
              w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
              gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &           +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
              gyp=gyc+gweno
	      gradDotCfl(i,j,k)=gradDotCfl(i,j,k)
     &		 +distfac*dep(nt,2)*gyp

	   endif

c23456789012345678901234567890123456789012345678901234567890123456789012
	
           gzc=8.0d0*(g(i,j,k+1)-g(i,j,k-1))-(g(i,j,k+2)-g(i,j,k-2))
           gzc=gzc/12.0d0
           am=g(i,j,k-3)-2.0d0*g(i,j,k-2)+g(i,j,k-1)
           bm=g(i,j,k-2)-2.0d0*g(i,j,k-1)+g(i,j,k)
           cm=g(i,j,k-1)-2.0d0*g(i,j,k)+g(i,j,k+1)
           dm=g(i,j,k)-2.0d0*g(i,j,k+1)+g(i,j,k+2)

	   if (dep(nt,3).ge.0.0d0) then

              s0=13.0d0*(am-bm)**2.0d0+3.0d0*(am-3.0d0*bm)**2.0d0
              s1=13.0d0*(bm-cm)**2.0d0+3.0d0*(bm+cm)**2.0d0
              s2=13.0d0*(cm-dm)**2.0d0+3.0d0*(3.0d0*cm-dm)**2.0d0
              a0=1.0d0/(eps+s0)**2.0d0
              a1=6.0d0/(eps+s1)**2.0d0
              a2=3.0d0/(eps+s2)**2.0d0
              w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
              gweno=w0*(am-2.0d0*bm+cm)/3.0d0
     &           +(w2-0.5d0)*(bm-2.0d0*cm+dm)/6.0d0
              gzm=gzc-gweno
	      gradDotCfl(i,j,k)=gradDotCfl(i,j,k)
     &		 +distfac*dep(nt,3)*gzm

	   else

              ap=g(i,j,k+1)-2.0d0*g(i,j,k+2)+g(i,j,k+3)
              bp=dm; dp=bm; cp=cm
              s0=13.0d0*(ap-bp)**2.0d0+3.0d0*(ap-3.0d0*bp)**2.0d0
              s1=13.0d0*(bp-cp)**2.0d0+3.0d0*(bp+cp)**2.0d0
              s2=13.0d0*(cp-dp)**2.0d0+3.0d0*(3.0d0*cp-dp)**2.0d0
              a0=1.0d0/(eps+s0)**2.0d0
              a1=6.0d0/(eps+s1)**2.0d0
              a2=3.0d0/(eps+s2)**2.0d0
              w0=a0/(a0+a1+a2); w2=a2/(a0+a1+a2); 
              gweno=w0*(ap-2.0d0*bp+cp)/3.0d0
     &           +(w2-0.5d0)*(bp-2.0d0*cp+dp)/6.0d0
              gzp=gzc+gweno
	      gradDotCfl(i,j,k)=gradDotCfl(i,j,k)
     &		 +distfac*dep(nt,3)*gzp

	   endif

           if (dabs(gradDotCfl(i,j,k)).lt.small) gradDotCfl(i,j,k)=0.0d0

  100   continue

        call bound(nx,ny,nz,gradDotCfl,ibcv)
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine reinwenorkts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &	   iscf,nreinit,nord,ibcg,ibcv,cfn,g)
        implicit double precision (a-h,o-z)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gtmp(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension gold(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension grad0(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension grad(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension f1(nx,ny,nz),f2(nx,ny,nz),f3(nx,ny,nz)
        dimension mask(nx,ny,nz),indx(nptot,4)
        dimension x(nx+1),y(ny+1),z(nz+1)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
        common /rparam/small,pi
	common /tubedat/twa,twb,gjump,nomask
	common /reindat/irk,igrd,cfl_re

	limit = 1
        dt=cfl_re*dx

c Need to define gtmp everywhere, not just on T, since wenogradts 
c needs adjacent points for calculation of its gradient

        do 50 k=-2,nz+3
        do 50 j=-2,ny+3
        do 50 i=-2,nx+3
           if (dabs(g(i,j,k)).lt.small) g(i,j,k)=0.0d0
           gtmp(i,j,k)=g(i,j,k)
 50     continue
 	call bound(nx,ny,nz,gtmp,ibcg)

c A HIGH-ORDER GRADIENT FOR THE ORIGINAL (PRE-REINITIALIZED) DISTRIBUTION
c (E.G., WITH WENOGRADTS), COMBINED WITH THE USE OF A FIRST-ORDER GRADIENT 
c IN THE SIGNED-DISTANCE EQUATION PROVIDES SUPERIOR RESULTS WITH FEWEST 
c GRID POINTS; I.E., SAME RESULTS WITH 32x32x32 AS WITH 64x64x64. FOR THE 
c HIGH-ORDER GRADIENT, WENO WORKS BETTER THAN ENO.

	if (igrd.eq.0) then
 	   call enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &        nord,g,g,grad0,ibcv,'rein',0)
	else
 	   call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &        nord,g,g,grad0,ibcv,'rein',0)
c	   call gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,g,grad0,ibcv,'rein',0)
	endif

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 300 n=1,nreinit
	
 	call gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &     nord,g,gtmp,grad,ibcv,'rein',20)

c	if (igrd.eq.0) then
c	   call enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',10*n)
c	else
c	   call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',20)
c	endif

 	do 100 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
c	   if ((indx(nt,4).ne.2).and.(indx(nt,4).ne.1)) goto 100
c	   if (indx(nt,4).ne.2) goto 100
	   gold(i,j,k)=gtmp(i,j,k)
 	   gm=gjump*dfloat(mask(i,j,k))
           gc=g(i,j,k)-gm
           gctmp=gtmp(i,j,k)-gm
	   if (dabs(gc).lt.small) gc=0.0d0

           sgnd=dsign(1.0d0,gc)
	   if (gc.eq.0.0d0) sgnd=0.0d0

c          if (dabs(gc).le.dx) then
c             sgn=gc/dsqrt(gc**2.0d0+dx**2.0d0)
c             sgn=gc/dsqrt(gc**2.0d0+(dx*grad(i,j,k))**2.0d0)
              sgn=gc/dsqrt(gc**2.0d0+(dx*grad0(i,j,k))**2.0d0)
c             sgn=gc/dx+dsin(pi*gc/dx)/pi
c          else
c             sgn=dsign(1.0d0,gc)
c          endif

cc         if (dabs(gctmp).le.dx) then
cc            sgn=gctmp/dsqrt(gctmp**2.0d0+dx**2.0d0)
c             sgn=gctmp/dx+dsin(pi*gctmp/dx)/pi
cc         else
cc            sgn=dsign(1.0d0,gctmp)
cc         endif

c23456789012345678901234567890123456789012345678901234567890123456789012

           f1(i,j,k)=sgnd*(1.0d0-grad(i,j,k))

	   iflg1=0
           if ((iscf.eq.1).and.
     &        ((gc*(g(i-1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i+1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j-1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j+1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k-1)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k+1)-gm).lt.0.0d0))) then

	      if (dabs(gc).gt.small) then
	         if (grad0(i,j,k).ge.(dabs(gc)/dx)) then
                    dist=gc/grad0(i,j,k)

c                   gxm1=(g(i,j,k)-g(i-1,j,k))/dx
c                   gxp1=(g(i+1,j,k)-g(i,j,k))/dx
c                   gym1=(g(i,j,k)-g(i,j-1,k))/dy
c                   gyp1=(g(i,j+1,k)-g(i,j,k))/dy
c                   gzm1=(g(i,j,k)-g(i,j,k-1))/dz
c                   gzp1=(g(i,j,k+1)-g(i,j,k))/dz
c	 	    if (gc.gt.0.0d0) then
c		       gx2=dmax1(dmax1(0.0d0,gxm1)**2.0d0,
c    &		          dmin1(0.0d0,gxp1)**2.0d0)
c		       gy2=dmax1(dmax1(0.0d0,gym1)**2.0d0,
c    &		          dmin1(0.0d0,gyp1)**2.0d0)
c		       gz2=dmax1(dmax1(0.0d0,gzm1)**2.0d0,
c    &		          dmin1(0.0d0,gzp1)**2.0d0)
c		    else
c                      gx2=dmax1(dmin1(0.0d0,gxm1)**2.0d0,
c    &                    dmax1(0.0d0,gxp1)**2.0d0)
c                      gy2=dmax1(dmin1(0.0d0,gym1)**2.0d0,
c    &                    dmax1(0.0d0,gyp1)**2.0d0)
c                      gz2=dmax1(dmin1(0.0d0,gzm1)**2.0d0,
c    &                    dmax1(0.0d0,gzp1)**2.0d0)
c		    endif
c		    gradg=dsqrt(gx2+gy2+gz2)
c                   dist=gc/gradg
	         else
                    gxm1=(g(i,j,k)-g(i-1,j,k))/dx
                    gxp1=(g(i+1,j,k)-g(i,j,k))/dx
                    gym1=(g(i,j,k)-g(i,j-1,k))/dy
                    gyp1=(g(i,j+1,k)-g(i,j,k))/dy
                    gzm1=(g(i,j,k)-g(i,j,k-1))/dz
                    gzp1=(g(i,j,k+1)-g(i,j,k))/dz
		    gxc1=0.5d0*(g(i+1,j,k)-g(i-1,j,k))/dx
		    gyc1=0.5d0*(g(i,j+1,k)-g(i,j-1,k))/dy
		    gzc1=0.5d0*(g(i,j,k+1)-g(i,j,k-1))/dz
	            gx2=dmax1(gxm1**2.0d0,gxp1**2.0d0,gxc1**2.0d0)
	            gy2=dmax1(gym1**2.0d0,gyp1**2.0d0,gyc1**2.0d0)
	            gz2=dmax1(gzm1**2.0d0,gzp1**2.0d0,gzc1**2.0d0)
	            gradg=dmax1(dabs(gc)/dx,dsqrt(gx2+gy2+gz2))
                    dist=gc/gradg
	         endif
	      else
	         dist=0.0d0
	      endif

c             gradg=dmax1(grad0(i,j,k),small)
c             dist=gc/gradg

              f1(i,j,k)=(dist-sgnd*dabs(gctmp))/dx
	      iflg1=1
	   endif

           if (irk.eq.1) then
              gtmp(i,j,k)=gold(i,j,k)+dt*f1(i,j,k)/2.0d0
           else
              gtmp(i,j,k)=gold(i,j,k)+dt*f1(i,j,k)

c	      gmin=dmin1(gtmp(i-1,j,k),gtmp(i+1,j,k),gtmp(i,j-1,k),
c    &		 gtmp(i,j+1,k),gtmp(i,j,k-1),gtmp(i,j,k+1))
c	      gmax=dmax1(gtmp(i-1,j,k),gtmp(i+1,j,k),gtmp(i,j-1,k),
c    &		 gtmp(i,j+1,k),gtmp(i,j,k-1),gtmp(i,j,k+1))
c	      gtmp(i,j,k)=dmin1(gtmp(i,j,k),gmax)
c	      gtmp(i,j,k)=dmax1(gtmp(i,j,k),gmin)
           endif

           if ((dabs(gtmp(i,j,k)-giso).gt.twb).or.
     &        (dabs(dabs(gtmp(i,j,k)-giso)-twb).lt.small)) then
              gtmp(i,j,k)=dsign(twb,(gtmp(i,j,k)-giso))+giso
	   endif

	   if (dabs(gtmp(i,j,k)-gm).lt.small) gtmp(i,j,k)=gm
           if (((gtmp(i,j,k)-gm)*gc).lt.0.0d0) then
              gtmp(i,j,k)=gm
              write (6,*) 'reinwenorkts(rk1):'
	      write (6,1000) i,j,k,nt,iflg1,mask(i,j,k),indx(nt,4),
     &		 g(i,j,k),gtmp(i,j,k),gc,grad0(i,j,k),grad(i,j,k)
           endif
 100    continue
 	call bound(nx,ny,nz,gtmp,ibcg)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,gtmp,ibcg)
	endif

	do 125 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)
 	   gmax=dmax1(gtmp(i-1,j,k),gtmp(i+1,j,k),gtmp(i,j-1,k),
     &        gtmp(i,j+1,k),gtmp(i,j,k-1),gtmp(i,j,k+1))
	   if (gtmp(i,j,k).gt.gmax) then
              write (6,2000) i,j,k,nx,g(i-1,j,k),g(i,j,k),
     &           g(i+1,j,k),g(i,j-1,k),g(i,j+1,k),g(i,j,k-1),g(i,j,k+1)
              write (6,2000) i,j,k,nx,gtmp(i-1,j,k),gtmp(i,j,k),
     &           gtmp(i+1,j,k),gtmp(i,j-1,k),gtmp(i,j+1,k),
     &		 gtmp(i,j,k-1),gtmp(i,j,k+1)
c	      write (6,*) ""
c	      write (6,2010) (g(idum,j,k),idum=-2,nx+3)
c	      write (6,2010) (gtmp(idum,j,k),idum=-2,nx+3)
c	      write (6,2010) (g(i,jdum,k),jdum=-2,ny+3)
c	      write (6,2010) (gtmp(i,jdum,k),jdum=-2,ny+3)
c	      write (6,2010) (g(i,j,kdum),kdum=-2,nz+3)
c	      write (6,2010) (gtmp(i,kdum,k),kdum=-2,nz+3)
c	      write (6,*) ""
	   endif
 125    continue
 2000   format (4(1x,i3),7(1x,e12.6))
 2010   format (40(1x,e12.6))

	if (irk.eq.0) goto 300
	
c23456789012345678901234567890123456789012345678901234567890123456789012

c	if (igrd.eq.0) then
c	   call enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',2)
c	else
c	   call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',2)
c	endif

 	call gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &     nord,g,gtmp,grad,ibcv,'rein',20)

 	do 150 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)
c	   if ((indx(nt,4).ne.2).or.(indx(nt,4).ne.1)) goto 150
 	   gm=gjump*dfloat(mask(i,j,k))
           gc=g(i,j,k)-gm
           gctmp=gtmp(i,j,k)-gm
	   if (dabs(gc).lt.small) gc=0.0d0

           sgnd=dsign(1.0d0,gc)
	   if (gc.eq.0.0d0) sgnd=0.0d0

           if (dabs(gc).le.dx) then
              sgn=gc/dsqrt(gc**2.0d0+dx**2.0d0)
c             sgn=gc/dsqrt(gc**2.0d0+(dx*grad(i,j,k))**2.0d0)
c             sgn=gc/dx+dsin(pi*gc/dx)/pi
           else
              sgn=dsign(1.0d0,gc)
           endif
cc         if (dabs(gctmp).le.dx) then
cc            sgn=gctmp/dsqrt(gctmp**2.0d0+dx**2.0d0)
c             sgn=gctmp/dx+dsin(pi*gctmp/dx)/pi
cc         else
cc            sgn=dsign(1.0d0,gctmp)
cc         endif

           f2(i,j,k)=sgnd*(1.0d0-grad(i,j,k))

	   iflg1=0
           if ((iscf.eq.1).and.
     &        ((gc*(g(i-1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i+1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j-1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j+1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k-1)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k+1)-gm).lt.0.0d0))) then

	      if (dabs(gc).gt.small) then
	         if (grad0(i,j,k).ge.(dabs(gc)/dx)) then
                    dist=gc/grad0(i,j,k)
	         else
                    gxm1=(g(i,j,k)-g(i-1,j,k))/dx
                    gxp1=(g(i+1,j,k)-g(i,j,k))/dx
                    gym1=(g(i,j,k)-g(i,j-1,k))/dy
                    gyp1=(g(i,j+1,k)-g(i,j,k))/dy
                    gzm1=(g(i,j,k)-g(i,j,k-1))/dz
                    gzp1=(g(i,j,k+1)-g(i,j,k))/dz
		    gxc1=0.5d0*(g(i+1,j,k)-g(i-1,j,k))/dx
		    gyc1=0.5d0*(g(i,j+1,k)-g(i,j-1,k))/dy
		    gzc1=0.5d0*(g(i,j,k+1)-g(i,j,k-1))/dz
	            gx2=dmax1(gxm1**2.0d0,gxp1**2.0d0,gxc1**2.0d0)
	            gy2=dmax1(gym1**2.0d0,gyp1**2.0d0,gyc1**2.0d0)
	            gz2=dmax1(gzm1**2.0d0,gzp1**2.0d0,gzc1**2.0d0)
	            gradg=dmax1(dabs(gc)/dx,dsqrt(gx2+gy2+gz2))
                    dist=gc/gradg
	         endif
	      else
	         dist=0.0d0
	      endif

c             gradg=dmax1(grad0(i,j,k),small)
c             dist=gc/gradg
              f2(i,j,k)=(dist-sgnd*dabs(gctmp))/dx
	      iflg1=1
	   endif

c23456789012345678901234567890123456789012345678901234567890123456789012

           gtmp(i,j,k)=gold(i,j,k)+dt*f2(i,j,k)/2.0d0

           if ((dabs(gtmp(i,j,k)-giso).gt.twb).or.
     &        (dabs(dabs(gtmp(i,j,k)-giso)-twb).lt.small)) then
              gtmp(i,j,k)=dsign(twb,(gtmp(i,j,k)-giso))+giso
	   endif

	   if (dabs(gtmp(i,j,k)-gm).lt.small) gtmp(i,j,k)=gm
           if (((gtmp(i,j,k)-gm)*gc).lt.0.0d0) then
              gtmp(i,j,k)=gm
              write (6,*) 'reinwenorkts(rk2):'
	      write (6,1000) i,j,k,nt,iflg1,mask(i,j,k),indx(nt,4),
     &		 g(i,j,k),gm,gc,gradg,grad(i,j,k)
           endif
 150    continue
 	call bound(nx,ny,nz,gtmp,ibcg)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,gtmp,ibcg)
	endif

c	if (igrd.eq.0) then
c	   call enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',3)
c	else
c	   call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',3)
c	endif
 
 	call gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &     nord,g,gtmp,grad,ibcv,'rein',20)

 	do 200 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
c	   if ((indx(nt,4).ne.2).or.(indx(nt,4).ne.1)) goto 200
 	   gm=gjump*dfloat(mask(i,j,k))
           gc=g(i,j,k)-gm
           gctmp=gtmp(i,j,k)-gm
	   if (dabs(gc).lt.small) gc=0.0d0

           sgnd=dsign(1.0d0,gc)
	   if (gc.eq.0.0d0) sgnd=0.0d0

           if (dabs(gc).le.dx) then
              sgn=gc/dsqrt(gc**2.0d0+dx**2.0d0)
c             sgn=gc/dsqrt(gc**2.0d0+(dx*grad(i,j,k))**2.0d0)
c             sgn=gc/dx+dsin(pi*gc/dx)/pi
           else
              sgn=dsign(1.0d0,gc)
           endif
cc         if (dabs(gctmp).le.dx) then
cc            sgn=gctmp/dsqrt(gctmp**2.0d0+dx**2.0d0)
c             sgn=gctmp/dx+dsin(pi*gctmp/dx)/pi
cc         else
cc            sgn=dsign(1.0d0,gctmp)
cc         endif

           f3(i,j,k)=sgnd*(1.0d0-grad(i,j,k))

	   iflg1=0
           if ((iscf.eq.1).and.
     &        ((gc*(g(i-1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i+1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j-1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j+1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k-1)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k+1)-gm).lt.0.0d0))) then

	      if (dabs(gc).gt.small) then
	         if (grad0(i,j,k).ge.(dabs(gc)/dx)) then
                    dist=gc/grad0(i,j,k)
	         else
                    gxm1=(g(i,j,k)-g(i-1,j,k))/dx
                    gxp1=(g(i+1,j,k)-g(i,j,k))/dx
                    gym1=(g(i,j,k)-g(i,j-1,k))/dy
                    gyp1=(g(i,j+1,k)-g(i,j,k))/dy
                    gzm1=(g(i,j,k)-g(i,j,k-1))/dz
                    gzp1=(g(i,j,k+1)-g(i,j,k))/dz
		    gxc1=0.5d0*(g(i+1,j,k)-g(i-1,j,k))/dx
		    gyc1=0.5d0*(g(i,j+1,k)-g(i,j-1,k))/dy
		    gzc1=0.5d0*(g(i,j,k+1)-g(i,j,k-1))/dz
	            gx2=dmax1(gxm1**2.0d0,gxp1**2.0d0,gxc1**2.0d0)
	            gy2=dmax1(gym1**2.0d0,gyp1**2.0d0,gyc1**2.0d0)
	            gz2=dmax1(gzm1**2.0d0,gzp1**2.0d0,gzc1**2.0d0)
	            gradg=dmax1(dabs(gc)/dx,dsqrt(gx2+gy2+gz2))
                    dist=gc/gradg
	         endif
	      else
	         dist=0.0d0
	      endif

c             gradg=dmax1(grad0(i,j,k),small)
c             dist=gc/gradg
              f3(i,j,k)=(dist-sgnd*dabs(gctmp))/dx
	      iflg1=1
	   endif

c23456789012345678901234567890123456789012345678901234567890123456789012

           gtmp(i,j,k)=gold(i,j,k)+dt*f3(i,j,k)

           if ((dabs(gtmp(i,j,k)-giso).gt.twb).or.
     &        (dabs(dabs(gtmp(i,j,k)-giso)-twb).lt.small)) then
              gtmp(i,j,k)=dsign(twb,(gtmp(i,j,k)-giso))+giso
	   endif

	   if (dabs(gtmp(i,j,k)-gm).lt.small) gtmp(i,j,k)=gm
           if (((gtmp(i,j,k)-gm)*gc).lt.0.0d0) then
              gtmp(i,j,k)=gm
              write (6,*) 'reinwenorkts(rk3):'
	      write (6,1000) i,j,k,nt,iflg1,mask(i,j,k),indx(nt,4),
     &		 g(i,j,k),gm,gc,gradg,grad(i,j,k)
           endif
 200    continue
 	call bound(nx,ny,nz,gtmp,ibcg)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,gtmp,ibcg)
	endif

c	if (igrd.eq.0) then
c	   call enogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',4)
c	else
c	   call wenogradts(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
c    &        nord,g,gtmp,grad,ibcv,'rein',4)
c	endif

 	call gradfirst(nx,ny,nz,dx,dy,dz,nptot,ntc,indx,mask,
     &     nord,g,gtmp,grad,ibcv,'rein',20)

 	do 250 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
c	   if ((indx(nt,4).ne.2).or.(indx(nt,4).ne.1)) goto 250
 	   gm=gjump*dfloat(mask(i,j,k))
           gc=g(i,j,k)-gm
           gctmp=gtmp(i,j,k)-gm
	   if (dabs(gc).lt.small) gc=0.0d0

           sgnd=dsign(1.0d0,gc)
	   if (gc.eq.0.0d0) sgnd=0.0d0

           if (dabs(gc).le.dx) then
              sgn=gc/dsqrt(gc**2.0d0+dx**2.0d0)
c             sgn=gc/dsqrt(gc**2.0d0+(dx*grad(i,j,k))**2.0d0)
c             sgn=gc/dx+dsin(pi*gc/dx)/pi
           else
              sgn=dsign(1.0d0,gc)
           endif
cc         if (dabs(gctmp).le.dx) then
cc            sgn=gctmp/dsqrt(gctmp**2.0d0+dx**2.0d0)
c             sgn=gctmp/dx+dsin(pi*gctmp/dx)/pi
cc         else
cc            sgn=dsign(1.0d0,gctmp)
cc         endif

           f4=sgnd*(1.0d0-grad(i,j,k))

	   iflg1=0
           if ((iscf.eq.1).and.
     &        ((gc*(g(i-1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i+1,j,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j-1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j+1,k)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k-1)-gm).lt.0.0d0).or.
     &         (gc*(g(i,j,k+1)-gm).lt.0.0d0))) then

c23456789012345678901234567890123456789012345678901234567890123456789012

	      if (dabs(gc).gt.small) then
	         if (grad0(i,j,k).ge.(dabs(gc)/dx)) then
                    dist=gc/grad0(i,j,k)
	         else
                    gxm1=(g(i,j,k)-g(i-1,j,k))/dx
                    gxp1=(g(i+1,j,k)-g(i,j,k))/dx
                    gym1=(g(i,j,k)-g(i,j-1,k))/dy
                    gyp1=(g(i,j+1,k)-g(i,j,k))/dy
                    gzm1=(g(i,j,k)-g(i,j,k-1))/dz
                    gzp1=(g(i,j,k+1)-g(i,j,k))/dz
		    gxc1=0.5d0*(g(i+1,j,k)-g(i-1,j,k))/dx
		    gyc1=0.5d0*(g(i,j+1,k)-g(i,j-1,k))/dy
		    gzc1=0.5d0*(g(i,j,k+1)-g(i,j,k-1))/dz
	            gx2=dmax1(gxm1**2.0d0,gxp1**2.0d0,gxc1**2.0d0)
	            gy2=dmax1(gym1**2.0d0,gyp1**2.0d0,gyc1**2.0d0)
	            gz2=dmax1(gzm1**2.0d0,gzp1**2.0d0,gzc1**2.0d0)
	            gradg=dmax1(dabs(gc)/dx,dsqrt(gx2+gy2+gz2))
                    dist=gc/gradg
	         endif
	      else
	         dist=0.0d0
	      endif

c             gradg=dmax1(grad0(i,j,k),small)
c             dist=gc/gradg
              f4=(dist-sgnd*dabs(gctmp))/dx
	      iflg1=1
	   endif

           gtmp(i,j,k)=gold(i,j,k)
     &        +dt*(f1(i,j,k)+2.0d0*(f2(i,j,k)+f3(i,j,k))+f4)/6.0d0

           if ((dabs(gtmp(i,j,k)-giso).gt.twb).or.
     &        (dabs(dabs(gtmp(i,j,k)-giso)-twb).lt.small)) then
              gtmp(i,j,k)=dsign(twb,(gtmp(i,j,k)-giso))+giso
	   endif

	   if (dabs(gtmp(i,j,k)-gm).lt.small) gtmp(i,j,k)=gm
           if (((gtmp(i,j,k)-gm)*gc).lt.0.0d0) then
              gtmp(i,j,k)=gm
              write (6,*) 'reinwenorkts(rk4):'
	      write (6,1000) i,j,k,nt,iflg1,mask(i,j,k),indx(nt,4),
     &		 g(i,j,k),gm,gc,gradg,grad(i,j,k)
           endif
 250    continue
 	call bound(nx,ny,nz,gtmp,ibcg)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,gtmp,ibcg)
	endif

 300    continue

 	do 350 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
c	   if ((indx(nt,4).ne.2).or.(indx(nt,4).ne.1)) goto 350
 	   if (dabs(gtmp(i,j,k)).lt.small) gtmp(i,j,k)=0.0d0
           g(i,j,k)=gtmp(i,j,k)
 350    continue
 	call bound(nx,ny,nz,g,ibcg)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,g,ibcg)
	endif

c       do 400 nt=1,ntc
c          i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
c          giso=gjump*dfloat(mask(i,j,k))
c          if ((dabs(g(i,j,k)-giso).gt.twb).or.
c    &        (dabs(dabs(g(i,j,k)-giso)-twb).lt.small)) then
c             g(i,j,k)=dsign(twb,(g(i,j,k)-giso))+giso
c          endif
c 400   continue

 1000   format ('reinwenorkts: sgn chg at ',7(1x,i5),5(2x,e18.12))
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine isosurf(nx,ny,nz,x,y,z,g,area,iflag)
        implicit double precision (a-h,o-z)
        dimension x(nx+1),y(ny+1),z(nz+1)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension matcon1(256,16),matcon2(12,2)
        dimension xc(8),yc(8),zc(8),fc(8),xt(16),yt(16),zt(16)
        LOGICAL :: logic1
        CHARACTER*15 filenam,isofile,entrada,str
        ALLOCATABLE :: vectiso(:),isofile(:),numvert(:),srfarea(:)
        common /rparam/small,pi
        common /isodat/matcon1,matcon2
	common /tubedat/twa,twb,gjump,nomask

        numiso=0
	gminL=minval(g(1:nx+1,1:ny+1,1:nz+1)/gjump)
	gmaxR=maxval(g(1:nx+1,1:ny+1,1:nz+1)/gjump)
c	gminL=minval(g(1:nx,1:ny,1:nz)/gjump)
c	gmaxR=maxval(g(1:nx,1:ny,1:nz)/gjump)
c       gminL=minval(g/gjump)
c       gmaxR=maxval(g/gjump)
        igstrt=ceiling(gminL)
        igstop=floor(gmaxR)
        do i=igstrt,igstop
           numiso=numiso+1
        end do

        ALLOCATE (vectiso(numiso),isofile(numiso),numvert(numiso),
     &     srfarea(numiso))

        do 100 i=1,numiso
           vectiso(i)=gjump*dfloat(igstrt)+dfloat(i-1)*gjump

           if (iflag.eq.1) then
              write (str,'(I10)') i
              str=adjustl(str)
              isofile(i)='patch'//trim(str)//'.txt'
              filenam = isofile(i)
              open (UNIT=i, FILE=filenam, STATUS='replace')
              write (i,*) ' isonum','     giso'
              write (i,FMT="(i6,7x,e18.12)") i,vectiso(i)
              write (i,*) 'Coordinates'
           end if
 100    continue
        
        write (6,*) (vectiso(i),i=1,numiso), gminL, gmaxR

c23456789012345678901234567890123456789012345678901234567890123456789012

        numvert = 1
        srfarea=0.0d0
        do 300 i=1,nx
           xrht=x(i+1); xlft=x(i)
           xc(1)=xrht; xc(2)=xrht; xc(3)=xrht; xc(4)=xrht
           xc(5)=xlft; xc(6)=xlft; xc(7)=xlft; xc(8)=xlft
        do 300 j=1,ny
           ybot=y(j); ytop=y(j+1)
           yc(1)=ybot; yc(2)=ybot; yc(5)=ybot; yc(6)=ybot
           yc(3)=ytop; yc(4)=ytop; yc(7)=ytop; yc(8)=ytop
        do 300 k=1,nz
           zbck=z(k); zfnt=z(k+1)
           zc(1)=zbck; zc(4)=zbck; zc(5)=zbck; zc(8)=zbck
           zc(2)=zfnt; zc(3)=zfnt; zc(6)=zfnt; zc(7)=zfnt

           fc(1)=g(i+1,j,k)
           fc(2)=g(i+1,j,k+1)
           fc(3)=g(i+1,j+1,k+1)
           fc(4)=g(i+1,j+1,k)
           fc(5)=g(i,j,k)
           fc(6)=g(i,j,k+1)
           fc(7)=g(i,j+1,k+1)
           fc(8)=g(i,j+1,k)

c23456789012345678901234567890123456789012345678901234567890123456789012

           do 200 isonum=1,numiso
              value = vectiso(isonum)
              indexcube = 1
              if (fc(1).lt.value) indexcube = indexcube + 1
              if (fc(2).lt.value) indexcube = indexcube + 2
              if (fc(3).lt.value) indexcube = indexcube + 4
              if (fc(4).lt.value) indexcube = indexcube + 8
              if (fc(5).lt.value) indexcube = indexcube + 16
              if (fc(6).lt.value) indexcube = indexcube + 32
              if (fc(7).lt.value) indexcube = indexcube + 64
              if (fc(8).lt.value) indexcube = indexcube + 128
              !
              logic1=.true.

c23456789012345678901234567890123456789012345678901234567890123456789012

              icol = 1
              do WHILE (logic1)
                 numedge = matcon1(indexcube,icol)
                 if (numedge.eq.0) then
                    logic1=.false.
                 else
                    num1 = matcon2(numedge,1)
                    num2 = matcon2(numedge,2)

                    xt(icol)=xc(num1)+(value-fc(num1))
     &                 *(xc(num2)-xc(num1))/(fc(num2)-fc(num1))
                    yt(icol)=yc(num1)+(value-fc(num1))
     &                 *(yc(num2)-yc(num1))/(fc(num2)-fc(num1))
                    zt(icol)=zc(num1)+(value-fc(num1))
     &                 *(zc(num2)-zc(num1))/(fc(num2)-fc(num1))

                    if (mod(icol,3).eq.0) then
                       s1=dsqrt((xt(icol-2)-xt(icol-1))**2.0d0
     &                    +(yt(icol-2)-yt(icol-1))**2.0d0
     &                    +(zt(icol-2)-zt(icol-1))**2.0d0)
                       s2=dsqrt((xt(icol-1)-xt(icol))**2.0d0
     &                    +(yt(icol-1)-yt(icol))**2.0d0
     &                    +(zt(icol-1)-zt(icol))**2.0d0)
                       s3=dsqrt((xt(icol-2)-xt(icol))**2.0d0
     &                    +(yt(icol-2)-yt(icol))**2.0d0
     &                    +(zt(icol-2)-zt(icol))**2.0d0)
                       st=0.5d0*(s1+s2+s3)
c                      srfarea(isonum)=srfarea(isonum)
c    &                    +dsqrt(st*(st-s1)*(st-s2)*(st-s3))
                       temp=st*(st-s1)*(st-s2)*(st-s3)
		       if (dabs(temp).lt.small) temp=0.0d0
                       srfarea(isonum)=srfarea(isonum)+dsqrt(temp)
                    end if
                    numvert(isonum) = numvert(isonum) + 1
                    icol = icol + 1

                    if (iflag.eq.1) then
                       write (isonum,FMT="(i6,3(1x,f10.2))") 
     &                    numvert(isonum),xt,yt,zt
                    end if

                 end if
              end do
 200       continue
 300    continue
        
        area=0.0d0
        do i=1,numiso
           area=area+srfarea(i)
        end do
        area=area/((y(ny+1)-y(1))*(z(nz+1)-z(1)))

c23456789012345678901234567890123456789012345678901234567890123456789012

        numvert = numvert - 1
        if (iflag.eq.1) then
           do 500 i=1,numiso
              write (i,*) 'end coordinates'
              write (i,*)
              write (i,*) 'Elements'
              nfaces = numvert(i)/3
              ia = 1 ; ib = 2 ; ic = 3
           
              do 400 j=1,nfaces
                 write (i,FMT="(4(1x,i6))") j,ia,ib,ic
                 ia = ia + 3 ; ib = ib + 3 ; ic = ic + 3
 400          continue

              write (i,*) 'end elements'
              close (i)
 500       continue
        end if
        deallocate (vectiso,isofile,numvert)
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine isosurfts(nt0,nx,ny,nz,x,y,z,nptot,ntc,indx,mask,
     &     igstrt,numiso,g,area,localmax,iflag,cx,cy,cz,igmin,igmax)
        implicit double precision (a-h,o-z)
        dimension x(nx+1),y(ny+1),z(nz+1)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension matcon1(256,16),matcon2(12,2)
        dimension xc(8),yc(8),zc(8),fc(8),xt(16),yt(16),zt(16)
        dimension mask(nx,ny,nz),indx(nptot,4)
     	dimension cx(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cy(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cz(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension igmin(ny,nz),igmax(ny,nz)
        LOGICAL :: logic1
        CHARACTER*15 filenam,isofile,entrada,str
        ALLOCATABLE :: vectiso(:),isofile(:),numvert(:),srfarea(:)
        common /rparam/small,pi
        common /isodat/matcon1,matcon2
	common /tubedat/twa,twb,gjump,nomask
	common /grid/icush,imin,imax,izone
        common /iparam/iquad

        numiso=0

c	gminL=minval(g(1:nx+1,1:ny+1,1:nz+1)/gjump)
c	gmaxR=maxval(g(1:nx+1,1:ny+1,1:nz+1)/gjump)

c       igstrt=ceiling(gminL)
c       igstop=floor(gmaxR)
c       do i=igstrt,igstop
c          numiso=numiso+1
c       end do

        igstrt=0; igstop=0; numiso=1

        ALLOCATE (vectiso(numiso),isofile(numiso),numvert(numiso),
     &     srfarea(numiso))

        do 100 i=1,numiso
           vectiso(i)=gjump*dfloat(igstrt)+dfloat(i-1)*gjump

           if (iflag.eq.1) then
              write (str,'(I10)') i
              str=adjustl(str)
              isofile(i)='patch'//trim(str)//'.txt'
              filenam = isofile(i)
              open (UNIT=i, FILE=filenam, STATUS='replace')
              write (i,*) ' isonum','     giso'
              write (i,FMT="(i6,7x,e18.12)") i,vectiso(i)
              write (i,*) 'Coordinates'
	   else if (iflag.eq.2) then
	      localmax=0
           end if
 100    continue
        
c23456789012345678901234567890123456789012345678901234567890123456789012

        numvert = 1
        srfarea=0.0d0

 	do 300 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
 	   if (indx(nt,4).ne.2) goto 300
 	   if ((iquad.eq.1).and.((j.eq.ny).or.(k.eq.nz))) goto 300

c       do 300 i=1,nx
           xrht=x(i+1); xlft=x(i)
           xc(1)=xrht; xc(2)=xrht; xc(3)=xrht; xc(4)=xrht
           xc(5)=xlft; xc(6)=xlft; xc(7)=xlft; xc(8)=xlft
c       do 300 j=1,ny
           ybot=y(j); ytop=y(j+1)
           yc(1)=ybot; yc(2)=ybot; yc(5)=ybot; yc(6)=ybot
           yc(3)=ytop; yc(4)=ytop; yc(7)=ytop; yc(8)=ytop
c       do 300 k=1,nz
           zbck=z(k); zfnt=z(k+1)
           zc(1)=zbck; zc(4)=zbck; zc(5)=zbck; zc(8)=zbck
           zc(2)=zfnt; zc(3)=zfnt; zc(6)=zfnt; zc(7)=zfnt

           fc(1)=g(i+1,j,k)
           fc(2)=g(i+1,j,k+1)
           fc(3)=g(i+1,j+1,k+1)
           fc(4)=g(i+1,j+1,k)
           fc(5)=g(i,j,k)
           fc(6)=g(i,j,k+1)
           fc(7)=g(i,j+1,k+1)
           fc(8)=g(i,j+1,k)

c23456789012345678901234567890123456789012345678901234567890123456789012

           do 200 isonum=1,numiso
              value = vectiso(isonum)
              indexcube = 1
              if ((fc(1)-value).lt.-small) indexcube = indexcube + 1
              if ((fc(2)-value).lt.-small) indexcube = indexcube + 2
              if ((fc(3)-value).lt.-small) indexcube = indexcube + 4
              if ((fc(4)-value).lt.-small) indexcube = indexcube + 8
              if ((fc(5)-value).lt.-small) indexcube = indexcube + 16
              if ((fc(6)-value).lt.-small) indexcube = indexcube + 32
              if ((fc(7)-value).lt.-small) indexcube = indexcube + 64
              if ((fc(8)-value).lt.-small) indexcube = indexcube + 128
              !
              logic1=.true.

c	      if ((indexcube.ne.1).and.(indexcube.ne.256)) then
c                write(6,*) ""
c                write(6,185) xlft,xrht,ybot,ytop,zbck,zfnt,indexcube
c                write(6,190) fc(5),fc(6),fc(7),fc(8),fc(1),fc(2),
c    &              fc(3),fc(4)
c 185            format (6(1x,e18.12),i5)
c 190            format (8(1x,e18.12))
c	      endif

c23456789012345678901234567890123456789012345678901234567890123456789012

              icol = 1
              do WHILE (logic1)
                 numedge = matcon1(indexcube,icol)
                 if (numedge.eq.0) then
                    logic1=.false.
                 else
                    num1 = matcon2(numedge,1)
                    num2 = matcon2(numedge,2)

                    xt(icol)=xc(num1)+(value-fc(num1))
     &                 *(xc(num2)-xc(num1))/(fc(num2)-fc(num1))
                    yt(icol)=yc(num1)+(value-fc(num1))
     &                 *(yc(num2)-yc(num1))/(fc(num2)-fc(num1))
                    zt(icol)=zc(num1)+(value-fc(num1))
     &                 *(zc(num2)-zc(num1))/(fc(num2)-fc(num1))

                    if (mod(icol,3).eq.0) then
                       s1=dsqrt((xt(icol-2)-xt(icol-1))**2.0d0
     &                    +(yt(icol-2)-yt(icol-1))**2.0d0
     &                    +(zt(icol-2)-zt(icol-1))**2.0d0)
                       s2=dsqrt((xt(icol-1)-xt(icol))**2.0d0
     &                    +(yt(icol-1)-yt(icol))**2.0d0
     &                    +(zt(icol-1)-zt(icol))**2.0d0)
                       s3=dsqrt((xt(icol-2)-xt(icol))**2.0d0
     &                    +(yt(icol-2)-yt(icol))**2.0d0
     &                    +(zt(icol-2)-zt(icol))**2.0d0)
                       st=0.5d0*(s1+s2+s3)
c                      srfarea(isonum)=srfarea(isonum)
c    &                    +dsqrt(st*(st-s1)*(st-s2)*(st-s3))
                       temp=st*(st-s1)*(st-s2)*(st-s3)
		       if (dabs(temp).lt.small) temp=0.0d0
                       srfarea(isonum)=srfarea(isonum)+dsqrt(temp)

c                      write(6,195) xt(icol-2),yt(icol-2),zt(icol-2),
c    &                    xt(icol-1),yt(icol-1),zt(icol-1),
c    &			  xt(icol),yt(icol),zt(icol),dsqrt(temp),
c    &			  numvert(isonum)/3
c 195                  format (3(/,3(1x,e18.12)),1x,e18.12,1x,i5)

                    end if
                    numvert(isonum) = numvert(isonum) + 1
                    icol = icol + 1

                    if (iflag.eq.1) then
                       write (isonum,FMT="(i6,3(1x,f10.2))") 
     &                    numvert(isonum),xt,yt,zt
                    end if

                 end if
              end do
 200       continue

 	gc = g(i,j,k)
 	if ((dabs(gc).gt.small).and.(iflag.eq.2).and.
     &	   ((gc-g(i+1,j,k)).gt.small).and.
     &	   ((gc-g(i-1,j,k)).gt.small).and.
     &	   ((gc-g(i,j+1,k)).gt.small).and.
     &	   ((gc-g(i,j-1,k)).gt.small).and.
     &	   ((gc-g(i,j,k+1)).gt.small).and.
     &	   ((gc-g(i,j,k-1)).gt.small)) then

c    &	   ((gc-g(i+1,j-1,k)).gt.small).and.
c    &	   ((gc-g(i-1,j-1,k)).gt.small).and.
c    &	   ((gc-g(i+1,j+1,k)).gt.small).and.
c    &	   ((gc-g(i-1,j+1,k)).gt.small).and.

c    &	   ((gc-g(i+1,j-1,k-1)).gt.small).and.
c    &	   ((gc-g(i-1,j-1,k-1)).gt.small).and.
c    &	   ((gc-g(i+1,j+1,k-1)).gt.small).and.
c    &	   ((gc-g(i-1,j+1,k-1)).gt.small).and.

c    &	   ((gc-g(i+1,j-1,k+1)).gt.small).and.
c    &	   ((gc-g(i-1,j-1,k+1)).gt.small).and.
c    &	   ((gc-g(i+1,j+1,k+1)).gt.small).and.
c    &	   ((gc-g(i-1,j+1,k+1)).gt.small)) then

 	   ishft=(nt0-nx)/2
           if ((i+ishft).gt.0) then
              iv=mod((i+ishft-1),nt0)+1
           else
              iv=nt0+mod((i+ishft),nt0)
           endif

      	   localmax=localmax+1
cc	   write (6,201) i,igmin(j,k),igmax(j,k),localmax,g(i-1,j,k),
c	   write (6,201) iv,j,k,localmax,g(i-1,j,k),
cc   &	      g(i,j,k),g(i+1,j,k),x(1),x(i),cx(iv-1,j,k),cx(iv,j,k),
cc   &	      cx(iv+1,j,k)
c	   write (6,201) iv,j,k,localmax,g(i-1,j,k),g(i,j,k),
c    &	      g(i+1,j,k),g(i,j-1,k),g(i,j+1,k),g(i,j,k-1),g(i,j,k+1)
 	endif
 201    format (4(1x,i3),8(1x,e18.12))

 300    continue
        
        area=0.0d0
        do i=1,numiso
           area=area+srfarea(i)
        end do
	if (iquad.eq.0) then
           area=area/((y(ny+1)-y(1))*(z(nz+1)-z(1)))
	else
           area=area/((y(ny)-y(1))*(z(nz)-z(1)))
	endif

c23456789012345678901234567890123456789012345678901234567890123456789012

        numvert = numvert - 1
        if (iflag.eq.1) then
           do 500 i=1,numiso
              write (i,*) 'end coordinates'
              write (i,*)
              write (i,*) 'Elements'
              nfaces = numvert(i)/3
              ia = 1 ; ib = 2 ; ic = 3
           
              do 400 j=1,nfaces
                 write (i,FMT="(4(1x,i6))") j,ia,ib,ic
                 ia = ia + 3 ; ib = ib + 3 ; ic = ic + 3
 400          continue

              write (i,*) 'end elements'
              close (i)
 500       continue
        end if
        deallocate (vectiso,isofile,numvert)
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine isoinit(matcon1,matcon2)
        implicit double precision (a-h,o-z)
        dimension matcon1(256,16),matcon2(12,2)
        dimension matvect(4096)

        matvect(1:512)=(/
     &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  2, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  9,  4, 10,  9,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  4,  2,  3, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  10,  3, 11,  1,  3, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   3,  9,  4,  3, 11,  9, 11, 10,  9,  0,  0,  0,  0,  0,  0,  0,
     &   4, 12,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1, 12,  3,  9, 12,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2, 10,  1,  3,  4, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2, 12,  3,  2, 10, 12, 10,  9, 12,  0,  0,  0,  0,  0,  0,  0,
     &   4, 11,  2, 12, 11,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1, 11,  2,  1,  9, 11,  9, 12, 11,  0,  0,  0,  0,  0,  0,  0,
     &   4, 10,  1,  4, 12, 10, 12, 11, 10,  0,  0,  0,  0,  0,  0,  0,
     &  10,  9, 11, 11,  9, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5,  8,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5,  4,  1,  8,  4,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  2, 10,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5,  2, 10,  5,  8,  2,  8,  4,  2,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 11,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   4,  5,  8,  4,  1,  5,  2,  3, 11,  0,  0,  0,  0,  0,  0,  0,
     &  10,  3, 11, 10,  1,  3,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,
     &   3, 11, 10,  3, 10,  8,  3,  8,  4,  8, 10,  5,  0,  0,  0,  0,
     &   9,  5,  8,  4, 12,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  12,  5,  8, 12,  3,  5,  3,  1,  5,  0,  0,  0,  0,  0,  0,  0,
     &  10,  1,  2,  9,  5,  8,  3,  4, 12,  0,  0,  0,  0,  0,  0,  0,
     &   5,  8, 12, 10,  5, 12, 10, 12,  3, 10,  3,  2,  0,  0,  0,  0,
     &   4, 11,  2,  4, 12, 11,  8,  9,  5,  0,  0,  0,  0,  0,  0,  0,
     &   2, 12, 11,  2,  5, 12,  2,  1,  5,  8, 12,  5,  0,  0,  0,  0,
     &   5,  8,  9, 10,  1, 12, 10, 12, 11, 12,  1,  4,  0,  0,  0,  0,
     &   5,  8, 12,  5, 12, 10, 10, 12, 11,  0,  0,  0,  0,  0,  0,  0/
     &  )
        matvect(513:1024)=(/
     &  10,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  10,  6,  5,  1,  9,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  6,  5,  2,  6,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   9,  6,  5,  9,  4,  6,  4,  2,  6,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 11, 10,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   4,  1,  9,  2,  3, 11,  5, 10,  6,  0,  0,  0,  0,  0,  0,  0,
     &   6,  3, 11,  6,  5,  3,  5,  1,  3,  0,  0,  0,  0,  0,  0,  0,
     &   3, 11,  6,  4,  3,  6,  4,  6,  5,  4,  5,  9,  0,  0,  0,  0,
     &  10,  6,  5,  3,  4, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1, 12,  3,  1,  9, 12,  5, 10,  6,  0,  0,  0,  0,  0,  0,  0,
     &   1,  6,  5,  1,  2,  6,  3,  4, 12,  0,  0,  0,  0,  0,  0,  0,
     &   3,  2,  6,  3,  6,  9,  3,  9, 12,  5,  9,  6,  0,  0,  0,  0,
     &  11,  4, 12, 11,  2,  4, 10,  6,  5,  0,  0,  0,  0,  0,  0,  0,
     &   5, 10,  6,  1,  9,  2,  9, 11,  2,  9, 12, 11,  0,  0,  0,  0,
     &   6,  5,  1,  6,  1, 12,  6, 12, 11, 12,  1,  4,  0,  0,  0,  0,
     &   6,  5,  9,  6,  9, 11, 11,  9, 12,  0,  0,  0,  0,  0,  0,  0,
     &  10,  8,  9,  6,  8, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  10,  4,  1, 10,  6,  4,  6,  8,  4,  0,  0,  0,  0,  0,  0,  0,
     &   1,  8,  9,  1,  2,  8,  2,  6,  8,  0,  0,  0,  0,  0,  0,  0,
     &   2,  6,  4,  4,  6,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  10,  8,  9, 10,  6,  8, 11,  2,  3,  0,  0,  0,  0,  0,  0,  0,
     &  11,  2,  3, 10,  6,  1,  6,  4,  1,  6,  8,  4,  0,  0,  0,  0,
     &   9,  1,  3,  9,  3,  6,  9,  6,  8, 11,  6,  3,  0,  0,  0,  0,
     &   3, 11,  6,  3,  6,  4,  4,  6,  8,  0,  0,  0,  0,  0,  0,  0,
     &   8, 10,  6,  8,  9, 10,  4, 12,  3,  0,  0,  0,  0,  0,  0,  0,
     &  10,  6,  8, 10,  8,  3, 10,  3,  1,  3,  8, 12,  0,  0,  0,  0,
     &   3,  4, 12,  1,  2,  9,  2,  8,  9,  2,  6,  8,  0,  0,  0,  0,
     &  12,  3,  2, 12,  2,  8,  8,  2,  6,  0,  0,  0,  0,  0,  0,  0,
     &  10,  6,  9,  9,  6,  8, 11,  2,  4, 11,  4, 12,  0,  0,  0,  0,
     &   6,  8,  1,  6,  1, 10,  8, 12,  1,  2,  1, 11, 12, 11,  1,  0,
     &  12, 11,  1, 12,  1,  4, 11,  6,  1,  9,  1,  8,  6,  8,  1,  0,
     &  12, 11,  6,  8, 12,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     &  )
        matvect(1025:1536)=(/
     &  11,  7,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  4,  6, 11,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  10,  1,  2,  6, 11,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  9,  4,  2, 10,  9,  6, 11,  7,  0,  0,  0,  0,  0,  0,  0,
     &   2,  7,  6,  3,  7,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  7,  6,  2,  3,  7,  4,  1,  9,  0,  0,  0,  0,  0,  0,  0,
     &  10,  7,  6, 10,  1,  7,  1,  3,  7,  0,  0,  0,  0,  0,  0,  0,
     &   6, 10,  9,  6,  9,  3,  6,  3,  7,  4,  3,  9,  0,  0,  0,  0,
     &   3,  4, 12, 11,  7,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  12,  1,  9, 12,  3,  1, 11,  7,  6,  0,  0,  0,  0,  0,  0,  0,
     &   1,  2, 10,  3,  4, 12,  6, 11,  7,  0,  0,  0,  0,  0,  0,  0,
     &   6, 11,  7,  2, 10,  3, 10, 12,  3, 10,  9, 12,  0,  0,  0,  0,
     &   7,  4, 12,  7,  6,  4,  6,  2,  4,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9, 12,  1, 12,  6,  1,  6,  2,  6, 12,  7,  0,  0,  0,  0,
     &   4, 12,  7,  1,  4,  7,  1,  7,  6,  1,  6, 10,  0,  0,  0,  0,
     &   7,  6, 10,  7, 10, 12, 12, 10,  9,  0,  0,  0,  0,  0,  0,  0,
     &   6, 11,  7,  5,  8,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5,  4,  1,  5,  8,  4,  7,  6, 11,  0,  0,  0,  0,  0,  0,  0,
     &   2, 10,  1,  6, 11,  7,  9,  5,  8,  0,  0,  0,  0,  0,  0,  0,
     &  11,  7,  6,  2, 10,  8,  2,  8,  4,  8, 10,  5,  0,  0,  0,  0,
     &   7,  2,  3,  7,  6,  2,  5,  8,  9,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3,  6,  6,  3,  7,  4,  1,  5,  4,  5,  8,  0,  0,  0,  0,
     &   9,  5,  8, 10,  1,  6,  1,  7,  6,  1,  3,  7,  0,  0,  0,  0,
     &   8,  4, 10,  8, 10,  5,  4,  3, 10,  6, 10,  7,  3,  7, 10,  0,
     &   4, 12,  3,  8,  9,  5, 11,  7,  6,  0,  0,  0,  0,  0,  0,  0,
     &   6, 11,  7,  5,  8,  3,  5,  3,  1,  3,  8, 12,  0,  0,  0,  0,
     &   1,  2, 10,  5,  8,  9,  3,  4, 12,  6, 11,  7,  0,  0,  0,  0,
     &  10,  3,  2, 10, 12,  3, 10,  5, 12,  8, 12,  5,  6, 11,  7,  0,
     &   9,  5,  8,  4, 12,  6,  4,  6,  2,  6, 12,  7,  0,  0,  0,  0,
     &   6,  2, 12,  6, 12,  7,  2,  1, 12,  8, 12,  5,  1,  5, 12,  0,
     &   1,  6, 10,  1,  7,  6,  1,  4,  7, 12,  7,  4,  9,  5,  8,  0,
     &   7,  6, 10,  7, 10, 12,  5,  8, 10,  8, 12, 10,  0,  0,  0,  0/
     &  )
        matvect(1537:2048)=(/
     &  11,  5, 10,  7,  5, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5, 11,  7,  5, 10, 11,  1,  9,  4,  0,  0,  0,  0,  0,  0,  0,
     &  11,  1,  2, 11,  7,  1,  7,  5,  1,  0,  0,  0,  0,  0,  0,  0,
     &   9,  4,  2,  9,  2,  7,  9,  7,  5,  7,  2, 11,  0,  0,  0,  0,
     &   2,  5, 10,  2,  3,  5,  3,  7,  5,  0,  0,  0,  0,  0,  0,  0,
     &   4,  1,  9,  2,  3, 10,  3,  5, 10,  3,  7,  5,  0,  0,  0,  0,
     &   1,  3,  5,  5,  3,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   9,  4,  3,  9,  3,  5,  5,  3,  7,  0,  0,  0,  0,  0,  0,  0,
     &  11,  5, 10, 11,  7,  5, 12,  3,  4,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  3,  3,  9, 12,  5, 10, 11,  5, 11,  7,  0,  0,  0,  0,
     &   4, 12,  3,  1,  2,  7,  1,  7,  5,  7,  2, 11,  0,  0,  0,  0,
     &   7,  5,  2,  7,  2, 11,  5,  9,  2,  3,  2, 12,  9, 12,  2,  0,
     &  10,  7,  5, 10,  4,  7, 10,  2,  4, 12,  7,  4,  0,  0,  0,  0,
     &   9, 12,  2,  9,  2,  1, 12,  7,  2, 10,  2,  5,  7,  5,  2,  0,
     &   4, 12,  7,  4,  7,  1,  1,  7,  5,  0,  0,  0,  0,  0,  0,  0,
     &   7,  5,  9, 12,  7,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   8, 11,  7,  8,  9, 11,  9, 10, 11,  0,  0,  0,  0,  0,  0,  0,
     &   1,  8,  4,  1, 11,  8,  1, 10, 11,  7,  8, 11,  0,  0,  0,  0,
     &  11,  7,  8,  2, 11,  8,  2,  8,  9,  2,  9,  1,  0,  0,  0,  0,
     &  11,  7,  8, 11,  8,  2,  2,  8,  4,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3,  7,  2,  7,  9,  2,  9, 10,  9,  7,  8,  0,  0,  0,  0,
     &   3,  7, 10,  3, 10,  2,  7,  8, 10,  1, 10,  4,  8,  4, 10,  0,
     &   8,  9,  1,  8,  1,  7,  7,  1,  3,  0,  0,  0,  0,  0,  0,  0,
     &   8,  4,  3,  7,  8,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   3,  4, 12, 11,  7,  9, 11,  9, 10,  9,  7,  8,  0,  0,  0,  0,
     &   3,  1,  8,  3,  8, 12,  1, 10,  8,  7,  8, 11, 10, 11,  8,  0,
     &   2,  9,  1,  2,  8,  9,  2, 11,  8,  7,  8, 11,  3,  4, 12,  0,
     &  12,  3,  2, 12,  2,  8, 11,  7,  2,  7,  8,  2,  0,  0,  0,  0,
     &   9, 10,  7,  9,  7,  8, 10,  2,  7, 12,  7,  4,  2,  4,  7,  0,
     &   1, 10,  2, 12,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   8,  9,  1,  8,  1,  7,  4, 12,  1, 12,  7,  1,  0,  0,  0,  0,
     &   8, 12,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     &  )
        matvect(2049:2560)=(/
     &   8,  7, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   4,  1,  9, 12,  8,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  2, 10, 12,  8,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   9,  2, 10,  9,  4,  2, 12,  8,  7,  0,  0,  0,  0,  0,  0,  0,
     &  11,  2,  3,  7, 12,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 11,  4,  1,  9,  7, 12,  8,  0,  0,  0,  0,  0,  0,  0,
     &   3, 10,  1,  3, 11, 10,  7, 12,  8,  0,  0,  0,  0,  0,  0,  0,
     &   7, 12,  8,  3, 11,  4, 11,  9,  4, 11, 10,  9,  0,  0,  0,  0,
     &   8,  3,  4,  7,  3,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   8,  1,  9,  8,  7,  1,  7,  3,  1,  0,  0,  0,  0,  0,  0,  0,
     &   3,  8,  7,  3,  4,  8,  1,  2, 10,  0,  0,  0,  0,  0,  0,  0,
     &   2,  7,  3,  2,  9,  7,  2, 10,  9,  9,  8,  7,  0,  0,  0,  0,
     &  11,  8,  7, 11,  2,  8,  2,  4,  8,  0,  0,  0,  0,  0,  0,  0,
     &  11,  8,  7,  2,  8, 11,  2,  9,  8,  2,  1,  9,  0,  0,  0,  0,
     &   1,  4,  8,  1,  8, 11,  1, 11, 10,  7, 11,  8,  0,  0,  0,  0,
     &   8,  7, 11,  8, 11,  9,  9, 11, 10,  0,  0,  0,  0,  0,  0,  0,
     &   7,  9,  5, 12,  9,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   4,  7, 12,  4,  1,  7,  1,  5,  7,  0,  0,  0,  0,  0,  0,  0,
     &   9,  7, 12,  9,  5,  7, 10,  1,  2,  0,  0,  0,  0,  0,  0,  0,
     &  10,  5,  7, 10,  7,  4, 10,  4,  2, 12,  4,  7,  0,  0,  0,  0,
     &   7,  9,  5,  7, 12,  9,  3, 11,  2,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 11,  4,  1, 12,  1,  7, 12,  1,  5,  7,  0,  0,  0,  0,
     &   5, 12,  9,  5,  7, 12,  1,  3, 10,  3, 11, 10,  0,  0,  0,  0,
     &  11, 10,  4, 11,  4,  3, 10,  5,  4, 12,  4,  7,  5,  7,  4,  0,
     &   9,  3,  4,  9,  5,  3,  5,  7,  3,  0,  0,  0,  0,  0,  0,  0,
     &   1,  5,  3,  5,  7,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2, 10,  1,  3,  4,  5,  3,  5,  7,  5,  4,  9,  0,  0,  0,  0,
     &   2, 10,  5,  2,  5,  3,  3,  5,  7,  0,  0,  0,  0,  0,  0,  0,
     &   9,  2,  4,  9,  7,  2,  9,  5,  7,  7, 11,  2,  0,  0,  0,  0,
     &  11,  2,  1, 11,  1,  7,  7,  1,  5,  0,  0,  0,  0,  0,  0,  0,
     &   5,  7,  4,  5,  4,  9,  7, 11,  4,  1,  4, 10, 11, 10,  4,  0,
     &  11, 10,  5,  7, 11,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     &  )
        matvect(2561:3072)=(/
     &   5, 10,  6,  8,  7, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  4,  5, 10,  6, 12,  8,  7,  0,  0,  0,  0,  0,  0,  0,
     &   6,  1,  2,  6,  5,  1,  8,  7, 12,  0,  0,  0,  0,  0,  0,  0,
     &  12,  8,  7,  9,  4,  5,  4,  6,  5,  4,  2,  6,  0,  0,  0,  0,
     &  10,  6,  5, 11,  2,  3,  8,  7, 12,  0,  0,  0,  0,  0,  0,  0,
     &   7, 12,  8,  2,  3, 11,  1,  9,  4,  5, 10,  6,  0,  0,  0,  0,
     &   8,  7, 12,  6,  5, 11,  5,  3, 11,  5,  1,  3,  0,  0,  0,  0,
     &   4,  5,  9,  4,  6,  5,  4,  3,  6, 11,  6,  3, 12,  8,  7,  0,
     &   8,  3,  4,  8,  7,  3,  6,  5, 10,  0,  0,  0,  0,  0,  0,  0,
     &  10,  6,  5,  1,  9,  7,  1,  7,  3,  7,  9,  8,  0,  0,  0,  0,
     &   4,  7,  3,  4,  8,  7,  2,  6,  1,  6,  5,  1,  0,  0,  0,  0,
     &   7,  3,  9,  7,  9,  8,  3,  2,  9,  5,  9,  6,  2,  6,  9,  0,
     &  10,  6,  5, 11,  2,  7,  2,  8,  7,  2,  4,  8,  0,  0,  0,  0,
     &   2,  7, 11,  2,  8,  7,  2,  1,  8,  9,  8,  1, 10,  6,  5,  0,
     &   5,  1, 11,  5, 11,  6,  1,  4, 11,  7, 11,  8,  4,  8, 11,  0,
     &   8,  7, 11,  8, 11,  9,  6,  5, 11,  5,  9, 11,  0,  0,  0,  0,
     &   7, 10,  6,  7, 12, 10, 12,  9, 10,  0,  0,  0,  0,  0,  0,  0,
     &   4,  7, 12,  1,  7,  4,  1,  6,  7,  1, 10,  6,  0,  0,  0,  0,
     &   1, 12,  9,  1,  6, 12,  1,  2,  6,  6,  7, 12,  0,  0,  0,  0,
     &   7, 12,  4,  7,  4,  6,  6,  4,  2,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 11, 10,  6, 12, 10, 12,  9, 12,  6,  7,  0,  0,  0,  0,
     &   1, 12,  4,  1,  7, 12,  1, 10,  7,  6,  7, 10,  2,  3, 11,  0,
     &  12,  9,  6, 12,  6,  7,  9,  1,  6, 11,  6,  3,  1,  3,  6,  0,
     &   7, 12,  4,  7,  4,  6,  3, 11,  4, 11,  6,  4,  0,  0,  0,  0,
     &   6,  9, 10,  6,  3,  9,  6,  7,  3,  4,  9,  3,  0,  0,  0,  0,
     &  10,  6,  7, 10,  7,  1,  1,  7,  3,  0,  0,  0,  0,  0,  0,  0,
     &   2,  6,  9,  2,  9,  1,  6,  7,  9,  4,  9,  3,  7,  3,  9,  0,
     &   2,  6,  7,  3,  2,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  4,  7,  2,  7, 11,  4,  9,  7,  6,  7, 10,  9, 10,  7,  0,
     &  11,  2,  1, 11,  1,  7, 10,  6,  1,  6,  7,  1,  0,  0,  0,  0,
     &   1,  4,  9,  6,  7, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  11,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     &  )
        matvect(3073:3584)=(/
     &  12,  6, 11,  8,  6, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  12,  6, 11, 12,  8,  6,  9,  4,  1,  0,  0,  0,  0,  0,  0,  0,
     &   6, 12,  8,  6, 11, 12,  2, 10,  1,  0,  0,  0,  0,  0,  0,  0,
     &  11,  8,  6, 11, 12,  8, 10,  9,  2,  9,  4,  2,  0,  0,  0,  0,
     &  12,  2,  3, 12,  8,  2,  8,  6,  2,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  4,  2,  3,  8,  2,  8,  6,  8,  3, 12,  0,  0,  0,  0,
     &  10,  8,  6, 10,  3,  8, 10,  1,  3,  3, 12,  8,  0,  0,  0,  0,
     &   8,  6,  3,  8,  3, 12,  6, 10,  3,  4,  3,  9, 10,  9,  3,  0,
     &   3,  6, 11,  3,  4,  6,  4,  8,  6,  0,  0,  0,  0,  0,  0,  0,
     &   9,  3,  1,  9,  6,  3,  9,  8,  6, 11,  3,  6,  0,  0,  0,  0,
     &  10,  1,  2,  6, 11,  4,  6,  4,  8,  4, 11,  3,  0,  0,  0,  0,
     &  10,  9,  3, 10,  3,  2,  9,  8,  3, 11,  3,  6,  8,  6,  3,  0,
     &   2,  4,  6,  4,  8,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  8,  1,  8,  2,  2,  8,  6,  0,  0,  0,  0,  0,  0,  0,
     &  10,  1,  4, 10,  4,  6,  6,  4,  8,  0,  0,  0,  0,  0,  0,  0,
     &  10,  9,  8,  6, 10,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   6,  9,  5,  6, 11,  9, 11, 12,  9,  0,  0,  0,  0,  0,  0,  0,
     &   6,  1,  5,  6, 12,  1,  6, 11, 12, 12,  4,  1,  0,  0,  0,  0,
     &   1,  2, 10,  9,  5, 11,  9, 11, 12, 11,  5,  6,  0,  0,  0,  0,
     &  11, 12,  5, 11,  5,  6, 12,  4,  5, 10,  5,  2,  4,  2,  5,  0,
     &   3,  6,  2,  3,  9,  6,  3, 12,  9,  5,  6,  9,  0,  0,  0,  0,
     &   1,  5, 12,  1, 12,  4,  5,  6, 12,  3, 12,  2,  6,  2, 12,  0,
     &   1,  3,  6,  1,  6, 10,  3, 12,  6,  5,  6,  9, 12,  9,  6,  0,
     &  10,  5,  6,  3, 12,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   3,  6, 11,  4,  6,  3,  4,  5,  6,  4,  9,  5,  0,  0,  0,  0,
     &   6, 11,  3,  6,  3,  5,  5,  3,  1,  0,  0,  0,  0,  0,  0,  0,
     &   4, 11,  3,  4,  6, 11,  4,  9,  6,  5,  6,  9,  1,  2, 10,  0,
     &   6, 11,  3,  6,  3,  5,  2, 10,  3, 10,  5,  3,  0,  0,  0,  0,
     &   9,  5,  6,  9,  6,  4,  4,  6,  2,  0,  0,  0,  0,  0,  0,  0,
     &   1,  5,  6,  2,  1,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   9,  5,  6,  9,  6,  4, 10,  1,  6,  1,  4,  6,  0,  0,  0,  0,
     &  10,  5,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     &  )
        matvect(3585:4096)=(/
     &   5, 12,  8,  5, 10, 12, 10, 11, 12,  0,  0,  0,  0,  0,  0,  0,
     &   1,  9,  4,  5, 10,  8, 10, 12,  8, 10, 11, 12,  0,  0,  0,  0,
     &   2, 11, 12,  2, 12,  5,  2,  5,  1,  8,  5, 12,  0,  0,  0,  0,
     &   4,  2,  5,  4,  5,  9,  2, 11,  5,  8,  5, 12, 11, 12,  5,  0,
     &   5, 12,  8, 10, 12,  5, 10,  3, 12, 10,  2,  3,  0,  0,  0,  0,
     &  10,  8,  5, 10, 12,  8, 10,  2, 12,  3, 12,  2,  1,  9,  4,  0,
     &  12,  8,  5, 12,  5,  3,  3,  5,  1,  0,  0,  0,  0,  0,  0,  0,
     &  12,  8,  5, 12,  5,  3,  9,  4,  5,  4,  3,  5,  0,  0,  0,  0,
     &   3, 10, 11,  3,  8, 10,  3,  4,  8,  8,  5, 10,  0,  0,  0,  0,
     &  10, 11,  8, 10,  8,  5, 11,  3,  8,  9,  8,  1,  3,  1,  8,  0,
     &   4,  8, 11,  4, 11,  3,  8,  5, 11,  2, 11,  1,  5,  1, 11,  0,
     &   2, 11,  3,  9,  8,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5, 10,  2,  5,  2,  8,  8,  2,  4,  0,  0,  0,  0,  0,  0,  0,
     &   5, 10,  2,  5,  2,  8,  1,  9,  2,  9,  8,  2,  0,  0,  0,  0,
     &   5,  1,  4,  8,  5,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   5,  9,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &  10, 11,  9, 11, 12,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   4,  1, 10,  4, 10, 12, 12, 10, 11,  0,  0,  0,  0,  0,  0,  0,
     &   1,  2, 11,  1, 11,  9,  9, 11, 12,  0,  0,  0,  0,  0,  0,  0,
     &   4,  2, 11, 12,  4, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  3, 12,  2, 12, 10, 10, 12,  9,  0,  0,  0,  0,  0,  0,  0,
     &   4,  1, 10,  4, 10, 12,  2,  3, 10,  3, 12, 10,  0,  0,  0,  0,
     &   1,  3, 12,  9,  1, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   4,  3, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   3,  4,  9,  3,  9, 11, 11,  9, 10,  0,  0,  0,  0,  0,  0,  0,
     &  10, 11,  3,  1, 10,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   3,  4,  9,  3,  9, 11,  1,  2,  9,  2, 11,  9,  0,  0,  0,  0,
     &   2, 11,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   2,  4,  9, 10,  2,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1, 10,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   1,  4,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
     &  )

        do i=1,256
           matcon1(i,1:16)=matvect(1+16*(i-1):16*i)
        end do

c Each row of matcon2 provides the two vertices of the edge, with the 
c row number representing the edge number

        matcon2=RESHAPE((/1,2,3,1,5,6,7,5,1,2,3,4,
     &                    2,3,4,4,6,7,8,8,5,6,7,8/),(/12,2/))

        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine tubemask(nx,ny,nz,nptot,ntc,mask,indx,cfn,gt,
     &	   ibcg,giso)
        implicit double precision (a-h,o-z)
        dimension mask(nx,ny,nz),indx(nptot,4)
	dimension gt(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	common /tubedat/twa,twb,gjump,nomask

	ntc=0
        do 100 k=1,nz
        do 100 j=1,ny
        do 100 i=1,nx

 	mask(i,j,k)=0
c	dabg=dabs(gt(i,j,k))
 	dabg=dabs(gt(i,j,k)-giso)

	if (dabg.le.twa) then
	   cfn(i,j,k)=1.0d0
	elseif (dabg.le.twb) then
	   cfn(i,j,k)=((2.0d0*dabg+twb-3.0d0*twa)*(dabg-twb)**2.0d0)
     &        /(twb-twa)**3.0d0
	else
	   cfn(i,j,k)=0.0d0
	endif

	if (dabg.lt.twb) then
	   mask(i,j,k)=2; ntc=ntc+1
	   indx(ntc,1)=i; indx(ntc,2)=j; indx(ntc,3)=k

	elseif ((dabs(gt(i-1,j,k)-giso).lt.twb)
     &	    .or.(dabs(gt(i+1,j,k)-giso).lt.twb)
     &	    .or.(dabs(gt(i,j-1,k)-giso).lt.twb)
     &	    .or.(dabs(gt(i,j+1,k)-giso).lt.twb)
     &	    .or.(dabs(gt(i,j,k-1)-giso).lt.twb)
     &	    .or.(dabs(gt(i,j,k+1)-giso).lt.twb)) then

	   mask(i,j,k)=1; ntc=ntc+1
	   indx(ntc,1)=i; indx(ntc,2)=j; indx(ntc,3)=k
	   gt(i,j,k)=dsign(twb,(gt(i,j,k)-giso))+giso
	   cfn(i,j,k)=0.0d0

	else
	   gt(i,j,k)=dsign(twb,(gt(i,j,k)-giso))+giso
	endif

  100	continue
 	call bound(nx,ny,nz,gt,ibcg)
 	call bound(nx,ny,nz,cfn,ibcv)
        return
        end


c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine tubes(nx,ny,nz,nptot,ntc,mask,indx,cfn,igstrt,
     &	   numiso,gt,ibcg,ibcv,ntube,mstep)
        implicit double precision (a-h,o-z)
        dimension mask(nx,ny,nz),indx(nptot,4)
	dimension gt(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cfn(-2:nx+3,-2:ny+3,-2:nz+3)
	common /tubedat/twa,twb,gjump,nomask
        common /rparam/small,pi

	limit = 1

	if (ntube.eq.0) goto 75
        do 50 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3);
	   if ((indx(nt,4).eq.1).or.(indx(nt,4).eq.2)) then
	      giso=gjump*dfloat(mask(i,j,k))
 	      if (dabs(gt(i,j,k)-giso).lt.small) gt(i,j,k)=giso
              if ((dabs(gt(i,j,k)-giso).gt.twb).or.
     &           (dabs(dabs(gt(i,j,k)-giso)-twb).lt.small)) then
		 gt(i,j,k)=dsign(twb,(gt(i,j,k)-giso))+giso
	      endif
	   endif
  50	continue
  75	continue
 	call bound(nx,ny,nz,gt,ibcg)

 	gmin=minval(gt/gjump)
 	gmax=maxval(gt/gjump)
 	gminL=minval(gt(1,1:ny,1:nz)/gjump)
 	gminR=minval(gt(nx+1,1:ny,1:nz)/gjump)
 	gmaxL=maxval(gt(1,1:ny,1:nz)/gjump)
 	gmaxR=maxval(gt(nx+1,1:ny,1:nz)/gjump)

        igstrt=ceiling(gmin)
        igstop=floor(gmax)
 	if (gmax.gt.(dfloat(igstop)+0.5d0)) igstop=igstop+1 
 	if (gmin.lt.(dfloat(igstrt)-0.5d0)) igstrt=igstrt-1 

        numiso=0; ntc=0; nomask=-10000
        do i=igstrt,igstop
           numiso=numiso+1
        end do

c       numiso=1; ntc=0; nomask=-10000
c	igstrt=0; igstop=0

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 300 k=1,nz
        do 300 j=1,ny
        do 300 i=1,nx

 	mask(i,j,k)=nomask; cfn(i,j,k)=0.0d0

        do 100 isonum=1,numiso
           iso=igstrt+isonum-1
           giso=gjump*dfloat(iso)
 	   dabg=dabs(gt(i,j,k)-giso)

 	   if (cfn(i,j,k).eq.0.0d0) then
 	
 	      if (dabg.le.twa) then
 	         cfn(i,j,k)=1.0d0
 	      elseif (dabg.le.twb) then
 	         cfn(i,j,k)=((2.0d0*dabg+twb-3.0d0*twa)
     &              *(dabg-twb)**2.0d0)/(twb-twa)**3.0d0
 	      else
 	         cfn(i,j,k)=0.0d0
 	      endif
 	
 	   endif

 	   if (mask(i,j,k).eq.nomask) then

	      if (dabg.lt.twb) then
	         mask(i,j,k)=iso; ntc=ntc+1; indx(ntc,4)=2
	         indx(ntc,1)=i; indx(ntc,2)=j; indx(ntc,3)=k
	      elseif ((dabs(gt(i-1,j,k)-giso).lt.twb-small)
     &	          .or.(dabs(gt(i+1,j,k)-giso).lt.twb-small)
     &	          .or.(dabs(gt(i,j-1,k)-giso).lt.twb-small)
     &	          .or.(dabs(gt(i,j+1,k)-giso).lt.twb-small)
     &	          .or.(dabs(gt(i,j,k-1)-giso).lt.twb-small)
     &	          .or.(dabs(gt(i,j,k+1)-giso).lt.twb-small)) then
	         mask(i,j,k)=iso; ntc=ntc+1; indx(ntc,4)=1
	         indx(ntc,1)=i; indx(ntc,2)=j; indx(ntc,3)=k
	         gt(i,j,k)=dsign(twb,(gt(i,j,k)-giso))+giso
	      endif

	   endif

  100	continue

c23456789012345678901234567890123456789012345678901234567890123456789012

	if (mask(i,j,k).eq.nomask) then

           gstrt=gjump*dfloat(igstrt)
	   if (numiso.eq.1) then
	      gt(i,j,k)=dsign(twb,(gt(i,j,k)-gstrt))+gstrt
	   elseif (gt(i,j,k).lt.gstrt) then
	      gt(i,j,k)=gstrt-twb
	   else
              iflaga=0; iflagb=0
              iso=igstrt
	      do while ((iflaga.eq.0).and.(iflagb.eq.0))
                 iso=iso+1
                 gisom=gjump*dfloat(iso-1)
                 gisop=gjump*dfloat(iso)
		 if ((gt(i,j,k).gt.gisom).and.(gt(i,j,k).lt.gisop)) then
	            gt(i,j,k)=0.5d0*(gisom+gisop)
		    iflaga=1
		 endif
		 if (iso.eq.igstop) iflagb=1
	      end do
	      if (iflaga.eq.0) then
                 gstop=gjump*dfloat(igstop)
		 if (gt(i,j,k).gt.gstop) then
		    gt(i,j,k)=gstop+twb
		    iflaga=1
		 else
		   write(6,*) 'tubes: level set not found'
		 endif
	      endif
	   endif

	endif

  300	continue

 1000   format ('tubes:',1x,2(1x,i3),8(2x,e18.12))
 1100   format (50(2x,e18.12))
 1200   format (1x,50(1x,i6))

 	call bound(nx,ny,nz,gt,ibcg)
 	call bound(nx,ny,nz,cfn,ibcv)
	if (limit.eq.1) then
           call ridlocalextrema(nx,ny,nz,nptot,ntc,mask,indx,gt,ibcg)
	endif
        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine cfield(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel,uconst,
     &	   cx_tmp,cy_tmp,cz_tmp,shift,nptot,ntc,indx,igetdt)
        implicit double precision (a-h,o-z)
	dimension x(nx+1),y(ny+1),z(nz+1)
     	dimension cx_tmp(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cy_tmp(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cz_tmp(-2:nx+3,-2:ny+3,-2:nz+3)
        dimension indx(nptot,4)
	common /veldat/umax,sl,umod,tke,kr,nmodes,iso,ifreq,cfl
        common /rparam/small,pi
        common /iparam/iquad
	common /kechk/uprime_rel

        ax=1.0d0; ay=1.0d0; az=1.0d0
	energy=0.0d0
	fm=2.0d0*pi
        umax=0.0d0


 	do 100 nt=1,ntc
 	   i=indx(nt,1); j=indx(nt,2); k=indx(nt,3); 
  	   if (indx(nt,4).ne.2) goto 100

 	   xvel=x(i)-shift*dx*cx_tmp(i,j,k)
  	   yvel=-0.5d0/dfloat(1+iquad)+dfloat(j-1)*dy
     &	      -shift*dy*cy_tmp(i,j,k)
  	   zvel=-0.5d0/dfloat(1+iquad)+dfloat(k-1)*dz
     &	      -shift*dz*cz_tmp(i,j,k)

           adiv=0.0d0
           utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0

           if (nmodes.gt.1) then
              do 50 n3=1,nmodes
                 rn3=dfloat(n3); wn3=fm*rn3
              do 50 n2=1,nmodes
                 rn2=dfloat(n2); wn2=fm*rn2
              do 50 n1=1,nmodes
                 if ((n1.eq.n2).and.(n2.eq.n3)) goto 50
                 rn1=dfloat(n1); wn1=fm*rn1
c                rnmag=dsqrt(rn1**2.0d0+rn2**2.0d0+rn3**2.0d0)
                 wnmag=dsqrt(wn1**2.0d0+wn2**2.0d0+wn3**2.0d0)
c                ai=rnmag**(-5.0d0/6.0d0)
                 ai=wnmag**(-5.0d0/6.0d0)
c                dum1=(1.0d0-3.0d0*(rn1/rnmag)**2.0d0)/rn1
                 dum1=(1.0d0-3.0d0*(wn1/wnmag)**2.0d0)/wn1
c                dum2=(1.0d0-3.0d0*(rn2/rnmag)**2.0d0)/rn2
                 dum2=(1.0d0-3.0d0*(wn2/wnmag)**2.0d0)/wn2
c                dum3=(1.0d0-3.0d0*(rn3/rnmag)**2.0d0)/rn3
                 dum3=(1.0d0-3.0d0*(wn3/wnmag)**2.0d0)/wn3
                 utmp=utmp+ai*dum1
     &              *dcos(wn1*(xvel-trel))*dsin(wn2*yvel)*dsin(wn3*zvel)
                 vtmp=vtmp+ai*dum2
     &              *dsin(wn1*(xvel-trel))*dcos(wn2*yvel)*dsin(wn3*zvel)
                 wtmp=wtmp+ai*dum3
     &              *dsin(wn1*(xvel-trel))*dsin(wn2*yvel)*dcos(wn3*zvel)
   50         continue
           else

c isotropic but non-solenoidal...
c             utmp =  dcos(fm*(xvel-trel))*dsin(fm*yvel)*dsin(fm*zvel)
c             vtmp = -dsin(fm*(xvel-trel))*dcos(fm*yvel)*dsin(fm*zvel)
c             wtmp = -dsin(fm*(xvel-trel))*dsin(fm*yvel)*dcos(fm*zvel)
c nonisotropic but solenoidal...
              utmp =  dcos(2.0d0*fm*(xvel-trel))*dsin(fm*yvel)
     &		 *dsin(fm*zvel)
              vtmp = -dsin(2.0d0*fm*(xvel-trel))*dcos(fm*yvel)
     &		 *dsin(fm*zvel)
              wtmp = -dsin(2.0d0*fm*(xvel-trel))*dsin(fm*yvel)
     &		 *dcos(fm*zvel)
           endif

c23456789012345678901234567890123456789012345678901234567890123456789012

	   umax=dmax1(umax,dsqrt(utmp**2.0d0+vtmp**2.0d0+wtmp**2.0d0))
	   energy=energy+(utmp**2.0d0+vtmp**2.0d0+wtmp**2.0d0)/2.0d0

 	   if (dabs(utmp).lt.small) utmp=0.0d0
 	   if (dabs(vtmp).lt.small) vtmp=0.0d0
 	   if (dabs(wtmp).lt.small) wtmp=0.0d0
	   cx_tmp(i,j,k)=utmp
	   cy_tmp(i,j,k)=vtmp
	   cz_tmp(i,j,k)=wtmp

  100	continue

        adiv=0.0d0
        do 150 n3=1,nmodes
           rn3=dfloat(n3)
        do 150 n2=1,nmodes
           rn2=dfloat(n2)
        do 150 n1=1,nmodes
           rn1=dfloat(n1)
           rnmag=dsqrt(rn1**2.0d0+rn2**2.0d0+rn3**2.0d0)
           ratio=rnmag/rn1
           adiv=adiv+((ratio-3.0d0/ratio)*rnmag**(-11.0d0/6.0d0))**2.0d0
  150   continue

	if (igetdt.eq.1) then
           if (ifreq.eq.0) then
              dt=cfl*dx/(sl+umax)
           else
              dt=cfl*dx/(uconst+sl+umax)
           endif
	endif


        sigma=dsqrt(adiv/8.0d0)
c	cx_tmp=dt*(uconst+cx_tmp/sigma)/dx
c	cy_tmp=dt*cy_tmp/(sigma*dy)
c	cz_tmp=dt*cz_tmp/(sigma*dz)
	cx_tmp=dt*(uconst+cx_tmp)/dx
	cy_tmp=dt*cy_tmp/dy
	cz_tmp=dt*cz_tmp/dz

c    Ratio of uprime calculated to the value set in the main program (uprime_rel)
c    should equal 1
        comp_vol=(1.0d0/dfloat(1+iquad))**2.0d0
c       energy=(energy*dx*dy*dz/comp_vol)/sigma**2.0d0
        energy=(energy*dx*dy*dz/comp_vol)
        uprime_rel=dsqrt(2.0d0*energy/3.0d0)

        return
        end

c23456789012345678901234567890123456789012345678901234567890123456789012

 	subroutine cfieldtot(nx,ny,nz,dx,dy,dz,dt,trel,uconst,
     &	   cx_tmp,cy_tmp,cz_tmp,shift,igetdt)
        implicit double precision (a-h,o-z)
     	dimension cx_tmp(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cy_tmp(-2:nx+3,-2:ny+3,-2:nz+3)
     	dimension cz_tmp(-2:nx+3,-2:ny+3,-2:nz+3)
	common /veldat/umax,sl,umod,tke,kr,nmodes,iso,ifreq,cfl
        common /rparam/small,pi
        common /iparam/iquad
	common /kechk/uprime_rel

        ax=1.0d0; ay=1.0d0; az=1.0d0
	energy=0.0d0
	fm=2.0d0*pi
        umax=0.0d0

c	utmpmax=0.0d0; vtmpmax=0.0d0; wtmpmax=0.0d0
c	ushftmax=0.0d0; vshftmax=0.0d0; wshftmax=0.0d0

        do 100 k=-2,nz+3
        do 100 j=-2,ny+3
        do 100 i=-2,nx+3

  	   xvel=-0.5d0+dfloat(i-1)*dx
     &	      -shift*dx*cx_tmp(i,j,k)
  	   yvel=-0.5d0/dfloat(1+iquad)+dfloat(j-1)*dy
     &	      -shift*dy*cy_tmp(i,j,k)
  	   zvel=-0.5d0/dfloat(1+iquad)+dfloat(k-1)*dz
     &	      -shift*dz*cz_tmp(i,j,k)

c          ushftmax=dmax1(ushftmax,dabs(shift*cx_tmp(i,j,k)))
c          vshftmax=dmax1(vshftmax,dabs(shift*cy_tmp(i,j,k)))
c          wshftmax=dmax1(wshftmax,dabs(shift*cz_tmp(i,j,k)))

	   adiv=0.0d0
	   utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0

	   goto 75

	   if (nmodes.gt.1) then
              do 50 n3=1,nmodes
                 rn3=dfloat(n3); wn3=fm*rn3
              do 50 n2=1,nmodes
                 rn2=dfloat(n2); wn2=fm*rn2
              do 50 n1=1,nmodes
                 if ((n1.eq.n2).and.(n2.eq.n3)) goto 50
                 rn1=dfloat(n1); wn1=fm*rn1
                 wnmag=dsqrt(wn1**2.0d0+wn2**2.0d0+wn3**2.0d0)

c		 wnmax=dmax1(wn1,wn2,wn3)
c		 wnfac=wnmax/dsqrt(0.5*(wnmag**2.0d0-wnmax**2.0d0))

c               epsfac = 1/8.0d0
c               a12 = (dabs(rn1-rn2) + epsfac)/
c    &		   dsqrt(rn1**2.0d0 + rn2**2.0d0)
c               a13 = (dabs(rn1-rn3) + epsfac)
c    &		   /dsqrt(rn1**2.0d0 + rn3**2.0d0)
c               a23 = (dabs(rn2-rn3)+ epsfac)
c    &		   /dsqrt(rn2**2.0d0 + rn3**2.0d0)
c               wnfac = a12*a13*a23*dsqrt(8.0d0)

                 ai=(wnmag**(-5.0d0/6.0d0))
c                ai=(wnmag**(-5.0d0/6.0d0))*wnfac
c                ai=(wnmag**(-5.0d0/6.0d0))/(wnfac**8.0d0)
c                ai=(wnmag**(-5.0d0/6.0d0))/(wnfac**15.0d0)

                 dum1=(1.0d0-3.0d0*(wn1/wnmag)**2.0d0)/wn1
                 dum2=(1.0d0-3.0d0*(wn2/wnmag)**2.0d0)/wn2
                 dum3=(1.0d0-3.0d0*(wn3/wnmag)**2.0d0)/wn3
                 utmp=utmp+ai*dum1
     &              *dcos(wn1*(xvel-trel))*dsin(wn2*yvel)*dsin(wn3*zvel)
                 vtmp=vtmp+ai*dum2
     &              *dsin(wn1*(xvel-trel))*dcos(wn2*yvel)*dsin(wn3*zvel)
                 wtmp=wtmp+ai*dum3
     &              *dsin(wn1*(xvel-trel))*dsin(wn2*yvel)*dcos(wn3*zvel)
   50         continue
	   else
c isotropic but non-solenoidal...
c             utmp =  dcos(fm*(xvel-trel))*dsin(fm*yvel)*dsin(fm*zvel)
c             vtmp = -dsin(fm*(xvel-trel))*dcos(fm*yvel)*dsin(fm*zvel)
c             wtmp = -dsin(fm*(xvel-trel))*dsin(fm*yvel)*dcos(fm*zvel)
c nonisotropic but solenoidal...
              utmp =  dcos(2.0d0*fm*(xvel-trel))*dsin(fm*yvel)
     &		 *dsin(fm*zvel)
              vtmp = -dsin(2.0d0*fm*(xvel-trel))*dcos(fm*yvel)
     &		 *dsin(fm*zvel)
              wtmp = -dsin(2.0d0*fm*(xvel-trel))*dsin(fm*yvel)
     &		 *dcos(fm*zvel)
	   endif

c23456789012345678901234567890123456789012345678901234567890123456789012

   75      continue
           do 85 n1=1,nmodes
 	      rnfac = 2.0d0
c	      rnfac = 1.0d0
              rn1=dfloat(2**(n1-1)); wn1=fm*rn1
              utmp = utmp + dcos(rnfac*fm*rn1*(xvel-trel))
     &	         *dsin(fm*rn1*yvel)*dsin(fm*rn1*zvel)
              vtmp = vtmp - dsin(rnfac*fm*rn1*(xvel-trel))
     &	         *dcos(fm*rn1*yvel)*dsin(fm*rn1*zvel)
              wtmp = wtmp - dsin(rnfac*fm*rn1*(xvel-trel))
     &           *dsin(fm*rn1*yvel)*dcos(fm*rn1*zvel)
   85      continue

	   if ((i.ge.1).and.(i.le.nx).and.(j.ge.1).and.(j.le.ny)
     &	      .and.(k.ge.1).and.(k.le.nz)) then
	      energy=energy+(utmp**2.0d0+vtmp**2.0d0+wtmp**2.0d0)/2.0d0
	   endif

c	   utmp = utmp/dfloat(nmodes)**3.0d0
c	   vtmp = vtmp/dfloat(nmodes)**3.0d0
c	   wtmp = wtmp/dfloat(nmodes)**3.0d0

           if (dabs(utmp).lt.small) utmp=0.0d0
           if (dabs(vtmp).lt.small) vtmp=0.0d0
           if (dabs(wtmp).lt.small) wtmp=0.0d0
           cx_tmp(i,j,k)=utmp
           cy_tmp(i,j,k)=vtmp
           cz_tmp(i,j,k)=wtmp

c          utmpmax=dmax1(utmpmax,dabs(utmp))
c          vtmpmax=dmax1(vtmpmax,dabs(vtmp))
c          wtmpmax=dmax1(wtmpmax,dabs(wtmp))

           umax=dmax1(umax,dsqrt(utmp**2.0d0+vtmp**2.0d0+wtmp**2.0d0))

  100   continue

        if (igetdt.eq.1) then
           if (ifreq.eq.0) then
              dt=cfl*dx/(sl+umax)
c             dt=0.25*cfl*dx/(sl+umax)
           else
              dt=cfl*dx/(uconst+sl+umax)
           endif
	endif
	write (6,*) 'dt: ',dt,' dx: ',dx,' uconst: ',uconst,' sl: ',sl,
     &     ' umax: ',umax,' cfl: ',cfl

c       adiv=0.0d0
c       do 150 n3=1,nmodes
c          rn3=dfloat(n3)
c       do 150 n2=1,nmodes
c          rn2=dfloat(n2)
c       do 150 n1=1,nmodes
c          rn1=dfloat(n1)
c          rnmag=dsqrt(rn1**2.0d0+rn2**2.0d0+rn3**2.0d0)
c          ratio=rnmag/rn1
c          adiv=adiv+((ratio-3.0d0/ratio)*rnmag**(-11.0d0/6.0d0))**2.0d0
c 150   continue

c       sigma=dsqrt(adiv/8.0d0)
c	cx_tmp=dt*(uconst+cx_tmp/sigma)/dx
c	cy_tmp=dt*cy_tmp/(sigma*dy)
c	cz_tmp=dt*cz_tmp/(sigma*dz)

	cx_tmp=dt*(uconst+cx_tmp)/dx
	cy_tmp=dt*cy_tmp/dy
	cz_tmp=dt*cz_tmp/dz

c    Ratio of uprime calculated to the value set in the main program (uprime_rel)
c    should equal 1

c	write (6,*) 'uprime_rel = ',uprime_rel,'   in cfieldtot()'
c	write (6,*) energy,comp_vol,sigma,adiv

c       write (6,*) ushftmax, vshftmax, wshftmax
c       write (6,*) utmpmax/sigma, vtmpmax/sigma, 
c    &	  wtmpmax/sigma, sigma
c       write (6,*) dt*utmpmax/dx, dt*vtmpmax/dy, 
c    &    dt*wtmpmax/dz, sigma

        comp_vol=(1.0d0/dfloat(1+iquad))**2.0d0
c       energy=(energy*dx*dy*dz/comp_vol)/sigma**2.0d0
        energy=(energy*dx*dy*dz/comp_vol)
        uprime_rel=dsqrt(2.0d0*energy/3.0d0)

        return
        end


c23456789012345678901234567890123456789012345678901234567890123456789012

        subroutine departure(nx,ny,nz,dx,dy,dz,dt,x,y,z,trel,uconst,
     &     dep,dephalf,nptot,ntc,indx)
        implicit double precision (a-h,o-z)
        dimension x(nx+1),y(ny+1),z(nz+1),indx(nptot,4)
        dimension dep(ntc,3),dephalf(ntc,3)
c       dimension dep(ntc,3),dephalf(ntc,3),deptmp(ntc,3),depold(ntc,3)
c       dimension f1(ntc,3),f2(ntc,3),f3(ntc,3)
        allocatable :: deptmp(:,:),depold(:,:),f1(:,:),f2(:,:),f3(:,:)
        common /veldat/umax,sl,umod,tke,kr,nmodes,iso,ifreq,cfl
        common /rparam/small,pi
        common /iparam/iquad
        common /reindat/irk,igrd,cfl_re

        allocate (deptmp(ntc,3)); deptmp=0.0d0
        allocate (depold(ntc,3)); depold=0.0d0
        allocate (f1(ntc,3),f2(ntc,3),f3(ntc,3))
	f1=0.0d0; f2=0.0d0; f3=0.0d0

        if (iso.ne.1) then ! For mass conservation (non-isotropic)
           ax=1.0d0; ay=-dfloat(kr)/2.0d0; az=ay
        else ! For isotropic fluctuation intensities (mass not conserved)
           ax=1.0d0; ay=-1.0d0; az=ay
        endif
        fm=2.0d0*pi; dtb=-0.5d0*dt ! backwards time integration

        do 50 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)
	   deptmp(nt,1)=x(i); deptmp(nt,2)=y(j); deptmp(nt,3)=z(k)
 50     continue

c23456789012345678901234567890123456789012345678901234567890123456789012

        do 300 n=1,2
	
	if (n.eq.2) then
	   dephalf=(depold-deptmp)/dx
	endif

	depold=deptmp

 	do 100 nt=1,ntc
 	   xdep=deptmp(nt,1); ydep=deptmp(nt,2); zdep=deptmp(nt,3)
	   utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0
           do m=1,nmodes
              fi=fm*dfloat(m)
              ai=dfloat(m)**(-5.0d0/6.0d0)
 	      utmp=utmp+ax*ai*dcos(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dsin(fi*zdep)
 	      vtmp=vtmp+ay*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dcos(fi*ydep)*dsin(fi*zdep)
 	      wtmp=wtmp+az*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dcos(fi*zdep)
   	   end do
 	   f1(nt,1)=utmp; f1(nt,2)=vtmp; f1(nt,3)=wtmp
           deptmp(nt,1)=depold(nt,1)+dtb*f1(nt,1)/2.0d0
           deptmp(nt,2)=depold(nt,2)+dtb*f1(nt,2)/2.0d0
           deptmp(nt,3)=depold(nt,3)+dtb*f1(nt,3)/2.0d0
 100    continue

c23456789012345678901234567890123456789012345678901234567890123456789012

 	do 150 nt=1,ntc
 	   xdep=deptmp(nt,1); ydep=deptmp(nt,2); zdep=deptmp(nt,3)
	   utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0
           do m=1,nmodes
              fi=fm*dfloat(m)
              ai=dfloat(m)**(-5.0d0/6.0d0)
 	      utmp=utmp+ax*ai*dcos(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dsin(fi*zdep)
 	      vtmp=vtmp+ay*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dcos(fi*ydep)*dsin(fi*zdep)
 	      wtmp=wtmp+az*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dcos(fi*zdep)
   	   end do
 	   f2(nt,1)=utmp; f2(nt,2)=vtmp; f2(nt,3)=wtmp
           deptmp(nt,1)=depold(nt,1)+dtb*f2(nt,1)/2.0d0
           deptmp(nt,2)=depold(nt,2)+dtb*f2(nt,2)/2.0d0
           deptmp(nt,3)=depold(nt,3)+dtb*f2(nt,3)/2.0d0
 150    continue

 	do 200 nt=1,ntc
 	   xdep=deptmp(nt,1); ydep=deptmp(nt,2); zdep=deptmp(nt,3)
	   utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0
           do m=1,nmodes
              fi=fm*dfloat(m)
              ai=dfloat(m)**(-5.0d0/6.0d0)
 	      utmp=utmp+ax*ai*dcos(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dsin(fi*zdep)
 	      vtmp=vtmp+ay*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dcos(fi*ydep)*dsin(fi*zdep)
 	      wtmp=wtmp+az*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dcos(fi*zdep)
   	   end do
 	   f3(nt,1)=utmp; f3(nt,2)=vtmp; f3(nt,3)=wtmp
           deptmp(nt,1)=depold(nt,1)+dtb*f3(nt,1)
           deptmp(nt,2)=depold(nt,2)+dtb*f3(nt,2)
           deptmp(nt,3)=depold(nt,3)+dtb*f3(nt,3)
 200    continue

c23456789012345678901234567890123456789012345678901234567890123456789012

 	do 250 nt=1,ntc
 	   xdep=deptmp(nt,1); ydep=deptmp(nt,2); zdep=deptmp(nt,3)
	   utmp=0.0d0; vtmp=0.0d0; wtmp=0.0d0
           do m=1,nmodes
              fi=fm*dfloat(m)
              ai=dfloat(m)**(-5.0d0/6.0d0)
 	      utmp=utmp+ax*ai*dcos(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dsin(fi*zdep)
 	      vtmp=vtmp+ay*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dcos(fi*ydep)*dsin(fi*zdep)
 	      wtmp=wtmp+az*ai*dsin(dfloat(kr)*fi*(xdep-trel))
     &	         *dsin(fi*ydep)*dcos(fi*zdep)
   	   end do
           deptmp(nt,1)=depold(nt,1)
     &        +dtb*(f1(nt,1)+2.0d0*(f2(nt,1)+f3(nt,1))+utmp)/6.0d0
           deptmp(nt,2)=depold(nt,2)
     &        +dtb*(f1(nt,2)+2.0d0*(f2(nt,2)+f3(nt,2))+vtmp)/6.0d0
           deptmp(nt,3)=depold(nt,3)
     &        +dtb*(f1(nt,3)+2.0d0*(f2(nt,3)+f3(nt,3))+wtmp)/6.0d0
 250    continue

 300    continue
 	do nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)
	   dep(nt,1)=(x(i)-deptmp(nt,1))/dx
	   dep(nt,2)=(y(j)-deptmp(nt,2))/dy
	   dep(nt,3)=(z(k)-deptmp(nt,3))/dz
	end do

        deallocate (deptmp,depold,f1,f2,f3)
	return
	end


c23456789012345678901234567890123456789012345678901234567890123456789012

	subroutine transf1(x,xtransf,nxold,nx,ishft)
        implicit double precision (a-h,o-z)
	dimension x(1:nxold+1),xtransf(1:nx+1)
	common /grid/icush,imin,imax,izone

	ishft=ishft+imin-icush-1
 	xtransf(icush+1:icush+1+izone)=x(imin:imax)
 	do 100 i=1,icush
  	   xtransf(i)=x(imin)
     &	      -(x(imax)-x(imin))*dfloat(icush+1-i)/dfloat(izone)
 	   xtransf(icush+1+izone+i)=x(imax)
     &	      +(x(imax)-x(imin))*dfloat(i)/dfloat(izone)
  100	continue
        return
        end


c23456789012345678901234567890123456789012345678901234567890123456789012

	subroutine transf2(g,gtransf,nxold,nx,ny,nz,ibnd)
        implicit double precision (a-h,o-z)
	dimension g(-2:nxold+3,-2:ny+3,-2:nz+3)
	dimension gtransf(-2:nx+3,-2:ny+3,-2:nz+3)
	common /grid/icush,imin,imax,izone
	common /tubedat/twa,twb,gjump,nomask

	gtransf(icush+1:icush+izone+1,-2:ny+3,-2:nz+3)
     &	   =g(imin:imax,-2:ny+3,-2:nz+3) 
	gtransf(-2:icush,-2:ny+3,-2:nz+3)=-twb*dfloat(ibnd); 
	gtransf(icush+izone+2:nx+3,-2:ny+3,-2:nz+3)=twb*dfloat(ibnd)
        return
        end

      	subroutine timestamp
      	integer*4 today(3), now(3)
      	call idate(today)   ! today(1)=day, (2)=month, (3)=year
      	call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
      	write ( *, 1000 )  today(2), today(1), today(3), now
 1000 	format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      	return
      	end


c23456789012345678901234567890123456789012345678901234567890123456789012

	subroutine chksuprema(nx,ny,nz,nt0,nptot,ntc,indx,localmax,g,
     &	   grad,cx,cy,cz)
        implicit double precision (a-h,o-z)
	dimension g(-2:nx+3,-2:ny+3,-2:nz+3),indx(nptot,4)
	dimension grad(-2:nx+3,-2:ny+3,-2:nz+3)
	dimension cx(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cy(-2:nt0+3,-2:ny+3,-2:nz+3)
	dimension cz(-2:nt0+3,-2:ny+3,-2:nz+3)
        common /rparam/small,pi

	localmax=0
	do 100 nt=1,ntc
           i=indx(nt,1); j=indx(nt,2); k=indx(nt,3)
 	   gc = g(i,j,k)
 	   if (((gc-g(i+1,j,k)).gt.small).and.
     &	       ((gc-g(i-1,j,k)).gt.small).and.
     &	       ((gc-g(i,j+1,k)).gt.small).and.
     &	       ((gc-g(i,j-1,k)).gt.small).and.
     &	       ((gc-g(i,j,k+1)).gt.small).and.
     &	       ((gc-g(i,j,k-1)).gt.small)) then

 	      localmax=localmax+1

              ishft=(nt0-nx)/2
              if ((i+ishft).gt.0) then
                 iv=mod((i+ishft-1),nt0)+1
              else
                 iv=nt0+mod((i+ishft),nt0)
              endif

cc            write (6,2000) iv,j,k,indx(nt,4),g(i-1,j,k),g(i,j,k),
c             write (6,2000) i,j,k,nx,g(i-1,j,k),g(i,j,k),
c    &           g(i+1,j,k),g(i,j-1,k),g(i,j+1,k),g(i,j,k-1),g(i,j,k+1)
c             write (6,2010) grad(i-1,j,k),grad(i,j,k),
c    &           grad(i+1,j,k),grad(i,j-1,k),grad(i,j+1,k),
c    &		 grad(i,j,k-1),grad(i,j,k+1)
c             write (6,2010) cx(iv-1,j,k),cx(iv,j,k),
c    &           cx(iv+1,j,k),cx(iv,j-1,k),cx(iv,j+1,k),
c    &		 cx(iv,j,k-1),cx(iv,j,k+1)
c             write (6,2010) cy(iv-1,j,k),cy(iv,j,k),
c    &           cy(iv+1,j,k),cy(iv,j-1,k),cy(iv,j+1,k),
c    &		 cy(iv,j,k-1),cy(iv,j,k+1)
c             write (6,2010) cz(iv-1,j,k),cz(iv,j,k),
c    &           cz(iv+1,j,k),cz(iv,j-1,k),cz(iv,j+1,k),
c    &		 cz(iv,j,k-1),cz(iv,j,k+1)
 	   endif

  100	continue

 2000   format (4(1x,i3),7(1x,e18.12))
 2010   format (7(1x,e18.12))
	return
	end

c **********************************************************************
c (C) Copyright 2009-2019 by Ralph C Aldredge, PhD, PE                 
c  All Rights Reserved.                           		       
c Licensed under the GNU Affero General Public License v3.0
c **********************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
