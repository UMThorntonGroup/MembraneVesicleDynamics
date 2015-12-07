!	Membrane_Vesicle.F90: A Membrane Vesicle Dynamics Simulation Code
!
!	Copyright (C) 2011 Chloe M. Funkhouser, Francisco J. Solis, Katsuyo Thornton
!	Materials Science and Engineering, University of Michigan
!
!
!	This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
program Membrane_Vesicle
implicit none
Include 'mpif.h'
!declare variables
!variable endings indicate which grid they apply to
!	'E' = Yin grid
!	'A' = Yang grid

!USER INPUT SECTION
!Azimuthal number of grid elements Npsi MUST be divisible by 6
!Polar number of grid elements Ntheta MUST be divisible by 2 (Ntheta=Npsi/3) 
!Note that dimension of a, b, ip_psi, ip_th, Corn_psi, Corn_th will change with Npsi and Ntheta 

!seedE and seed A : Random seeds for composition initialization

!endt : total loop iterations

!plotstep : initial frequency of data output (data output every plotstep iterations).
!If variable reduce is set to 1, plotstep will decrease with increasing loop iterations 
!(i.e., less frequent output a later times)

!cntt : initial loop iteration count - allows user to start from 0 (cntt=0) or pick up
!a previous run by setting cntt to the output plot number to begin with

!M : Constant in Cahn-Hilliard dynamics

!gam : Constant in surface morphology dynamics (lower-case gamma)

!sigma_init : Initial value of surface tension sigma

!areafrac : Area fraction of the two phases

!lambda : Bending rigidity (upper-case Lambda)

!C1, C2 : Spontaneous curvatures of the two phases

!rinit : Initial radius of vesicle

!noise : Amplitude of random noise to be used for compositional initialization

!lambda_line : Line tension (lower-case lambda)


integer(kind=4), parameter :: Npsi=270, Ntheta=90
integer(kind=4) :: seedE=4, seedA=5
integer(kind=4) :: endt=500000000, plotstep=1, cntt=0, reduce=1
real(kind=8) :: M=1.d0, gam=0.05d0, sigma_init=5.d0, areafrac=0.4d0, lambda=1.d0
real(kind=8) :: C1=1.d0, C2=-1.d0, rinit=3.d0, noise=0.1d0, lambda_line=1.5d0

!END OF USER INPUT SECTION

real(kind=8) :: dpsi, dtheta, dt, dx, W, eps, P, u, v, s, cutoff, off=2.d0
real(kind=4) :: r4_uniform_01
integer(kind=4) :: Center_neg, Center_pos, Min_psi, Min_th, Max_psi, Max_th, cntt_orig, realt
real(kind=8) :: x, y, z, xp, yp, zp, psip, thp, th_prime, psi_prime, dist, Rint
real(kind=8), dimension(0:4000) :: a, b
integer(kind=4), dimension(0:4000) :: ip_psi, ip_th, Corn_psi, Corn_th
real(kind=8) :: sigma, dsig, area0, area1, area0diff, weight_elem
real(kind=8), dimension(1:7) :: sum_arr, sums
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)) :: phiE, rE, phiA, rA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)) :: d1rEall, d2rEall, sqgEall
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)) :: d1rAall, d2rAall, sqgAall
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)) :: cos_sq_th_all
real(kind=8), dimension((-Npsi/2-5):(Npsi/2+5),(-Ntheta/2-4):(Ntheta/2+4)):: phiEIN, phiAIN, rEIN, rAIN

integer(kind=4) :: i, j, j_all, t, count, endcount, plott, plotnum
character :: step*6, namePE*35, namePA*35, namerE*12, namerA*12
integer(kind=4) :: myrank, nprocs, err, psi_start, psi_end
integer(kind=4) :: disp, size_total, size_noghost, size_noghost_psi
integer(kind=4), dimension(0:3) :: displ, size_in, size_out
integer(kind=4) :: status(MPI_STATUS_SIZE)


integer(kind=4), dimension(0:26500) :: area_psi, area_th
real(kind=8), dimension(0:26500) :: weight
integer(kind=4), dimension(0:23500) :: area_psi_mod, area_th_mod
real(kind=8), dimension(0:23500) :: weight_mod
real(kind=8) :: dist_near, dist_far, psi_a_dist, psi_c_dist, th_b_dist, th_d_dist, weighting
real(kind=8) :: psi_alpha, psi_gamma, th_beta, th_delta, area_cutoff, area_cutoff_sq, cos_area_cutoff
integer(kind=4) :: alpha, beta, gamma, delta, count_a, count_a_mod, endcount_a, endcount_a_mod, index
real(kind=8) :: x_a, y_a, z_a, x_p, y_p, z_p, xC, yC, zC


real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: cos_th, cos_sq_th, sin_th
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: pE, p_primeE, phi_adE, dphitE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: pA, p_primeA, phi_adA, dphitA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: d1phiE, d2phiE, d11phiE, d22phiE, d12phiE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: d1phiA, d2phiA, d11phiA, d22phiA, d12phiA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: muE, d1muE, d2muE, d11muE, d22muE, d12muE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: muA, d1muA, d2muA, d11muA, d22muA, d12muA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: d1rE, d2rE, d11rE, d22rE, d12rE, drtE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: d1rA, d2rA, d11rA, d22rA, d12rA, drtA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: rsqE, d1rsqE, d2rsqE, sqgE, TmE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: rsqA, d1rsqA, d2rsqA, sqgA, TmA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: gE, g11E, g12E, g22E, Gam1E, Gam2E
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: gA, g11A, g12A, g22A, Gam1A, Gam2A
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: KE, K11E, K12E, K22E, KabE, SE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: KA, K11A, K12A, K22A, KabA, SA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: Lap_phiE, Lap_muE, Lap_QE, grad_phiE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: Lap_phiA, Lap_muA, Lap_QA, grad_phiA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: QE, d1QE, d2QE, d11QE, d22QE, d12QE
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)):: QA, d1QA, d2QA, d11QA, d22QA, d12QA
	

call MPI_INIT(err)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, err)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, err)

Rint=real(Ntheta,8)/2.d0 + 2.d0

plott=plotstep*1000

!Calculate mesh step size
dpsi=(1.5d0*atan(1.d0)*4.d0)/real(Npsi,8)
dtheta=(0.5d0*atan(1.d0)*4.d0)/real(Ntheta,8)
dx=rinit*dpsi
dt=0.005d0*dx*dx*dx*dx


W=12.d0*lambda_line/(2.0d0*rinit*dpsi)
eps=lambda_line*6.d0*SQRT(2.d0)/SQRT(W)

P=2.d0*sigma_init/rinit

Center_neg=-Npsi/3
Center_pos=Npsi/3

Min_psi=-Npsi/2-5
Max_psi=Npsi/2+5

Min_th=-Ntheta/2-4
Max_th=Ntheta/2+4

cntt_orig=cntt

do i=Min_th,Max_th,1
	sin_th(i,:)=sin(real(i,8)*dtheta)
	cos_th(i,:)=cos(real(i,8)*dtheta)
	cos_sq_th(i,:)=cos_th(i,:)*cos_th(i,:)
	cos_sq_th_all(i,:)=cos_sq_th(i,1)
end do


!The following is all for initialization setup and output, and also setup for the interpolation,
!all of which need only be done on the master processor
if(myrank==0) then

	!***********************************************************************
	!!INITIALIZATIONS
	if(cntt==0) then
		!initialize concentration matrix with random noise
		do i=Min_th,Max_th,1
			do j=Min_psi,Max_psi,1
	
				phiE(i,j)=(areafrac-noise) + (2.d0*noise)*r4_uniform_01(seedE)
				phiA(i,j)=(areafrac-noise) + (2.d0*noise)*r4_uniform_01(seedA)

			end do
		end do

		
		rE(:,:)=rinit
		rA(:,:)=rinit
		
!Offcenter sphere	
!		do i=Min_th,Max_th,1
!			rE(i,:)=-off*sin(real(i,8)*dtheta) + SQRT(rinit*rinit - off*off*cos(real(i,8)*dtheta)*cos(real(i,8)*dtheta))
!			
!		end do
!	
!	
!		do i=Min_th,Max_th,1
!			do j=Min_psi,Max_psi,1
!	
!				th_prime=asin(cos(real(i,8)*dtheta)*sin(real(j,8)*dpsi))
!			
!				rA(i,j)=-off*sin(th_prime) + SQRT(rinit*rinit - off*off*cos(th_prime)*cos(th_prime))
!		
!			end do
!		end do
		
	else	
	
		plotnum=(100000) + cntt
		
		!If plot labels go over 99999, start over at 10000
		if (plotnum.GT.999999) then
			plotnum=plotnum-900000
		end if
	
		!change numerical file names to strings and label as .txt files
		write(step, "(I6)") plotnum
		
		namePE=''//step //'pE.txt'
		namePA=''//step //'pA.txt'
		
		namerE=''//step //'rE.txt'
		namerA=''//step //'rA.txt'
		

		open(unit=10, file=namePE, status='old', blank='null')
			read(10,992) ((phiE(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
			992 format(99(E16.9, X))
		close(unit=10, status='keep')
		
		open(unit=10, file=namePA, status='old', blank='null')
			read(10,993) ((phiA(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
			993 format(99(E16.9, X))
		close(unit=10, status='keep')
		
		
		open(unit=10, file=namerE, status='old', blank='null')
			read(10,994) ((rE(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
			994 format(99(E16.9, X))
		close(unit=10, status='keep')
		
		open(unit=10, file=namerA, status='old', blank='null')
			read(10,995) ((rA(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
			995 format(99(E16.9, X))
		close(unit=10, status='keep')

	
	end if
	
	call expdata(cntt, phiE, phiA, rE, rA, Npsi, Ntheta, Min_psi, Max_psi, Min_th, Max_th)
	
	
!*****************************************
!SETUP FOR INTERPOLATION BETWEEN YIN AND YANG GRIDS
	!Find coordinates of corresponding points
	
	count=0
	
	!Negative psi corners
	do i=Min_th,Max_th,1
		do j=Min_psi,Center_neg,1
			
			!Identify points in corner that need interpolated values (outside a specified radius)
			dist=SQRT(real((Center_neg-j)*(Center_neg-j),8) + real(i*i,8))
	
			if (dist.GT.Rint) then
				
				!Insert indices of point needing interpolation into arrays
				Corn_psi(count)=j
				Corn_th(count)=i
				
				!Cartesian coordinates of point (i,j)
				x=rinit*cos(real(i,8)*dtheta)*cos(real(j,8)*dpsi)
				y=rinit*cos(real(i,8)*dtheta)*sin(real(j,8)*dpsi)			
				z=-rinit*sin(real(i,8)*dtheta)
	
				!Cartesian coordinates of point (i,j) in other grid
				xp=-x
				yp=-z
				zp=-y
				
				!psip and thp are the polar coordinates of point (i,j) in the other grid
				
				!When i=Center_neg (90 degrees), xp=0 and psip cannot be calculated normally (div. by 0)
				!Set ip_psi(count) for these points manually
				if(j==Center_neg .AND. i.GT.0) then
					thp=asin(-zp/rinit)
					
					!Find nearest points, in polar coordinates, to point (i,j) in the other grid
					!Insert the indices of this interpolation point into arrays
					ip_psi(count)=Center_pos
					ip_th(count)=nint(thp/dtheta)
					
					!Calculate offset distances of interpolation point, insert into arrays
					a(count)=0.d0
					b(count)=thp-real(ip_th(count),8)*dtheta
				
				else if(j==Center_neg .AND. i.LT.0) then
					thp=asin(-zp/rinit)
					
					!Find nearest points, in polar coordinates, to point (i,j) in the other grid
					!Insert the indices of this interpolation point into arrays
					ip_psi(count)=Center_neg
					ip_th(count)=nint(thp/dtheta)
	
					!Calculate offset distances of interpolation point, insert into arrays
					a(count)=0.d0
					b(count)=thp-real(ip_th(count),8)*dtheta
					
				else
					
					thp=asin(-zp/rinit)
					psip=acos(min(xp/(rinit*cos(thp)),1.0d0))*SIGN(1.d0,yp)
			
					!Find nearest points, in polar coordinates, to point (i,j) in the other grid
					!Insert the indices of this interpolation point into arrays
					ip_psi(count)=nint(psip/dpsi)
					ip_th(count)=nint(thp/dtheta)
			
					!Calculate offset distances of interpolation point, insert into arrays
					a(count)=psip-real(ip_psi(count),8)*dpsi
					b(count)=thp-real(ip_th(count),8)*dtheta
	
				end if
				
				!Increment array counter for next point needing interpolation
				count=count + 1
		
			end if
	
		end do
	end do
	
	!*****************************
	!Positive psi corners
	
	do i=Min_th,Max_th,1
		do j=Center_pos,Max_psi,1
			
			!Identify points in corner that need interpolated values
			dist=SQRT(real((Center_pos-j)*(Center_pos-j),8) + real(i*i,8))
	
			if (dist.GT.Rint) then
				
				!Insert indices of point needing interpolation into arrays
				Corn_psi(count)=j
				Corn_th(count)=i
				
				!Cartesian coordinates of point (i,j)
				x=rinit*cos(real(i,8)*dtheta)*cos(real(j,8)*dpsi)
				y=rinit*cos(real(i,8)*dtheta)*sin(real(j,8)*dpsi)			
				z=-rinit*sin(real(i,8)*dtheta)
	
				!Cartesian coordinates of point (i,j) in other grid
				xp=-x
				yp=-z
				zp=-y
				
				!psip and thp are the polar coordinates of point (i,j) in the other grid
				
				!When i=Center_pos (90 degrees), xp=0 and psip cannot be calculated normally (div. by 0)
				!Set ip_psi(count) for these points manually
				if(j==Center_pos .AND. i.GT.0) then
					thp=asin(-zp/rinit)
	
					!Find nearest points, in polar coordinates, to point (i,j) in the other grid
					!Insert the indices of this interpolation point into arrays
					ip_psi(count)=Center_pos
					ip_th(count)=nint(thp/dtheta)
					
					!Calculate offset distances of interpolation point, insert into arrays
					a(count)=0.d0
					b(count)=thp-real(ip_th(count),8)*dtheta
				
				else if (j==Center_pos .AND. i.LT.0) then
					thp=asin(-zp/rinit)
					
					!Find nearest points, in polar coordinates, to point (i,j) in the other grid
					!Insert the indices of this interpolation point into arrays
					ip_psi(count)=Center_neg
					ip_th(count)=nint(thp/dtheta)
					
					!Calculate offset distances of interpolation point, insert into arrays
					a(count)=0.d0
					b(count)=thp-real(ip_th(count),8)*dtheta
					
				else
					thp=asin(-zp/rinit)
					psip=acos(min(xp/(rinit*cos(thp)),1.0d0))*SIGN(1.d0,yp)
	
					!Find nearest points, in polar coordinates, to point (i,j) in the other grid
					!Insert the indices of this interpolation point into arrays
					ip_psi(count)=nint(psip/dpsi)
					ip_th(count)=nint(thp/dtheta)
					
					!Calculate offset distances of interpolation point, insert into arrays
					a(count)=psip-real(ip_psi(count),8)*dpsi
					b(count)=thp-real(ip_th(count),8)*dtheta
	
				end if
				!Increment array counter for next point needing interpolation
				count=count + 1
				
			end if
	
		end do
	end do
	
	!Boundary conditions for non-interpolated parts of grid (bottom edge)
	do i=Min_th,Min_th+1,1
		do j=Center_neg+1,Center_pos-1,1
		
		Corn_psi(count)=j
		Corn_th(count)=i
		
		x=rinit*cos(real(i,8)*dtheta)*cos(real(j,8)*dpsi)
		y=rinit*cos(real(i,8)*dtheta)*sin(real(j,8)*dpsi)			
		z=-rinit*sin(real(i,8)*dtheta)
	
		xp=-x
		yp=-z
		zp=-y
				
		thp=asin(-zp/rinit)
		psip=acos(min(xp/(rinit*cos(thp)),1.0d0))*SIGN(1.d0,yp)
				
		ip_psi(count)=nint(psip/dpsi)
		ip_th(count)=nint(thp/dtheta)
				
		a(count)=psip-real(ip_psi(count),8)*dpsi
		b(count)=thp-real(ip_th(count),8)*dtheta
		
		count=count+1
		end do
	end do
		
	!Boundary conditions for non-interpolated parts of grid (top edge)	
	do i=Max_th-1,Max_th,1
		do j=Center_neg+1,Center_pos-1,1
	
		Corn_psi(count)=j
		Corn_th(count)=i
		
		x=rinit*cos(real(i,8)*dtheta)*cos(real(j,8)*dpsi)
		y=rinit*cos(real(i,8)*dtheta)*sin(real(j,8)*dpsi)			
		z=-rinit*sin(real(i,8)*dtheta)
	
		xp=-x
		yp=-z
		zp=-y
			
		thp=asin(-zp/rinit)
		psip=acos(min(xp/(rinit*cos(thp)),1.0d0))*SIGN(1.d0,yp)
				
		ip_psi(count)=nint(psip/dpsi)
		ip_th(count)=nint(thp/dtheta)
				
		a(count)=psip-real(ip_psi(count),8)*dpsi
		b(count)=thp-real(ip_th(count),8)*dtheta
		
		count=count+1
		end do
	end do
	
	endcount=count-1
		
	!END SETUP FOR INTERPOLATION BETWEEN YIN AND YANG GRIDS
	!*****************************************
	
	
	!*****************************************
	!SETUP FOR WEIGHTING IN AREA CONSERVATION SCHEME
	
	!Note that the area calc setup has psi as the first index and theta as the second
	!This is accounted for at the end when the indices are switched to match with the data arrays
	
	count_a=0
	
	area_cutoff=45.0d0
	area_cutoff_sq=area_cutoff*area_cutoff
	cos_area_cutoff=cos(area_cutoff*dpsi)
	
	!Quadrant 1 (negative theta and psi) for area calc
	
	!Cartesian coordinates of Center_neg point on unit sphere (psi,theta)=(-90,0) is degrees
	xC=cos(0.d0)*cos(real(Center_neg,8)*dpsi)
	yC=cos(0.d0)*sin(real(Center_neg,8)*dpsi)			
	zC=0.d0
	
	do i=Min_psi,Center_neg-1,1
		do j=Min_th,-1,1
		
			area_psi(count_a)=i
			area_th(count_a)=j
			
			!Cartesian coordinates of NEAR cell corner associated with point (i,j) on unit sphere
			x=cos((real(j,8)+0.5d0)*dtheta)*cos((real(i,8)+0.5d0)*dpsi)
			y=cos((real(j,8)+0.5d0)*dtheta)*sin((real(i,8)+0.5d0)*dpsi)			
			z=-sin((real(j,8)+0.5d0)*dtheta)
	
			!Check if this point should be used for area calculation (nonzero weight)
			!Formula from 'Spherical distance' on Mathworld
			dist_near=acos(x*xC + y*yC + z*zC)/dpsi
	
			if(dist_near .LT. area_cutoff) then
			
				!Cartesian coordinates of FAR cell corner associated with point (i,j) on unit sphere
				x=cos((real(j,8)-0.5d0)*dtheta)*cos((real(i,8)-0.5d0)*dpsi)
				y=cos((real(j,8)-0.5d0)*dtheta)*sin((real(i,8)-0.5d0)*dpsi)			
				z=-sin((real(j,8)-0.5d0)*dtheta)
			
				!Check if any part of this point's cell overlaps with other grid (distance from Center_neg to far corner of cell)
				dist_far=acos(x*xC + y*yC + z*zC)/dpsi
	
				if(dist_far .LE. area_cutoff) then
				
					weight(count_a) = 1.d0
				
				else
				
					!Find whether the cell has configuration A,B,C,or D
					!See if alpha point exists
					psi_alpha = -180.d0 - asin(-cos_area_cutoff/cos((real(j,8)-0.5d0)*dtheta))/dpsi
	
					if(psi_alpha .GE. (real(i,8) - 0.5d0) .AND. psi_alpha .LE. (real(i,8) + 0.5d0) .AND. j.NE.-45) then
						alpha=1
						!Location of (psi_alpha,theta_j-0.5)
						x_a=cos((real(j,8)-0.5d0)*dtheta)*cos(psi_alpha*dpsi)
						y_a=cos((real(j,8)-0.5d0)*dtheta)*sin(psi_alpha*dpsi)
						z_a=-sin((real(j,8)-0.5d0)*dtheta)
						
						!Location of (psi_i+0.5,theta_j-0.5)
						x_p=cos((real(j,8)-0.5d0)*dtheta)*cos((real(i,8)+0.5d0)*dpsi)
						y_p=cos((real(j,8)-0.5d0)*dtheta)*sin((real(i,8)+0.5d0)*dpsi)
						z_p=-sin((real(j,8)-0.5d0)*dtheta)
						
						psi_a_dist=acos(x_a*x_p + y_a*y_p + z_a*z_p)
	
					else
						alpha=0
					end if
					
					!See if gamma point exists
					psi_gamma = -180.d0 - asin(-cos_area_cutoff/cos((real(j,8)+0.5d0)*dtheta))/dpsi
	
					if(psi_gamma .GE. (real(i,8) - 0.5d0) .AND. psi_gamma .LE. (real(i,8) + 0.5d0)) then
						gamma=1
						!Location of (psi_gamma,theta_j+0.5)
						x_a=cos((real(j,8)+0.5d0)*dtheta)*cos(psi_gamma*dpsi)
						y_a=cos((real(j,8)+0.5d0)*dtheta)*sin(psi_gamma*dpsi)
						z_a=-sin((real(j,8)+0.5d0)*dtheta)
						
						!Location of (psi_i+0.5,theta_j+0.5)
						x_p=cos((real(j,8)+0.5d0)*dtheta)*cos((real(i,8)+0.5d0)*dpsi)
						y_p=cos((real(j,8)+0.5d0)*dtheta)*sin((real(i,8)+0.5d0)*dpsi)
						z_p=-sin((real(j,8)+0.5d0)*dtheta)					
						
						psi_c_dist=acos(x_a*x_p + y_a*y_p + z_a*z_p)
	
					else
						gamma=0
					end if
					
					!See if beta point exists
					th_beta = -acos(-cos_area_cutoff/sin((real(i,8)-0.5d0)*dpsi))/dtheta
	
					if(th_beta .GE. (real(j,8) - 0.5d0) .AND. th_beta .LE. (real(j,8) + 0.5d0) .AND. i.NE.-135) then
						beta=1
						!Location of (psi_i-0.5,theta_beta)
						x_a=cos(th_beta*dtheta)*cos((real(i,8)-0.5d0)*dpsi)
						y_a=cos(th_beta*dtheta)*sin((real(i,8)-0.5d0)*dpsi)
						z_a=-sin(th_beta*dtheta)
						
						!Location of (psi_i-0.5,theta_j+0.5)
						x_p=cos((real(j,8)+0.5d0)*dtheta)*cos((real(i,8)-0.5d0)*dpsi)
						y_p=cos((real(j,8)+0.5d0)*dtheta)*sin((real(i,8)-0.5d0)*dpsi)
						z_p=-sin((real(j,8)+0.5d0)*dtheta)
	
						th_b_dist=acos(x_a*x_p + y_a*y_p + z_a*z_p)
	
					else
						beta=0
					end if
					
					!See if delta point exists
					th_delta = -acos(-cos_area_cutoff/sin((real(i,8)+0.5d0)*dpsi))/dtheta
	
					if(th_delta .GE. (real(j,8) - 0.5d0) .AND. th_delta .LE. (real(j,8) + 0.5d0)) then
						delta=1
						!Location of (psi_i+0.5,theta_delta)
						x_a=cos(th_delta*dtheta)*cos((real(i,8)+0.5d0)*dpsi)
						y_a=cos(th_delta*dtheta)*sin((real(i,8)+0.5d0)*dpsi)
						z_a=-sin(th_delta*dtheta)
						
						!Location of (psi_i+0.5,theta_j+0.5)
						x_p=cos((real(j,8)+0.5d0)*dtheta)*cos((real(i,8)+0.5d0)*dpsi)
						y_p=cos((real(j,8)+0.5d0)*dtheta)*sin((real(i,8)+0.5d0)*dpsi)
						z_p=-sin((real(j,8)+0.5d0)*dtheta)
	
						th_d_dist=acos(x_a*x_p + y_a*y_p + z_a*z_p)					
	
					else
						delta=0
					end if
	
						
					!Now use the alpha,beta,gamma,delta values to determine A,B,C, or D
					
					if(alpha==1 .AND. beta==1) then !this is case A
						weight(count_a) = (psi_a_dist*dtheta + 0.5d0*(dtheta + th_b_dist) &
						
							*(dpsi - psi_a_dist))/(dpsi*dtheta)
				
					else if(alpha==1 .AND. gamma==1) then !this is case B
						weight(count_a) = 0.5d0*dtheta*(psi_a_dist + psi_c_dist)/(dpsi*dtheta)
				
					else if(beta==1 .AND. delta==1) then !this is case C
						weight(count_a) = 0.5d0*dpsi*(th_b_dist + th_d_dist)/(dpsi*dtheta)
				
					else if(gamma==1 .AND. delta==1) then !this is case D
						weight(count_a) = 0.5d0*psi_c_dist*th_d_dist/(dpsi*dtheta)
					end if
	
				end if
				
			else
				weight(count_a)=0.d0
	
			end if
			
			count_a=count_a+1
		end do	
	end do
	
	endcount_a=count_a-1
	
	do index=0,endcount_a
		i=area_psi(index)
		j=area_th(index)
		
		!Quadrant 2
		area_psi(count_a)=-i
		area_th(count_a)=j
		weight(count_a)=weight(index)
		
		count_a=count_a+1
		
		!Quadrant 4
		area_psi(count_a)=i
		area_th(count_a)=-j
		weight(count_a)=weight(index)
		
		count_a=count_a+1
		
		!Quadrant 3
		area_psi(count_a)=-i
		area_th(count_a)=-j
		weight(count_a)=weight(index)
		
		count_a=count_a+1
	end do
	
	do i=Min_psi,Center_neg,1 !for j=0 line, between Quadrants 1 and 4
		area_psi(count_a)=i
		area_th(count_a)=0
		
		if(real(i,8)+0.5d0 .LT. (real(Center_neg,8) - area_cutoff)) then
			weight(count_a) = 0.d0
		else if(real(i,8)-0.5d0 .GE. (real(Center_neg,8) - area_cutoff)) then
			weight(count_a) = 1.d0
		else
			weight(count_a) = 0.5d0
		end if
	
		count_a=count_a+1
	end do
	
	
	do j=Min_th,-1,1 !for i=Center_neg line in Quadrant 1
		area_psi(count_a)=Center_neg
		area_th(count_a)=j
		
		if(real(j,8)+0.5d0 .LT. -area_cutoff) then
			weight(count_a) = 0.d0
		else if (real(j,8)-0.5d0 .GE. -area_cutoff) then
			weight(count_a) = 1.d0
		else
			weight(count_a) = 0.5d0
		end if
	
		count_a=count_a+1
	end do
	
	do j=1,Max_th,1 !for i=Center_neg line in Quadrant 4
		area_psi(count_a)=Center_neg
		area_th(count_a)=j
		
		if(real(j,8)-0.5d0 .GT. area_cutoff) then
			weight(count_a) = 0.d0
		else if (real(j,8)+0.5d0 .LE. area_cutoff) then
			weight(count_a) = 1.d0
		else
			weight(count_a) = 0.5d0
		end if
	
		count_a=count_a+1
	end do
	
	
	do i=Center_pos,Max_psi,1 !for j=0 line between Quadrants 2 and 3
		area_psi(count_a)=i
		area_th(count_a)=0
		
		if(real(i,8)-0.5d0 .GT. (real(Center_pos,8) + area_cutoff)) then
			weight(count_a) = 0.d0
		else if(real(i,8)+0.5d0 .LE. (real(Center_pos,8) + area_cutoff)) then
			weight(count_a) = 1.d0
		else
			weight(count_a) = 0.5d0
		end if
	
		count_a=count_a+1
	end do
	
	do j=Min_th,-1,1 !for i=Center_pos line in Quadrant 2
		area_psi(count_a)=Center_pos
		area_th(count_a)=j
		
		if(real(j,8)+0.5d0 .LT. -area_cutoff) then
			weight(count_a) = 0.d0
		else if(real(j,8)-0.5d0 .GE. -area_cutoff) then
			weight(count_a) = 1.d0
		else
			weight(count_a) = 0.5d0
		end if
	
		count_a=count_a+1
	end do
	
	do j=1,Max_th,1 !for i=Center_pos line in Quadrant 3
		area_psi(count_a)=Center_pos
		area_th(count_a)=j
		
		if(real(j,8)-0.5d0 .GT. area_cutoff) then
			weight(count_a) = 0.d0
		else if(real(j,8)+0.5d0 .LE. area_cutoff) then
			weight(count_a) = 1.d0
		else
			weight(count_a) = 0.5d0
		end if
	
		count_a=count_a+1
	end do
			
	!***************
	
	do i=Center_neg+1,Center_pos-1,1 !Entire center part of grid, between psi=Center_neg and Center_pos
		do j=Min_th+5,Max_th-5,1
			area_psi(count_a)=i
			area_th(count_a)=j
			weight(count_a)=1.d0
			count_a=count_a+1
		end do
	end do
	
	do i=Center_neg+1,Center_pos-1,1
		j=Min_th+4
		area_psi(count_a)=i
		area_th(count_a)=j
		weight(count_a)=0.5d0
		count_a=count_a+1
	end do
	
	do i=Center_neg+1,Center_pos-1,1 
		j=Max_th-4
		area_psi(count_a)=i
		area_th(count_a)=j
		weight(count_a)=0.5d0
		count_a=count_a+1
	end do
	
	endcount_a=count_a-1
	
	!Remove all points with weight of 0 from weighting arrays
	count_a_mod=0
	do count_a=0,endcount_a,1
		if (weight(count_a) .NE. 0.d0) then
			area_psi_mod(count_a_mod)=area_psi(count_a)
			area_th_mod(count_a_mod)=area_th(count_a)
			weight_mod(count_a_mod)=weight(count_a)
			count_a_mod=count_a_mod + 1
		end if
	end do
	
	endcount_a_mod=count_a_mod-1
	
	!END SETUP FOR WEIGHTING IN AREA CONSERVATION SCHEME
	!*****************************************
	
	!Actually perform interpolation
	
	do count=0,endcount
		
		phiE(Corn_th(count),Corn_psi(count))=phiA(ip_th(count),ip_psi(count)) &
		
				+ a(count)*(phiA(ip_th(count),ip_psi(count)+1)-phiA(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
	
				+ b(count)*(phiA(ip_th(count)+1,ip_psi(count))-phiA(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
		
				+ 0.5d0*a(count)*a(count)*(phiA(ip_th(count),ip_psi(count)+1) &
				
						- 2.d0*phiA(ip_th(count),ip_psi(count)) &
				
						+ phiA(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
					
				+ 0.5d0*b(count)*b(count)*(phiA(ip_th(count)+1,ip_psi(count)) &
				
						- 2.d0*phiA(ip_th(count),ip_psi(count)) &
				
						+ phiA(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
					
				+ a(count)*b(count)*(phiA(ip_th(count)+1,ip_psi(count)+1) - phiA(ip_th(count)+1,ip_psi(count)-1) &
	
					- phiA(ip_th(count)-1,ip_psi(count)+1) + phiA(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
					
	
		phiA(Corn_th(count),Corn_psi(count))=phiE(ip_th(count),ip_psi(count)) &
		
				+ a(count)*(phiE(ip_th(count),ip_psi(count)+1)-phiE(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
	
				+ b(count)*(phiE(ip_th(count)+1,ip_psi(count))-phiE(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
		
				+ 0.5d0*a(count)*a(count)*(phiE(ip_th(count),ip_psi(count)+1) &
				
						- 2.d0*phiE(ip_th(count),ip_psi(count)) &
				
						+ phiE(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
					
				+ 0.5d0*b(count)*b(count)*(phiE(ip_th(count)+1,ip_psi(count)) &
				
						- 2.d0*phiE(ip_th(count),ip_psi(count)) &
				
						+ phiE(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
					
				+ a(count)*b(count)*(phiE(ip_th(count)+1,ip_psi(count)+1) - phiE(ip_th(count)+1,ip_psi(count)-1) &
	
					- phiE(ip_th(count)-1,ip_psi(count)+1) + phiE(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
					
					
		rE(Corn_th(count),Corn_psi(count))=rA(ip_th(count),ip_psi(count)) &
		
				+ a(count)*(rA(ip_th(count),ip_psi(count)+1)-rA(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
	
				+ b(count)*(rA(ip_th(count)+1,ip_psi(count))-rA(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
		
				+ 0.5d0*a(count)*a(count)*(rA(ip_th(count),ip_psi(count)+1) &
				
						- 2.d0*rA(ip_th(count),ip_psi(count)) &
				
						+ rA(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
					
				+ 0.5d0*b(count)*b(count)*(rA(ip_th(count)+1,ip_psi(count)) &
				
						- 2.d0*rA(ip_th(count),ip_psi(count)) &
				
						+ rA(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
					
				+ a(count)*b(count)*(rA(ip_th(count)+1,ip_psi(count)+1) - rA(ip_th(count)+1,ip_psi(count)-1) &
	
					- rA(ip_th(count)-1,ip_psi(count)+1) + rA(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
					
	
		rA(Corn_th(count),Corn_psi(count))=rE(ip_th(count),ip_psi(count)) &
		
				+ a(count)*(rE(ip_th(count),ip_psi(count)+1)-rE(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
	
				+ b(count)*(rE(ip_th(count)+1,ip_psi(count))-rE(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
		
				+ 0.5d0*a(count)*a(count)*(rE(ip_th(count),ip_psi(count)+1) &
				
						- 2.d0*rE(ip_th(count),ip_psi(count)) &
				
						+ rE(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
					
				+ 0.5d0*b(count)*b(count)*(rE(ip_th(count)+1,ip_psi(count)) &
				
						- 2.d0*rE(ip_th(count),ip_psi(count)) &
				
						+ rE(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
					
				+ a(count)*b(count)*(rE(ip_th(count)+1,ip_psi(count)+1) - rE(ip_th(count)+1,ip_psi(count)-1) &
	
					- rE(ip_th(count)-1,ip_psi(count)+1) + rE(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
	
	
	end do	

	!********************************************************************************************************************

	!For initial area calculations
	call Derivative_all(rE, rA, Npsi, Ntheta, dpsi, dtheta, d1rEall, d2rEall, d1rAall, d2rAall)
	
	sqgEall=SQRT(rE*rE*(rE*rE*cos_sq_th_all + cos_sq_th_all*d1rEall*d1rEall + d2rEall*d2rEall))
	
	sqgAall=SQRT(rA*rA*(rA*rA*cos_sq_th_all + cos_sq_th_all*d1rAall*d2rAall + d2rAall*d2rAall))

	area0=0.d0
	area1=0.d0
	
	do count_a_mod=0,endcount_a_mod,1
		j=area_psi_mod(count_a_mod)
		i=area_th_mod(count_a_mod)
		weight_elem=weight_mod(count_a_mod)*dpsi*dtheta
		
		!Total area
		area0=area0 + (sqgEall(i,j) + sqgAall(i,j))*weight_elem
			
		!Area of phase 1 (alpha phase, phi=1)
		area1=area1 + (phiE(i,j)*sqgEall(i,j) + phiA(i,j)*sqgAall(i,j))*weight_elem					   			   
	
	end do
	
	area0diff=area1 - (area0 - area1)
	!********************************************************************************************************************


end if

call MPI_BCAST(area_psi_mod(0), 23117, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
call MPI_BCAST(area_th_mod(0), 23117, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
call MPI_BCAST(weight_mod(0), 23117, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
call MPI_BCAST(area0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
call MPI_BCAST(area0diff, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)

endcount_a_mod=23116
!psi_start and psi_end set up so that each processor works with a domain
!encompassing the entire theta range (-49,49) => 99 elements, and 74 elements in psi (approx. 1/4 of psi range)

!Size of entire block in each processor, including all ghost zones (sent from processor 0)
size_total=(Ntheta+9)*(Max_psi/2+4)!99*74
size_out(:)=size_total

!Size of block in each processor sent back to processor 0 (block minus ghost zones in psi direction)
size_noghost=(Ntheta+9)*(Max_psi/2-1)!99*69
size_in(:)=size_noghost

size_noghost_psi=size_noghost/(Ntheta+9)

!Amount of displacement between starting points of where the good data from processors is sent back to proc. 0
!Also the amount of displacement between starting points of where the updated data is sent to worker processors
!Similar to size_noghost, but needs to account for the extra displacement needed to skip over the ghost zones
disp=(Ntheta+9)*(Max_psi/2-1)!99*69

do i=0,nprocs-1
	displ(i)=disp*i
end do

!psi_start and psi_end set up so that each processor works with a domain
!encompassing the entire theta range (-49,49) => 99 elements, and 74 elements in psi (approx. 1/4 of psi range)
!Each of these processor domains overlap by 5 points to account for ghost zone needs

psi_start=Min_psi + myrank*size_noghost_psi
psi_end=myrank*size_noghost_psi - (size_noghost_psi - 2)

call MPI_Barrier(MPI_COMM_WORLD, err)
!********************************************************************************************************************

do t=1,endt

	call MPI_Scatterv(phiE(Min_th,Min_psi), size_out, displ, MPI_DOUBLE_PRECISION, phiE(Min_th,psi_start), &
		size_total, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)	
		
	call MPI_Scatterv(phiA(Min_th,Min_psi), size_out, displ, MPI_DOUBLE_PRECISION, phiA(Min_th,psi_start), &
		size_total, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)	
		
	call MPI_Scatterv(rE(Min_th,Min_psi), size_out, displ, MPI_DOUBLE_PRECISION, rE(Min_th,psi_start), &
		size_total, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)	
		
	call MPI_Scatterv(rA(Min_th,Min_psi), size_out, displ, MPI_DOUBLE_PRECISION, rA(Min_th,psi_start), &
		size_total, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)	
		
	call Derivative(phiE(:,psi_start:psi_end), phiA(:,psi_start:psi_end), Max_psi, Npsi, Ntheta, dpsi, dtheta, &
					d1phiE, d2phiE, d12phiE, d11phiE, d22phiE, d1phiA, d2phiA, d12phiA, d11phiA, d22phiA)
					
	call Derivative(rE(:,psi_start:psi_end), rA(:,psi_start:psi_end), Max_psi, Npsi, Ntheta, dpsi, dtheta, &
					d1rE, d2rE, d12rE, d11rE, d22rE, d1rA, d2rA, d12rA, d11rA, d22rA)
	
	pE=phiE(:,psi_start:psi_end)*phiE(:,psi_start:psi_end)*phiE(:,psi_start:psi_end)*(10.d0 - 15.d0*phiE(:,psi_start:psi_end) &
	
		+ 6.d0*phiE(:,psi_start:psi_end)*phiE(:,psi_start:psi_end))
	
	pA=phiA(:,psi_start:psi_end)*phiA(:,psi_start:psi_end)*phiA(:,psi_start:psi_end)*(10.d0 - 15.d0*phiA(:,psi_start:psi_end) &
	
		+ 6.d0*phiA(:,psi_start:psi_end)*phiA(:,psi_start:psi_end))
	
	p_primeE=30.d0*phiE(:,psi_start:psi_end)*phiE(:,psi_start:psi_end)*(phiE(:,psi_start:psi_end)*phiE(:,psi_start:psi_end) &
	
		- 2.d0*phiE(:,psi_start:psi_end) + 1.d0)

	p_primeA=30.d0*phiA(:,psi_start:psi_end)*phiA(:,psi_start:psi_end)*(phiA(:,psi_start:psi_end)*phiA(:,psi_start:psi_end) &
	
		- 2.d0*phiA(:,psi_start:psi_end) + 1.d0)
	
	rsqE=rE(:,psi_start:psi_end)*rE(:,psi_start:psi_end)
	
	rsqA=rA(:,psi_start:psi_end)*rA(:,psi_start:psi_end)
	
	d1rsqE=d1rE*d1rE
	
	d1rsqA=d1rA*d1rA
	
	d2rsqE=d2rE*d2rE
	
	d2rsqA=d2rA*d2rA
	
	gE=rsqE*(rsqE*cos_sq_th + cos_sq_th*d1rsqE + d2rsqE)

	gA=rsqA*(rsqA*cos_sq_th + cos_sq_th*d1rsqA + d2rsqA)
	
	sqgE=SQRT(gE)
	sqgA=SQRT(gA)
	
	!Here, g11, g12, and g22 represent the inverse metric (numbers are superscripts)
	g11E=(rsqE*cos_sq_th + d2rsqE)/gE
	g11A=(rsqA*cos_sq_th + d2rsqA)/gA
	
	g12E=-d1rE*d2rE/gE
	g12A=-d1rA*d2rA/gA
	
	g22E=(rsqE + d1rsqE)/gE
	g22A=(rsqA + d1rsqA)/gA
	
	!K11, K12, and K22 represent elements in the curvature tensor.  Numbers are subscripts.
	K11E=rE(:,psi_start:psi_end)/SQRT(gE)*cos_th*(rsqE + 2.d0*d1rsqE - rE(:,psi_start:psi_end)*d11rE)

	K11A=rA(:,psi_start:psi_end)/SQRT(gA)*cos_th*(rsqA + 2.d0*d1rsqA - rA(:,psi_start:psi_end)*d11rA)
	
	K12E=rE(:,psi_start:psi_end)/SQRT(gE)*(cos_th*(2.d0*d1rE*d2rE - rE(:,psi_start:psi_end)*d12rE) &
	
		- sin_th*rE(:,psi_start:psi_end)*d2rE)

	K12A=rA(:,psi_start:psi_end)/SQRT(gA)*(cos_th*(2.d0*d1rA*d2rA - rA(:,psi_start:psi_end)*d12rA) &
	
		- sin_th*rA(:,psi_start:psi_end)*d2rA)
	
	K22E=rE(:,psi_start:psi_end)/SQRT(gE)*cos_th*(rsqE*cos_sq_th + 2.d0*d2rsqE &
	
		+ sin_th*cos_th*rE(:,psi_start:psi_end)*d1rE - rE(:,psi_start:psi_end)*d22rE)

	K22A=rA(:,psi_start:psi_end)/SQRT(gA)*cos_th*(rsqA*cos_sq_th + 2.d0*d2rsqA &
	
		+ sin_th*cos_th*rA(:,psi_start:psi_end)*d1rA - rA(:,psi_start:psi_end)*d22rA)
	
	KE=g11E*K11E + 2.d0*g12E*K12E + g22E*K22E
	KA=g11A*K11A + 2.d0*g12A*K12A + g22A*K22A
	
	SE =  -1.d0*(K11E*(K11E*g11E*g11E + 2.d0*K12E*g11E*g12E + K22E*g12E*g12E) &
	
		+ 2.d0*K12E*(K11E*g11E*g12E + K12E*(g11E*g22E + g12E*g12E) + K22E*g12E*g22E) &
		
		+ K22E*(K22E*g22E*g22E + 2.d0*K12E*g12E*g22E + K11E*g12E*g12E))
		
	SA =  -1.d0*(K11A*(K11A*g11A*g11A + 2.d0*K12A*g11A*g12A + K22A*g12A*g12A) &
	
		+ 2.d0*K12A*(K11A*g11A*g12A + K12A*(g11A*g22A + g12A*g12A) + K22A*g12A*g22A) &
		
		+ K22A*(K22A*g22A*g22A + 2.d0*K12A*g12A*g22A + K11A*g12A*g12A))
		
	Gam1E=-1.d0/(gE*gE)*rsqE*cos_th*(-rsqE*rsqE*sin_th*cos_sq_th &
		
		+ rE(:,psi_start:psi_end)*cos_th*d1rE*(cos_sq_th*d1rsqE + d2rsqE) &
				
		- d1rE*(d2rsqE*(-sin_th*d1rE + cos_th*d11rE) &
				
			- 2.d0*cos_th*d1rE*d2rE*d12rE + cos_th*d1rsqE*d22rE) &
					
		+ rsqE*(-sin_th*d2rsqE - cos_th*d1rE &
		
		*(sin_th*cos_th*d1rE + cos_sq_th*d11rE + d22rE)))
	
	Gam1A=-1.d0/(gA*gA)*rsqA*cos_th*(-rsqA*rsqA*sin_th*cos_sq_th &
		
		+ rA(:,psi_start:psi_end)*cos_th*d1rA*(cos_sq_th*d1rsqA + d2rsqA) &
				
		- d1rA*(d2rsqA*(-sin_th*d1rA + cos_th*d11rA) &
				
			- 2.d0*cos_th*d1rA*d2rA*d12rA + cos_th*d1rsqA*d22rA) &
					
		+ rsqA*(-sin_th*d2rsqA - cos_th*d1rA &
		
		*(sin_th*cos_th*d1rA + cos_sq_th*d11rA + d22rA)))
				
	Gam2E=1.d0/(gE*gE)*rsqE*d2rE*(cos_th*sin_th*d1rsqE*d1rE &
	
		- rE(:,psi_start:psi_end)*(cos_sq_th*d1rsqE + d2rsqE) + d2rsqE*d11rE &
		
		- 2.d0*d1rE*d2rE*d12rE + d1rsqE*d22rE &
		
		+ rsqE*(cos_th*(sin_th*d1rE + cos_th*d11rE) + d22rE))
		
	Gam2A=1.d0/(gA*gA)*rsqA*d2rA*(cos_th*sin_th*d1rsqA*d1rA &
	
		- rA(:,psi_start:psi_end)*(cos_sq_th*d1rsqA + d2rsqA) + d2rsqA*d11rA &
		
		- 2.d0*d1rA*d2rA*d12rA + d1rsqA*d22rA &
		
		+ rsqA*(cos_th*(sin_th*d1rA + cos_th*d11rA) + d22rA))
		
	
	Lap_phiE=g11E*d11phiE + 2.d0*g12E*d12phiE + g22E*d22phiE - d1phiE*Gam1E - d2phiE*Gam2E
				
	Lap_phiA=g11A*d11phiA + 2.d0*g12A*d12phiA + g22A*d22phiA - d1phiA*Gam1A - d2phiA*Gam2A
		
	grad_phiE=g11E*d1phiE*d1phiE + 2.d0*g12E*d1phiE*d2phiE + g22E*d2phiE*d2phiE
	
	grad_phiA=g11A*d1phiA*d1phiA + 2.d0*g12A*d1phiA*d2phiA + g22A*d2phiA*d2phiA
	
	QE=pE*(KE - C1) + (1.d0 - pE)*(KE - C2)
	
	QA=pA*(KA - C1) + (1.d0 - pA)*(KA - C2)
	
	call Derivative(QE, QA, Max_psi, Npsi, Ntheta, dpsi, dtheta, &
					d1QE, d2QE, d12QE, d11QE, d22QE, d1QA, d2QA, d12QA, d11QA, d22QA)	
	
	Lap_QE=g11E*d11QE + 2.d0*g12E*d12QE + g22E*d22QE - d1QE*Gam1E - d2QE*Gam2E
				
	Lap_QA=g11A*d11QA + 2.d0*g12A*d12QA + g22A*d22QA - d1QA*Gam1A - d2QA*Gam2A
	
	!"Kab" actually represents Kab*del_sub_a(phi)*del_sub_b(phi), the gradient curvature contraction of phi
	
	KabE = d1phiE*d1phiE*(K11E*g11E*g11E + 2.d0*K12E*g11E*g12E + K22E*g12E*g12E) &
	
		+ 2.d0*d1phiE*d2phiE*(K11E*g11E*g12E + K12E*(g11E*g22E + g12E*g12E) + K22E*g12E*g22E) &
		
		+ d2phiE*d2phiE*(K22E*g22E*g22E + 2.d0*K12E*g12E*g22E + K11E*g12E*g12E)
		
	KabA = d1phiA*d1phiA*(K11A*g11A*g11A + 2.d0*K12A*g11A*g12A + K22A*g12A*g12A) &
	
		+ 2.d0*d1phiA*d2phiA*(K11A*g11A*g12A + K12A*(g11A*g22A + g12A*g12A) + K22A*g12A*g22A) &
		
		+ d2phiA*d2phiA*(K22A*g22A*g22A + 2.d0*K12A*g12A*g22A + K11A*g12A*g12A)	
		
		
	TmE=-(KE*(0.25d0*W*phiE(:,psi_start:psi_end)*phiE(:,psi_start:psi_end) &
	
		*(1.d0-phiE(:,psi_start:psi_end))*(1.d0-phiE(:,psi_start:psi_end)) &
	
		+ 0.5d0*eps*eps*grad_phiE &
		
		+ 0.5d0*lambda*(pE*(KE-C1)*(KE-C1) + (1.d0-pE)*(KE-C2)*(KE-C2))) &
		
		+ lambda*QE*SE - lambda*Lap_QE - eps*eps*KabE - P)
		
	TmA=-(KA*(0.25d0*W*phiA(:,psi_start:psi_end)*phiA(:,psi_start:psi_end)&
	
		*(1.d0-phiA(:,psi_start:psi_end))*(1.d0-phiA(:,psi_start:psi_end)) &
	
		+ 0.5d0*eps*eps*grad_phiA &
		
		+ 0.5d0*lambda*(pA*(KA-C1)*(KA-C1) + (1.d0-pA)*(KA-C2)*(KA-C2))) &
		
		+ lambda*QA*SA - lambda*Lap_QA - eps*eps*KabA - P)

	sum_arr(:)=0.d0

	do count_a_mod=0,endcount_a_mod,1
		i=area_th_mod(count_a_mod)
		j_all=area_psi_mod(count_a_mod)
		weight_elem=weight_mod(count_a_mod)*dpsi*dtheta

		if(j_all.GE.psi_start+3 .AND. j_all.LE.psi_end-2) then

			j=j_all + (1 - psi_start)

			!Element A in notes
			sum_arr(1)=sum_arr(1) + (KE(i,j)*KE(i,j)*sqgE(i,j) + KA(i,j)*KA(i,j)*sqgA(i,j))*weight_elem

			!Element B (and C) in notes			
			sum_arr(2)=sum_arr(2) + (KE(i,j)*KE(i,j)*(2.d0*phiE(i,j_all) - 1.d0)*sqgE(i,j) &
			
								  +  KA(i,j)*KA(i,j)*(2.d0*phiA(i,j_all) - 1.d0)*sqgA(i,j))*weight_elem

			!Element D in notes	
			sum_arr(3)=sum_arr(3) + (KE(i,j)*KE(i,j) &
			
								*(2.d0*phiE(i,j_all) - 1.d0)*(2.d0*phiE(i,j_all) - 1.d0)*sqgE(i,j) &
			
					   			  + KA(i,j)*KA(i,j) &
					   			  
					   			*(2.d0*phiA(i,j_all) - 1.d0)*(2.d0*phiA(i,j_all) - 1.d0)*sqgA(i,j))*weight_elem					 
							 
			!Element E in notes	(first part of it - second part is subtracted after loop here)
			sum_arr(4)=sum_arr(4) + (TmE(i,j)*KE(i,j)*sqgE(i,j) + TmA(i,j)*KA(i,j)*sqgA(i,j))*weight_elem			 
					
			!Element F in notes	(first part of it - second part is subtracted after loop here)				
			sum_arr(5)=sum_arr(5) + (TmE(i,j)*KE(i,j)*(2.d0*phiE(i,j_all) - 1.d0)*sqgE(i,j) &
			
					   			   + TmA(i,j)*KA(i,j)*(2.d0*phiA(i,j_all) - 1.d0)*sqgA(i,j))*weight_elem
			!Total area
			sum_arr(6)=sum_arr(6) + (sqgE(i,j) + sqgA(i,j))*weight_elem
			
			!Area of phase 1 (alpha phase, phi=1)
			sum_arr(7)=sum_arr(7) + (phiE(i,j_all)*sqgE(i,j) + phiA(i,j_all)*sqgA(i,j))*weight_elem
		end if

	end do
	
	call MPI_Allreduce(sum_arr(1), sums(1), 7, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)

	!Subtract second part of element E
	sums(4)=sums(4) - 1.d0/(gam*dt)*(area0 - sums(6))

	!Subtract second part of element F
	sums(5)=sums(5) - 1.d0/(gam*dt)*(area0diff - (sums(7) - (sums(6) - sums(7))))

	sigma=(sums(5)*sums(2) - sums(3)*sums(4))/(sums(2)*sums(2) - sums(3)*sums(1))
!	sigma=(sums(5) - sums(3)*sums(4)/sums(2))/(sums(2) - sums(3)*sums(1)/sums(2))
	
	dsig=(sums(2)*sums(4)/sums(1) - sums(5))/(sums(2)*sums(2)/sums(1) - sums(3))

	muE=0.5d0*W*phiE(:,psi_start:psi_end)*(1.d0 - phiE(:,psi_start:psi_end)) &
	
		*(1.d0 - 2.d0*phiE(:,psi_start:psi_end)) - eps*eps*Lap_phiE &
	
		+ lambda*p_primeE*(KE*(C2 - C1) + 0.5d0*(C1*C1 - C2*C2)) + 2.d0*dsig

	muA=0.5d0*W*phiA(:,psi_start:psi_end)*(1.d0 - phiA(:,psi_start:psi_end)) &
	
		*(1.d0 - 2.d0*phiA(:,psi_start:psi_end)) - eps*eps*Lap_phiA &
	
		+ lambda*p_primeA*(KA*(C2 - C1) + 0.5d0*(C1*C1 - C2*C2)) + 2.d0*dsig
	
	call Derivative(muE, muA, Max_psi, Npsi, Ntheta, dpsi, dtheta, &
					d1muE, d2muE, d12muE, d11muE, d22muE, d1muA, d2muA, d12muA, d11muA, d22muA)
	
	Lap_muE=g11E*d11muE + 2.d0*g12E*d12muE + g22E*d22muE - d1muE*Gam1E - d2muE*Gam2E
				
	Lap_muA=g11A*d11muA + 2.d0*g12A*d12muA + g22A*d22muA - d1muA*Gam1A - d2muA*Gam2A

	drtE=dt*gam*SQRT(gE)/(rsqE*cos_th)*(TmE - KE*sigma - KE*(2.d0*phiE(:,psi_start:psi_end) - 1.d0)*dsig)
		
	drtA=dt*gam*SQRT(gA)/(rsqA*cos_th)*(TmA - KA*sigma - KA*(2.d0*phiA(:,psi_start:psi_end) - 1.d0)*dsig)

	phi_adE=drtE*(d1phiE*(d1rE*g11E + d2rE*g12E) + d2phiE*(d1rE*g12E + d2rE*g22E))
	
	phi_adA=drtA*(d1phiA*(d1rA*g11A + d2rA*g12A) + d2phiA*(d1rA*g12A + d2rA*g22A))
	
	dphitE=dt*(M*Lap_muE + phi_adE)
	dphitA=dt*(M*Lap_muA + phi_adA)
	
	phiE(:,psi_start:psi_end)=phiE(:,psi_start:psi_end) + dphitE
	phiA(:,psi_start:psi_end)=phiA(:,psi_start:psi_end) + dphitA
	
	rE(:,psi_start:psi_end)=rE(:,psi_start:psi_end) + drtE
	rA(:,psi_start:psi_end)=rA(:,psi_start:psi_end) + drtA
	
	call MPI_Gatherv(phiE(Min_th,psi_start+3), size_noghost, MPI_DOUBLE_PRECISION, phiE(Min_th,psi_start+3), &
		size_in, displ, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
		
	call MPI_Gatherv(phiA(Min_th,psi_start+3), size_noghost, MPI_DOUBLE_PRECISION, phiA(Min_th,psi_start+3), &
		size_in, displ, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
		
	call MPI_Gatherv(rE(Min_th,psi_start+3), size_noghost, MPI_DOUBLE_PRECISION, rE(Min_th,psi_start+3), &
		size_in, displ, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
		
	call MPI_Gatherv(rA(Min_th,psi_start+3), size_noghost, MPI_DOUBLE_PRECISION, rA(Min_th,psi_start+3), &
		size_in, displ, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)


	!Send updated arrays back to processor 0
	if(myrank==0) then
	
		do count=0,endcount
			
			phiE(Corn_th(count),Corn_psi(count))=phiA(ip_th(count),ip_psi(count)) &
			
					+ a(count)*(phiA(ip_th(count),ip_psi(count)+1)-phiA(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
		
					+ b(count)*(phiA(ip_th(count)+1,ip_psi(count))-phiA(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
			
					+ 0.5d0*a(count)*a(count)*(phiA(ip_th(count),ip_psi(count)+1) &
					
							- 2.d0*phiA(ip_th(count),ip_psi(count)) &
					
							+ phiA(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
						
					+ 0.5d0*b(count)*b(count)*(phiA(ip_th(count)+1,ip_psi(count)) &
					
							- 2.d0*phiA(ip_th(count),ip_psi(count)) &
					
							+ phiA(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
						
					+ a(count)*b(count)*(phiA(ip_th(count)+1,ip_psi(count)+1) - phiA(ip_th(count)+1,ip_psi(count)-1) &
		
						- phiA(ip_th(count)-1,ip_psi(count)+1) + phiA(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
						
		
			phiA(Corn_th(count),Corn_psi(count))=phiE(ip_th(count),ip_psi(count)) &
			
					+ a(count)*(phiE(ip_th(count),ip_psi(count)+1)-phiE(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
		
					+ b(count)*(phiE(ip_th(count)+1,ip_psi(count))-phiE(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
			
					+ 0.5d0*a(count)*a(count)*(phiE(ip_th(count),ip_psi(count)+1) &
					
							- 2.d0*phiE(ip_th(count),ip_psi(count)) &
					
							+ phiE(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
						
					+ 0.5d0*b(count)*b(count)*(phiE(ip_th(count)+1,ip_psi(count)) &
					
							- 2.d0*phiE(ip_th(count),ip_psi(count)) &
					
							+ phiE(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
						
					+ a(count)*b(count)*(phiE(ip_th(count)+1,ip_psi(count)+1) - phiE(ip_th(count)+1,ip_psi(count)-1) &
		
						- phiE(ip_th(count)-1,ip_psi(count)+1) + phiE(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
						
						
			rE(Corn_th(count),Corn_psi(count))=rA(ip_th(count),ip_psi(count)) &
			
					+ a(count)*(rA(ip_th(count),ip_psi(count)+1)-rA(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
		
					+ b(count)*(rA(ip_th(count)+1,ip_psi(count))-rA(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
			
					+ 0.5d0*a(count)*a(count)*(rA(ip_th(count),ip_psi(count)+1) &
					
							- 2.d0*rA(ip_th(count),ip_psi(count)) &
					
							+ rA(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
						
					+ 0.5d0*b(count)*b(count)*(rA(ip_th(count)+1,ip_psi(count)) &
					
							- 2.d0*rA(ip_th(count),ip_psi(count)) &
					
							+ rA(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
						
					+ a(count)*b(count)*(rA(ip_th(count)+1,ip_psi(count)+1) - rA(ip_th(count)+1,ip_psi(count)-1) &
		
						- rA(ip_th(count)-1,ip_psi(count)+1) + rA(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
						
		
			rA(Corn_th(count),Corn_psi(count))=rE(ip_th(count),ip_psi(count)) &
			
					+ a(count)*(rE(ip_th(count),ip_psi(count)+1)-rE(ip_th(count),ip_psi(count)-1))/(2.d0*dpsi) &
		
					+ b(count)*(rE(ip_th(count)+1,ip_psi(count))-rE(ip_th(count)-1,ip_psi(count)))/(2.d0*dtheta) &
			
					+ 0.5d0*a(count)*a(count)*(rE(ip_th(count),ip_psi(count)+1) &
					
							- 2.d0*rE(ip_th(count),ip_psi(count)) &
					
							+ rE(ip_th(count),ip_psi(count)-1))/(dpsi*dpsi) &
						
					+ 0.5d0*b(count)*b(count)*(rE(ip_th(count)+1,ip_psi(count)) &
					
							- 2.d0*rE(ip_th(count),ip_psi(count)) &
					
							+ rE(ip_th(count)-1,ip_psi(count)))/(dtheta*dtheta) &
						
					+ a(count)*b(count)*(rE(ip_th(count)+1,ip_psi(count)+1) - rE(ip_th(count)+1,ip_psi(count)-1) &
		
						- rE(ip_th(count)-1,ip_psi(count)+1) + rE(ip_th(count)-1,ip_psi(count)-1))/(4.d0*dpsi*dtheta) 
		
		
		end do	
	
		!export results to text file
		if (mod(real(t,8),real(plott,8))==0) then
			cntt=cntt+plotstep
			call expdata(cntt, phiE, phiA, rE, rA, Npsi, Ntheta, Min_psi, Max_psi, Min_th, Max_th)
		end if
		
		!To output data less often at later times when evolution is slow
		!'plot#1' means 1000 times through time loop
		!this equality does not change, though not every plot# is output
		!For plot#0-10, output after every plot
		!For plot#10-50, output after every 5 plots
		!For plot#50-100, output after every 10 plots	
		!For plot#100-1000, output after every 50 plots	
		!For plot#1000-5000, output after every 250 plots
		!For plot#5000-10,000, output after every 500 plots
		!For plot#10,000-50,000, output after every 1000 plots
		!For plot#50,000-, output after every 5000 plots
		
!		realt=t + 1000*cntt_orig
!		
!		if(reduce==1) then
!			if(realt==50000 .OR. realt==5000000 .OR. realt==10000000) then		
!				plotstep=plotstep*2
!				plott=1000*plotstep
!			end if
!			
!			if(realt==10000 .OR. realt==100000 .OR. realt==1000000 .OR. realt==50000000) then	
!				plotstep=plotstep*5
!				plott=1000*plotstep
!			end if
!		end if

		
	end if
end do

call MPI_FINALIZE(err)

end program

!********************************
!********************************

!subroutine for calculating Derivatives

subroutine Derivative(MeshE, MeshA, Max_psi, Npsi, Ntheta, dpsi, dtheta, &					
					  d1E, d2E, d12E, d11E, d22E, d1A, d2A, d12A, d11A, d22A)

integer(kind=4), intent(in) :: Max_psi, Npsi, Ntheta
real(kind=8), intent(in) :: dpsi, dtheta
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)), intent(in) :: MeshE, MeshA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)), intent(out) :: d1E, d2E, d12E, d11E, d22E
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),1:((Npsi+26)/4)), intent(out) :: d1A, d2A, d12A, d11A, d22A
integer(kind=4) :: Min_p, Min_t, Max_p, Max_t

Min_p=2
Max_p=Max_psi/2+3

Min_t=-Ntheta/2-3
Max_t=Ntheta/2+3


d1E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t+1:Max_t+1,Min_p:Max_p)-MeshE(Min_t-1:Max_t-1,Min_p:Max_p))/(2.d0*dtheta)

d1A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t+1:Max_t+1,Min_p:Max_p)-MeshA(Min_t-1:Max_t-1,Min_p:Max_p))/(2.d0*dtheta)

d2E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t:Max_t,Min_p+1:Max_p+1)-MeshE(Min_t:Max_t,Min_p-1:Max_p-1))/(2.d0*dpsi)

d2A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t:Max_t,Min_p+1:Max_p+1)-MeshA(Min_t:Max_t,Min_p-1:Max_p-1))/(2.d0*dpsi)


d11E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t+1:Max_t+1,Min_p:Max_p) + MeshE(Min_t-1:Max_t-1,Min_p:Max_p) &

								- 2.d0*MeshE(Min_t:Max_t,Min_p:Max_p))/(dtheta*dtheta)
								
d11A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t+1:Max_t+1,Min_p:Max_p) + MeshA(Min_t-1:Max_t-1,Min_p:Max_p) &

								- 2.d0*MeshA(Min_t:Max_t,Min_p:Max_p))/(dtheta*dtheta)

d22E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t:Max_t,Min_p+1:Max_p+1) + MeshE(Min_t:Max_t,Min_p-1:Max_p-1) &

								- 2.d0*MeshE(Min_t:Max_t,Min_p:Max_p))/(dpsi*dpsi)
								
d22A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t:Max_t,Min_p+1:Max_p+1) + MeshA(Min_t:Max_t,Min_p-1:Max_p-1) &

								- 2.d0*MeshA(Min_t:Max_t,Min_p:Max_p))/(dpsi*dpsi)

d12E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t+1:Max_t+1,Min_p+1:Max_p+1) - MeshE(Min_t+1:Max_t+1,Min_p-1:Max_p-1) &

				- MeshE(Min_t-1:Max_t-1,Min_p+1:Max_p+1) + MeshE(Min_t-1:Max_t-1,Min_p-1:Max_p-1))/(4.d0*dtheta*dpsi)
								
d12A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t+1:Max_t+1,Min_p+1:Max_p+1) - MeshA(Min_t+1:Max_t+1,Min_p-1:Max_p-1) &

				- MeshA(Min_t-1:Max_t-1,Min_p+1:Max_p+1) + MeshA(Min_t-1:Max_t-1,Min_p-1:Max_p-1))/(4.d0*dtheta*dpsi)


end subroutine Derivative
!********************************

!subroutine for calculating Derivatives on entire mesh, not split up MPI pieces
!This is necessary for the area0,area1 calculations done on the master processor for initialization only

subroutine Derivative_all(MeshE, MeshA, Npsi, Ntheta, dpsi, dtheta, &					
					  d1E, d2E, d1A, d2A)

integer(kind=4), intent(in) :: Npsi, Ntheta
real(kind=8), intent(in) :: dpsi, dtheta
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)), intent(in) :: MeshE, MeshA
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)), intent(out) :: d1E, d2E
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)), intent(out) :: d1A, d2A
integer(kind=4) :: Min_p, Min_t, Max_p, Max_t

Min_p=-Npsi/2-4
Max_p=Npsi/2+4

Min_t=-Ntheta/2-3
Max_t=Ntheta/2+3


d1E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t+1:Max_t+1,Min_p:Max_p)-MeshE(Min_t-1:Max_t-1,Min_p:Max_p))/(2.d0*dtheta)

d1A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t+1:Max_t+1,Min_p:Max_p)-MeshA(Min_t-1:Max_t-1,Min_p:Max_p))/(2.d0*dtheta)

d2E(Min_t:Max_t,Min_p:Max_p)=(MeshE(Min_t:Max_t,Min_p+1:Max_p+1)-MeshE(Min_t:Max_t,Min_p-1:Max_p-1))/(2.d0*dpsi)

d2A(Min_t:Max_t,Min_p:Max_p)=(MeshA(Min_t:Max_t,Min_p+1:Max_p+1)-MeshA(Min_t:Max_t,Min_p-1:Max_p-1))/(2.d0*dpsi)


end subroutine Derivative_all

!********************************


!********************************

!subroutine to export data to a text file
subroutine expdata(cntt, phiE, phiA, rE, rA, Npsi, Ntheta, Min_psi, Max_psi, Min_th, Max_th)

integer(kind=4), intent(in) :: Npsi, Ntheta
integer(kind=4), intent(in) :: cntt, Min_psi, Min_th, Max_psi, Max_th
real(kind=8), dimension((-Ntheta/2-4):(Ntheta/2+4),(-Npsi/2-5):(Npsi/2+5)), intent(in) :: phiE, rE, phiA, rA
integer(kind=4) :: plotnum
integer(kind=4) :: i,j
character :: namePE*12, step*6, namePA*12, namerE*12, namerA*12

plotnum=100000+cntt

!If plot labels go over 99999, start over at 10000
if (plotnum.GT.999999) then
	plotnum=plotnum-900000
end if

!change numerical file names to strings and label as .txt files
write(step, "(I6)") plotnum

namePE=''//step //'pE.txt'
namePA=''//step //'pA.txt'

namerE=''//step //'rE.txt'
namerA=''//step //'rA.txt'

open(unit=10, file=namePE, status='replace')
	write(10,996) ((phiE(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
	996 format(99(E16.9, X))
close(unit=10, status='keep')

open(unit=10, file=namePA, status='replace')
	write(10,997) ((phiA(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
	997 format(99(E16.9, X))
close(unit=10, status='keep')


open(unit=10, file=namerE, status='replace')
	write(10,998) ((rE(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
	998 format(99(E16.9, X))
close(unit=10, status='keep')

open(unit=10, file=namerA, status='replace')
	write(10,999) ((rA(i,j), i=Min_th,Max_th,1), j=Min_psi,Max_psi,1)
	999 format(99(E16.9, X))
close(unit=10, status='keep')



end subroutine expdata

!*************************



function r4_uniform_01 ( seed )


!! R4_UNIFORM_01 returns a unit real ( kind = 4 ) pseudorandom number.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
