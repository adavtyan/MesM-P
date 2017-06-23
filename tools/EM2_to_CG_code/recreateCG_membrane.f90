program recreateCG_membrane

   !implicit none
   Implicit Double Precision(a-h,o-z)
   
   Integer, Parameter                                :: DP = Selected_real_kind(15,307)
   integer											 :: m,additional_tri
   integer, Dimension(:), Allocatable				 :: atom_one,atom_two,atom_three  !This is to distinguish top and bottom layer
   Integer, Dimension(:), Allocatable				 :: atom_s_one,atom_s_two,atom_s_three
   Real(kind=DP), Dimension(:) , Allocatable         :: rx,ry,rz
   Real(kind=DP), Dimension(:) , Allocatable		 :: rx_one,ry_one,rz_one,rx_two,ry_two,rz_two,rx_three,ry_three,rz_three
   Integer, Dimension(:), Allocatable                :: print_flag,NGB_in_TRI,gb_patch_flag, keep_em_flag
   Integer, Dimension(:,:), Allocatable              :: TRI_to_GB_map,triangles,edges
   Real(kind=DP), Dimension(:) , Allocatable         :: rx_em,ry_em,rz_em
   Real(kind=DP), Dimension(:,:) , Allocatable       :: EM_distance
   Real(kind=DP), Dimension(:) , Allocatable         :: rx_emm,ry_emm,rz_emm,rx_tri,ry_tri,rz_tri
   Real(kind=DP), Dimension(:) , Allocatable         :: rx_s,ry_s,rz_s
   Real(kind=DP), Dimension(:) , Allocatable		 :: rx_s_one,ry_s_one,rz_s_one,rx_s_two,ry_s_two,&
   														rz_s_two,rx_s_three,ry_s_three,rz_s_three
   Real(kind=DP), Dimension(:) , Allocatable         :: sorntx,sornty,sorntz
   Real(kind=DP), Dimension(:) , Allocatable         :: qw,qx,qy,qz
   Real(kind=DP), Dimension(:) , Allocatable         :: qw_em,qx_em,qy_em,qz_em,qw_emm,qx_emm,qy_emm,qz_emm
   Real(kind=DP), Dimension(:) , Allocatable         :: qw_s,qx_s,qy_s,qz_s
   
   character(len=31)								 :: config, patch, output_xyz, output_em2
   real												 :: cut_1, cut_2
   Real(kind=DP)									 :: normx,normy,normz,norm_mod,phi
   Real(kind=DP)									 :: qw_av,qx_av,qy_av,qz_av

   real*8 :: xlo, xhi, ylo, yhi, zlo, zhi

   ncfile = 9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!															!!!
!!		Read in the preequilibrated GB patch	!!!
!!															!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,*) "Membrane patch file in Angstroms? (The atom desig should be integers)"
   read(*,*) patch
   write(*,*) "Name of the EM2 file (membrane only). Max 31 chars please."
   read(*,*) config  
   !write(*,*) "Reading configuration file. Number of EM sites?"    !Membrane particles only
   !read(*,*) n_em
   write(*,*) "Output xyz file? The output will be in nm"
   read(*,*) output_xyz
   write(*,*) "Output em2 file?"
   read(*,*) output_em2

! read in the pre-equilibrated patch of GB membrane
   open(unit=ncfile, file=patch)
   read(ncfile,*) ngb
   read(ncfile,*)
   
   ngb=ngb/3
   
   allocate(rx_one(ngb),rx_two(ngb),rx_three(ngb))
   allocate(ry_one(ngb),ry_two(ngb),ry_three(ngb))
   allocate(rz_one(ngb),rz_two(ngb),rz_three(ngb))
   allocate(atom_one(ngb),atom_two(ngb),atom_three(ngb))

   allocate(sorntx(ngb))
   allocate(sornty(ngb))
   allocate(sorntz(ngb))

   !allocate(qw(ngb))
   !allocate(qx(ngb))
   !allocate(qy(ngb))
   !allocate(qz(ngb))

   !read(ncfile,*) boxx,boxy,boxz
   do igb = 1,ngb
	  read(ncfile,*) atom_one(igb), rx_one(igb),ry_one(igb),rz_one(igb)
	  read(ncfile,*) atom_two(igb), rx_two(igb),ry_two(igb),rz_two(igb)
	  read(ncfile,*) atom_three(igb), rx_three(igb),ry_three(igb),rz_three(igb)
   enddo

!do igb = 1,ngb
!	read(ncfile,*) qw(igb),qx(igb),qy(igb),qz(igb)
!   qq = qw(igb)**2 + qx(igb)**2 +qy(igb)**2 +qz(igb)**2
!   qq = sqrt(qq)
!	qw(igb) = qw(igb)/qq
!	qx(igb) = qx(igb)/qq
!	qy(igb) = qy(igb)/qq
!	qz(igb) = qz(igb)/qq
!enddo

close(ncfile)

rx_one = 0.10d0*rx_one
ry_one = 0.10d0*ry_one
rz_one = 0.10d0*rz_one
rx_two = 0.10d0*rx_two
ry_two = 0.10d0*ry_two
rz_two = 0.10d0*rz_two
rx_three = 0.10d0*rx_three
ry_three = 0.10d0*ry_three
rz_three = 0.10d0*rz_three

! EM2 sigma is 6.6
! peak in g(r) is about 7.5 nm
! GB sigma is 0.76

! shift CoM of GB patch to the origin

cmx = 0.d0
cmy = 0.d0
cmz = 0.d0

do igb=1,ngb
	cmx = cmx + rx_one(igb)
	cmx = cmx + rx_two(igb)
	cmx = cmx + rx_three(igb)
	cmy = cmy + ry_one(igb)
	cmy = cmy + ry_two(igb)
	cmy = cmy + ry_three(igb)
	cmz = cmz + rz_one(igb)
	cmz = cmz + rz_two(igb)
	cmz = cmz + rz_three(igb)
enddo

cmx = cmx/dble(3*ngb) 
cmy = cmy/dble(3*ngb) 
cmz = cmz/dble(3*ngb)

! shift the origin into the 3rd quadrant, since the triangle will be drawn with the lower left
! vertex at 0,0. don't want to run out of GB particle before the triangle is filled! 
!do igb=1,ngb
	!rx(igb) = rx(igb) - cmx
	!ry(igb) = ry(igb) - cmy
!	rx(igb) = rx(igb) - cmx + 0.30*boxx
!	ry(igb) = ry(igb) - cmy + 0.30*boxy
!	rz(igb) = rz(igb) - cmz
!enddo

sorntx_sum = 0.0
sornty_sum = 0.0
sorntz_sum = 0.0

sbudix = 0.d0
sbudiy = 0.d0
sbudiz = 1.d0

!do i=1,ngb
!	axx = qw(i)**2 + qx(i)**2 - qy(i)**2 - qz(i)**2
!	axy = 2.d0 * ( qx(i) * qy(i) - qw(i) * qz(i) )
!	axz = 2.d0 * ( qx(i) * qz(i) + qw(i) * qy(i) )
!	ayx = 2.d0 * ( qx(i) * qy(i) + qw(i) * qz(i) )
!	ayy = qw(i) **2 - qx(i)**2 + qy(i)**2 -qz(i)**2
!	ayz = 2.d0 * ( qy(i) * qz(i) - qw(i) * qx(i) )
!	azx = 2.d0 * ( qx(i) * qz(i) - qw(i) * qy(i) )
!	azy = 2.d0 * ( qy(i) * qz(i) + qw(i) * qx(i) )
!	azz = qw(i)**2 - qx(i)**2 - qy(i)**2 + qz(i)**2
!	sorntx(i) = ( axx * sbudix + axy * sbudiy + axz * sbudiz )
!	sornty(i) = ( ayx * sbudix + ayy * sbudiy + ayz * sbudiz )
!	sorntz(i) = ( azx * sbudix + azy * sbudiy + azz * sbudiz )

!	sorntx_sum = sorntx_sum + abs(sorntx(i))
!	sornty_sum = sornty_sum + abs(sornty(i))
!	sorntz_sum = sorntz_sum + abs(sorntz(i))
!enddo

!sorntx_avg = sorntx_sum/dble(ngb)
!sornty_avg = sornty_sum/dble(ngb)
!sorntz_avg = sorntz_sum/dble(ngb)

!write(*,*) "oritentation of GB patch: ", sorntx_avg, sornty_avg, sorntz_avg
   
!open(unit=ncfile,file='rasmol-gb.xyz')
!	Write(ncfile,'(i5)') 3*ngb
!	Write(ncfile,'("          ")')
!	Do i=1,ngb
!		Write(ncfile,'(I10,3f10.4)') atom_one(i), rx_one(i)-cmx,ry_one(i)-cmy,rz_one(i)-cmz
!		Write(ncfile,'(I10,3f10.4)') atom_two(i), rx_two(i)-cmx,ry_two(i)-cmy,rz_two(i)-cmz
!		Write(ncfile,'(I10,3f10.4)') atom_three(i), rx_three(i)-cmx,ry_three(i)-cmy,rz_three(i)-cmz
!		!Write(ncfile,'(a10,3f10.4)') 'P', (rx(i)-2.d-1*sorntx(i))*5.d0,(ry(i)-2.d-1*sornty(i))*5.d0,(rz(i)-2.d-1*sorntz(i))*5.d0
!	Enddo
!close(ncfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																		!!!
!!   done reading/centering preequilibrated GB patch	!!!
!!																		!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!									!!					
!!		read in EM2 config	!!
!!									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



Open(unit=ncfile, file=config)
! 3 lines at top of dump
read(ncfile,*)
read(ncfile,*)
read(ncfile,*)
read(ncfile,*) n_em ! read number of particles
!allocate arrays with number of particles
allocate(rx_emm(n_em))
allocate(ry_emm(n_em))
allocate(rz_emm(n_em))
allocate(qw_emm(n_em))
allocate(qx_emm(n_em))
allocate(qy_emm(n_em))
allocate(qz_emm(n_em))
allocate(keep_em_flag(n_em))
read(ncfile,*)
read(ncfile,*) xlo, xhi; boxx1 = xhi - xlo ! read in upper lower bounds, calc box size
read(ncfile,*) ylo, yhi; boxy1 = yhi - ylo ! read in upper lower bounds, calc box size
read(ncfile,*) zlo, zhi; boxz1 = zhi - zlo ! read in upper lower bounds, calc box size
read(ncfile,*) ! can add check here later
do iem = 1,n_em
	read(ncfile,*) jem,jtype,rx_emm(iem),ry_emm(iem),rz_emm(iem), qw_emm(iem),qx_emm(iem),qy_emm(iem),qz_emm(iem)
enddo
close(ncfile)
write(*,*) 'done reading config file'

!! cut out a small, preselected piece of the EM2 config.
!! I use vmd to examine the EM2 config and choose an interesting 
!! piece to reverse map

!! cut out a donut junction defined by all EM2 particles within a radius
!! emm(3302).

keep_em_flag = 0
ikeep = 0
icenter = 3302                    !What is this? M. This is probably redundant.
write(*,*) "Center for the membrane selection: x y z?"
read(*,*) rx_center, ry_center, rz_center             
rx_center=real(rx_center)
ry_center=real(ry_center)
rz_center=real(rz_center)

write(*,*) "Cut-off for the selection? Provide a range from-to: e.g. 0.0 12.0"
read(*,*) cut_1, cut_2

cut_1=real(cut_1)
cut_2=real(cut_2)
do iem = 1,n_em
	!d = sqrt((rx_emm(iem) - rx_emm(icenter))**2 + (ry_emm(iem) - ry_emm(icenter))**2 + (rz_emm(iem) - rz_emm(icenter))**2)
	d = sqrt((rx_emm(iem) - rx_center)**2 + (ry_emm(iem) - ry_center)**2 + (rz_emm(iem) - rz_center)**2)
	!if(d .lt. cut) then
	if((d .lt. cut_2) .and. (d .gt. cut_1)) then
		ikeep = ikeep + 1
		keep_em_flag(iem) = 1
	endif
enddo

nkeep = ikeep
allocate(rx_em(nkeep))
allocate(ry_em(nkeep))
allocate(rz_em(nkeep))
allocate(qw_em(nkeep))
allocate(qx_em(nkeep))
allocate(qy_em(nkeep))
allocate(qz_em(nkeep))

ikeep = 0
do iem = 1,n_em
	if(keep_em_flag(iem) .eq. 1) then
		ikeep = ikeep + 1
		rx_em(ikeep) = rx_emm(iem)
		ry_em(ikeep) = ry_emm(iem)
		rz_em(ikeep) = rz_emm(iem)
		qw_em(ikeep) = qw_emm(iem)
		qx_em(ikeep) = qx_emm(iem)
		qy_em(ikeep) = qy_emm(iem)
		qz_em(ikeep) = qz_emm(iem)
		write(*,*) 'C',  rx_emm(iem),ry_emm(iem),rz_emm(iem) 
	endif
enddo

n_em = nkeep

! now find center of mass of this thing
cmx = 0.d0
cmy = 0.d0
cmz = 0.d0
            
do i=1,n_em   
	cmx = cmx + rx_em(i)
	cmy = cmy + ry_em(i)
	cmz = cmz + rz_em(i)
enddo

cmx = cmx / dble(n_em)
cmy = cmy / dble(n_em)
cmz = cmz / dble(n_em)

open(unit=ncfile, file=output_em2)
Write(ncfile,'(3f14.4)'), boxx1, boxy1, boxz1
do i=1,n_em
	rx_em(i) = rx_em(i) + (boxx1*5.d-1 - cmx)
	ry_em(i) = ry_em(i) + (boxy1*5.d-1 - cmy)
	rz_em(i) = rz_em(i) + (boxz1*5.d-1 - cmz)
	Write(ncfile,'(3f14.4)') rx_em(i),ry_em(i),rz_em(i)	
enddo
do i=1,n_em
	Write(ncfile,'(4f14.4)') qw_em(i),qx_em(i),qy_em(i),qz_em(i)
enddo
close(ncfile)

! now the EM2 patch is stored in rx_em, etc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																					!!
!!		Find triangles made of nearest neighbor EM2 particles		!!
!!																					!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! need to allocate an array to store the vertices...
! not sure how big it will be, let's oversize it:

max_tri = 15000
max_hex_points = 1000
max_gb = max_tri*max_hex_points 
write(*,*) 'size of triangle array: ',max_tri
write(*,*) 'size of total gb array: ',max_gb

allocate(NGB_in_TRI(max_tri))
allocate(TRI_to_GB_map(max_tri,max_gb))
allocate(rx_tri(max_tri))
allocate(ry_tri(max_tri))
allocate(rz_tri(max_tri))
allocate(triangles(max_tri,3))
allocate(EM_distance(n_em,n_em))
allocate(rx_s_one(max_gb),rx_s_two(max_gb),rx_s_three(max_gb))
allocate(ry_s_one(max_gb),ry_s_two(max_gb),ry_s_three(max_gb))
allocate(rz_s_one(max_gb),rz_s_two(max_gb),rz_s_three(max_gb))
allocate(atom_s_one(max_gb),atom_s_two(max_gb),atom_s_three(max_gb))
!allocate(qw_s(max_gb))
!allocate(qx_s(max_gb))
!allocate(qy_s(max_gb))
!allocate(qz_s(max_gb))

! this param decides who are neighbors for the purpose of forming triangles
EM_cutoff = 11.9                	!Originally, 11.9. Larger number ensures that there are no holes
!brute force way...triple loop checking for pair dist below cutoff
do i=1,n_em
	do j=i+1,n_em
		EM_distance(i,j) = sqrt((rx_em(i) - rx_em(j))**2 + (ry_em(i) - ry_em(j))**2 + (rz_em(i) - rz_em(j))**2)
		EM_distance(j,i) = EM_distance(i,j)
	enddo
enddo

do i=1,n_em
   EM_distance(i,i) = -1.0
enddo

allocate(edges(n_em,n_em)) ! this will track how many triangle vertices each EM2 is a part of. 
! six is the max...
edges = 0

! i think we can use edges(i,j) together with the half plane check to figure out whether the next tri should be drawn
! the trick is to find the two k's that go with i and j...well, the other k. we already have one k in hand, need
! to find the other k and check whether they are on the same side of ij. 

!open(unit=19, file='triangles.txt')

i_tri = 0
do i=1,n_em
	do j=i+1,n_em
		do k=j+1,n_em
			if(EM_distance(i,j).lt.EM_cutoff .and. EM_distance(j,k).lt.EM_cutoff .and. EM_distance(k,i).lt.EM_cutoff) then
				if(edges(i,j) .lt. 2 .and. edges(i,k) .lt. 2 .and. edges(j,k) .lt. 2) then
					dot = -1.0 ! this flag will determine whether a new triangle can be drawn. it will flip to positive in the following
									! if statement is the new triangle shares an edge with an already formed triangle and the k and kprime
									! vertices are one the same side of the shared edge
					if(edges(i,j) .eq. 1) then ! i,j already form one edge...find the third vertex from the already 
						do j_tri = 1,i_tri		! formed triangle to check whether the new vertex is in the same half plane
							if(triangles(i_tri,1) .eq. i .and. triangles(i_tri,2) .eq. j) then ! found the first tri with edge i,j
								k_prime = triangles(i_tri,3)
								rx_ij = rx_em(i) - rx_em(j)
								ry_ij = ry_em(i) - ry_em(j)
								rz_ij = rz_em(i) - rz_em(j)
								rx_kj = rx_em(k) - rx_em(j)
								ry_kj = ry_em(k) - ry_em(j)
								rz_kj = rz_em(k) - rz_em(j)
								rx_kpj = rx_em(k_prime) - rx_em(j)
								ry_kpj = ry_em(k_prime) - ry_em(j)
								rz_kpj = rz_em(k_prime) - rz_em(j)

								rx_cross1 = ry_kj*rz_ij - rz_kj*ry_ij
								ry_cross1 = -1.0*(rx_kj*rz_ij - rz_kj*rx_ij)
								rz_cross1 = rx_kj*ry_ij - ry_kj*rx_ij

								rx_cross2 = ry_kpj*rz_ij - rz_kpj*ry_ij
                                ry_cross2 = -1.0*(rx_kpj*rz_ij - rz_kpj*rx_ij)
                                rz_cross2 = rx_kpj*ry_ij - ry_kpj*rx_ij

								dot = rx_cross1*rx_cross2 + ry_cross1*ry_cross2 + rz_cross1*rz_cross2
							endif
						enddo
					endif ! here ends the if statement checking the k and k_prime vertices for triangles with shared edge
					
					if(dot .lt. 0.0) then ! k and k_prime on opposite sides of line ij, or no shared edge       
						i_tri = i_tri+1
						edges(i,j) = edges(i,j) + 1
						edges(i,k) = edges(i,k) + 1
						edges(j,k) = edges(j,k) + 1
						triangles(i_tri,1) = i
						triangles(i_tri,2) = j
						triangles(i_tri,3) = k
						write(*,*), i, j, k
						!write(19,*) i, j, k
					endif
				
				endif
			endif
		enddo
	enddo
enddo

!Now adding additional triangles from standard input. According to Euler, you cannot mesh a sphere
!with triangles without defects - there should be 12 I believe

write(*,*) "How many additional triangles?"
read(*,*) additional_tri

do m = 1,additional_tri
   i_tri = i_tri + 1
   read(*,*) i,j,k
   triangles(i_tri,1) = i
   triangles(i_tri,2) = j
   triangles(i_tri,3) = k
   write(*,*) "Add. triangles", i, j, k
   write(*,*) i, j, k
enddo

n_tri = i_tri
write(*,*) 'counted ', n_tri, ' triangles' 
!do i=1,n_em
!	do j=1,n_em
!		write(*,*) 'edge', i,j,edges(i,j)
!	enddo
!enddo

!close(19)

allocate(gb_patch_flag(ngb))

! now here comes the action...loop over the list of triangles. in the lab frame pick out the particles in the 
! prequilibrated GB patch that are inside the triangle. then transform back to the local frame of the triangle

!open(unit=ncfile,file='temp.xyz')
!write(ncfile,*)
!write(ncfile,*)
! is will index the final GB particles
is = 0
do i=1,n_tri
   ! get the vertices
   write(*,*) 'looping over triangles. on tri ', i, triangles(i,1), triangles(i,2), triangles(i,3)
   i_em = triangles(i,1)
   j_em = triangles(i,2)
   k_em = triangles(i,3)
   ! the center of the triangle
   rx_tri(i) = (rx_em(i_em) + rx_em(j_em) + rx_em(k_em))/3.0  
   ry_tri(i) = (ry_em(i_em) + ry_em(j_em) + ry_em(k_em))/3.0  
   rz_tri(i) = (rz_em(i_em) + rz_em(j_em) + rz_em(k_em))/3.0  
   
   !This part fixes the misorientation of the bilayer - as a consequence of the chirality of the
   !triangles (vertices are counted sometimes clock-wise, sometimes counter-clockwise), the patched
   !bilayer sometimes faces upward, and sometimes downward. As a fix, I calculate a cross product of the
   !vectors formed by vertices (12 cross 23), and if the resulting vector is (closely) parallel to the
   !norm of the patch (determined from the average quaternion of the three vertices) then the configuration
   !is R, otherwise, it is S. (Jan 2012, M.S.)
      
   !vertices of the triangle in local frame
   r1_prime_x = rx_em(i_em)
   r1_prime_y = ry_em(i_em)
   r1_prime_z = rz_em(i_em)
   r2_prime_x = rx_em(j_em)
   r2_prime_y = ry_em(j_em)
   r2_prime_z = rz_em(j_em)
   r3_prime_x = rx_em(k_em)
   r3_prime_y = ry_em(k_em)
   r3_prime_z = rz_em(k_em)

   r13_x = r3_prime_x - r1_prime_x
   r13_y = r3_prime_y - r1_prime_y
   r13_z = r3_prime_z - r1_prime_z
   r12_x = r2_prime_x - r1_prime_x
   r12_y = r2_prime_y - r1_prime_y
   r12_z = r2_prime_z - r1_prime_z

   !the following four lines were added Jan 2012, M.S.
   r23_x = r3_prime_x - r2_prime_x
   r23_y = r3_prime_y - r2_prime_y
   r23_z = r3_prime_z - r2_prime_z
   r23_mod = sqrt(r23_x**2 + r23_y**2 + r23_z**2)

   r13_mod = sqrt(r13_x**2 + r13_y**2 + r13_z**2)
   r12_mod = sqrt(r12_x**2 + r12_y**2 + r12_z**2)
   cos_theta = (r13_x*r12_x + r13_y*r12_y + r13_z*r12_z)/(r13_mod*r12_mod)
   theta = acos(cos_theta)
   
   !the following added Jan 2012, M.S. This calculates the cross product.
   cross_x = r12_y*r23_z - r12_z*r23_y
   cross_y = r12_z*r23_x - r12_x*r23_z
   cross_z = r12_x*r23_y - r12_y*r23_x
   cross_mod = sqrt(cross_x**2 + cross_y**2 + cross_z**2)
   
   qw_av = qw_em(i_em) !(qw_em(i_em) + qw_em(j_em) + qw_em(k_em))/3.0d0 !
   qx_av = qx_em(i_em) !(qx_em(i_em) + qx_em(j_em) + qx_em(k_em))/3.0d0 !
   qy_av = qy_em(i_em) !(qy_em(i_em) + qy_em(j_em) + qy_em(k_em))/3.0d0 !
   qz_av = qz_em(i_em) !(qz_em(i_em) + qz_em(j_em) + qz_em(k_em))/3.0d0 !
   !qq = sqrt(qw_av**2 + qx_av**2 + qy_av**2 + qz_av**2)
   !qw_av = qw_av/qq
   !qx_av = qx_av/qq
   !qy_av = qy_av/qq
   !qz_av = qz_av/qq
   
    axx = qw_av**2 + qx_av**2 - qy_av**2 - qz_av**2
    axy = 2.d0 * ( qx_av * qy_av - qw_av * qz_av )
    axz = 2.d0 * ( qx_av * qz_av + qw_av * qy_av )
	ayx = 2.d0 * ( qx_av * qy_av + qw_av * qz_av )
	ayy = qw_av**2 - qx_av**2 + qy_av**2 -qz_av**2
	ayz = 2.d0 * ( qy_av * qz_av - qw_av * qx_av )
	azx = 2.d0 * ( qx_av * qz_av - qw_av * qy_av )
	azy = 2.d0 * ( qy_av * qz_av + qw_av * qx_av )
	azz = qw_av**2 - qx_av**2 - qy_av**2 + qz_av**2

	sbudixi = 0.d0
	sbudiyi = 0.d0
	sbudizi = 1.d0
	normx = ( axx * sbudixi + axy * sbudiyi + axz * sbudizi )
	normy = ( ayx * sbudixi + ayy * sbudiyi + ayz * sbudizi )
	normz = ( azx * sbudixi + azy * sbudiyi + azz * sbudizi )
   
   !normx = 0.0d0
   !normy = 0.0d0
   !normz = 1.0d0
   norm_mod = sqrt(normx**2+normy**2+normz**2)
   phi = (normx*cross_x + normy*cross_y + normz*cross_z)/(cross_mod*norm_mod)
   phi = ACOS(phi)   

   !vertices in the lab frame   !This is local, not lab - the one above is lab (global)
   r1_x = 0.0
   r1_y = 0.0
   r1_z = 0.0
   r2_x = r12_mod*cos_theta
   r2_y = r12_mod*sin(theta)
   r2_z = 0.0
   r3_x = r13_mod
   r3_y = 0.0
   r3_z = 0.0

   !initialize the array which track hexagonal gb particles that are inside the triangle
   gb_patch_flag = 0
    
	! checking which GB particles are inside the triangle formed by three EM2 particles
   slope1 = r2_y/r2_x
   slope2 = (r3_y - r2_y)/(r3_x - r2_x)
   b2 = r2_y - slope2*r2_x
	do igb=1,ngb
		ymax1 = slope1*rx_one(igb) 
		ymax2 = slope2*rx_one(igb) + b2 
		x = rx_one(igb)
		y = ry_one(igb)
		if(x.gt.r1_x .and. x.lt.r2_x .and. y.gt.r1_y .and. y.lt.ymax1) then
			gb_patch_flag(igb) = 1
		elseif(x.gt.r2_x .and. x.lt.r3_x .and. y.gt.r1_y .and. y.lt.ymax2) then
			gb_patch_flag(igb) = 1
 		endif
	enddo

! now define the local coord system that will determine the trans from lab to local frame
   x_prime_x = r13_x/r13_mod
   x_prime_y = r13_y/r13_mod
   x_prime_z = r13_z/r13_mod
   z_prime_x =       r13_y*r12_z - r13_z*r12_y
   z_prime_y = -1.0*(r13_x*r12_z - r13_z*r12_x)
   z_prime_z =       r13_x*r12_y - r13_y*r12_x
   zprime_mod = sqrt(z_prime_x**2 + z_prime_y**2 + z_prime_z**2)
   z_prime_x = z_prime_x/zprime_mod
   z_prime_y = z_prime_y/zprime_mod
   z_prime_z = z_prime_z/zprime_mod
   y_prime_x =       z_prime_y*x_prime_z - z_prime_z*x_prime_y
   y_prime_y = -1.0*(z_prime_x*x_prime_z - z_prime_z*x_prime_x)
   y_prime_z =       z_prime_x*x_prime_y - z_prime_y*x_prime_x

   trace = 1.0 + x_prime_x + y_prime_y + z_prime_z
   if(trace .le. 0.0d0) then
     write(*,*) 'trace less than zero...R -> Q trans unstable', trace
   else
      write(*,*) "trace:", trace
   endif

	!! the rotation to go from the lab frame to the frame of the EM2 triangle/local frame
	qw_t = 0.5*sqrt(trace)
    qx_t = (y_prime_z - z_prime_y)/(4.0*qw_t)
    qy_t = (z_prime_x - x_prime_z)/(4.0*qw_t)
    qz_t = (x_prime_y - y_prime_x)/(4.0*qw_t)

	qq = qw_t**2 + qx_t**2 + qy_t**2 + qz_t**2
	qw_t = qw_t/qq 
	qx_t = qx_t/qq 
	qy_t = qy_t/qq 
	qz_t = qz_t/qq 
 

! transform each gb site that is inside the triangle back to the local frame
   n_triGB = 0 ! this will count GB particles in the current triangle
   do igb=1,ngb
    	if (gb_patch_flag(igb) .eq. 1) then

			is = is + 1
			! rotate the coords of the GB particle 
			rx_s_one(is) = (rx_one(igb)*x_prime_x + ry_one(igb)*y_prime_x + rz_one(igb)*z_prime_x + r1_prime_x)
			rx_s_two(is) = (rx_two(igb)*x_prime_x + ry_two(igb)*y_prime_x + rz_two(igb)*z_prime_x + r1_prime_x)
			rx_s_three(is) = (rx_three(igb)*x_prime_x + ry_three(igb)*y_prime_x + rz_three(igb)*z_prime_x + r1_prime_x)
			ry_s_one(is) = (rx_one(igb)*x_prime_y + ry_one(igb)*y_prime_y + rz_one(igb)*z_prime_y + r1_prime_y)
			ry_s_two(is) = (rx_two(igb)*x_prime_y + ry_two(igb)*y_prime_y + rz_two(igb)*z_prime_y + r1_prime_y)
			ry_s_three(is) = (rx_three(igb)*x_prime_y + ry_three(igb)*y_prime_y + rz_three(igb)*z_prime_y + r1_prime_y)
			rz_s_one(is) = (rx_one(igb)*x_prime_z + ry_one(igb)*y_prime_z + rz_one(igb)*z_prime_z + r1_prime_z)
			rz_s_two(is) = (rx_two(igb)*x_prime_z + ry_two(igb)*y_prime_z + rz_two(igb)*z_prime_z + r1_prime_z)
			rz_s_three(is) = (rx_three(igb)*x_prime_z + ry_three(igb)*y_prime_z + rz_three(igb)*z_prime_z + r1_prime_z)
			atom_s_one(is) = atom_one(igb)
			atom_s_two(is) = atom_two(igb)
			atom_s_three(is) = atom_three(igb)
			!write(*,*) is, rx(igb),ry(igb),rz(igb), rx_s(is), ry_s(is), rz_s(is) 
			! rotate the rotation to get the q's for the GB particle in the new, rotated frame
			!qw_s(is) = qw_t*qw(igb) - qx_t*qx(igb) - qy_t*qy(igb) - qz_t*qz(igb) 
			!qx_s(is) = qw_t*qx(igb) + qx_t*qw(igb) + qy_t*qz(igb) - qz_t*qy(igb)
			!qy_s(is) = qw_t*qy(igb) - qx_t*qz(igb) + qy_t*qw(igb) + qz_t*qx(igb)
			!qz_s(is) = qw_t*qz(igb) + qx_t*qy(igb) - qy_t*qx(igb) + qz_t*qw(igb)
			!write(*,*) is, qw(igb), qx(igb), qy(igb), qz(igb)
			! these next two lines track to which triangle the GB particle belongs. 
			n_triGB = n_triGB + 1
			TRI_to_GB_map(i,n_triGB) = is 
    
            if (phi > 3.14/2) then
               if (atom_s_one(is) == 1) then
                  atom_s_one(is) = 4
               elseif (atom_s_one(is) == 10) then
                  atom_s_one(is) = 40
               elseif (atom_s_one(is) == 4) then
                  atom_s_one(is) = 1
               elseif (atom_s_one(is) == 40) then
                  atom_s_one(is) = 10
               endif
               if (atom_s_two(is) == 2) then
                  atom_s_two(is) = 5
               elseif (atom_s_two(is) == 5) then
                  atom_s_two(is) = 2
               endif
               if (atom_s_three(is) == 3) then
                  atom_s_three(is) = 6
               elseif (atom_s_three(is) == 6) then
                  atom_s_three(is) = 3
               endif
            endif               
            !write(ncfile,*) atom_s_one(is),rx_s_one(is),ry_s_one(is),rz_s_one(is)
            !write(ncfile,*) atom_s_two(is),rx_s_two(is),ry_s_two(is),rz_s_two(is)
            !write(ncfile,*) atom_s_three(is),rx_s_three(is),ry_s_three(is),rz_s_three(is)   
        endif
   enddo
   NGB_in_TRI(i) = n_triGB ! this array stores the number of GB particles in each triangle

enddo ! end loop over triangles
close(ncfile)


! total GB particles, currently
n_s = is

! now use TRI_to_GB_map to search GB particles belonging to neighboring triangles for overlaps
write(*,*) 'starting overlaps search...'
allocate(print_flag(n_s))
print_flag = 1
TRI_cut_sqr = 120.0 ! notice the cutoffs are the squared distance to eliminate the sqrt
upper_GB_cut = 12.0	! originally 12.0


! stuff for the orientation dependent overlap check:
aspect_ratio = 3.0
chi = (aspect_ratio**2 - 1.0)/(aspect_ratio**2 + 1.0)
sigma_0 = 0.76 ! GB core radius, originally 0.76

do i_tri=1,n_tri
  write(*,*) 'searching neighbors of triangle ',i_tri, '  for overlaps'
  do j_tri=i_tri+1,n_tri
     if(j_tri .ne. i_tri) then
     d_tri = (rx_tri(i_tri) - rx_tri(j_tri))**2 + (ry_tri(i_tri) - ry_tri(j_tri))**2 + (rz_tri(i_tri) - rz_tri(j_tri))**2  
     !write(*,*) 'sqr dist between tri ', i_tri, ' and ',j_tri, ' : ',d_tri
     !if(d_tri .lt. TRI_cut_sqr) then
         n_i = NGB_in_TRI(i_tri)
         n_j = NGB_in_TRI(j_tri)
         do i_gb=1,n_i
            do j_gb=1,n_j
               is_i = TRI_to_GB_map(i_tri,i_gb)
               is_j = TRI_to_GB_map(j_tri,j_gb)
               d_gb = sqrt((rx_s_one(is_i) - rx_s_one(is_j))**2 + (ry_s_one(is_i) - ry_s_one(is_j))**2 + &
               		  (rz_s_one(is_i) - rz_s_one(is_j))**2)
               if(d_gb .lt. upper_GB_cut) then
                 ! need to compute orientation dependent sigma in order to check properly for overlaps...
                 ! first for particle is_i...
                 !axx = qw_s(is_i)**2 + qx_s(is_i)**2 - qy_s(is_i)**2 - qz_s(is_i)**2
                 !axy = 2.d0 * ( qx_s(is_i) * qy_s(is_i) - qw_s(is_i) * qz_s(is_i) )
                 !axz = 2.d0 * ( qx_s(is_i) * qz_s(is_i) + qw_s(is_i) * qy_s(is_i) )
                 !ayx = 2.d0 * ( qx_s(is_i) * qy_s(is_i) + qw_s(is_i) * qz_s(is_i) )
                 !ayy = qw_s(is_i) **2 - qx_s(is_i)**2 + qy_s(is_i)**2 -qz_s(is_i)**2
                 !ayz = 2.d0 * ( qy_s(is_i) * qz_s(is_i) - qw_s(is_i) * qx_s(is_i) )
                 !azx = 2.d0 * ( qx_s(is_i) * qz_s(is_i) - qw_s(is_i) * qy_s(is_i) )
                 !azy = 2.d0 * ( qy_s(is_i) * qz_s(is_i) + qw_s(is_i) * qx_s(is_i) )
                 !azz = qw_s(is_i)**2 - qx_s(is_i)**2 - qy_s(is_i)**2 + qz_s(is_i)**2
                 !sbudixi = 0.d0
                 !sbudiyi = 0.d0
                 !sbudizi = 1.d0
                 !sorntxi = ( axx * sbudixi + axy * sbudiyi + axz * sbudizi )
                 !sorntyi = ( ayx * sbudixi + ayy * sbudiyi + ayz * sbudizi )
                 !sorntzi = ( azx * sbudixi + azy * sbudiyi + azz * sbudizi )
                 
                 ! now for particle is_j...
                 !axx = qw_s(is_j)**2 + qx_s(is_j)**2 - qy_s(is_j)**2 - qz_s(is_j)**2
                 !axy = 2.d0 * ( qx_s(is_j) * qy_s(is_j) - qw_s(is_j) * qz_s(is_j) )
                 !axz = 2.d0 * ( qx_s(is_j) * qz_s(is_j) + qw_s(is_j) * qy_s(is_j) )
                 !ayx = 2.d0 * ( qx_s(is_j) * qy_s(is_j) + qw_s(is_j) * qz_s(is_j) )
                 !ayy = qw_s(is_j) **2 - qx_s(is_j)**2 + qy_s(is_j)**2 -qz_s(is_j)**2
                 !ayz = 2.d0 * ( qy_s(is_j) * qz_s(is_j) - qw_s(is_j) * qx_s(is_j) )
                 !azx = 2.d0 * ( qx_s(is_j) * qz_s(is_j) - qw_s(is_j) * qy_s(is_j) )
                 !azy = 2.d0 * ( qy_s(is_j) * qz_s(is_j) + qw_s(is_j) * qx_s(is_j) )
                 !azz = qw_s(is_j)**2 - qx_s(is_j)**2 - qy_s(is_j)**2 + qz_s(is_j)**2
                 !sbudixj = 0.d0
                 !sbudiyj = 0.d0
                 !sbudizj = 1.d0
                 !sorntxj = ( axx * sbudixi + axy * sbudiyi + axz * sbudizi )
                 !sorntyj = ( ayx * sbudixi + ayy * sbudiyi + ayz * sbudizi )
                 !sorntzj = ( azx * sbudixi + azy * sbudiyi + azz * sbudizi )
                 
                 ! interparticle unit vector:
                 !s_x = (rx_s(is_i) - rx_s(is_j))/d_gb
                 !s_y = (ry_s(is_i) - ry_s(is_j))/d_gb
                 !s_z = (rz_s(is_i) - rz_s(is_j))/d_gb
                 ! orientation vectors dot prod
                 !sbud_dot = (sorntxi*sorntxj + sorntyi*sorntyj + sorntzi*sorntzj)
                 !tmp_1 = s_x*sorntxi + s_y*sorntyi + s_z*sorntzi + s_x*sorntxj + s_y*sorntyj + s_z*sorntzj
                 !tmp_2 = s_x*sorntxi + s_y*sorntyi + s_z*sorntzi - s_x*sorntxj - s_y*sorntyj - s_z*sorntzj
                 !term_1 = tmp_1/(1.0 + chi*sbud_dot)
                 !term_2 = tmp_2/(1.0 - chi*sbud_dot)
                 !sigma = sigma_0/(sqrt(1.0 - 0.5*chi*(term_1 + term_2)))
                 cut = 0.50*sigma_0                 
                 if ((d_gb .lt. cut) .and. (print_flag(is_j) .ne. 0)) then
                     !write(*,*) 'eliminating GB particle ',is_i,' due to overlap'
                     print_flag(is_i) = 0
                 endif
               endif
            enddo
         enddo
     endif
     !endif
  enddo
enddo

ixyz = 0
! write TANTALUS config
open(unit=ncfile,file=output_xyz)
	! just picking a huge boxsize for now...
   !write(ncfile,'(4f15.5)') 500.0, 500.0, 500.0, 0.d0
   write(ncfile,*) 3*n_s
   write(ncfile,*)
   do igb = 1,n_s
     if(print_flag(igb) .eq. 1)then
		ixyz = ixyz + 1
        write(ncfile,'(I2,3f15.5)') atom_s_one(igb), rx_s_one(igb),ry_s_one(igb),rz_s_one(igb)
        write(ncfile,'(I2,3f15.5)') atom_s_two(igb), rx_s_two(igb),ry_s_two(igb),rz_s_two(igb)
        write(ncfile,'(I2,3f15.5)') atom_s_three(igb), rx_s_three(igb),ry_s_three(igb),rz_s_three(igb)
     endif
   enddo
   !write(ncfile,'(3f15.5,i15,f15.5)') 0.d0,0.d0,0.d0,10,0.d0
   !do igb = 1,n_s
   !   if(print_flag(igb) .eq. 1)then
   !      write(ncfile,'(4f15.5)') qw_s(igb),qx_s(igb),qy_s(igb),qz_s(igb)
   !   endif
   !enddo
   !write(ncfile,'(4f15.5)') qw_s(1),qx_s(1),qy_s(1),qz_s(1)
close(ncfile)



!open(unit=ncfile,file='config-gb-em.xyz')
!write(ncfile,'(i8)') ixyz+n_em
!write(ncfile,'("          ")')

!do i=1,n_em
!  tmp_x = 10.0*rx_em(i)
!  tmp_y = 10.0*ry_em(i)
!  tmp_z = 10.0*rz_em(i)
!   write(ncfile,'(a10,3f14.4)') 'C', tmp_x,tmp_y,tmp_z
!enddo

do i=1,n_s
 if(print_flag(i) .eq. 10)then
     	axx = qw_s(i)**2 + qx_s(i)**2 - qy_s(i)**2 - qz_s(i)**2
     	axy = 2.d0 * ( qx_s(i) * qy_s(i) - qw_s(i) * qz_s(i) )
     	axz = 2.d0 * ( qx_s(i) * qz_s(i) + qw_s(i) * qy_s(i) )
     	ayx = 2.d0 * ( qx_s(i) * qy_s(i) + qw_s(i) * qz_s(i) )
     	ayy = qw_s(i) **2 - qx_s(i)**2 + qy_s(i)**2 -qz_s(i)**2
     	ayz = 2.d0 * ( qy_s(i) * qz_s(i) - qw_s(i) * qx_s(i) )
     	azx = 2.d0 * ( qx_s(i) * qz_s(i) - qw_s(i) * qy_s(i) )
     	azy = 2.d0 * ( qy_s(i) * qz_s(i) + qw_s(i) * qx_s(i) )
     	azz = qw_s(i)**2 - qx_s(i)**2 - qy_s(i)**2 + qz_s(i)**2

		sbudixi = 0.d0
		sbudiyi = 0.d0
        sbudizi = 1.d0
     	sorntxi = ( axx * sbudixi + axy * sbudiyi + axz * sbudizi )
     	sorntyi = ( ayx * sbudixi + ayy * sbudiyi + ayz * sbudizi )
     	sorntzi = ( azx * sbudixi + azy * sbudiyi + azz * sbudizi )
     	tmp_x = 10.0*rx_s(i)
     	tmp_y = 10.0*ry_s(i)
     	tmp_z = 10.0*rz_s(i)
     	!Write(ncfile,'(a10,3f14.4)') 'N', (tmp_x+8.d-1*sorntxi),(tmp_y+8.d-1*sorntyi),(tmp_z+8.d-1*sorntzi)
     	!Write(ncfile,'(a10,3f14.4)') 'P', (tmp_x-8.d-1*sorntxi),(tmp_y-8.d-1*sorntyi),(tmp_z-8.d-1*sorntzi)
     	!Write(ncfile,'(a10,3f14.4)') 'N', tmp_x,tmp_y,tmp_z

   endif
enddo
close(ncfile)

!end program reverse_map_em_to_GB


end program recreateCG_membrane
