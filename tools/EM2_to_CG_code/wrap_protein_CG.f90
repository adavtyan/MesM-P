Program wrap_prot_CG

	implicit none
	
	Integer, Parameter :: DP = Selected_real_kind(15,307)
	integer :: i,j,N,NBAR,iem,ied,n_total,n_edbar,Ncut
	integer :: is,itmp,jed,jem,js,ked,ks,ls
	integer :: junk,shall_I_rotate
	Integer, Dimension(:,:), Allocatable				  :: bar_map	
	
	Real(kind=DP), Dimension(:) , Allocatable :: rx_mem1,ry_mem1,rz_mem1,rx_mem2,ry_mem2,rz_mem2
	Real(kind=DP), Dimension(:) , Allocatable :: rx_mem3,ry_mem3,rz_mem3,qw_mem,qx_mem,qy_mem,qz_mem
	Real(kind=DP) :: xnorm,ynorm,znorm,xcom,ycom,zcom
	Real(kind=DP) :: r,x0,y0,z0,crossx,crossy,crossz,qq
	real(kind=DP) :: axx,axy,axz,ayx,ayy,ayz,azx,azy,azz
	real :: r_xy, rcut_enh, rcut_reg, rij_x, rij_y, rij_z,rx,ry,rz
	real :: c_pi,boxx,boxy,boxz
	real :: d,d_cut,delta,delta_d
	real :: dip_x,dip_y,dip_z
	real :: e_cut, eps_enh, eps_reg, sig_enh, sig_reg
	real :: phi,theta
	real :: pow, sbudixi, sbudiyi, sbudizi, sorntxi,sorntyi,sorntzi
	real :: tmpx,tmpy,tmpz,tmp1,tmp2
	real :: cutoff,dist
	real :: rotate,comx,comy,comz
	
	
	Real(kind=DP), Dimension(:) , Allocatable :: rx_EDBAR, ry_EDBAR, rz_EDBAR
	Real(kind=DP), Dimension(:) , Allocatable :: qw_EDBAR, qx_EDBAR, qy_EDBAR, qz_EDBAR
	Real(kind=DP), Dimension(:) , Allocatable :: rx_s, ry_s, rz_s, qw_s, qx_s, qy_s, qz_s
	
    character :: atom
    Character, Dimension(:) , Allocatable             :: atom_name,atom_name_s

	character(len=41)								  :: input,output
	real											  :: rot

	c_pi = 3.1415927
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! read in EM2 config !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) "Distance from BAR?"
    read(*,*) delta
    write(*,*) "Output xyz BAR file?"
    read(*,*) output
    
    
    write(*,*) "Box dimensions? x y z"
    read(*,*) boxx,boxy,boxz
    
	open(unit=12, file="Selection.xyz")
	read(12,*) N
	read(12,*)
	
	NBAR = N
	
	allocate(rx_mem1(N),ry_mem1(N),rz_mem1(N))
	!allocate(rx_mem2(N),ry_mem2(N),rz_mem2(N))
	!allocate(rx_mem3(N),ry_mem3(N),rz_mem3(N))
	allocate(qw_mem(N),qx_mem(N),qy_mem(N),qz_mem(N))
	
    xcom=0.000000
    ycom=0.000000
    zcom=0.000000	
	do i = 1,N
	   read(12,*) atom,rx_mem1(i),ry_mem1(i),rz_mem1(i)
	   !read(12,*) atom,rx_mem2(i),ry_mem2(i),rz_mem2(i)
	   !read(12,*) atom,rx_mem3(i),ry_mem3(i),rz_mem3(i)
       xcom=xcom+rx_mem1(i)
       ycom=ycom+ry_mem1(i)
       zcom=zcom+rz_mem1(i)
	enddo
    xcom=xcom/N
    ycom=ycom/N
    zcom=zcom/N
    
    rx_mem1 = rx_mem1/10.0
    ry_mem1 = ry_mem1/10.0
    rz_mem1 = rz_mem1/10.0
	!rx_mem2 = rx_mem2/10.0
    !ry_mem2 = ry_mem2/10.0
    !rz_mem2 = rz_mem2/10.0
    !rx_mem3 = rx_mem3/10.0
    !ry_mem3 = ry_mem3/10.0
    !rz_mem3 = rz_mem3/10.0
	
	xcom = xcom/10.0
	ycom = ycom/10.0
	zcom = zcom/10.0
	
	open(unit=15,file="Mem_normals_quaternions.txt")
	
    do i=1,N
       read(15,*) qw_mem(i),qx_mem(i),qy_mem(i),qz_mem(i)
    enddo
	
	close(12)	
	
	x0=0.00000
    y0=0.00000
    z0=1.00000

	! lets move the CoM of the membrane coords to the center of the box
	! this way the config lines up with the mapped config

   !do i=1,NBAR
   !   rx_mem1(i) = rx_mem1(i) + (boxx*5.d-1 - xcom)
   !   ry_mem1(i) = ry_mem1(i) + (boxy*5.d-1 - ycom)
   !	  rz_mem1(i) = rz_mem1(i) + (boxz*5.d-1 - zcom)
      !rx_mem2(i) = rx_mem2(i) + (boxx*5.d-1 - xcom)
      !ry_mem2(i) = ry_mem2(i) + (boxy*5.d-1 - ycom)
   	  !rz_mem2(i) = rz_mem2(i) + (boxz*5.d-1 - zcom)   	  
      !rx_mem3(i) = rx_mem3(i) + (boxx*5.d-1 - xcom)
      !ry_mem3(i) = ry_mem3(i) + (boxy*5.d-1 - ycom)
   	  !rz_mem3(i) = rz_mem3(i) + (boxz*5.d-1 - zcom)   	  
   !enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  read in coords of EDBAR and orientations of dipoles, 			 !!
!!  in Ang and aligned to xyz in lab frame 								 !!
!!  (long axis of BAR along y, membrane normal along z) 				 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!n_edbar = 26

   open(unit=13, file='BARforTube.xyz')
   read(13,*) n_edbar
   read(13,*)
   
   allocate(rx_EDBAR(n_edbar))
   allocate(ry_EDBAR(n_edbar))
   allocate(rz_EDBAR(n_edbar))
   allocate(atom_name(n_edbar))
   allocate(qw_EDBAR(n_edbar))
   allocate(qx_EDBAR(n_edbar))
   allocate(qy_EDBAR(n_edbar))
   allocate(qz_EDBAR(n_edbar))

   do ied = 1,n_edbar
      read(13,*) atom_name(ied), rx_EDBAR(ied),ry_EDBAR(ied),rz_EDBAR(ied),junk,tmp1,tmp2,dip_x,dip_y,dip_z
		!rx_EDBAR(ied) = rx_EDBAR(ied)/10.0; ! BARforTube.xyz is in Angstroms
		!ry_EDBAR(ied) = ry_EDBAR(ied)/10.0;
		!rz_EDBAR(ied) = rz_EDBAR(ied)/10.0;
		! convert the dipole axis into a quaternion representation. note there is an undetermined DoF since
		! we only care which direction the dipole points
	  r = sqrt(dip_x**2 + dip_y**2 + dip_z**2)
	  r_xy = sqrt(dip_x**2 + dip_y**2)
	  theta = acos(dip_z/r)
	  phi = asin(dip_y/r_xy)
	  qw_EDBAR(ied) = cos(0.5*theta)*cos(0.5*(phi + c_pi)) 
	  qx_EDBAR(ied) = sin(0.5*theta)*cos(0.5*(phi + c_pi))
	  qy_EDBAR(ied) = sin(0.5*theta)*sin(0.5*(phi + c_pi))
	  qz_EDBAR(ied) = cos(0.5*theta)*sin(0.5*(phi + c_pi)) 
   enddo
   close(13)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  place the EDBARs...aligned to bottom-top vector (approximate normal)  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	n_total = n_edbar*NBAR
	allocate(rx_s(n_total))
	allocate(ry_s(n_total))
	allocate(rz_s(n_total))
	allocate(qw_s(n_total))
	allocate(qx_s(n_total))
	allocate(qy_s(n_total))
	allocate(qz_s(n_total))
	allocate(atom_name_s(n_total))
	allocate(bar_map(NBAR,n_edbar))

	! might as well write the xyz file while calculating the rotations
	open(unit=16, file='EM-BARs.xyz')
	write(16,'(i8)') n_total
	write(16,'("          ")')



!!!!!!!!!!!!!!!!!!!!!!!!!!
	is = 0
	do iem = 1,NBAR
	    	
	
		axx = qw_mem(iem)**2 + qx_mem(iem)**2 - qy_mem(iem)**2 - qz_mem(iem)**2
		axy = 2.d0 * ( qx_mem(iem) * qy_mem(iem) - qw_mem(iem) * qz_mem(iem) )
		axz = 2.d0 * ( qx_mem(iem) * qz_mem(iem) + qw_mem(iem) * qy_mem(iem) )
		ayx = 2.d0 * ( qx_mem(iem) * qy_mem(iem) + qw_mem(iem) * qz_mem(iem) )
		ayy = qw_mem(iem) **2 - qx_mem(iem)**2 + qy_mem(iem)**2 -qz_mem(iem)**2
		ayz = 2.d0 * ( qy_mem(iem) * qz_mem(iem) - qw_mem(iem) * qx_mem(iem) )
		azx = 2.d0 * ( qx_mem(iem) * qz_mem(iem) - qw_mem(iem) * qy_mem(iem) )
		azy = 2.d0 * ( qy_mem(iem) * qz_mem(iem) + qw_mem(iem) * qx_mem(iem) )
		azz = qw_mem(iem)**2 - qx_mem(iem)**2 - qy_mem(iem)**2 + qz_mem(iem)**2

	!! an adjustable shift along the direction normal to the membrane...this can be used to 
	!! shift the BARs in or out to get a config that does not clash horribly with the membrane
	!! delta controls the amount of the shift
	!delta = 4.5				!originally 3.0
	sbudixi = 0.d0
	sbudiyi = 0.d0
	sbudizi = 1.d0
	sorntxi = ( axx * sbudixi + axy * sbudiyi + axz * sbudizi )
	sorntyi = ( ayx * sbudixi + ayy * sbudiyi + ayz * sbudizi )
	sorntzi = ( azx * sbudixi + azy * sbudiyi + azz * sbudizi )
	
	
		do ied=1,n_edbar
			rx = rx_mem1(iem) + axx*rx_EDBAR(ied) + axy*ry_EDBAR(ied) + axz*rz_EDBAR(ied) + delta*sorntxi
			ry = ry_mem1(iem) + ayx*rx_EDBAR(ied) + ayy*ry_EDBAR(ied) + ayz*rz_EDBAR(ied) + delta*sorntyi
			rz = rz_mem1(iem) + azx*rx_EDBAR(ied) + azy*ry_EDBAR(ied) + azz*rz_EDBAR(ied) + delta*sorntzi
			atom = atom_name(ied)
			tmpx = rx*10.0
			tmpy = ry*10.0
			tmpz = rz*10.0
			write(16,'(a10,3f14.4)') atom_name(ied), tmpx, tmpy, tmpz
			! store the coordinates of the BAR sites for the config file
			is = is+1
			bar_map(iem,ied) = is
			rx_s(is) = rx
			ry_s(is) = ry
			rz_s(is) = rz
			atom_name_s(is) = atom
			! now compute the quaternions for the dipole orientation. note that 
			! there is an unspecified DoF since we only care about the direction 
			! of the dipole and not the rotation about the dipole axis
			qw_s(is) = qw_mem(iem)*qw_EDBAR(ied) - qx_mem(iem)*qx_EDBAR(ied) - qy_mem(iem)*qy_EDBAR(ied) - qz_mem(iem)*qz_EDBAR(ied)
			qx_s(is) = qw_mem(iem)*qx_EDBAR(ied) + qx_mem(iem)*qw_EDBAR(ied) + qy_mem(iem)*qz_EDBAR(ied) - qz_mem(iem)*qy_EDBAR(ied)
			qy_s(is) = qw_mem(iem)*qy_EDBAR(ied) - qx_mem(iem)*qz_EDBAR(ied) + qy_mem(iem)*qw_EDBAR(ied) + qz_mem(iem)*qx_EDBAR(ied)
			qz_s(is) = qw_mem(iem)*qz_EDBAR(ied) + qx_mem(iem)*qy_EDBAR(ied) - qy_mem(iem)*qx_EDBAR(ied) + qz_mem(iem)*qw_EDBAR(ied)
			qq = qw_s(is)**2 + qx_s(is)**2 + qy_s(is)**2 + qz_s(is)**2
			qw_s(is) = qw_s(is)/qq
			qx_s(is) = qx_s(is)/qq
			qy_s(is) = qy_s(is)/qq
			qz_s(is) = qz_s(is)/qq
		enddo
	enddo
	close(16)

	!! precalculate the cutoff in r based on a set energy cutoff. 
	e_cut = 10000.0
	eps_enh = 45.0
	eps_reg = 5.0
	sig_enh = 1.0
	sig_reg = 1.5
	pow = -1.0/6.0
	rcut_enh = sig_enh*((sqrt(e_cut/(4.0*eps_enh) + 0.25) + 0.5)**(pow)) ! rcut ~= 0.591 for e_cut = 10^5
	rcut_reg = sig_reg*((sqrt(e_cut/(4.0*eps_reg) + 0.25) + 0.5)**(pow)) ! rcut ~= 0.738 for e_cut = 10^5

	
	!! BAR-BAR overlap check
	do iem=1,NBAR
		do jem=iem+1,NBAR
			!if(iem .ne. jem)then
				do ied=1,n_edbar
					do jed=1,n_edbar
						is = bar_map(iem,ied)
						js = bar_map(jem,jed)
						d = sqrt((rx_s(is) - rx_s(js))**2 + (ry_s(is) - ry_s(js))**2 + (rz_s(is) - rz_s(js))**2)
						if(((ied .eq. 1) .or. (ied .eq. 2)) .and. ((jed .eq. 14) .or. (jed .eq. 15)))then
							d_cut = rcut_enh
						elseif(((jed .eq. 1) .or. (jed .eq. 2)) .and. ((ied .eq. 14) .or. (ied .eq. 15)))then
							d_cut = rcut_enh
						else
							d_cut = rcut_reg
						endif

						if(d .lt. d_cut)then
            			write(*,*) 'relieving clash between ', is, js, 'd_cut = ', d_cut, 'd = ', d
						!	write(*,*) 'd before move: ',d
							! relieve that clash by moving each BAR away from each other along the is ij vector
							! by an amount 0.5*(d_cut - |r_ij|) 
							rij_x = (rx_s(js) - rx_s(is))/d
							rij_y = (ry_s(js) - ry_s(is))/d
							rij_z = (rz_s(js) - rz_s(is))/d
							delta_d = 0.5*(d_cut - d)
							!delta_d = 1.0
							!write(*,*) "delta_d: ",delta_d
							do ked=1,n_edbar
								ks = bar_map(iem,ked)
								ls = bar_map(jem,ked)
								!write(*,*) 'moving sites ', ks,ls, 'from ', rx_s(ks),ry_s(ks),rz_s(ks),rx_s(ls),ry_s(ls),rz_s(ls)
								rx_s(ks) = rx_s(ks) - 1.0*delta_d*rij_x ! move BAR i along -rij by 1/2 d-d_cut
								ry_s(ks) = ry_s(ks) - 1.0*delta_d*rij_y 
								rz_s(ks) = rz_s(ks) - 1.0*delta_d*rij_z
								rx_s(ls) = rx_s(ls) + 1.0*delta_d*rij_x ! move BAR i along rij by 1/2 d-d_cut
								ry_s(ls) = ry_s(ls) + 1.0*delta_d*rij_y 
								rz_s(ls) = rz_s(ls) + 1.0*delta_d*rij_z
								!write(*,*) 'to', rx_s(ks),ry_s(ks),rz_s(ks),rx_s(ls),ry_s(ls),rz_s(ls) 
							enddo
							!dnew = sqrt((rx_s(is) - rx_s(js))**2 + (ry_s(is) - ry_s(js))**2 + (rz_s(is) - rz_s(js))**2)
							!write(*,*) 'd after move: ',dnew
         			endif
					enddo
				enddo
			!endif
		enddo
	enddo

	write(*,*) 'FINISHED FIRST OVERLAP CHECK'

	!! BAR-BAR overlap check
	do iem=1,NBAR
		do jem=iem+1,NBAR
			!if(iem .ne. jem)then
				do ied=1,n_edbar
					do jed=1,n_edbar
						is = bar_map(iem,ied)
						js = bar_map(jem,jed)
						d = sqrt((rx_s(is) - rx_s(js))**2 + (ry_s(is) - ry_s(js))**2 + (rz_s(is) - rz_s(js))**2)
						if(((ied .eq. 1) .or. (ied .eq. 2)) .and. ((jed .eq. 14) .or. (jed .eq. 15)))then
							d_cut = rcut_enh
						elseif(((jed .eq. 1) .or. (jed .eq. 2)) .and. ((ied .eq. 14) .or. (ied .eq. 15)))then
							d_cut = rcut_enh
						else
							d_cut = rcut_reg
						endif

						if(d .lt. d_cut)then
            			write(*,*) 'relieving clash between ', is, js, 'd_cut = ', d_cut, 'd = ', d
							!write(*,*) 'd before move: ',d
							! relieve that clash by moving each BAR away from each other along the is ij vector
							! by an amount 0.5*(d_cut - |r_ij|) 
							rij_x = (rx_s(js) - rx_s(is))/d
							rij_y = (ry_s(js) - ry_s(is))/d
							rij_z = (rz_s(js) - rz_s(is))/d
							delta_d = 0.5*(d_cut - d)
							!delta_d = 1.0
							!write(*,*) "delta_d: ",delta_d
							do ked=1,n_edbar
								ks = bar_map(iem,ked)
								ls = bar_map(jem,ked)
								!write(*,*) 'moving sites ', ks,ls, 'from ', rx_s(ks),ry_s(ks),rz_s(ks),rx_s(ls),ry_s(ls),rz_s(ls)
								rx_s(ks) = rx_s(ks) - 1.0*delta_d*rij_x ! move BAR i along -rij by 1/2 d-d_cut
								ry_s(ks) = ry_s(ks) - 1.0*delta_d*rij_y 
								rz_s(ks) = rz_s(ks) - 1.0*delta_d*rij_z
								rx_s(ls) = rx_s(ls) + 1.0*delta_d*rij_x ! move BAR i along rij by 1/2 d-d_cut
								ry_s(ls) = ry_s(ls) + 1.0*delta_d*rij_y 
								rz_s(ls) = rz_s(ls) + 1.0*delta_d*rij_z
								!write(*,*) 'to', rx_s(ks),ry_s(ks),rz_s(ks),rx_s(ls),ry_s(ls),rz_s(ls) 
							enddo
							!dnew = sqrt((rx_s(is) - rx_s(js))**2 + (ry_s(is) - ry_s(js))**2 + (rz_s(is) - rz_s(js))**2)
							!write(*,*) 'd after move: ',dnew
         			endif
					enddo
				enddo
			!endif
		enddo
	enddo

	write(*,*) 'FINISHED SECOND OVERLAP CHECK'

	!! BAR-BAR overlap check
	do iem=1,NBAR
		do jem=iem+1,NBAR
			!if(iem .ne. jem)then
				do ied=1,n_edbar
					do jed=1,n_edbar
						is = bar_map(iem,ied)
						js = bar_map(jem,jed)
						d = sqrt((rx_s(is) - rx_s(js))**2 + (ry_s(is) - ry_s(js))**2 + (rz_s(is) - rz_s(js))**2)
						if(((ied .eq. 1) .or. (ied .eq. 2)) .and. ((jed .eq. 14) .or. (jed .eq. 15)))then
							d_cut = rcut_enh
						elseif(((jed .eq. 1) .or. (jed .eq. 2)) .and. ((ied .eq. 14) .or. (ied .eq. 15)))then
							d_cut = rcut_enh
						else
							d_cut = rcut_reg
						endif

						if(d .lt. d_cut)then
            			write(*,*) 'relieving clash between ', is, js, 'd_cut = ', d_cut, 'd = ', d
							!write(*,*) 'd before move: ',d
							! relieve that clash by moving each BAR away from each other along the is ij vector
							! by an amount 0.5*(d_cut - |r_ij|) 
							rij_x = (rx_s(js) - rx_s(is))/d
							rij_y = (ry_s(js) - ry_s(is))/d
							rij_z = (rz_s(js) - rz_s(is))/d
							delta_d = 0.5*(d_cut - d)
							!delta_d = 1.0
							!write(*,*) "delta_d: ",delta_d
							do ked=1,n_edbar
								ks = bar_map(iem,ked)
								ls = bar_map(jem,ked)
								!write(*,*) 'moving sites ', ks,ls, 'from ', rx_s(ks),ry_s(ks),rz_s(ks),rx_s(ls),ry_s(ls),rz_s(ls)
								rx_s(ks) = rx_s(ks) - 1.0*delta_d*rij_x ! move BAR i along -rij by 1/2 d-d_cut
								ry_s(ks) = ry_s(ks) - 1.0*delta_d*rij_y 
								rz_s(ks) = rz_s(ks) - 1.0*delta_d*rij_z
								rx_s(ls) = rx_s(ls) + 1.0*delta_d*rij_x ! move BAR i along rij by 1/2 d-d_cut
								ry_s(ls) = ry_s(ls) + 1.0*delta_d*rij_y 
								rz_s(ls) = rz_s(ls) + 1.0*delta_d*rij_z
								!write(*,*) 'to', rx_s(ks),ry_s(ks),rz_s(ks),rx_s(ls),ry_s(ls),rz_s(ls) 
							enddo
							!dnew = sqrt((rx_s(is) - rx_s(js))**2 + (ry_s(is) - ry_s(js))**2 + (rz_s(is) - rz_s(js))**2)
							!write(*,*) 'd after move: ',dnew
         			endif
					enddo
				enddo
			!endif
		enddo
	enddo


   write(*,*) "Shall I rotate the BARs? 1 for yes, 0 for no"
   read(*,*) shall_I_rotate
   
   open(unit=17, file=output)
   write(17,'(i8)') n_total
   write(17,'("          ")')

  
   if (shall_I_rotate == 1) then
      is = 1;js = 1 
      do i=1,NBAR
         rotate = RAND(0)
         rotate = 3.14*rotate
      
         comx = 0.0d0
         comy = 0.0d0
         do j=1,n_edbar
            comx = comx + rx_s(js)
            comy = comy + ry_s(js)
            comz = comz + rz_s(js)
            js = js + 1      
         enddo
     
         comx = comx/n_edbar
         comy = comy/n_edbar
         comz = comz/n_edbar

      
        do j=1,n_edbar
	       itmp = is-1
           ied = MOD(itmp,26) + 1
		   atom = atom_name_s(is)
		   tmpx = (rx_s(is)-comx)*COS(rotate)-(ry_s(is)-comy)*SIN(rotate)
		   tmpy = (rx_s(is)-comx)*SIN(rotate)+(ry_s(is)-comy)*COS(rotate)
		   tmpx = tmpx + comx
		   tmpy = tmpy + comy
		   tmpz = rz_s(is) 
		   write(17,*) atom, tmpx, tmpy, tmpz
		   is = is + 1
	     enddo  
      enddo
   
   else
	  do is=1,n_total
		itmp = is-1
		ied = MOD(itmp,26) + 1
		!ied = ((i-1) % 26) + 1
		tmpx = rx_s(is)
		tmpy = ry_s(is)
		tmpz = rz_s(is)
		!if (MOD(is,3) .ne. 0) then				!Only write 2/3 of BARs
		write(17,*) atom_name(ied), tmpx, tmpy, tmpz
		!endif
	  enddo	
   endif
   close(17)

end Program
