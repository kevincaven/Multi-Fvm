module domainmesh
	use toggles
 	use constants
	implicit none
	
	private
	public get_domainmesh, ntran, &
	mx, my, nel, hx, hy, h, h1, h2, area, xleft, xright, ybottom, ytop

	real(kind = dp) :: xleft, xright, ybottom, ytop
	integer(kind = ki) :: mx, my
	integer(kind = ki) :: nel
	real(kind = dp) :: hx, hy, h, h1, h2, area

contains

	subroutine get_domainmesh()
		implicit none

		print*, '----- step1: get parameters -----'
		print*
		call get_domain()
		call get_mesh()
	 	call check_info()
! 	 	call test()

	 	return
	end subroutine get_domainmesh

	integer(kind=ki) &
	function ntran(nid)
		implicit none
		integer(kind=ki), intent(in) :: nid
		ntran = mod(nid + 2, 3) + 1
		return
	end function ntran

	subroutine get_domain()
		implicit none

		open(11, file='data_in/domain.dat', status='old')
		read(11, *) xleft, xright, ybottom, ytop
		close(11)

		return
	end subroutine get_domain

	subroutine get_mesh()
		implicit none

! ! 		way 1
! 		print*, "please input mx & my! "
! 		read(*,*)  mx, my
! ! 		way 2	
! 		open(11, file='data_in/mesh.dat', status='old')
! 			read(11, *) mx, my
! 		close(11)
! 		way3
		mx = coarsemesh_index
		my = mx

		nel = mx*my*2
		hx = (xright - xleft)/mx 
		hy = (ytop - ybottom)/my
		h = sqrt(hx**2 + hy**2)
		h1 = sqrt(hx**2 + 4*hy**2)
		h2 = sqrt(hy**2 + 4*hx**2)
		area = hx * hy /2.0_dp

		return
	end subroutine get_mesh

! 	subroutine get_pcoef_inpes()
! 		implicit none

! 		open(11, file='data_in/pcoef.dat', status='old')
! 		read(11, *) pcoef
! 		close(11)

! 		open(11, file='data_in/ineps.dat', status='old')
! 		read(11, *) ineps
! 		close(11)

! 		return
! 	end subroutine get_pcoef_inpes

	subroutine check_info()
		implicit none
		
		print*, 'domain information: '
		print*, 'xleft/xright is:', xleft, xright
		print*, 'ybottom/ytop is:', ybottom, ytop
		print*

		print*, 'mesh information: '
		print*, 'mx/my is:', mx, my
		print*, 'total number of elements is: ', nel
		print*, 'hx/hy is:', hx, hy
		print*, 'h is :', h
		print*, 'element area is:', area
		print*

		print*, 'other information: '
		print*, 'penalty coef is: ', pcoef
		print*, 'ineps is: ', ineps 
		print*

		return		
	end subroutine check_info

	subroutine test()
		implicit none
		
		print*, ntran(4), ntran(0)

		return
	end subroutine test

end module domainmesh