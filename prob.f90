module prob
 	use constants, only: ki, dp, pi
	use domainmesh, only: mx, my, nel, hx, hy
	use toggles, only: eg_index, ineps, finestmesh_index
	use getlgmap, only: kidof, kbdof, lgmap, cor
	implicit none

	private
	public get_prob, coef, force, ntruesol, &
	truex, ub, eg_index
! 	true_sol, true_sol_gradx, true_sol_grady, &

	integer(kind=ki) :: ntrue
	real(kind=dp), allocatable :: truex(:, :)
	real(kind=dp), allocatable :: ub(:)

contains

	subroutine get_prob()
		implicit none
		
		print*, '----- step3: get prob -----'
		print*
		call ntruesol()
		call get_ub()
		call check_info()
! 		call test()

		return		
	end subroutine get_prob

	real(kind=dp) &
	function coef(x, y)
		implicit none
		real(kind=dp) :: x, y
		real(kind=dp) :: pp, eps_1, eps_2, eps 

		pp = 1.8_dp
		eps = 1._dp/ineps
		eps_1 = eps 
		eps_2 = eps

		select case (eg_index)
			case (1)
				coef = 1.0_dp
			case (2)		
				coef = 1.0_dp
			case (3)
				coef = &
				(2._dp + pp*sin(2._dp*pi*x/eps_1)) / ( 2._dp + pp*cos(2._dp*pi*y/eps_1)) + &
				(2._dp + pp*sin(2._dp*pi*y/eps_2)) / (2._dp + pp*sin(2._dp*pi*x/eps_2))
			case (4)
				coef = &
				(2._dp + pp*sin(2._dp*pi*x/eps_1)) / ( 2._dp + pp*cos(2._dp*pi*y/eps_1)) + &
				(2._dp + sin(2._dp*pi*y/eps_2)) / (2._dp + pp*cos(2._dp*pi*x/eps_2))

			case default
				print*, 'error: coef, eg_index'
		end select

		return
	end function coef

! source item: f(x)

	real(kind=dp) &
	function force(x, y)
		implicit none
		real(kind=dp) :: x, y

		select case (eg_index)
			case (1)
				force = -1.0_dp
			case (2)		
				force = 2*pi*pi * cos((x + 1/2._dp)*pi)*sin(y*pi)
			case (3)
				force = -1.0_dp
			case (4)
				force = -1.0_dp
			case default
				print*, 'problem error: force item f(x, y)! '
		end select

		return
	end function force

!	boundary truesol
	
	real(kind=dp)  &
	function boundvalue(x, y)
		implicit none
		real(kind=dp) :: x, y

		select case (eg_index)
			case (1)
				boundvalue = 1/4._dp*(x**2+y**2) + 2._dp
			case (2)		
! 				boundvalue = cos((x + 1/2._dp)*pi)*sin(y*pi)
				boundvalue = 0_dp
			case (3)
				boundvalue = 0._dp
			case (4)
				boundvalue = 0._dp
			case default
				print*, 'problem error: true solution u(x, y)! '
		end select
		
		return
	end function boundvalue

! 	get ub

	subroutine get_ub()
		implicit none
		integer(kind=ki) :: i, j, nid

		allocate(ub(kbdof))
		open(11, file='data_out/ub.dat')
		do i = 1, nel
		do j = 1, 3
			nid = lgmap(i, j)
			if (nid .le. kbdof) then
				ub(nid) = boundvalue(cor(nid ,1), cor(nid ,2))
				write(11,*) i, j, nid, ub(nid)
			end if
		enddo
		end do
		close(11)
		
		return
	end subroutine get_ub

! 	true solution, get true solution over \Omega

	subroutine ntruesol()
		implicit none	
		real(kind=dp) :: x, y
		integer(kind=ki) :: i, j, k, localid, kdofid
		integer(kind=ki) :: mi, mj
		real(kind=dp) :: value
		integer(kind=ki) :: ktemp = 0
		character(len=100) :: datafile, npath

		ntrue = finestmesh_index

		write(npath, *) ntrue

		select case (eg_index)
			case (1)
			datafile = '../truesol/eg1_linear_'//trim(adjustl(npath))//'.dat'
			case (2)
			datafile = '../truesol/eg2_linear_'//trim(adjustl(npath))//'.dat'
			case (3)
			datafile = '../truesol/eg3_linear_'//trim(adjustl(npath))//'.dat'
			case (4)
			datafile = '../truesol/eg4_linear_'//trim(adjustl(npath))//'.dat'			
		end select

		print*, 'truesol is from: '
		print*, datafile

		allocate(truex(ntrue + 1, ntrue + 1))

		open(12, file=datafile, status='old')
			do mj = 1, ntrue + 1
			do mi = 1, ntrue + 1
				read(12, *) truex(mi, mj)
			end do
			end do
		close(12)

		return
	end subroutine ntruesol

	subroutine check_info()
		implicit none
		
		print*, 'true solution infomation: '
		print*
				
		return		
	end subroutine check_info

	subroutine test()
		implicit none
		integer(kind=ki) :: i, j

		print*, coef(0.1_dp, 0.2_dp)
			
		do i = 1, 1025
		do j = 1, 1025
			print*, coef(1/1024._dp*(i - 1), 1/1024._dp*(j - 1))
		enddo
		enddo

		return
	end subroutine test

end module prob