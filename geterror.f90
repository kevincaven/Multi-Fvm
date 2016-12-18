module geterror
 	use constants, only: ki, dp, pi
 	use vectors
 	use domainmesh, only: mx, my, nel, h, hx, hy, area, ntran
 	use openbasis, only: p, n, hhx, hhy
	use getlgmap, only: kdof, kidof, kbdof, lgmap, cor, get_elid
	use prob, only: truex, ub, coef
	use localid, only: bf_nid2
	use solveequation, only: xspan, yspan
	implicit none

	private
	public get_error

	real(kind=dp) :: truel2, trueh1, errorl2, errorh1
	
contains

	subroutine get_error()
		implicit none

		print*, '----- step8: caculate the error -----'
		print*
		call pointerr()
		call globalerr()
		call check_info()
! 		call test()

		return
	end subroutine get_error

	subroutine pointerr()
		implicit none
		integer(kind=ki) :: k, mi, mj, index, elid
		integer(kind=ki) :: node, itrue(3), jtrue(3)
		open(11, file='err/err_absolute.dat')
		open(12, file='err/err_relative.dat')
		do mj = 1, my
		do mi = 1, mx
			index = 1
			elid = get_elid(mi, mj, index)
			itrue = [(mi - 1)*n + 1, mi*n + 1, (mi - 1)*n + 1]
			jtrue = [(mj - 1)*n + 1, mj*n + 1, mj*n + 1]
			do node = 1, 3
				k = lgmap(elid, node)
				if ( k .gt. kbdof ) then
				write(11, *) elid, k, mi, mj, yspan(k) - truex(itrue(node), jtrue(node))
				write(12, *) elid, k, mi, mj, (yspan(k) - truex(itrue(node), jtrue(node)))/truex(itrue(node), jtrue(node))
				endif
			end do

			index = 2
			elid = get_elid(mi, mj, index)
			itrue = [(mi - 1)*n + 1, mi*n + 1, mi*n + 1]
			jtrue = [(mj - 1)*n + 1, (mj - 1)*n + 1, mj*n + 1]
			do node = 1, 3
				k = lgmap(elid, node)
				if ( k .gt. kbdof ) then
				write(11, *) elid, k, mi, mj, yspan(k) - truex(itrue(node), jtrue(node))
				write(12, *) elid, k, mi, mj, (yspan(k) - truex(itrue(node), jtrue(node)))/truex(itrue(node), jtrue(node))
				endif
			end do
		end do
		end do
		close(11)		
		close(12)

		return
	end subroutine pointerr

	subroutine globalerr()
		implicit none
		integer(kind=ki) :: mi, mj, ni, nj, index1, index2
		real(kind=dp) :: area_fine
		real(kind=dp) :: utrue_l2, utrue_h1, err_l2, err_h1
		real(kind=dp) :: usol(3), utrue(3) 
		real(kind=dp) :: perm
		type(vector) :: grad_nbf(3)

		area_fine = hhx*hhy/2

		truel2 = 0._dp
		trueh1 = 0._dp
		errorl2 = 0._dp
		errorh1 = 0._dp

		do mj = 1, my
		do mi = 1, mx
			do nj = 1, n
			do ni = 1, n
				do index2 = 1, 2
	! 				get grad_nbf
					call grad_nodebasis(grad_nbf, index2, hhx, hhy)
	! 				get usol, utrue, perm, index1
					call getsol_local(usol, utrue, perm, index1, mi, mj, ni, nj, index2)
	! 				get utrue_l2, utrue_h1, err_l2, err_h1
					call geterror_local(utrue_l2, err_l2, utrue_h1, err_h1, usol, utrue, grad_nbf, perm)

					truel2 = truel2 + utrue_l2
					trueh1 = trueh1 + utrue_h1
					errorl2 = errorl2 + err_l2
					errorh1 = errorh1 + err_h1
				end do
			end do			
			end do		
		end do
		end do
		
		truel2 = sqrt(truel2*area_fine)
		errorl2 = sqrt(errorl2*area_fine)
		trueh1 = sqrt(trueh1*area_fine)
		errorh1 = sqrt(errorh1*area_fine)

		open(11, file = 'err/err_l2')
			write(11, *) mx
			write(11, *) n
			write(11, *) errorl2, truel2, errorl2/truel2
		close(11)

		open(12, file = 'err/err_h1')
			write(12, *) mx
			write(12, *) n
			write(12, *) errorh1, trueh1, errorh1/trueh1
		close(12)

		return
	end subroutine globalerr

	subroutine getsol_local(usol, utrue, perm, index1, mi, mj, ni, nj, index2)
		implicit none
		real(kind=dp), intent(out) :: usol(3), utrue(3)
		integer(kind=ki), intent(out) :: index1
		real(kind=dp), intent(out) :: perm
		integer(kind=ki), intent(in) :: mi, mj, ni, nj, index2

		integer(kind=ki) :: elid
		integer(kind=ki) :: nid
		integer(kind=ki) :: i, j, itrue(3), jtrue(3)


		if (ni .lt. nj) then
			index1 = 1
		elseif (ni .gt. nj) then
			index1 = 2
		elseif (ni .eq. nj) then
			index1 = index2
		else
			print*, 'error: getsol_local, index1'
		endif

! 		get usol(3)
		usol = 0._dp
		elid = get_elid(mi, mj, index1)
		do i = 1, 3	! fine element node index
			nid = bf_nid2(index1, ni, nj, index2, i)
			do j = 1, 3	! basis function index
				usol(i) = usol(i) + xspan(lgmap(elid, j))*p(index1, j, mi, mj, nid)
			end do
		end do

! 		get utrue(3)
		utrue = 0._dp
		i = (mi - 1)*n + (ni - 1) + 1
		j = (mj - 1)*n + (nj - 1) + 1
		if ( index2 .eq. 1 ) then
			itrue = [i, i + 1, i]
			jtrue = [j, j + 1, j + 1]
			utrue(1) = truex(itrue(1), jtrue(1))			
			utrue(2) = truex(itrue(2), jtrue(2))			
			utrue(3) = truex(itrue(3), jtrue(3))

			perm = 1/3._dp*( &
			coef((itrue(1) - 1)*hhx, (jtrue(1) - 1)*hhy) + &
			coef((itrue(2) - 1)*hhx, (jtrue(2) - 1)*hhy) + &
			coef((itrue(3) - 1)*hhx, (jtrue(3) - 1)*hhy))
		elseif ( index2 .eq. 2 ) then
			itrue = [i, i + 1, i + 1]
			jtrue = [j, j, j + 1]
			utrue(1) = truex(itrue(1), jtrue(1))			
			utrue(2) = truex(itrue(2), jtrue(2))			
			utrue(3) = truex(itrue(3), jtrue(3))

			perm = 1/3._dp*( &
			coef((itrue(1) - 1)*hhx, (jtrue(1) - 1)*hhy) + &
			coef((itrue(2) - 1)*hhx, (jtrue(2) - 1)*hhy) + &
			coef((itrue(3) - 1)*hhx, (jtrue(3) - 1)*hhy))
		else
			print*, 'error: getsol_local, utrue' 			
		end if

! 		perm = 1._dp

		return		
	end subroutine getsol_local

! 	local fine element integral, three middle point method
! 	get local norm^2 / area

	subroutine geterror_local(utrue_l2, err_l2, utrue_h1, err_h1, usol, utrue, grad_nbf, perm)
		implicit none
		real(kind=dp), intent(out) :: utrue_l2, utrue_h1, err_l2, err_h1
		real(kind=dp), intent(in) :: usol(3), utrue(3)
		real(kind=dp), intent(in) :: perm
		type(vector), intent(in) :: grad_nbf(3)

		real(kind=dp) :: utrue_l2tmp, err_l2tmp
		type(vector) :: grad_utrue, grad_err
		integer(kind=ki) :: node, node1, node2

		utrue_l2tmp = 0._dp
		err_l2tmp = 0._dp
		do node = 1, 3
			node1 = ntran(node + 1)
			node2 = ntran(node - 1)
			utrue_l2tmp = utrue_l2tmp + ( ( utrue(node1) + utrue(node2) )/2 )**2
			err_l2tmp = err_l2tmp + ( ( usol(node1) + usol(node2) - utrue(node1) - utrue(node2) )/2 )**2
		end do
		utrue_l2 = utrue_l2tmp/3
		err_l2 = err_l2tmp/3

		grad_utrue = utrue(1)*grad_nbf(1) + utrue(2)*grad_nbf(2) + utrue(3)*grad_nbf(3)
		grad_err = (usol(1) - utrue(1))*grad_nbf(1) + (usol(2) - utrue(2))*grad_nbf(2) + (usol(3) - utrue(3))*grad_nbf(3)

		utrue_h1 = grad_utrue*grad_utrue*perm
		err_h1 = grad_err*grad_err*perm

		return		
	end subroutine geterror_local

	subroutine grad_nodebasis(grad_nbf, index2, hxx, hyy)
		implicit none
		type(vector), intent(out) :: grad_nbf(3)
		integer(kind=ki), intent(in) :: index2
		real(kind=dp), intent(in) :: hxx, hyy

! 		if (index2 .eq. 1) then 
! 			grad_nbf%x = [0._dp, 1._dp/hxx, -1._dp/hxx]
! 			grad_nbf%y = [-1._dp/hyy, 0._dp, 1._dp/hyy]
! 		elseif (index2 .eq. 2) then
! 			grad_nbf%x = [-1._dp/hxx, 1._dp/hxx, 0._dp]
! 			grad_nbf%y = [0._dp, -1._dp/hyy, 1._dp/hyy]
! 		endif

		if (index2 .eq. 1) then 
			grad_nbf%x = [0._dp, 1/hxx, -1/hxx]
			grad_nbf%y = [-1/hyy, 0._dp, 1/hyy]
		elseif (index2 .eq. 2) then
			grad_nbf%x = [-1/hxx, 1/hxx, 0._dp]
			grad_nbf%y = [0._dp, -1/hyy, 1/hyy]
		endif


		return		
	end subroutine grad_nodebasis

	subroutine check_info()
		implicit none
		
		print*, 'point error infomation: '
		print*, 'L^2 error is:', mx, errorl2, errorl2/truel2
		print*, 'H^1 error is:', mx, errorh1, errorh1/trueh1
		print*
		print*, 'L^2, H^1norm of true solution: '
		print*, mx, truel2, trueh1
		print*

		print*, 'please check: '
		print*, 'err/err_absolute.dat'
		print*, 'err/err_relative.dat'
		print*, 'err/err_l2.dat'
		print*, 'err/err_h1.dat'
		print*

		return
	end subroutine check_info

	subroutine test()
		implicit none
		real(kind=dp) :: utrue_l2, utrue_h1, err_l2, err_h1
		real(kind=dp) :: usol(3), utrue(3), bk1(3), bk2(3)
		real(kind=dp) :: perm
		integer(kind=ki) :: mi, mj, ni, nj, index1, index2
		type(vector) :: grad_nbf(3)

		usol = [1._dp, 1._dp, 1._dp]
		utrue = [0.5_dp, 0.5_dp, 0.5_dp]
		bk1 = [0._dp, 1._dp, -1._dp]
		bk2 = [-1._dp, 1._dp, 0._dp]
		grad_nbf%x = bk1
		grad_nbf%y = bk2

		print*, '----- geterror test -----'
		print*

!!		test geterror_local
! 		call get_integral(err_l2, utrue_l2, err_h1, utrue_h1, usol, utrue, bk1, bk2)
! 		print*, err_l2, utrue_l2, err_h1, utrue_h1

! 		call geterror_local(utrue_l2, utrue_h1, err_l2, err_h1, usol, utrue, grad_nbf)
! 		print*, err_l2, utrue_l2, err_h1, utrue_h1

! 		test getsol_local
		mi = 5; mj = 5; ni = 2; nj = 2; index2 = 1
		call getsol_local(usol, utrue, perm, index1, mi, mj, ni, nj, index2)
		print*, usol
		print*, utrue

		print*

		return
	end subroutine test

end module geterror
