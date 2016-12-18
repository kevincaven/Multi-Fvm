module eleinfo
	use vectors
	use constants, only: ki, dp, pi
	use toggles, only: gamma_index
	use domainmesh, only: mx, my, hx, hy, h, h1, h2
	use openbasis, only: p, &
	n, nx, ny, nz, hhx, hhy
	use localid, only: bf_nid1, bf_nid2
	implicit none

	private
	public get_ele, element, test_eleinfo

	type element
		integer(kind=ki) :: mi_me, mj_me, index_me
		type(vector) :: node(3)
		type(vector) :: barycenter
		real(kind=dp) :: pa_lenth(3)
		type(vector) :: nvec_matrix(3, 3)

		real(kind=dp), allocatable :: value(:, :) ! value(3, nz)
		real(kind=dp) :: gamma_matrix(3, 3)
		type(vector), allocatable :: grad(:, :) ! grad(3, n**2)
	end type element

contains

	type(element) &
	function get_ele(mi, mj, index1) result( ele )
		implicit none
		integer(kind=ki), intent(in) :: mi, mj, index1
		integer(kind=ki) :: node, nid

! 		mi_me, mj_me, index_me
		ele%mi_me = mi
		ele%mj_me = mj
		ele%index_me = index1		
! 		node(3)
		if ( index1 == 1 ) then
			ele%node%x = [ hx*(mi - 1), hx*mi, hx*(mi - 1) ]
			ele%node%y = [ hy*(mj - 1), hy*mj, hy*mj ]
		elseif (index1 == 2) then
			ele%node%x = [ hx*(mi - 1), hx*mi, hx*mi ]
			ele%node%y = [ hy*(mj - 1), hy*(mj - 1), hy*mj ]
		else
			print*, 'error: get_ele, node '
		end if
! 		barycenter
		ele%barycenter = 1/3._dp*(ele%node(1) + ele%node(2) + ele%node(3))
! 		distance of barycenter(p) and node(ai)
		if ( index1 == 1 ) then
			ele%pa_lenth = [ 1/3*h1, 1/3*h2, 1/3*h ]
		elseif (index1 == 2) then
			ele%pa_lenth = [ 1/3*h2, 1/3*h, 1/3*h1 ]
		else
			print*, 'error: get_ele, pa_lenth '
		end if
! 		nvec_matrix(3, 3)
		call get_nvec(index1, ele%nvec_matrix)
! 		value(nz)
		allocate(ele%value(3, nz))
		do node = 1, 3
		do nid = 1, nz
			ele%value(node, nid) = p(index1, node, mi, mj, nid)
		end do
		end do
! 		gamma_matrix(3, 3)
		call get_gamma(mi, mj, index1, ele%gamma_matrix)
! 		grad(3, n**2)
		allocate(ele%grad(3, n**2))
		do node = 1, 3
			call get_grad(index1, ele%value(node, :), ele%grad(node, :))
		end do

		return		
	end function get_ele

	subroutine get_nvec(index, nmat)
		implicit none
		integer(kind=ki), intent(in) :: index
		type(vector), intent(out) :: nmat(3, 3)

		if ( index == 1 ) then
			nmat(1, 1) = [0, 1]
			nmat(2, 2) = [-1, 0]
			nmat(3, 3) = [hy/h, -hx/h]

			nmat(1, 2) = [hy/h, hx/h]
			nmat(2, 1) = -nmat(1, 2)

			nmat(2, 3) = [-2*hy/h1, hx/h1]
			nmat(3, 2) = -nmat(2, 3)

			nmat(1, 3) = [-hy/h2, 2*hx/h2]
			nmat(3, 1) = -nmat(1, 3)
		elseif (index == 2) then
			nmat(1, 1) = [1, 0]
			nmat(2, 2) = [-hy/h, hx/h]
			nmat(3, 3) = [0, -1]

			nmat(1, 3) = [hy/h, hx/h]
			nmat(3, 1) = -nmat(1, 3)

			nmat(2, 1) = [-2*hy/h1, hx/h1]
			nmat(1, 2) = -nmat(2, 1)

			nmat(2, 3) = [-hy/h2, 2*hx/h2]
			nmat(3, 2) = -nmat(2, 3)
		else
			print*, 'error: get_ele, get_nvec '
		end if

		return
	end subroutine get_nvec

	subroutine get_gamma(mi, mj, index1, matg)
		implicit none
		integer(kind=ki), intent(in) :: mi, mj, index1
		real(kind=dp), intent(out) :: matg(3, 3)
		integer(kind=ki) :: phi_i, phi_j, i, j

		matg = 0._dp

		if ( gamma_index .eq. 1 ) then
			if ( index1 == 1 ) then
				do phi_i = 1, 3
					do i = 1, n
						matg(phi_i, 1) = matg(phi_i, 1) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i, n + 1)) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i + 1, n + 1))
					end do

					do i = 1, n
						matg(phi_i, 2) = matg(phi_i, 2) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, 1, i)) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, 1, i + 1))
					end do

					do i = 1, n
						matg(phi_i, 3) = matg(phi_i, 3) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i, i)) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i + 1, i + 1))
					end do
				end do
			elseif ( index1 == 2 ) then
				do phi_i = 1, 3
					do i = 1, n
						matg(phi_i, 1) = matg(phi_i, 1) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, n + 1, i)) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, n + 1, i + 1))
					end do

					do i = 1, n
						matg(phi_i, 2) = matg(phi_i, 2) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i, i)) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i + 1, i + 1))
					end do

					do i = 1, n
						matg(phi_i, 3) = matg(phi_i, 3) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i, 1)) + &
						0.5_dp/n*p(index1, phi_i, mi, mj, bf_nid1(index1, i + 1, 1))
					end do
				end do
			else
				print*, 'error: get_gamma, index1'
			end if			

		elseif (gamma_index .eq. 2) then
			matg = 0.5_dp
			matg(1, 1) = 0._dp
			matg(2, 2) = 0._dp
			matg(3, 3) = 0._dp
		else
			print*, 'error: get_gamma, gamma_index'
		end if

		return
	end subroutine get_gamma

	subroutine get_grad(index1, varray, grad)
		implicit none
		integer(kind=ki), intent(in) :: index1
		real(kind=dp), intent(in) :: varray(nz)
		type(vector), intent(out) :: grad(n**2)
		integer(kind=ki) :: ni, nj, index2, eleid, nid(3)
		integer(kind=ki) :: i

		grad%x = 0._dp
		grad%y = 0._dp

		if ( index1 == 1 ) then
			eleid = 0
			do nj = 1, n
			do ni = 1, nj
				index2 = 1
				eleid = eleid + 1
				nid(1) = bf_nid2(index1, ni, nj, index2, 1)
				nid(2) = bf_nid2(index1, ni, nj, index2, 2)
				nid(3) = bf_nid2(index1, ni, nj, index2, 3)
				grad(eleid)%x = ( varray(nid(2)) - varray(nid(3)) )/hhx
				grad(eleid)%y = ( varray(nid(3)) - varray(nid(1)) )/hhy

				if ( ni .ne. nj )  then
				index2 = 2
				eleid = eleid + 1
				nid(1) = bf_nid2(index1, ni, nj, index2, 1)
				nid(2) = bf_nid2(index1, ni, nj, index2, 2)
				nid(3) = bf_nid2(index1, ni, nj, index2, 3)
				grad(eleid)%x = ( varray(nid(2)) - varray(nid(1)) )/hhx
				grad(eleid)%y = ( varray(nid(3)) - varray(nid(2)) )/hhy
				end if
			end do
			end do
		elseif ( index1 == 2 ) then
			eleid = 0
			do nj = 1, n
			do ni = nj, n
				if ( ni .ne. nj )  then
				index2 = 1
				eleid = eleid + 1
				nid(1) = bf_nid2(index1, ni, nj, index2, 1)
				nid(2) = bf_nid2(index1, ni, nj, index2, 2)
				nid(3) = bf_nid2(index1, ni, nj, index2, 3)
				grad(eleid)%x = ( varray(nid(2)) - varray(nid(3)) )/hhx
				grad(eleid)%y = ( varray(nid(3)) - varray(nid(1)) )/hhy
				end if

				index2 = 2
				eleid = eleid + 1
				nid(1) = bf_nid2(index1, ni, nj, index2, 1)
				nid(2) = bf_nid2(index1, ni, nj, index2, 2)
				nid(3) = bf_nid2(index1, ni, nj, index2, 3)
				grad(eleid)%x = ( varray(nid(2)) - varray(nid(1)) )/hhx
				grad(eleid)%y = ( varray(nid(3)) - varray(nid(2)) )/hhy
			end do
			end do
		else
			print*, 'error: get_grad, index1'
		end if

! 		do i = 1, n**2
! 			print*,index1, i, grad(i)
! 		end do

		return
	end subroutine get_grad

	subroutine test_eleinfo()
		implicit none
		integer(kind=ki) :: i, j, node, nid, eleid
		type(vector) :: v1, v2
		type(element) :: ele

		ele = get_ele(1, 1, 2)
		v1%x = 1; v1%y = 2; v2%x = 3; v2%y = 4

		print*, '----- test_eleinfo -----'
		print*, 3._dp*v1, v2*4._dp
! 		print*, v1*v2, v1+v2, v1-v2
		
! 		print*, ele%node(2)

! 		print*, ele%barycenter

! 		do i = 1, 3
! 			do j = 1, 3
! 				print*, i, j, ele%nvec_matrix(i, j)
! 			end do
! 		end do

! 		print*, ele%pa_lenth

! 		do eleid = 1, n**2
! 			print*, eleid, ele%cor(eleid)
! 		end do

! 		do node = 1, 3
! 			do nid = 1, nz
! 				print*, node, nid, ele%value(node, nid)
! 			end do
! 		end do

! 		do i = 1, 3
! 			do j = 1, 3
! 				print*, i, j, ele%gamma_matrix(i, j)
! 			end do
! 		end do

! 		do i = 1, 3
! 			do j = 1, n**2
! 				print*, i, j, ele%grad(i, j)%x, ele%grad(i, j)%y
! 			end do
! 		end do

		return
	end subroutine test_eleinfo

end module eleinfo
