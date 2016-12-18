module integfunction
	
 	use constants, only: ki, dp, pi 
 	use vectors
	use domainmesh, only: mx, my, nel, hx, hy, h, area, ntran
	use toggles, only: penalty_index, edgeint_index, pcoef
	use prob, only: coef, force
	use openbasis, only: hhx, hhy, n, p
	use localid, only: segn, edge_nid, edge_eleid, pa_eleid, corn, corm
	use eleinfo
	implicit none

	public test_integfunction

contains

! ------- right hands item -------
	real(kind=dp) &
	function righthand_integ(mi, mj, index1, node)
		implicit none
		integer(kind=ki), intent(in) :: index1, mi, mj, node
		integer(kind=ki) :: i
		real(kind=dp) :: forcevalue(3)
		type(element) :: ele
		type(vector) :: midnode(3), midp(3)

		ele = get_ele(mi, mj, index1)

		do i = 1, 3
			midnode(i) = 1/2._dp*(ele%node(ntran(i + 1)) + ele%node(ntran(i - 1)))
			midp(i) = 1/2._dp*(ele%node(i) + ele%barycenter)
		end do

		do i = 1, 3
			forcevalue(i) = 1/3._dp*( &
				force(midnode(i)%x, midnode(i)%y) + &
				force(midp(ntran(i + 1))%x, midp(ntran(i + 1))%y) + &
				force(midp(ntran(i - 1))%x, midp(ntran(i - 1))%y) )
		end do

		righthand_integ = 1/3._dp*area*(&
		ele%gamma_matrix(node, 1)*forcevalue(1) + &
		ele%gamma_matrix(node, 2)*forcevalue(2) + &
		ele%gamma_matrix(node, 3)*forcevalue(3) )

		return
	end function righthand_integ

! ------- element_integ---------

	subroutine element_integ(stif, mi, mj, index1)
		implicit none
		integer(kind=ki), intent(in) :: index1, mi, mj
		real(kind=dp), intent(out) :: stif(3, 3)
		
		real(kind=dp) :: integ(3)
		real(kind=dp) :: coeffa

		integer(kind=ki) :: solb, testb, kt ! 1, 2, 3
		integer(kind=ki) :: ni, nj, kpa, eid, part
		type(element) :: ele

		ele = get_ele(mi, mj, index1)

		stif = 0._dp
		do testb = 1, 3
		do solb = 1, 3
! 			pa_nt, pa_{nt + 1}, pa_{nt - 1}
			do kt = 1, 3
				integ(kt) = 0._dp
				do kpa = 1, segn + 1
					eid = pa_eleid(index1, kt, kpa)

					coeffa = coef_int(corn(mi, mj, index1, kt, kpa)%x, corn(mi, mj, index1, kt, kpa)%y, &
						corm(mi, mj, index1, kt, kpa)%x, corm(mi, mj, index1, kt, kpa)%y )

					integ(kt) = integ(kt) + &
					ele%grad(solb, eid)*ele%nvec_matrix(ntran(kt - 1), ntran(kt + 1))*coeffa
				end do
				stif(testb, solb) = stif(testb, solb) + integ(kt)*&
				(ele%gamma_matrix(testb, ntran(kt - 1)) - ele%gamma_matrix(testb, ntran(kt + 1)))
			end do
		end do
		end do

		return
	end subroutine element_integ

! ---------------edge_integ-------------------
	
	subroutine get_elestiff(stif, mi, mj)
		implicit none
		integer(kind=ki), intent(in) :: mi, mj
! 		stif(index1, edge, neigh_or_self, i, j)
		real(kind=dp), intent(out) :: stif(2, 3, 2, 3, 3)

		real(kind=dp) :: stif1(3, 3), stif2(3, 3), &
		int1self(3, 3), int2self(3, 3), penself(3, 3), int1neigh(3, 3), int2neigh(3, 3), penneigh(3, 3)
		real(kind=dp) :: he
		real(kind=dp) :: coeffa

		integer(kind=ki) :: index1, index2
		integer(kind=ki) :: isolateid, isolateid_n
		integer(kind=ki) :: i, j	! 1, 2, 3
		integer(kind=ki) :: kpa, eid, ni, nj
		type(vector) :: edgenode(2)
		type(element) :: eleself, eleneigh, ele1, ele2

		ele1 = get_ele(mi, mj, 1)
		ele2 = get_ele(mi, mj, 2)

		stif = 0._dp

! 		1_up
		eleself = ele1
		he = hx
		index1 = 1
		isolateid = 1
		edgenode(1) = eleself%node(3)
		edgenode(2) = eleself%node(2)
		if ( mj .ne. my ) then
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			isolateid_n = 3
			eleneigh = get_ele(mi, mj + 1, 2)
			call get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
			stif(index1, isolateid, 1, :, :) = int1self + int2self + penself
			stif(index1, isolateid, 2, :, :) = int1neigh + int2neigh + penneigh
		else
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			stif(index1, isolateid, 1, :, :) = 2*(int1self + int2self) + penself
			stif(index1, isolateid, 2, :, :) = 0._dp
		endif
		
! 		1_left
		eleself = ele1
		he = hy
		index1 = 1
		isolateid = 2
		edgenode(1) = eleself%node(1)
		edgenode(2) = eleself%node(3)
		if ( mi .ne. 1 ) then
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			isolateid_n = 1
			eleneigh = get_ele(mi - 1, mj , 2)
			call get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
			stif(index1, isolateid, 1, :, :) = int1self + int2self + penself
			stif(index1, isolateid, 2, :, :) = int1neigh + int2neigh + penneigh
		else
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			stif(index1, isolateid, 1, :, :) = 2*(int1self + int2self) + penself
			stif(index1, isolateid, 2, :, :) = 0._dp		
		endif

! 		1_right
		eleself = ele1
		he = h
		index1 = 1
		isolateid = 3
		edgenode(1) = eleself%node(1)
		edgenode(2) = eleself%node(2)
		call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
		eleneigh = ele2
		isolateid_n = 2
		call get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
		stif(index1, isolateid, 1, :, :) = int1self + int2self + penself
		stif(index1, isolateid, 2, :, :) = int1neigh + int2neigh + penneigh

! 		2_right
		eleself = ele2
		he = hy
		index1 = 2
		isolateid = 1
		edgenode(1) = eleself%node(2)
		edgenode(2) = eleself%node(3)
		if ( mi .ne. mx ) then
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			isolateid_n = 2
			eleneigh = get_ele(mi + 1, mj, 1)
			call get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
			stif(index1, isolateid, 1, :, :) = int1self + int2self + penself
			stif(index1, isolateid, 2, :, :) = int1neigh + int2neigh + penneigh
		else
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			stif(index1, isolateid, 1, :, :) = 2*(int1self + int2self) + penself
			stif(index1, isolateid, 2, :, :) = 0._dp
		endif

! 		2_left
		eleself = ele2
		he = h
		index1 = 2
		isolateid = 2
		edgenode(1) = eleself%node(1)
		edgenode(2) = eleself%node(3)
		call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
		isolateid_n = 3
		eleneigh = ele1
		call get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
		stif(index1, isolateid, 1, :, :) = int1self + int2self + penself
		stif(index1, isolateid, 2, :, :) = int1neigh + int2neigh + penneigh

! 		2_bottom
		eleself = ele2
		he = hx
		index1 = 2
		isolateid = 3
		edgenode(1) = eleself%node(1)
		edgenode(2) = eleself%node(2)
		if ( mj .ne. 1 ) then
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			eleneigh = get_ele(mi, mj - 1, 1)
			isolateid_n = 1
			call get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
			stif(index1, isolateid, 1, :, :) = int1self + int2self + penself
			stif(index1, isolateid, 2, :, :) = int1neigh + int2neigh + penneigh
		else
			call get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
			stif(index1, isolateid, 1, :, :) = 2*(int1self + int2self) + penself
			stif(index1, isolateid, 2, :, :) = 0._dp
! 			stif(index1, isolateid, :, :, :) = 0._dp
		endif

		return
	end subroutine get_elestiff

	subroutine get_self(int1self, int2self, penself, he, eleself, index1, isolateid, edgenode)
		implicit none
		real(kind=dp), intent(out) :: int1self(3, 3), int2self(3, 3), penself(3, 3)
		real(kind=dp), intent(in) :: he
		type(element), intent(in) :: eleself
		integer(kind=ki), intent(in) :: index1, isolateid
		type(vector), intent(in) :: edgenode(2)

		real(kind=dp) :: integ(2, 3)
		real(kind=dp) :: coeffa, xh, yh
		real(kind=dp) :: vtemp1(2), vtemp2(2)
		integer(kind=ki) :: i, j, kpa, eid

		xh = (edgenode(2)%x - edgenode(1)%x)/n
		yh = (edgenode(2)%y - edgenode(1)%y)/n
		integ = 0._dp
		do i = 1, 3
			do kpa = 1, n
				coeffa = coef_int( &
				edgenode(1)%x + xh*(kpa - 1), edgenode(1)%y + yh*(kpa - 1), &
				edgenode(1)%x + xh*(kpa - 0), edgenode(1)%y + yh*(kpa - 0) )

				eid = edge_eleid(index1, isolateid, kpa)

				integ(1, i) = integ(1, i) + &
				eleself%grad(i, eid)*eleself%nvec_matrix(isolateid, isolateid)*coeffa
			end do
		end do
		int1self = - 1/2._dp*matmul( &
			reshape(eleself%gamma_matrix(:, isolateid), [3, 1]), &
			reshape(integ(1, :), [1, 3]) )
		int2self = transpose(int1self)
		if ( penalty_index .eq. 1 ) then
			penself = pcoef/h*he*matmul( &
				reshape(eleself%gamma_matrix(:, isolateid), [3, 1]), &
				reshape(eleself%gamma_matrix(:, isolateid), [1, 3]) )
		elseif ( penalty_index .eq. 2) then
			penself = 0._dp
			do i = 1, 3
			do j = 1, 3
				if ( i .eq. isolateid ) then
					vtemp1 = [0._dp, 0._dp]
				elseif ( i .eq. ntran(isolateid + 1) ) then
					vtemp1 = [1._dp, 0._dp]
				elseif ( i .eq. ntran(isolateid - 1) ) then
					vtemp1 = [0._dp, 1._dp]
				else
					print*, 'error: vtemp1'
				end if
				if ( j .eq. isolateid ) then
					vtemp2 = [0._dp, 0._dp]
				elseif ( j .eq. ntran(isolateid + 1) ) then
					vtemp2 = [1._dp, 0._dp]
				elseif ( j .eq. ntran(isolateid - 1) ) then
					vtemp2 = [0._dp, 1._dp]
				else
					print*, 'error: vtemp2'
				end if
				penself(i, j) = pcoef/h*edge_integ(he, vtemp1, vtemp2, 1)
			end do
			end do
		else	
			print*, 'error: get_elestiff, penalty_index'			
		end if

		return
	end subroutine get_self

	subroutine get_neigh(int1neigh, int2neigh, penneigh, he, eleself, eleneigh, index1, isolateid, isolateid_n, edgenode)
		implicit none
		real(kind=dp), intent(out) :: int1neigh(3, 3), int2neigh(3, 3), penneigh(3, 3)
		real(kind=dp), intent(in) :: he
		type(element), intent(in) :: eleself, eleneigh
		integer(kind=ki), intent(in) :: index1, isolateid, isolateid_n
		type(vector), intent(in) :: edgenode(2)

		real(kind=dp) :: integ(2, 3)
		real(kind=dp) :: coeffa, xh, yh
		real(kind=dp) :: vtemp1(2), vtemp2(2)
		integer(kind=ki) :: i, j, kpa, eid

		xh = (edgenode(2)%x - edgenode(1)%x)/n
		yh = (edgenode(2)%y - edgenode(1)%y)/n
		integ = 0._dp
		do i = 1, 3
			do kpa = 1, n
				coeffa = coef_int( &
				edgenode(1)%x + xh*(kpa - 1), edgenode(1)%y + yh*(kpa - 1), &
				edgenode(1)%x + xh*(kpa - 0), edgenode(1)%y + yh*(kpa - 0) )

				eid = edge_eleid(index1, isolateid, kpa)
				integ(1, i) = integ(1, i) + &
				eleself%grad(i, eid)*eleself%nvec_matrix(isolateid, isolateid)*coeffa

				eid = edge_eleid(3 - index1, isolateid_n, kpa)
				integ(2, i) = integ(2, i) + &
				eleneigh%grad(i, eid)*eleself%nvec_matrix(isolateid, isolateid)*coeffa
			end do
		end do
		int1neigh = - 1/2._dp*matmul( &
			reshape(eleself%gamma_matrix(:, isolateid), [3, 1]), &
			reshape(integ(2, :), [1, 3]) )
		int2neigh = 1/2._dp*matmul( &
			reshape(integ(1, :), [3, 1]), &
			reshape(eleneigh%gamma_matrix(:, isolateid_n), [1, 3]) )
		if ( penalty_index == 1 ) then
			penneigh = - pcoef/h*he*matmul( &
				reshape(eleself%gamma_matrix(:, isolateid), [3, 1]), &
				reshape(eleneigh%gamma_matrix(:, isolateid_n), [1, 3]) )
		elseif ( penalty_index == 2) then
			penneigh = 0._dp
			do i = 1, 3
			do j = 1, 3
				if ( i .eq. isolateid ) then
					vtemp1 = [0._dp, 0._dp]
				elseif ( i .eq. ntran(isolateid + 1) ) then
					vtemp1 = [1._dp, 0._dp]
				elseif ( i .eq. ntran(isolateid - 1) ) then
					vtemp1 = [0._dp, 1._dp]
				else
					print*, 'error: vtemp'
				end if
				if ( j .eq. isolateid_n ) then
					vtemp2 = [0._dp, 0._dp]
				elseif ( j .eq. ntran(isolateid_n + 1) ) then
					vtemp2 = [0._dp, 1._dp]
				elseif ( j .eq. ntran(isolateid_n - 1) ) then
					vtemp2 = [1._dp, 0._dp]
				else
					print*, 'error: vtemp'
				end if
				penneigh(i, j) = - pcoef/h*edge_integ(he, vtemp1, vtemp2, 1)
			end do
			end do
		else	
			print*, 'error: get_elestiff, penalty_index'			
		end if

		return		
	end subroutine get_neigh

	real(kind=dp) &
	function edge_integ(hlength, value1, value2, nn)
		implicit none
		real(kind=dp), intent(in) :: hlength, value1(nn + 1), value2(nn + 1)
		integer(kind=ki), intent(in) :: nn
		integer(kind=ki) :: i

		edge_integ = 0._dp
		do i = 1, nn
			edge_integ = edge_integ + 1/6._dp*hlength/nn*(&
			2*value1(i)*value2(i) + 2*value1(i + 1)*value2(i + 1) +&
			value1(i)*value2(i + 1) + value1(i + 1)*value2(i) )
		end do
		return
	end function edge_integ

	real(kind=dp) &
	function coef_int(x1, y1, x2, y2)
		implicit none
		real(kind=dp), intent(in) :: x1, y1, x2, y2
		integer(kind=ki) :: i
		real(kind=dp) :: xh, yh, htemp

		xh = x2 - x1
		yh = y2 - y1
		htemp = sqrt(xh**2 + yh**2)

		if ( edgeint_index .eq. 1 ) then
			coef_int = htemp/2*(coef(x1, y1) + coef(x2, y2))
		elseif ( edgeint_index .eq. 2) then
			coef_int = htemp/6*(coef(x1, y1) + coef(x2, y2) + 4*coef((x1 + x2)/2, (y1 + y2)/2))
		end if

		return
	end function coef_int

	subroutine test_integfunction
		implicit none
			
		integer(kind=ki) :: mi, mj, index1
		real(kind=dp) :: stif(3, 3)
		real(kind=dp) :: stife(2, 3, 2, 3, 3)

		print*, '----- test integfunction -----'

		mi = 2
		mj = 2
		index1 = 2


		call element_integ(stif, mi, mj, index1)
		call get_elestiff(stife, mi, mj)

		print*
		print*, 'mi, mj, index1: ', mi, mj, index1
		print*

		print*, 'element_integ'
		print*, stif
		print*

		print*, '1_up_self'
		print*, stife(1, 1, 1, :, :)
		print*

		print*, '1_up_neigh'
		print*, stife(1, 1, 2, :, :)
		print*

		print*, '1_left_self'
		print*, stife(1, 2, 1, :, :)
		print*

		print*, '1_left_neigh'
		print*, stife(1, 2, 2, :, :)
		print*

		print*, '1_right_self'
		print*, stife(1, 3, 1, :, :)
		print*

		print*, '1_right_neigh'
		print*, stife(1, 3, 2, :, :)
		print*

		print*, '2_right_self'
		print*, stife(2, 1, 1, :, :)
		print*

		print*, '2_right_neigh'
		print*, stife(2, 1, 2, :, :)
		print*

		print*, '2_left_self'
		print*, stife(2, 2, 1, :, :)
		print*

		print*, '2_left_neigh'
		print*, stife(2, 2, 2, :, :)
		print*

		print*, '2_bottom_self'
		print*, stife(2, 3, 1, :, :)
		print*

		print*, '2_bottom_neigh'
		print*, stife(2, 3, 2, :, :)
		print*

		return		
	end subroutine test_integfunction

end module integfunction