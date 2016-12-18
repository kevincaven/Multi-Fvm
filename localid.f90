module localid
	use vectors
	use constants, only: ki, dp, pi
	use domainmesh, only: mx, my, hx, hy, xleft, xright, ybottom, ytop
	use openbasis, only: n, hhx, hhy
	implicit none
	
	private
	public local_eleid, bf_nid1, bf_eleid, bf_nid2, &
	segn, edge_nid, edge_eleid, pa_eleid, corn, corm

	integer(kind=ki) :: segn

	integer(kind=ki), allocatable :: edge_nid(:, :, :) ! edge_nid(2, 3, n + 1)

	integer(kind=ki), allocatable :: edge_eleid(:, :, :) ! edge_eleid(2, 3, n)
	integer(kind=ki), allocatable :: pa_eleid(:, :, :) ! pa_eleid(2, 3, segn + 2)

	type(vector), allocatable :: corn(:, :, :, :, :) ! corn(mx, my, 2, 3, segn + 1)
	type(vector), allocatable :: corm(:, :, :, :, :) ! corm(mx, my, 2, 3, segn + 1)

contains

	subroutine local_eleid()
		implicit none
		
		print*, '----- step5: get edge element local id -----'
		print*
		call getedge_nid()
		call getedge_eleid()
		call getpa_eleid()
		call getpa_cor()
		call check_info()
		call test()

		return		
	end subroutine local_eleid

	subroutine check_info()
		implicit none
		
		print*, 'each edge has ', n, 'refined elements'
		print*, 'each pa has ', segn , 'complete refined elements', 1, 'incomplete refined element '
		print*

		return		
	end subroutine check_info

	subroutine test()
		implicit none
		integer(kind=ki) :: loop, index1, node, ni, nj, mi, mj

		index1 = 2
		node = 3

! !		test pa_eleid
! 		do loop = 1, segn + 1
! 			print*, loop, pa_eleid(1, 1, loop), pa_eleid(1, 2, loop), pa_eleid(1, 3, loop)
! 		end do
! 		print*
! 		do loop = 1, segn + 1
! 			print*, loop, pa_eleid(2, 1, loop), pa_eleid(2, 2, loop), pa_eleid(2, 3, loop)
! 		end do
! 		print*

! 		test corn, corm
		mi = 32
		mj = 32
		open(11, file='figure/testcorm.dat')
		do index1 = 1, 2
			do node = 1, 3
				do loop = 1, segn + 1
					write(11, *) loop, corm(mi, mj, index1, node, loop)
					write(11, *) loop, corn(mi, mj, index1, node, loop)
				end do
			end do		
		end do
		close(11)

		open(12, file='figure/testcorn.dat')
		do index1 = 1, 2
			do node = 1, 3
				do loop = 1, segn + 1
					write(12, *) loop, corn(mi, mj, index1, node, loop)
				end do
			end do		
		end do
		close(12)

		return
	end subroutine test

	integer(kind=ki)  &
	function bf_nid1(index1, ni, nj) result(nid)
!      	ni = 1, n + 1; nj = 1, n + 1
		implicit none
		integer(kind=ki), intent(in) :: index1, ni, nj

		if ( index1 == 1 ) then
			if ( ni .gt. nj ) print*, 'error: bf_nid1, ij'
			nid = nj*(nj - 1)/2 + ni
		elseif (index1 == 2) then
			if ( ni .lt. nj ) print*, 'error: bf_nid1, ij'
			nid = (nj - 1)*(n + 1 + n + 1 - nj + 2)/2 + ni - nj + 1
		else
			print*, 'error: bf_nid1, index1' 	
		end if 

		return
	end function bf_nid1

	integer(kind=ki) &
	function bf_eleid(index1, ni, nj, index2) result(eleid)
! 		ni = 1, n; nj = 1, n
		implicit none
		integer(kind=ki), intent(in) :: index1, index2, ni, nj

		select case (index1)
			case (1)
				if (ni .gt. nj) print*, 'error: bf_eleid, ni, nj'
				if ( (ni .eq. nj) .and. (index2 .eq. 2) ) print*, 'error: bf_eleid, index2' 
				eleid = (nj - 1)**2 + (ni - 1)*2 + index2
			case (2)
				if (ni .lt. nj) print*, 'error: bf_eleid, ni, nj'
				if ( (ni .eq. nj) .and. (index2 .eq. 1) ) print*, 'error: bf_eleid, index2' 
				eleid = n**2 - (n - nj + 1)**2 + (ni - nj)*2 + index2 - 1
			case default
				print*, 'error: bf_eleid, index1'
		end select

		return		
	end function bf_eleid

	integer(kind=ki) &
	function bf_nid2(index1, ni, nj, index2, node) result(nid)
! 		ni = 1, n; nj = 1, n
		implicit none
		integer(kind=ki), intent(in) :: ni, nj, index1, index2, node

		if ( index1 == 1 ) then
			if ( ni .gt. nj ) print*, 'error: bf_nid2, index1'
			if ( index2 == 1 ) then
				select case (node)
					case (1)
						nid=(nj - 1)*nj/2 + ni
					case (2)
						nid=nj*(nj + 1)/2 + ni + 1
					case (3)
						nid=nj*(nj + 1)/2 + ni
					case default
						print*, 'error: bf_nid2, node'
				end select
			elseif (index2 == 2) then
				if ( ni .ge. nj ) print*, 'error: bf_nid2, ni-nj'
				select case (node)
					case (1)
						nid=(nj - 1)*nj/2 + ni
					case (2)
						nid=(nj - 1)*nj/2 + ni + 1
					case (3)
						nid=nj*(nj + 1)/2 + ni + 1
					case default
						print*, 'error: bf_nid2, node'
				end select
			else
				print*, 'error: bf_nid2, index2'
			end if
		elseif (index1 == 2) then
			if ( ni .lt. nj ) print*, 'error: bf_nid2, ni-nj'
			if ( index2 == 1 ) then
				if ( ni .le. nj ) print*, 'error: bf_nid2, ni-nj' 
				select case (node)
					case (1)
						nid=(nj - 1)*(n + 1 + n + 1 - nj + 2)/2 + ni - nj + 1
					case (2)
						nid=nj*(n + 1 + n + 1 - nj + 1)/2 + ni - nj + 1
					case (3)
						nid=nj*(n + 1 + n + 1 - nj + 1)/2 + ni - nj
					case default
						print*, 'error: bf_nid2, node'
				end select
			elseif (index2 == 2) then
				select case (node)
					case (1)
						nid=(nj - 1)*(n + 1 + n + 1 - nj + 2)/2 + ni - nj + 1
					case (2)
						nid=(nj - 1)*(n + 1 + n + 1 - nj + 2)/2 + ni - nj + 2
					case (3)
						nid=nj*(n + 1 + n + 1 - nj + 1)/2 + ni - nj + 1
					case default
						print*, 'error: bf_nid2, node'
				end select
			else
				print*, 'error: bf_nid2, index2'
			end if
		else
			print*, 'error: bf_nid2, index1 '
		end if

		return
	end function bf_nid2

	subroutine getedge_nid()
		implicit none
		integer(kind=ki) :: index1, index2, i, node

		allocate(edge_nid(2, 3, n + 1))
		
		do i = 1, n + 1
			index1 = 1
			edge_nid(index1, 1, i) = bf_nid1(index1, i, n + 1)
			edge_nid(index1, 2, i) = bf_nid1(index1, 1, i) 
			edge_nid(index1, 3, i) = bf_nid1(index1, i, i) 

			index1 = 2
			edge_nid(index1, 1, i) = bf_nid1(index1, n + 1, i) 
			edge_nid(index1, 2, i) = bf_nid1(index1, i, i) 
			edge_nid(index1, 3, i) = bf_nid1(index1, i, 1) 
		end do

		return		
	end subroutine getedge_nid

	subroutine getedge_eleid()
		implicit none
		integer(kind=ki) :: index1, index2, i, node

		allocate(edge_eleid(2, 3, n))

		do i = 1, n
			index1 = 1
			index2 = index1
			edge_eleid(index1, 1, i) = bf_eleid(index1, i, n, index2)
			edge_eleid(index1, 2, i) = bf_eleid(index1, 1, i, index2)
			edge_eleid(index1, 3, i) = bf_eleid(index1, i, i, index2)

			index1 = 2
			index2 = index1
			edge_eleid(index1, 1, i) = bf_eleid(index1, n, i, index2)
			edge_eleid(index1, 2, i) = bf_eleid(index1, i, i, index2)
			edge_eleid(index1, 3, i) = bf_eleid(index1, i, 1, index2)
		end do

		return		
	end subroutine getedge_eleid

	subroutine getpa_eleid()
		implicit none
		integer(kind=ki) :: i, index1, ni, nj, index2, node

		segn = floor(2._dp*n/3)

		allocate(pa_eleid(2, 3, segn + 2))

! 		id1, pa1
		index1 = 1
		node = 1
		index2 = 3 - index1
		i = 0
		ni = 0
		nj = 0
		do
			if ( i .ge. (segn + 1) ) exit
			i = i + 1; ni = ni + 1; nj = nj + 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
			i = i + 1; ni = ni; nj = nj + 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
		end do
! 		id1, pa2
		index1 = 1
		node = 2
		index2 = 3 - index1
		i = 0
		ni = n + 1
		nj = n + 1
		do
			if ( i .ge. (segn + 1) ) exit
			i = i + 1; ni = ni - 1; nj = nj - 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
			i = i + 1; ni = ni - 1; nj = nj; index2 = 3 -index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
		end do
! 		id1, pa3
		index1 = 1
		node = 3
		index2 = 3 - index1
		i = 0
		ni = 0
		nj = n + 1
		do
			if ( i .ge. (segn + 1) ) exit
			i = i + 1; ni = ni + 1; nj = nj - 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
			i = i + 1; ni = ni; nj = nj; index2 = 3 -index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
		end do

! 		id2, pa1
		index1 = 2
		node = 1
		index2 = 3 - index1
		i = 0
		ni = 0
		nj = 0
		do
			if ( i .ge. (segn + 1) ) exit
			i = i + 1; ni = ni + 1; nj = nj + 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
			i = i + 1; ni = ni + 1; nj = nj; index2 = 3 -index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
		end do
! 		id2, pa2
		index1 = 2
		node = 2
		index2 = 3 - index1
		i = 0
		ni = n + 1
		nj = 0
		do
			if ( i .ge. (segn + 1) ) exit
			i = i + 1; ni = ni - 1; nj = nj + 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
			i = i + 1; ni = ni; nj = nj; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
		end do
! 		id, pa3
		index1 = 2
		node = 3
		index2 = 3 - index1
		i = 0
		ni = n + 1
		nj = n + 1
		do
			if ( i .ge. (segn + 1) ) exit
			i = i + 1; ni = ni - 1; nj = nj - 1; index2 = 3 - index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
			i = i + 1; ni = ni; nj = nj - 1; index2 = 3 -index2
			pa_eleid(index1, node, i) = bf_eleid(index1, ni, nj, index2)
		end do

		return	
	end subroutine getpa_eleid

	subroutine getpa_cor()
		implicit none
		integer(kind=ki) :: mi, mj, index1, ni, nj, index2, node
		real(kind=dp) :: xx, yy, xtemp, ytemp
		integer(kind=ki) :: loop

		allocate(corn(mx, my, 2, 3, segn + 1))
		allocate(corm(mx, my, 2, 3, segn + 1))

		do mj = 1, my
		do mi = 1, mx
			xx = xleft + (mi - 1)*hx
			yy = ybottom + (mj - 1)*hy


! 			id1, pa1
			index1 = 1
			node = 1
			index2 = 3 - index1
			loop = 0
			ni = 0
			nj = 0
			do
				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni + 1; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*nj

				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
			enddo
			! segn + 1
			if (index2 .ne. index1 ) then
				loop = loop + 1; ni = ni + 1; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 2/3._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/3._dp)
			else
				loop = loop + 1; ni = ni; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/3._dp)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 2/3._dp)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
			endif

! 			id1, pa2
			index1 = 1
			node = 2
			index2 = 3 - index1
			loop = 0
			ni = n + 1
			nj = n + 1
			do
				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni - 1; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)

				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni - 1; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)
			enddo
			! segn + 1
			if (index2 .ne. index1 ) then
				loop = loop + 1; ni = ni - 1; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 2/3._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/3._dp)
			else
				loop = loop + 1; ni = ni - 1; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/3._dp)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 2/3._dp)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)
			endif

! 			id1, pa3
			index1 = 1
			node = 3
			index2 = 3 - index1
			loop = 0
			ni = 0
			nj = n + 1
			do
				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni + 1; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/2._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/2._dp)

				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/2._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/2._dp)
			enddo
			! segn + 1
			if (index2 .ne. index1 ) then
				loop = loop + 1; ni = ni + 1; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 2/3._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/3._dp)
			else
				loop = loop + 1; ni = ni; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/3._dp)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 2/3._dp)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/2._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/2._dp)
			endif

! 			id2, pa1
			index1 = 2
			node = 1
			index2 = 3 - index1
			loop = 0
			ni = 0
			nj = 0
			do
				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni + 1; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)

				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni + 1; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)
			enddo
			! segn + 1
			if (index2 .ne. index1 ) then
				loop = loop + 1; ni = ni + 1; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/3._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 2/3._dp)
			else
				loop = loop + 1; ni = ni + 1; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 2/3._dp)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/3._dp)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)
			endif

! 			id2, pa2
			index1 = 2
			node = 2
			index2 = 3 - index1
			loop = 0
			ni = n + 1
			nj = 0
			do
				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni - 1; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)

				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)
			enddo
			! segn + 1
			if (index2 .ne. index1 ) then
				loop = loop + 1; ni = ni - 1; nj = nj + 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/3._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 2/3._dp)
			else
				loop = loop + 1; ni = ni; nj = nj; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 2/3._dp)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/3._dp)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 0.5_dp)
			endif			

! 			id2, pa3
			index1 = 2
			node = 3
			index2 = 3 - index1
			loop = 0
			ni = n + 1
			nj = n + 1
			do
				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni - 1; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)

				if ( loop .ge. segn) exit
				loop = loop + 1; ni = ni; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*nj
			enddo
			! segn + 1
			if (index2 .ne. index1 ) then
				loop = loop + 1; ni = ni - 1; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*ni
				corn(mi, mj, index1, node, loop)%y = yy + hhy*nj
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 1/3._dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 2/3._dp)
			else
				loop = loop + 1; ni = ni; nj = nj - 1; index2 = 3 - index2
				corn(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 2/3._dp)
				corn(mi, mj, index1, node, loop)%y = yy + hhy*(nj - 1/3._dp)
				corm(mi, mj, index1, node, loop)%x = xx + hhx*(ni - 0.5_dp)
				corm(mi, mj, index1, node, loop)%y = yy + hhy*nj
			endif		

		end do
		end do

		return
	end subroutine getpa_cor

end module localid
