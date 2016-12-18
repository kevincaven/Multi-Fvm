module assembstiff

 	use constants, only: ki, dp, pi 
	use domainmesh, only: mx, my, nel, hx, hy
	use getlgmap, only: kdof, kidof, kbdof, get_elid, lgmap
	use eleinfo, only: test_eleinfo
	use integfunction
	
	private
	public assemb_stiff, &
	row_el, pointerarray, &
	a, drhs0

	type row_el
		integer :: col
		real(kind=dp) :: value
		type(row_el), pointer :: next
	end type row_el
	
	type pointerarray
		type(row_el), pointer :: el
	end type pointerarray

	type(pointerarray), allocatable :: a(:)
	real(kind=dp), allocatable :: drhs0(:)

contains

	subroutine assemb_stiff()
		implicit none
		
		print*, '----- step6: assemble stiff matrix -----'
		print*
		call init_a()
		call init_drhs()
		call update_a()
		call update_drhs()
		call check_info()
! 		call test()

		return
	end subroutine assemb_stiff

	subroutine init_a()
		implicit none
		integer(kind=ki) :: k

		allocate(a(kdof))
		do k=1, kdof
			allocate(a(k)%el)
			a(k)%el%col = k
			a(k)%el%value = 0._dp
			a(k)%el%next=>null()
		end do

		return
	end subroutine init_a

	subroutine init_drhs()
		implicit none
		integer(kind=ki) :: k
		
		allocate(drhs0(kdof))
		do k=1, kdof
			drhs0(k) = 0._dp
		enddo

		return
	end subroutine init_drhs

	subroutine addvalue(head, coln, value)
		implicit none
		type(row_el), pointer :: head, p, q
		integer(kind=ki) :: coln, is_exsit, allstat
		real(kind=dp) :: value
		
		p=>head
		q=>head    
		is_exsit=0

		do while (associated(p))
			if(p%col .eq. coln) then
				p%value = p%value + value
				is_exsit = 1
				exit
			else
				q => p
				p => p%next
			endif
		end do

		if(is_exsit .eq. 0) then
			allocate(p, stat = allstat)
			if(allstat .ne. 0) print*, 'Wrong allocate in arrange'
			p%col = coln
			p%value = value
			p%next => null()
			q%next => p
		endif
		
		return
	end subroutine addvalue

	subroutine assemb_elestiff(elid1, elid2 , stif)
		implicit none
		integer(kind=ki), intent(in) :: elid1, elid2
		real(kind=dp), intent(in) :: stif(3, 3)
		integer(kind=ki) :: testb, solb, mrow, mcol
		real(kind=dp) :: value

		do testb = 1, 3
			mrow = lgmap(elid1, testb)
			do solb = 1, 3
				mcol = lgmap(elid2, solb)
				call addvalue(a(mrow)%el, mcol, stif(testb, solb))
			end do 
		end do 
				
		return		
	end subroutine assemb_elestiff

	subroutine update_a()
		implicit none
		integer(kind=ki) :: mi,mj, elid, elid1, elid2
		real(kind=dp) :: stif(2, 3, 2, 3, 3)
		real(kind=dp) :: stif1(3, 3), stif2(3, 3)
		integer(kind=ki) :: testb, solb, mrow, mcol
		real(kind=dp) :: value
		real(kind=dp) :: start, finish
		
		print*, 'assemble stiff matrix a, please wait ... ... '
		call cpu_time(start)

		do mj = 1, my
		do mi = 1, mx
! 			element integ
			elid1 = get_elid(mi, mj, 1)
			elid2 = get_elid(mi, mj, 2)
			call element_integ(stif1, mi, mj, 1)
			call assemb_elestiff(elid1, elid1, stif1)
			call element_integ(stif2, mi, mj, 2)
			call assemb_elestiff(elid2, elid2, stif2)

! 			edge integ
			call get_elestiff(stif, mi, mj)

! 			upside
			elid1 = get_elid(mi, mj, 1) ! self
			stif1 = stif(1, 1, 1, :, :)
			call assemb_elestiff(elid1, elid1, stif1)

			if(mj .ne. my) then
			elid2 = get_elid(mi, mj + 1, 2) ! neigh
			stif2 = stif(1, 1, 2, :, :)
			call assemb_elestiff(elid1, elid2, stif2)
			end if 

! 			leftside
			elid1 = get_elid(mi, mj, 1) ! self
			stif1 = stif(1, 2, 1, :, :)
			call assemb_elestiff(elid1, elid1, stif1)

			if(mi .ne. 1) then
			elid2 = get_elid(mi - 1, mj, 2) ! neigh
			stif2 = stif(1, 2, 2, :, :)
			call assemb_elestiff(elid1, elid2, stif2)
			end if 

!			1self to 2neigh
			elid1 = get_elid(mi, mj, 1) ! self
			stif1 = stif(1, 3, 1, :, :)
			call assemb_elestiff(elid1, elid1, stif1)

			elid2 = get_elid(mi, mj, 2) ! neigh
			stif2 = stif(1, 3, 2, :, :)
			call assemb_elestiff(elid1, elid2, stif2)

! 			rightside
			elid1 = get_elid(mi, mj, 2) ! self
			stif1 = stif(2, 1, 1, :, :)
			call assemb_elestiff(elid1, elid1, stif1)

			if(mi .ne. mx) then
			elid2 = get_elid(mi + 1, mj, 1) ! neigh
			stif2 = stif(2, 1, 2, :, :)
			call assemb_elestiff(elid1, elid2, stif2)
			end if          

! 			2self to 1neigh
			elid1 = get_elid(mi, mj, 2) ! self
			stif1 = stif(2, 2, 1, :, :)
			call assemb_elestiff(elid1, elid1, stif1)

			elid2 = get_elid(mi, mj, 1) ! neigh
			stif2 = stif(2, 2, 2, :, :)
			call assemb_elestiff(elid1, elid2, stif2)

! 			bottomside
			elid1 = get_elid(mi, mj, 2) ! self
			stif1 = stif(2, 3, 1, :, :)
			call assemb_elestiff(elid1, elid1, stif1)

			if(mj .ne. 1) then
			elid2 = get_elid(mi, mj - 1, 1) ! neigh
			stif2 = stif(2, 3, 2, :, :)
			call assemb_elestiff(elid1, elid2, stif2)
			end if          
		end do
		end do

		call cpu_time(finish)
		print*, 'time of assemble stiff matrix: '
		print'(f6.3," seconds")', finish - start
		print*

		return
	end subroutine update_a

	subroutine update_drhs()
		implicit none
		integer(kind=ki) :: mi, mj, elid1, elid2, k
		integer(kind=ki) :: testb
		real(kind=dp) :: start, finish

		print*, 'assemble right hand item, please wait ... ... '
		call cpu_time(start)
! 		right hand item
		do mj = 1, my
		do mi = 1, mx
			elid1 = get_elid(mi, mj, 1)
			elid2 = get_elid(mi, mj, 2)
			do testb = 1, 3
				drhs0(lgmap(elid1, testb)) = drhs0(lgmap(elid1, testb)) + &
				righthand_integ(mi, mj, 1, testb)
				drhs0(lgmap(elid2, testb)) = drhs0(lgmap(elid2, testb)) + &
				righthand_integ(mi, mj, 2, testb)
			end do
		enddo
		enddo

		open(11, file='data_out/drhs0.dat')
		do k=1,kdof
			write(11,*) k, drhs0(k)
		enddo
		close(11)

		call cpu_time(finish)
		print*, 'time of assemble right hand item: '
		print'(f6.3," seconds")', finish - start
		print*

		return
	end subroutine update_drhs

	subroutine check_info()
		implicit none
		
		print*, 'right hand item infomation: '
		print*, 'please check: '
		print*, 'data_out/drhs0.dat (drhs0(kdof))'
		print*

		return		
	end subroutine check_info

	subroutine test()
		implicit none

		call test_integfunction

! 		print*, righthand_integ(1, 1, 1, 1)

		return
	end subroutine test

end module assembstiff
	