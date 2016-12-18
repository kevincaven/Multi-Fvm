module getlgmap
 	
 	use constants, only: ki, dp, pi 
	use domainmesh, only: mx, my, nel, hx, hy
	implicit none
	
	private
	public get_elid, get_lgmap, &
	kdof, kidof, kbdof, lgmap, cor

	integer(kind=ki), allocatable :: pointcls(:, :)
	integer(kind=ki), allocatable :: lgmap(:, :)
	real(kind=dp), allocatable :: cor(:, :)
	integer(kind=ki) :: kidof, kbdof, kdof

contains

	subroutine get_lgmap()
		implicit none

		print*, '----- step2: get lgmap -----'
		print*
		call create_pointcls()
		call create_lgmap()
		call check_info()
! 		call test()

		return
	end subroutine get_lgmap

	integer(kind=ki) &
	function get_elid(mi, mj, index1) result(elid)
		implicit none
		integer(kind=ki), intent(in) :: mi, mj, index1

		elid = ( mj - 1 )*2*mx + 2*(mi - 1) + index1

		return
	end function get_elid

	subroutine create_pointcls()
		implicit none
		
		allocate(pointcls(mx + 1, my + 1))
		pointcls = 0
		pointcls(1, 1) = 1
		pointcls(mx + 1, 1) = 2
		pointcls(1, my+1) = 3
		pointcls(mx + 1, my + 1) = 4
		pointcls(2:mx, 1) = 5
		pointcls(2:mx, my + 1) = 6
		pointcls(1, 2:my) = 7
		pointcls(mx + 1, 2:my) = 8
		
		return
	end subroutine create_pointcls

	subroutine create_lgmap()
		implicit none

! 		useful output: 
! 		lgmap(nel, 2), cor(nel*3,2), kidof, kbdof, kdof
		integer(kind=ki) :: i, j
		integer(kind=ki) :: ktemp, elid
		integer(kind=ki) :: kbmax 
		real(kind=dp) :: xcor, ycor

		allocate(lgmap(nel, 3))
		lgmap = 0
		allocate(cor(nel*3, 2))
		cor = 0.0_dp

		kbdof = 0
		kbmax = 3*(mx + my - 1) * 2
		kidof = kbmax

		do j = 1, mx + 1
			do i = 1, my + 1
				xcor = hx*(i - 1)
				ycor = hy*(j - 1)

				select case (pointcls(i,j))
					case (1)
! 						left-bottom: i = 1 & j = 1
! 						2 freedoms
						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(1, 1) = ktemp
						cor(ktemp, 1) = 0.0_dp
						cor(ktemp, 2) = 0.0_dp

						kbdof  = kbdof + 1
						ktemp = kbdof
						lgmap(2, 1) = ktemp
						cor(ktemp, 1) = 0.0_dp
						cor(ktemp, 2) = 0.0_dp
						
					case (2)
! 						right-bottom: i = mx + 1 & j = 1
! 						1 freedom
						elid = get_elid(i - 1, j, 2)

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 2) = ktemp					
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case(3)
! 						left-up: j = mx + 1 & i = 1
! 						1 freedom
						elid = get_elid(i, j - 1,1)

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case (4)
! 						right-up: i = mx+1 & j = my+1
!						2 freedoms
						elid = get_elid(i - 1, j - 1, 1)
						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 2) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid + 1, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case(5)
! 						bottom: i = 2:mx, j = 1
! 						3 freedoms
						elid = get_elid(i, j, 1)

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid - 1, 2) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
				 		lgmap(elid, 1) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid + 1, 1) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case(6)
! 						up: i = 2:mx, j = my+1
! 						3 freedoms
						elid = get_elid(i, j - 1, 1)

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid - 2, 2) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
				 		lgmap(elid - 1, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case(7)
! 						left: i = 1, j = 2:my
! 						3 freedoms
						elid = get_elid(i, j - 1, 1)
						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						elid = get_elid(i, j, 1)
						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 1) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid + 1, 1) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case(8)
! 						right: i = mx+1, j = 2:my
! 						3 freedoms
						elid = get_elid(i - 1, j - 1, 1)
						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 2) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid + 1, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						elid = get_elid(i - 1, j, 2)
						kbdof = kbdof + 1
						ktemp = kbdof
						lgmap(elid, 2) = ktemp 
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case (0)
! 						center: i = 2:mx, j = 2:my
! 						6 freedoms
						elid = get_elid(i - 1,j - 1, 1)
						kidof = kidof + 1
						ktemp = kidof
						lgmap(elid, 2) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor
						
						kidof = kidof + 1
						ktemp = kidof
						lgmap(elid + 1, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kidof = kidof + 1
						ktemp = kidof
						lgmap(elid + 2, 3) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						elid = get_elid(i, j, 1)
						kidof = kidof + 1
						ktemp = kidof
						lgmap(elid - 1, 2) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kidof = kidof + 1
						ktemp = kidof
						lgmap(elid, 1) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

						kidof = kidof + 1
						ktemp = kidof
						lgmap(elid + 1, 1) = ktemp
						cor(ktemp, 1) = xcor
						cor(ktemp, 2) = ycor

					case default
						print*, "lgmap error"
						exit
				end select
			end do
		end do

		kdof = kidof
		kidof = kdof - kbdof 

		return
	end subroutine create_lgmap

	subroutine check_info()
		implicit none
		
		print*, 'freedom infomation: '
		print*, 'total number of the freedom is:', kdof
		print*, 'total number of the internal freedom is: ', kidof
		print*, 'total number of the border freedom is:', kbdof
		print*
		
		return
	end subroutine check_info

	subroutine test()
		implicit none
		integer(ki) :: i, j

		print*
		print*, 'lgmap'
		do i = 1, nel
			print*, i, (lgmap(i, j), j=1, 3)
		enddo

		print*, 'cor'
		do i = 1, nel
			do j = 1, 3
				print*, i, j, cor(lgmap(i, j),1), cor(lgmap(i,j), 2)					
			enddo
		enddo

		return
	end subroutine test

end module getlgmap