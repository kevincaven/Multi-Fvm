module solveequation
 	use constants, only: ki, dp, pi 
 	use domainmesh, only: mx, my, nel
	use getlgmap, only: kdof, kidof, kbdof, lgmap, cor, get_elid
	use prob, only: ub
	use openbasis, only: p, n, nz
 	use assembstiff, only: a, drhs0, row_el, pointerarray

	implicit none
	
	private
	public solve_equation, &
	x, y, xspan, yspan

	real(kind=dp), allocatable :: drhs(:)

	integer(kind=ki), allocatable :: ida(:), jda(:), idb(:), jdb(:)
	real(kind=dp), allocatable :: vda(:), vdb(:)
	real(kind=dp), allocatable :: ax(:)
	integer(kind=ki) :: limax, lbmax

	real(kind=dp), allocatable :: x(:), y(:)
	real(kind=dp), allocatable :: xspan(:), yspan(:)

	integer(kind=ki), allocatable :: nai(:), nap(:),index1(:)

contains

	subroutine solve_equation()
		implicit none
		
		print*, '----- step7: solve the equation -----'
		print*
		call pretreat()
		call umfpart()
		call check_info()

		return
	end subroutine solve_equation

	subroutine pretreat()
		implicit none
		integer(kind=ki) :: kimax, kbmax
		integer(kind=ki) :: mrow, mcol, k
		type(row_el), pointer :: head

		kimax = kidof*12
		kbmax = kbdof*12

		allocate(vda(kimax), ida(kimax), jda(kimax))
		allocate(vdb(kbmax), idb(kbmax), jdb(kbmax))
		
		limax = 0
		lbmax = 0

		do mrow=kbdof + 1, kdof
			head=>a(mrow)%el
			do while (associated(head))
				mcol = head%col
				if(mcol .gt. kbdof) then
					limax = limax + 1
					ida(limax) = mrow - kbdof
					jda(limax) = mcol - kbdof
					vda(limax) = head%value
				else
					lbmax = lbmax + 1
					idb(lbmax) = mrow - kbdof
					jdb(lbmax) = mcol 
					vdb(lbmax) = head%value
				end if
				head => head%next
			end do
		end do

		allocate(drhs(kidof))
		do k=1, kidof
			drhs(k) = drhs0(k + kbdof)
		enddo

		do k=1, lbmax
			mrow = idb(k)
			mcol = jdb(k)
			drhs(mrow) = drhs(mrow) - vdb(k)*ub(mcol)
		enddo

		open(11, file='data_out/ijvda.dat')
		do k=1, limax
			write(11,*) ida(k), jda(k), vda(k)
		enddo
		close(11)

		open(12, file='data_out/ijvdb.dat')
		do k=1, lbmax
			write(12,*) idb(k), jdb(k), vdb(k)
		enddo
		close(12)

		open(13, file='data_out/drhs.dat')  
		do k=1, kidof
			write(13,*) k, drhs(k)
		enddo
		close(13)

		return
	end subroutine pretreat

	subroutine umfpart()
		implicit none
		integer(kind=ki) :: k
		integer(kind=ki) :: mi, mj, index, elid, node, nid
		integer(kind=ki) :: np(3)
		
! 		x(:) umf solution, y(:) numerical solution

		allocate(x(kidof))
		allocate(ax(kidof*12))
		allocate(nai(kidof*12))
		allocate(nap(kidof + 1))
		allocate(index1(kidof))

		print*, 'go to umfpack'
		call umfdglr(kidof, drhs, ida, jda, vda, x, ax, nap, nai, index1, limax, 0)
		print*, 'out of umfpack'
		print*

		open(13, file='data_out/umfsolution.dat')
		do k=1, kidof
			write(13,*) k, x(k)
		end do 
		close(13)

		allocate(xspan(kdof))

		do k = 1, kbdof
			xspan(k) = ub(k)
		end do
		do k = 1, kidof
			xspan(k + kbdof) = x(k)
		end do

		allocate(yspan(kdof))

		do mj = 1, my
		do mi = 1, mx
			index = 1
			elid = get_elid(mi, mj, index)
			np = [1, nz, nz - n]
			do node = 1, 3
				yspan(lgmap(elid, node)) = &
				xspan(lgmap(elid, 1))*p(index, 1, mi, mj, np(node)) + &
				xspan(lgmap(elid, 2))*p(index, 2, mi, mj, np(node)) + &
				xspan(lgmap(elid, 3))*p(index, 3, mi, mj, np(node))
			end do

			index = 2
			elid = get_elid(mi, mj, index)
			np = [1, n + 1, nz]
			do node = 1, 3
				yspan(lgmap(elid, node)) = &
				xspan(lgmap(elid, 1))*p(index, 1, mi, mj, np(node)) + &
				xspan(lgmap(elid, 2))*p(index, 2, mi, mj, np(node)) + &
				xspan(lgmap(elid, 3))*p(index, 3, mi, mj, np(node))
			end do
		end do
		end do
		
		allocate(y(kidof))
 		y = yspan(kbdof + 1: kdof)

!  		save the solution data of the equation set
 		open(15, file='data_out/numsolution.dat')
 		do k=1, kidof
 			write(15, *) k, y(k)
 		end do 
 		close(15)
!  		save the solution data as  a plot form
 		open(14, file='data_out/numsol_plotdata.dat')
 		do k = 1, nel
 			do node = 1, 3
 				nid = lgmap(k, node)
 				if (nid .le. kbdof) then
 					write(14,*) cor(nid, 1), cor(nid, 2), ub(nid)
 				else
 					write(14,*) cor(nid, 1), cor(nid, 2), y(nid - kbdof)
 				end if
 			end do 
 		end do  
 		close(14)

 		return
	end subroutine umfpart

	subroutine check_info()
		implicit none
		
		print*, 'ub & drhs: '
		print*, 'data_out/ub.dat (ub)'
		print*, 'data_out/drhs.dat (drhs(kidof))'
		print*

		print*, 'equation infomation: '
 		print*,'the number of the stiff matrix non-zero elements:', limax
 		print*, 'please check: '
		print*, 'data_out/ijvda (ida, jda, vda)'
 		print*

		return		
	end subroutine check_info

end module solveequation