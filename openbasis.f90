module openbasis
	use toggles
	use constants
	use domainmesh, only: mx, my, hx, hy
	use toggles, only: bf_index
	implicit none
	
	private
	public open_basis, &
	p, n, nx, ny, nz, hhx, hhy

	integer(kind=ki) :: n, nx, ny, nz
	real(kind=dp) :: hhx, hhy
	real(kind=dp), allocatable :: p(:, :, :, :, :)

contains

	subroutine open_basis()
		implicit none

		print*, '----- step4: open multi-basis function -----'
		print*
		call read_basis()
		call check_info()
! 		call test()
		
		return		
	end subroutine open_basis

	subroutine read_basis()
		implicit none
		integer(kind=ki) :: mi, mj, k, nxyf
		real(kind=dp) :: start, finish
		character(len=100) :: datapath, basepath, npath, egpath

		write(npath, *) finemesh_index
		write(basepath, *) finestmesh_index 
		write(egpath, *) eg_index 

		datapath = '../basisfunction/base'//trim(adjustl(basepath))//'/'// &
		trim(adjustl(egpath))//'mulbase'//trim(adjustl(npath))//'/'


! 		if ( bf_index == 1 ) then
! 			datapath = '../basisfunction/base'//trim(adjustl(basepath))//'/linearbase'//trim(adjustl(npath))//'/'
! 		elseif ( bf_index == 2) then
! 			datapath = '../basisfunction/base'//trim(adjustl(basepath))//'/linearbase'//trim(adjustl(npath))//'/'
! 		elseif ( bf_index == 3) then
! 			datapath = '../basisfunction/base'//trim(adjustl(basepath))//'/1mulbase'//trim(adjustl(npath))//'/'
! 		elseif ( bf_index == 4) then
! 			datapath = '../basisfunction/base'//trim(adjustl(basepath))//'/2mulbase'//trim(adjustl(npath))//'/'
! 		else
! 			print*, 'error: read_basis, bf_index'
! 		end if

		print*, 'basis function files are from: '
		print*, trim(datapath)//'*.dat'
		print*, 'read multi-basis function, please wait ... ... '
		print*

		n = finemesh_index
		nx = n
		ny = n
		hhx = hx/nx
		hhy = hy/ny
		nz = (n + 2)*(n + 1)/2
		allocate(p(2, 3, mx, my, nz))
		p = 0._dp

		open(11, file=trim(datapath)//'nbase11.dat', status='old')
		open(12, file=trim(datapath)//'nbase12.dat', status='old')
		open(13, file=trim(datapath)//'nbase13.dat', status='old')
		open(14, file=trim(datapath)//'nbase21.dat', status='old')
		open(15, file=trim(datapath)//'nbase22.dat', status='old')
		open(16, file=trim(datapath)//'nbase23.dat', status='old')

		call cpu_time(start)

		do mj=1,my
		do mi=1,mx
			read(11,*) nxyf, (p(1, 1, mi, mj, k), k = 1, nxyf)
			read(12,*) nxyf, (p(1, 2, mi, mj, k), k = 1, nxyf)
			read(13,*) nxyf, (p(1, 3, mi, mj, k), k = 1, nxyf) 
			read(14,*) nxyf, (p(2, 1, mi, mj, k), k = 1, nxyf)
			read(15,*) nxyf, (p(2, 2, mi, mj, k), k = 1, nxyf)
			read(16,*) nxyf, (p(2, 3, mi, mj, k), k = 1, nxyf) 
		enddo
		enddo

		call cpu_time(finish)
		print'(" time of loading basis files: ")'
		print'(f6.3," seconds")', finish - start
		print*

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

		return
	end subroutine read_basis

	subroutine check_info()
		implicit none
		
		print*, 'basis function infomation: '
		print*, 'nx/ny/n: ', n
		print*, 'nz: ', nz
		print*
		
		return
	end subroutine check_info

	subroutine test()
		implicit none
		integer(kind=ki) :: nid
! 		print*, p(1, 2, 4, 4, nz), p(1, 2, 4, 4, nz-1)

		return
	end subroutine test

end module openbasis
