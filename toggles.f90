module toggles
	use constants, only: ki, dp, pi
	implicit none
	
	integer(kind=ki) :: eg_index, ineps
	integer(kind=ki) :: bf_index
	integer(kind=ki) :: gamma_index
	integer(kind=ki) :: penalty_index, pcoef
	integer(kind=ki) :: edgeint_index
	integer(kind=ki) :: finestmesh_index, coarsemesh_index, finemesh_index

contains

	subroutine set_toggles()
		implicit none
		
		ineps = 100

! 		print*, 'please input the problem index! ( 1 or 2 or 3)'
! 		read(*, *) eg_index
		eg_index = 3

! 		bf_index = 3
		bf_index = eg_index

		pcoef = 20

		penalty_index = 2

		gamma_index = 2

		edgeint_index = 2

! 		mesh

		finestmesh_index = 1024

		coarsemesh_index = 32

		finemesh_index = finestmesh_index/coarsemesh_index


		return
	end subroutine set_toggles

end module toggles