program main
	use vectors
	use constants, only: ki, dp, pi
 	use domainmesh, only: get_domainmesh
 	use getlgmap, only: get_lgmap
 	use toggles
 	use prob, only: get_prob
 	use openbasis, only: open_basis
 	use localid, only: local_eleid
 	use assembstiff, only: assemb_stiff
 	use solveequation, only: solve_equation
 	use geterror, only: get_error
	implicit none

	print*
	print*, "----- ----- THIS IS MULTI-DFVM METHOD ----- -----"
	print*

! 	call vec_test()

! 	step0: set toggles
	call set_toggles()
! 	step1: get parameters
	call get_domainmesh()
! 	step2: get lgmap
	call get_lgmap()
! 	step3: get problem
	call get_prob()
! 	step4: open multi-basis function
	call open_basis()
! 	step5: get edge element local id
	call local_eleid()
! 	step6: assemble stiff matrix
	call assemb_stiff()
! 	step7: solve the equation
	call solve_equation()
! 	step8: get error
	call get_error()

end program main
