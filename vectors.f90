module vectors
	implicit none

	private ki, dp

 	integer, parameter :: ki = selected_int_kind(9)
	integer(kind=ki), parameter :: dp = selected_real_kind(r=50, p=14)

	type vector 
		real(kind=dp) :: x, y
	end type vector

	interface assignment ( = )
		module procedure :: realarray_to_vector
		module procedure :: intarray_to_vector
	end interface

	interface operator( * )
		module procedure :: dotproduct
		module procedure :: vecvalmulti
		module procedure :: valvecmulti
	end interface

	interface operator( + )
		module procedure :: vecadd
		module procedure :: vec_add_val
		module procedure :: val_add_vec
	end interface

	interface operator( - )
		module procedure :: vecsub
		module procedure :: vec_sub_val
		module procedure :: val_sub_vec
		module procedure :: minus
	end interface

! 	interface assignment ( - )
! 	end interface

contains

! 	=

	subroutine realarray_to_vector(v1, array)
		implicit none
		type(vector), intent(out) :: v1
		real(kind=dp), intent(in) :: array(2)

		v1%x = array(1)
		v1%y = array(2)

		return
	end subroutine realarray_to_vector

	subroutine intarray_to_vector(v1, array)
		implicit none
		type(vector), intent(out) :: v1
		integer(kind=ki), intent(in) :: array(2)

		v1%x = array(1)
		v1%y = array(2)

		return
	end subroutine intarray_to_vector

! 	*

	real(kind=dp) &
	function dotproduct(v1, v2)
		implicit none
		type(vector), intent(in) :: v1, v2
		
		dotproduct = v1%x*v2%x + v1%y*v2%y

		return		
	end function dotproduct

	type(vector) &
	function vecvalmulti(vec, value)
		implicit none
		type(vector), intent(in) :: vec
		real(kind=dp), intent(in) :: value
		
		vecvalmulti%x = vec%x*value
		vecvalmulti%y = vec%y*value

		return		
	end function vecvalmulti

	type(vector) &
	function valvecmulti(value, vec)
		implicit none
		real(kind=dp), intent(in) :: value
		type(vector), intent(in) :: vec
		
		valvecmulti%x = value*vec%x
		valvecmulti%y = value*vec%y

		return		
	end function valvecmulti

! 	+

	type(vector) &
	function vecadd(v1, v2)
		implicit none
		type(vector), intent(in) :: v1, v2
		
		vecadd%x = v1%x + v2%x
		vecadd%y = v1%y + v2%y

		return		
	end function vecadd

	type(vector) &
	function vec_add_val(v1, value)
		implicit none
		type(vector), intent(in) :: v1
		real(kind=dp), intent(in) :: value
		
		vec_add_val%x = v1%x + value
		vec_add_val%y = v1%y + value

		return		
	end function vec_add_val

	type(vector) &
	function val_add_vec(value, v1)
		implicit none
		real(kind=dp), intent(in) :: value
		type(vector), intent(in) :: v1
		
		val_add_vec%x = value + v1%x
		val_add_vec%y = value + v1%y

		return		
	end function val_add_vec

! 	-

	type(vector) &
	function vecsub(v1, v2)
		implicit none
		type(vector), intent(in) :: v1, v2
		
		vecsub%x = v1%x - v2%x
		vecsub%y = v1%y - v2%y

		return		
	end function vecsub

	type(vector) &
	function vec_sub_val(v1, value)
		implicit none
		type(vector), intent(in) :: v1
		real(kind=dp), intent(in) :: value
		
		vec_sub_val%x = v1%x - value
		vec_sub_val%y = v1%y - value

		return		
	end function vec_sub_val

	type(vector) &
	function val_sub_vec(value, v1)
		implicit none
		real(kind=dp), intent(in) :: value
		type(vector), intent(in) :: v1
		
		val_sub_vec%x = value - v1%x
		val_sub_vec%y = value - v1%y

		return		
	end function val_sub_vec

	type(vector) &
	function minus(v2)
		implicit none
		type(vector), intent(in) :: v2

		minus%x = -v2%x
		minus%y = -v2%y
		
		return
	end function minus

	subroutine vec_test()
		implicit none
		type(vector) :: v1, v2
		real(kind=dp) :: value

		v1 = [1._dp, 2._dp]
		v2 = [3, 4]
		value = 10._dp
		print*, v1, v2 
		print*, v1 + value, value + v2
		print*, v1 - value, value - v2
		print*, -v1, -v2
				
		return		
	end subroutine vec_test

end module vectors
