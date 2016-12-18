module constants

	implicit none

  	private
  	public pi, dp, ki

 	integer, parameter :: ki = selected_int_kind(9)
	integer(kind=ki), parameter :: dp = selected_real_kind(r=50, p=14)

	real(kind=dp), parameter :: pi = acos(-1.d0)

end module constants