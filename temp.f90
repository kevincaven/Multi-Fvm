program main

	implicit none
	
	integer :: a(3, 1), b(1, 3)
	a(1, 1) = 1; a(2, 1) = 2; a(3, 1) = 3
	b(1, 1) = 1; b(1, 2) = 2; b(1, 3) = 3

	print*, matmul(a, b)


end program main
