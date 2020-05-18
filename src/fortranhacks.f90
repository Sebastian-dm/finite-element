module fortranhacks

use types
implicit none
save
    
private
public :: identity, increment_filename
    
contains
	subroutine identity(mat)
    	! This subroutine sets the matrix m to the identity matrix

        integer :: n, s(2)
        real(wp), intent(out) :: mat(:,:)

		mat = 0.0_wp
        s = shape(mat)
        do n = 1, s(1)
          	mat(n,n) = 1.0_wp
       	end do
    end subroutine identity

    subroutine increment_filename(outfile_path)
    	!This subroutine returns the number as
        !the lowest unused number for the given file name base
        ! eg. filebasename = base and base_1 and base_2 exists
        ! call increment_filename will then return base_3

        character(len=256), intent(out) :: outfile_path
        character(len=3) :: number
        integer :: n
        logical :: file_exists

        file_exists = .true.
        number = " 0"
        n = 0
        
        do while(file_exists)
          	inquire(file=trim(outfile_path)//"_"//trim(adjustl(number))//".m", exist=file_exists)
            if (.not. file_exists) then
              	outfile_path = trim(outfile_path)//"_"//trim(adjustl(number))//".m"
            else
              	!print *, "File number",n,"existed"
            end if
            n = n+1
            write(number,'(I2)') n
        end do
        
    end subroutine increment_filename

        

end module fortranhacks

! PRINTS
!$$$$$$ print *,"me"                            ! ME PRINT BLOCK
!$$$$$$ do n = 1,8                              ! ME PRINT BLOCK
!$$$$$$     print "(24(f4.2,tr1))",me(n,1:8)    ! ME PRINT BLOCK
!$$$$$$ end do                                  ! ME PRINT BLOCK