program main

    !! The "main" routine of the program  
    !!   
    !! This is where execution starts when running the program
    
    use processor
    use fea

    implicit none

    ! Read model data
    call input

    ! Initialize problem
    call initial
	
	!Choose what kind of FEA
	select case (antype)
        case ('static')
            ! Calculate displacements
    		call displ
        case ('static_nl')
            write (*, *) 'static_nl analysis is not implemented'
            stop
        case ('modal')
            ! Calculate eigen values
    		call eigen
        case ('angle')
            write (*, *) 'angle analysis is not implemented'
            stop
        case ('trans')
            ! Calculate transient response history
            call stopwatch('star')
            call trans
            call stopwatch('stop')
        case default
            write (*, *) 'No type was given for analysis, Using default static analysis'
            call displ
    end select
    
    ! Close plot window(s)
    call plot( done )
   
end program main
