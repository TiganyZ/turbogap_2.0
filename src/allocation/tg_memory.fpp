#:def ranksuffix(RANK)
   $:'' if RANK == 0 else '(' + ':' + ',:' * (RANK - 1) + ')'
#:enddef ranksuffix

#:set PRECISIONS = ['int', 'dp', 'logical']
#:set PREC_NAME = ['integer', 'real(dp)', 'logical']

#:set RANKS = range(1, 6)

#:def dimensions(RANK)
   $: ''.join(  ['dim_' + str(r) + ', ' for r in range( 1, RANK+1 ) ] )
#:enddef dimensions

#:def dimensions_def(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['integer, intent(in) :: dim_' + str(r) for r in range( 1, RANK+1 ) ] )
#:enddef dimensions_def

#:def dimensions_def_de(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['integer :: dim_' + str(r) for r in range( 1, RANK+1 ) ] )
#:enddef dimensions_def_de

#:def dimensions_size(RANK)
   $:'dim_1' + ' '.join(  [' * dim_' + str(r)   for r in range( 2, RANK + 1 ) ] )
#:enddef dimensions_size

#:def dimensions_array_size(RANK)
   $:'1:dim_1' + ' '.join(  [', 1:dim_' + str(r)   for r in range( 2, RANK + 1 ) ] )
#:enddef dimensions_array_size

module tg_memory
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   interface tg_allocate
      #:for PREC in PRECISIONS
         #:for RANK in RANKS
            module procedure tg_allocate_${RANK}$_${PREC}$
         #:endfor
      #:endfor
   end interface tg_allocate

   interface tg_deallocate
      #:for PREC in PRECISIONS
         #:for RANK in RANKS
            module procedure tg_deallocate_${RANK}$_${PREC}$
         #:endfor
      #:endfor
   end interface tg_deallocate
contains

   #:for PREC in PRECISIONS
      #:for RANK in RANKS

         subroutine tg_allocate_${RANK}$_${PREC}$ (array, ${dimensions(RANK)}$memory_total, memory_max, rank, name)
            character(len=*), intent(in) :: name
            integer, intent(in) :: rank
            #:if PREC == 'int'
               integer, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            #:if PREC == 'dp'
               real(${PREC}$), intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               logical, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            integer, intent(inout) :: memory_total
            integer, intent(inout) :: memory_max
            ${dimensions_def(RANK) }$

            integer :: total_array_size

            total_array_size = ${dimensions_size(RANK) }$*size_type

#ifdef _CHECK_ALLOCATE

            write (*, *) 'TGAllOC ', name, ' on rank ', rank, ' size ', total_array_size, ' bytes, allocated? ', allocated(array), ' current memory ', memory_total / 1024**2, ' Mb'

            if (allocated(array)) then
               write (*, *) 'array ', name, ' has already been allocated'
            end if

#endif

            memory_total = memory_total + total_array_size
            memory_max = max(memory_total, memory_max)

            allocate (array(${dimensions_array_size(RANK)}$))

         end subroutine tg_allocate_${RANK}$_${PREC}$

      #:endfor
   #:endfor

   #:for PREC in PRECISIONS
      #:for RANK in RANKS

         subroutine tg_deallocate_${RANK}$_${PREC}$ (array, memory_total, memory_max, rank, name)
            character(len=*), intent(in) :: name
            integer, intent(in) :: rank
            #:if PREC == 'int'
               integer, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            #:if PREC == 'dp'
               real(${PREC}$), intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               logical, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            integer, intent(inout) :: memory_total
            integer, intent(inout) :: memory_max
            integer :: number_of_elements
            integer :: i
            integer :: total_array_size

            number_of_elements = 1
            do i = 1, ${RANK }$
               number_of_elements = number_of_elements*size(array, i)
            end do

            total_array_size = number_of_elements*size_type

#ifdef _CHECK_ALLOCATE

            write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, ' size ', total_array_size, ' bytes, allocated? ', allocated(array), ' current memory ', memory_total / 1024**2, ' Mb'

            if (.nor.allocated(array)) then
               write (*, *) 'array ', name, ' has already been deallocated'
            end if

#endif

            memory_total = memory_total - total_array_size
            memory_max = max(memory_total, memory_max)

            deallocate (array)

         end subroutine tg_deallocate_${RANK}$_${PREC}$

      #:endfor
   #:endfor

end module tg_memory
