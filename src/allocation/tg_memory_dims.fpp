#:def ranksuffix(RANK)
   $:'' if RANK == 0 else '(' + ':' + ',:' * (RANK - 1) + ')'
#:enddef ranksuffix

#:set PRECISIONS = ['int', 'dp', 'logical']
#:set PREC_NAME = ['integer', 'real(dp)', 'logical']

#:set RANKS = range(1, 6)

#:def dimensions(RANK)
   $: ' dims, '
#:enddef dimensions

#:def dimensions_def(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['integer, intent(in) :: dims(' + str(RANK) + ')'  ] )
#:enddef dimensions_def

#:def size_def(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['integer :: sizes(' + str(RANK) + ')'  ] )
#:enddef size_def

#:def size_assign(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['sizes(' + str(RANKI) + ') = size( array, ' + str( RANKI ) + " )" for RANKI in range( 1, RANK+1 ) ] )
#:enddef size_assign

#:def size_assigntg(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['sizes(' + str(RANKI) + ') = size( tg_array%array, ' + str( RANKI ) + " )" for RANKI in range( 1, RANK+1 ) ] )
#:enddef size_assigntg

#:def dimensions_def_de(RANK)
   $:'' if RANK == 0 else '\n'.join(  ['integer :: dim_' + str(r) for r in range( 1, RANK+1 ) ] )
#:enddef dimensions_def_de

#:def dimensions_size(RANK)
   $:'dims(1)' + ' '.join(  [' * dims(' + str(r) + ')'   for r in range( 2, RANK + 1 ) ] )
#:enddef dimensions_size

#:def sizes_size(RANK)
   $:'sizes(1)' + ' '.join(  [' * sizes(' + str(r) + ')'   for r in range( 2, RANK + 1 ) ] )
#:enddef sizes_size

#:def dimensions_array_size(RANK)
   $:'1:dims(1)' + ' '.join(  [', 1:dims(' + str(r) + ')'   for r in range( 2, RANK + 1 ) ] )
#:enddef dimensions_array_size

module tg_memory
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   interface tg_alloc
      #:for PREC in PRECISIONS
         #:for RANK in RANKS
            module procedure tg_alloc_${RANK}$_${PREC}$
         #:endfor
      #:endfor
   end interface tg_alloc

   interface tg_dealloc
      #:for PREC in PRECISIONS
         #:for RANK in RANKS
            module procedure tg_dealloc_${RANK}$_${PREC}$
         #:endfor
      #:endfor
   end interface tg_dealloc

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

   #:for PREC in PRECISIONS
      #:for RANK in RANKS
         type tg_array_${RANK}$_${PREC}$
            #:if PREC == 'int'
               integer, allocatable :: array${ranksuffix(RANK) }$
               integer :: size_type = 4
            #:endif
            #:if PREC == 'dp'
               real(${PREC}$), allocatable :: array${ranksuffix(RANK) }$
               integer :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               logical, allocatable :: array${ranksuffix(RANK) }$
               integer :: size_type = 4
            #:endif
            integer :: dims(${RANK}$) = 0
            integer :: used_dims(${RANK}$) = 0
            logical :: allocated = .false.
            character(64) :: name = "none"
         end type tg_array_${RANK}$_${PREC}$

      #:endfor
   #:endfor
contains

   #:for PREC in PRECISIONS
      #:for RANK in RANKS

         subroutine tg_alloc_${RANK}$_${PREC}$ (tg_array, ${dimensions(RANK)}$memory_total, memory_max, rank, name)
            character(len=*), intent(in) :: name
            integer, intent(in) :: rank
            type(tg_array_${RANK}$_${PREC}$), intent(inout) :: tg_array
            integer, intent(inout) :: memory_total
            integer, intent(inout) :: memory_max
            ${dimensions_def(RANK) }$
            ${size_def(RANK) }$
            #:if PREC == 'int'
               integer, parameter :: size_type = 4
            #:endif
            #:if PREC == 'dp'
               integer, parameter :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               integer, parameter :: size_type = 4
            #:endif

            integer :: total_array_size
            integer :: prev_array_size
            integer :: i

            total_array_size = ${dimensions_size(RANK) }$*size_type

#ifdef _CHECK_ALLOCATE

            write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
               ' size ', total_array_size, ' bytes, allocated? ', &
               allocated(tg_array%array), ' current memory ', &
               dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
               ' memory max ', memory_max, " bytes"

            if (allocated(tg_array%array)) then
               write (*, *) 'array ', name, ' has already been allocated'
               write (*, *) 'input dims   ', dims
               write (*, *) 'dims         ', tg_array%dims
               write (*, *) 'used_dims    ', tg_array%used_dims
            end if

#endif

            if (allocated(tg_array%array)) then
               ! Check the sizes currently to modify the memory
               ${size_assigntg(RANK)}$

#ifdef _CHECK_ALLOCATE
               write (*, *) 'found dims   ', sizes
               write (*, *) 'already dims ', tg_array%dims
#endif
               prev_array_size = ${sizes_size(RANK) }$*size_type
               if (any((dims - tg_array%dims) > 0)) then
                  ! We must reallocate
                  ! If this condition is not true then the array is of sufficient size
                  deallocate (tg_array%array)
                  memory_total = memory_total - prev_array_size

#ifdef _CHECK_ALLOCATE
                  write (*, *) 'reallocating due to dim mismatch from request'
#endif
               else
                  tg_array%used_dims = dims
#ifdef _CHECK_ALLOCATE
                  write (*, *) 'not reallocating as size is larger than request'
#endif
                  return
               end if
            end if

            memory_total = memory_total + total_array_size
            memory_max = max(memory_total, memory_max)

            allocate (tg_array%array(${dimensions_array_size(RANK)}$))
            tg_array%used_dims = dims
            tg_array%dims = dims

         end subroutine tg_alloc_${RANK}$_${PREC}$

      #:endfor
   #:endfor

   #:for PREC in PRECISIONS
      #:for RANK in RANKS

         subroutine tg_dealloc_${RANK}$_${PREC}$ (tg_array, memory_total, rank, name)
            character(len=*), intent(in) :: name
            integer, intent(in) :: rank
            type(tg_array_${RANK}$_${PREC}$), intent(inout) :: tg_array
            #:if PREC == 'int'
               integer, parameter :: size_type = 4
            #:endif
            #:if PREC == 'dp'
               integer, parameter :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               integer, parameter :: size_type = 4
            #:endif
            integer, intent(inout) :: memory_total
            integer :: number_of_elements
            integer :: i
            integer :: total_array_size

            number_of_elements = 1
            do i = 1, ${RANK }$
               number_of_elements = number_of_elements*size(tg_array%array, i)
            end do

            total_array_size = number_of_elements*size_type

#ifdef _CHECK_ALLOCATE

            write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
               ' size ', total_array_size, ' bytes, allocated? ', &
               allocated(tg_array%array), ' current memory ', &
               dfloat(memory_total)/dfloat(1024)**2, ' Mb'

            if (.not. allocated(tg_array%array)) then
               write (*, *) 'array ', name, ' has already been deallocated'
            end if

#endif
            if (.not. allocated(tg_array%array)) then
               return
            end if

            memory_total = memory_total - total_array_size

            deallocate (tg_array%array)
            tg_array%used_dims = 0
            tg_array%dims = 0

         end subroutine tg_dealloc_${RANK}$_${PREC}$

      #:endfor
   #:endfor

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
               real(${PREC}$), allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               logical, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            integer, intent(inout) :: memory_total
            integer, intent(inout) :: memory_max
            ${dimensions_def(RANK) }$
            ${size_def(RANK) }$

            integer :: total_array_size
            integer :: prev_array_size
            integer :: i

            total_array_size = ${dimensions_size(RANK) }$*size_type

#ifdef _CHECK_ALLOCATE

            write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
               ' size ', total_array_size, ' bytes, allocated? ', &
               allocated(array), ' current memory ', &
               dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
               ' memory max ', memory_max, " bytes"

            if (allocated(array)) then
               write (*, *) 'array ', name, ' has already been allocated'
            end if

#endif

            if (allocated(array)) then
               ! Check the sizes currently to modify the memory
               ${size_assign(RANK)}$
               prev_array_size = ${sizes_size(RANK) }$*size_type
               if (any((dims - sizes) > 0)) then
                  ! We must reallocate
                  ! If this condition is not true then the array is of sufficient size
                  deallocate (array)
                  memory_total = memory_total - prev_array_size
               else
                  return
               end if
            end if

            memory_total = memory_total + total_array_size
            memory_max = max(memory_total, memory_max)

            allocate (array(${dimensions_array_size(RANK)}$))

         end subroutine tg_allocate_${RANK}$_${PREC}$

      #:endfor
   #:endfor

   #:for PREC in PRECISIONS
      #:for RANK in RANKS

         subroutine tg_deallocate_${RANK}$_${PREC}$ (array, memory_total, rank, name)
            character(len=*), intent(in) :: name
            integer, intent(in) :: rank
            #:if PREC == 'int'
               integer, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            #:if PREC == 'dp'
               real(${PREC}$), allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 8
            #:endif
            #:if PREC == 'logical'
               logical, allocatable, intent(inout) :: array${ranksuffix(RANK) }$
               integer, parameter :: size_type = 4
            #:endif
            integer, intent(inout) :: memory_total
            integer :: number_of_elements
            integer :: i
            integer :: total_array_size

            number_of_elements = 1
            do i = 1, ${RANK }$
               number_of_elements = number_of_elements*size(array, i)
            end do

            total_array_size = number_of_elements*size_type

#ifdef _CHECK_ALLOCATE

            write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
               ' size ', total_array_size, ' bytes, allocated? ', &
               allocated(array), ' current memory ', &
               dfloat(memory_total)/dfloat(1024)**2, ' Mb'

            if (.not. allocated(array)) then
               write (*, *) 'array ', name, ' has already been deallocated'
            end if

#endif

            if (.not. allocated(array)) then
               return
            end if

            memory_total = memory_total - total_array_size

            deallocate (array)

         end subroutine tg_deallocate_${RANK}$_${PREC}$

      #:endfor
   #:endfor

end module tg_memory
