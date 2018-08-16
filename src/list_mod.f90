module list_mod

  use iso_fortran_env

  implicit none

  type :: item
    real(real64), dimension(3) :: r
    type(item), pointer :: next
  end type

  type :: list
    type(item), pointer :: first, last, current
    integer(int32) :: size, index
  contains
    procedure :: create
    procedure :: append
    procedure :: prepend
    procedure :: delete_index
    procedure :: delete_indices
    procedure :: destroy
    procedure :: to_array
    procedure :: print
    procedure :: reverse
  end type

  contains

    subroutine create(self)

      class(list) :: self
      
      self%size = 0

    end subroutine

    !********************************************************************************

    subroutine create_from_array(self, veclist)

      class(list) :: self
      real(real64), dimension(:, :) :: veclist
      integer :: ivec

      call self%create()

      do ivec = 1, size(veclist, 2)
        call self%append(veclist(:, ivec))
      enddo

    end subroutine

    !********************************************************************************

    subroutine append(self, vec)

      class(list) :: self
      real(real64), dimension(3) :: vec
      type(item), pointer :: newitem

      allocate(newitem)
      newitem%r = vec
      if (self%size == 0) then
        self%index = 1
        self%first => newitem
        self%current => newitem
      else
        self%last%next => newitem
      endif
      self%last => newitem
      self%size = self%size + 1

    end subroutine

    !********************************************************************************

    subroutine prepend(self, vec)

      class(list) :: self
      real(real64), dimension(3) :: vec
      type(item), pointer :: newitem

      allocate(newitem)
      newitem%r = vec
      if (self%size == 0) then
        self%index = 1
        self%last => newitem
        self%current => newitem
      else
        newitem%next => self%first
      endif
      self%first => newitem
      self%size = self%size + 1

    end subroutine

    !********************************************************************************

    subroutine delete_index(self, index)

      class(list) :: self
      integer :: index
      integer :: ilist
      type(item), pointer :: delete

      if (index == 1) then
        delete => self%first
        self%first => self%first%next
      else
        if (index >= self%index) then
          self%current => self%first
          self%index = 1
        endif
        do ilist = 1, index - self%index - 1
          self%current => self%current%next
          self%index = self%index + 1
        enddo
        delete => self%current%next
        if (index /= self%size) self%current%next => self%current%next%next
      endif
      deallocate(delete)
      nullify(delete)
      self%size = self%size - 1

    end subroutine

    !********************************************************************************

    subroutine delete_indices(self, indices)

      class(list) :: self
      integer, dimension(:) :: indices
      integer :: ilist, iind
      type(item), pointer :: delete

      print*, indices

      do iind = 1, size(indices)
        call self%delete_index(indices(iind))
        indices = indices - 1
      enddo

    end subroutine

    !********************************************************************************

    subroutine destroy(self)

      class(list) :: self

      do while (self%size > 0)
        self%current => self%first
        if (self%size > 1) self%first => self%first%next
        deallocate(self%current)
        self%size = self%size - 1
      enddo
      nullify(self%first)
      nullify(self%current)
      nullify(self%last)
      
    end subroutine

    !********************************************************************************

    function to_array(self)

      class(list) :: self
      real(real64), dimension(:, :), allocatable :: to_array
      integer :: i

      allocate(to_array(3, self%size))

      self%current => self%first
      do i = 1, self%size
        to_array(:, i) = self%current%r
        self%current => self%current%next
      enddo
      self%current => self%first
      self%index = 1

    end function

    !********************************************************************************

    subroutine print(self)

      class(list) :: self
      integer :: i

      self%current => self%first
      do i = 1, self%size
        print*, self%current%r
        self%current => self%current%next
      enddo

    end subroutine

    !********************************************************************************

    subroutine reverse(self)

      class(list) :: self
      integer :: i

      self%last => self%first
      do i = 3, self%size
        self%current => self%last%next
        self%last%next => self%current%next
        self%current%next => self%first
        self%first => self%current
      enddo
      self%current => self%last%next
      nullify(self%last%next)
      self%current%next => self%first
      self%first => self%current

    end subroutine

end module