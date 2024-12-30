
module class_linkedlist
  implicit none

  private
  !!!!!!!!!!!!!!!!!!!
  type :: Link
     integer ::i
     complex(8) :: a=(0.0d0,0.0d0),b=(0.0d0,0.0d0)
   !   double precision :: a=0.0d0,b=0.0d0
     type(Link),pointer :: next => null() 
  end type Link
  !!!!!!!!!!!!!!!!!!!!

  type,public :: LinkedList ! double ended Link-list for integer
     type(Link), pointer :: first => null()
     type(Link), pointer :: last => null()
     integer :: size=0
     
   contains
     procedure::insertAB,multiplyAB,isEmpty,find,insert,insertFirst,insertLast,deleteFirst,display,get,getAB,findindex,overwriteAB
  end type LinkedList

contains

subroutine overwriteAB(this,i,x,y)
   class(LinkedList),intent(inout) :: this
   complex(8) :: x,y
   ! double precision :: x,y
   integer ,intent(in):: i
   
   type(Link),pointer :: current

   

   if (this%isEmpty()) return
   
   current=>this%first
   do while (associated(current))
      if (current%i==i) then
         current%a=x
         current%b=y
         return
      end if
   current => current%next
   end do

end subroutine overwriteAB




   subroutine insertAB(this,i,x,y)
      class(LinkedList),intent(inout) :: this
      complex(8) :: x,y
      ! double precision :: x,y
      integer ,intent(in):: i
      
      type(Link),pointer :: current

      

      if (this%isEmpty()) return
      
      current=>this%first
      do while (associated(current))
         if (current%i==i) then
            current%a=current%a+x
            current%b=current%b+y
            return
         end if
      current => current%next
      end do

   end subroutine insertAB


   subroutine multiplyAB(this,i,x,y)
      class(LinkedList),intent(inout) :: this
      complex(8) :: x,y
      integer ,intent(in):: i
      
      type(Link),pointer :: current

      

      if (this%isEmpty()) return
      
      current=>this%first
      do while (associated(current))
         if (current%i==i) then
            current%a=current%a*x
            current%b=current%b*y
            return
         end if
      current => current%next
      end do

   end subroutine multiplyAB



  function isEmpty(this)
    class(LinkedList),intent(inout) :: this
    logical :: isEmpty
    isEmpty=.not.associated(this%first)
  end function isEmpty

      subroutine insertFirst(this,i)
        class(LinkedList),intent(inout) :: this
        integer,intent(in) :: i
        type(Link),pointer :: newlink
        
        ! create the new link
        allocate(newlink)
        newlink%i=i
        
        ! special case
        if (this%isEmpty()) this%last=> newLink
        !linked it
        newlink%next=>this%first
        this%first=>newlink
        this%size=this%size+1
      end subroutine insertFirst


      subroutine insertLast(this,i)
        class(LinkedList),intent(inout) :: this
        integer,intent(in) :: i
        type(Link),pointer :: newlink
        
        ! create the new link
        allocate(newlink)
        newlink%i=i
       
        if (this%isEmpty()) then
           this%first=> newLink  ! special case
        else
           this%last%next=>newLink
        end if
        !linked it
        this%last=>newlink
            this%size=this%size+1
      end subroutine insertLast

      
      
   !   function deleteFirst(this)
      subroutine deleteFirst(this)
        class(LinkedList),intent(inout) :: this
        integer :: deleteFirst0
        deleteFirst0=0
        if (.not.(this%isEmpty())) then ! list not empty
           deleteFirst0=this%first%i
           this%first=>this%first%next
           this%size=this%size-1
        end if
      ! end function deleteFirst
      end subroutine deleteFirst


 subroutine insert(this,i) ! insert in ordered list
        class(LinkedList),intent(inout) :: this
        integer,intent(in) :: i
        type(Link),pointer :: newlink,current,previous
        logical :: test

        this%size=this%size+1
        
 ! create the new link
        allocate(newlink)
        newlink%i=i
           
     
if (this%isEmpty()) then ! list is empty
   this%first => newLink
return
end if

if (this%first%i >= i) then ! item needs to go at first spot of list
   newlink%next => this%first;
   this%first => newlink;
return
end if

previous => this%first
!do while (associated(previous%next).and.(previous%next%i<i))
!      ! traverse until end of list or spot is found
!  previous => previous%next
!end do

test=.true.
do while (test)
   if (.not.(associated(previous%next))) exit
   if (.not.(previous%next%i<i)) exit
      ! traverse until end of list or spot is found
  previous => previous%next
   end do
 
      if (.not.(associated(previous%next))) then ! check if at end of list
      previous%next => newLink
else !// implement insertion
newLink%next => previous%next
previous%next => newLink
end if

this%last=>newLink

       end subroutine insert


      
 subroutine display(this)
   class(LinkedList),intent(inout) :: this
   type(Link),pointer :: current
   current=>this%first
   do while (associated(current))
      print *,current%i
      current => current%next
   end do
 end subroutine display

 function get(this,m) result(j)
   class(LinkedList),intent(inout) :: this
   integer ,intent(in):: m
   integer :: e,j
   
   type(Link),pointer :: current
 
   

   if (this%isEmpty()) return
   

   current=>this%first
   if (m==1) then
      j=current%i
      return
   end if
   do e=1,m-1
   current => current%next
   end do
   j=current%i
   
   end function get

   function getAB(this,m) result(xy)
      class(LinkedList),intent(inout) :: this
      integer ,intent(in):: m
      integer :: e
      complex(8),dimension(1:2) :: xy
      ! double precision,dimension(1:2) :: xy
      logical:: find
      
      type(Link),pointer :: current
    
      
      find=.false. ! default
      if (this%isEmpty()) return
      
   
      current=>this%first
      if (m==1) then
         xy=(/current%a,current%b/)
         return
      end if
      do e=1,m-1
      current => current%next
      end do
      xy=(/current%a,current%b/)
      
      end function getAB
 

function find(this,i) 
  class(LinkedList),intent(inout) :: this
  integer ,intent(in):: i
  logical:: find
  
  type(Link),pointer :: current

  
  find=.false. ! default
  if (this%isEmpty()) return
  
  current=>this%first
  do while (associated(current))
     if (current%i==i) then
        find=.true.
        return
     end if
current => current%next
end do
end function find


function findindex(this,i) result(index)
   class(LinkedList),intent(inout) :: this
   integer ,intent(in):: i
   logical:: find
   integer :: index,j
   
   type(Link),pointer :: current
 
   
   if (this%isEmpty()) return
   

   j=1
   current=>this%first
   do while (associated(current))
      if (current%i==i) then
         index=j
         return
      end if
      j=j+1
 current => current%next
 end do
 end function findindex


end module class_linkedlist

