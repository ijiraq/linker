module utils
  implicit none

  public
  
contains
  function replace(string, c1, c2) result(modifiedString)
    implicit none
    character(len=*), intent(in) :: string, c1, c2
    character(len=:), allocatable :: modifiedString
    integer :: k
    
    modifiedString=string
    do while(INDEX(modifiedString,c1).ne.0)
       k = INDEX(modifiedString,c1)
       modifiedString = modifiedString(:k-1)//c2//modifiedString(k+1:)
    enddo
  end function replace

  function argsort(r) result(args)
    real(kind=8), intent(in), dimension(:) :: r
    integer, dimension(size(r)) :: d
    integer, allocatable :: args(:)
    
    integer, dimension(size(r)) :: il
    
    integer :: stepsize
    integer :: i,j,n,left,k,ksize
    
    n = size(r)
    
    do i=1,n
       d(i)=i
    end do
    
    if ( n==1 ) return
    
    stepsize = 1
    do while (stepsize<n)
       do left=1,n-stepsize,stepsize*2
          i = left
          j = left+stepsize
          ksize = min(stepsize*2,n-left+1)
          k=1
          
          do while ( i<left+stepsize .and. j<left+ksize )
             if ( r(d(i))>r(d(j)) ) then
                il(k)=d(i)
                i=i+1
                k=k+1
             else
                il(k)=d(j)
                j=j+1
                k=k+1
             endif
          enddo
          
          if ( i<left+stepsize ) then
             ! fill up remaining from left
             il(k:ksize) = d(i:left+stepsize-1)
          else
             ! fill up remaining from right
             il(k:ksize) = d(j:left+ksize-1)
          endif
          d(left:left+ksize-1) = il(1:ksize)
       end do
       stepsize=stepsize*2
    end do
    args = d 
    return
  end function argsort


end module utils
