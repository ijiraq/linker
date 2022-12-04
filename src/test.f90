program test
  real a(4,100000), b(4), c(100000), d(100000,4)
  integer i, j
  real start, finish

  CALL RANDOM_NUMBER(a)
  CALL RANDOM_NUMBER(b)
  CALL RANDOM_NUMBER(d)
  
  call cpu_time(start)
  do i=1,4
     do j=1,100000
        c(j) = a(i,j)*b(i)
     enddo
  enddo
  call cpu_time(finish)
  write(0,*) "(4x1000000) matrix multiple in loop took: ",finish-start," seconds"
  call cpu_time(start)
  do j=1,100000
     do i=1,4
        c(j) = d(j,i)*b(i)
     enddo
  enddo
  call cpu_time(finish)
  write(0,*) "(100000x4) matrix multiple in loop took: ",finish-start," seconds"
  call cpu_time(start)
  c = MATMUL(d,b)
  call cpu_time(finish)
  write(0,*) "(100000x4) matmul took: ",finish-start," seconds"
end program test
  
