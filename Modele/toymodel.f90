module param
real, parameter ::l=1,h=0.3 !longueur, hauteur
real ::dx, dy
integer, parameter :: nl=100, nh=300, nstep=1000 !nb pts 
end module param


program Stage
use param
implicit none
real, dimension(nh,nl) :: x,y
integer :: i,j,k
real :: x0,y0,x1,y1
  
!Discrétisation de l'espace

open(10, file='Espace')

dx=l/nl
dy=h/nh

do i=1,nh
   do j=1,nl
      x(i,j)=j*dx
      y(i,j)=i*dy
      write(10,*) x(i,j), y(i,j)
   end do
end do

close(10)

! Marche aléatoire et probabilités
call srand(3)
open(11, file='position')
write(*,*) 'Position initiale : x0, y0'
read(*,*) x0, y0
write(11,*) x0, y0

do k=1,nstep
   call marchealea(x0,y0,x1,y1)
   write(11,*) x1, y1
   x0=x1
   y0=y1
end do

close(11)

end program Stage

subroutine marchealea(x0,y0,x1,y1)
use param
implicit none
real :: a, b, x0, y0, x1, y1, dlim, dmin
real, dimension(nh,nl):: xa, ya,d
real, dimension(2):: mi

call attracteurs(xa,ya,x0,y0)

dlim=sqrt(25*(dx**2+dy**2))

if (dmin<dlim .and. dmin>0) then !si attracteur
   x1=x0+sign(dx,xa(mi(1),mi(2))-x0)
   y1=y0+sign(dy,ya(mi(1),mi(2))-y0)
else !sinon ecoulement
   a=0.5
   b=rand()
   x1=x0

   if (b>a) then
      x1=x0+dx
   else if (b<2*a/3) then
      x1=x0-dx
   end if

   y1=y0

   if (b>0.65) then
      y1=y0+dy
   else if (b<0.35) then
      y1=y0-dy
   end if
end if    

end subroutine marchealea

subroutine attracteurs(xa,ya,x0,y0)
use param
implicit none
integer :: i,j,m
real, dimension(nh,nl):: xa, ya,d
real, dimension(2):: mi
real :: dmin, x0, y0

xa=0
ya=0

open(12, file='attracteurs')
m=int(nh/2)

do i=m-10,m+10
   do j=1,nl
      xa(i,j)=j*dx
      ya(i,j)=i*dy
      d(i,j)=sqrt((xa(i,j)-x0)**2+(ya(i,j)-y0)**2) !distances
      write(12,*) xa(i,j), ya(i,j)
   end do
end do
close(12)

dmin= minval(d)
mi=minloc(d)

end subroutine attracteurs



  
  
  
