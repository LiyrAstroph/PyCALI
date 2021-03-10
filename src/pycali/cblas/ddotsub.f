*     ddotsub.f
 
*     The program is a fortran wrapper for ddot.
*     Witten by Keita Teranishi.  2/11/1998
*
       subroutine ddotsub(n,x,incx,y,incy,dot)
*
       external ddot
       double precision ddot
       integer n,incx,incy
       double precision x(*),y(*),dot
*
       dot=ddot(n,x,incx,y,incy)
       return
       end
