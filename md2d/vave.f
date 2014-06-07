      subroutine vave(v,va,nmol,ndim)
c     precedure to compute CM drift velocity

c     declaraction of variables
      implicit real*4(a-h,o-z)
      dimension v(64,3),va(3)

      do j=1,ndim
            va(j)=0.0
            do i=1,nmol
                  va(j)=va(j)+v(i,j)
            enddo
            va(j)=va(j)/float(nmol)
      enddo

      return
      end
