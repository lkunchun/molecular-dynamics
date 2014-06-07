      subroutine dqint(q,boxx,boxy,nside,nmol,v,f)
c     procedure to initialize positions
c     in a 2-d hexagonal fashion

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension q(64,3),v(64,3),f(64,3)

c     initialize constants

      dqx=boxx/float(nside)
      dqy=boxy/float(nside)

c     main loop for position initialization

      if (nside.gt.0) then

      do j=1,nside
            if (mod(j,2).eq.0) then
                  dq2=dqx/2.0
            else
                  dq2=0.0
            endif
            do k=1,nside
                  indx=float(k+nside*(j-1))
                  q(indx,1)=dqx*float(k-1)+dq2
                  q(indx,2)=dqy*float(j-1)
            enddo
      enddo

      else
      open (1,file='md2d.cfg')
      read (1,*) (q(i,1),q(i,2),v(i,1),v(i,2),f(i,1),f(i,2),i=1,nmol)
      close(1)
      endif

      return
      end
