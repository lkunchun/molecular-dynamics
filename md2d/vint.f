      subroutine vint(v,vbar,twopi,i1,i2,nmol,ndim,flg,nside)
c     precedure to initialize velocities

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension v(64,3),va(3)

c     set flag to initial velocities

      flg=1.0

c     main loop for randomizing velocities

      if (nside.gt.0) then
      do j=1,ndim
            do i=1,nmol

c                 pick random velocity from gaussian distribution

                  x=rand(0)
                  y=rand(0)
                  vrand=vbar*sqrt(-2.0*log(x))*cos(twopi*y)
                  v(i,j)=vrand
            enddo
      enddo
      endif

c     determine center of mass drift velocities

      call vave(v,va,nmol,ndim)

c     loop to set CM velocity of system to zero

      do j=1,ndim
            do i=1,nmol
                  v(i,j)=v(i,j)-va(j)
            enddo
      enddo

      return
      end
