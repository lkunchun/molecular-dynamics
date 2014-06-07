      subroutine advance(q,v,oq,f,dt,flg,nmol,ndim)
c     procedure to advance by one time step

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension q(64,3),oq(64,3),v(64,3),f(64,3)

c     initialize constants

      dt2=dt*dt
      dt2i=1.0/(2.0*dt)
      dti=1.0/dt

c     check if starting from velocity and position

      if (nint(flg).eq.0) then

c           verlet with positions
c     write(*,*) 'verlet with positions'

            do j=1,ndim
                  do i=1,nmol
                        temp=2.0*q(i,j)-oq(i,j)+f(i,j)*dt2
                        v(i,j)=(temp-oq(i,j))*dt2i
                        oq(i,j)=q(i,j)
                        q(i,j)=temp
                  enddo
            enddo
      else

c           verlet with position and velocity
c     write(*,*) 'verlet with positions and veloticy'

            do j=1,ndim
                  do i=1,nmol
                        temp=q(i,j)+v(i,j)*dt+.50*f(i,j)*dt2
                        v(i,j)=(temp-q(i,j))*dti
                        oq(i,j)=q(i,j)
                        q(i,j)=temp
                  enddo
            enddo
            flg=0.0
      endif

      return
      end
