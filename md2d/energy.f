      subroutine energy(v,pe,e,nmol,ndim)
c     precedure to compute total energy

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension v(64,3)

c     initialize constants

      c=.50
      ekin=0.0

c     loop to compute kinetic energies

      do j=1,ndim
            do i=1,nmol
                  ekin=ekin+c*v(i,j)*v(i,j)
            enddo
      enddo

c     compute total energy

      e=pe+ekin

      return
      end
