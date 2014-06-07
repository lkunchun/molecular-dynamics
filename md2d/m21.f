c     program to simulation a system of 2-d Lennard-Jones spheres
c     Ref: David Chandler, Lecture Notes, Stat. Mech. of Liquids.

c     declare variables and types

      implicit real*4(a-h,o-z)
      integer wtype,errind,dcunit,lx,ly,rand
      dimension q(64,3),oq(64,3),v(64,3),f(64,3),va(3),xv(64),yv(64)
      dimension oxv(64),oyv(64)


c     input parameters

      write(*,*)'time step, temperature, density, atoms per side'
      write(*,*)'no. of atoms, no. of time steps, potential cutoff'
      write(*,*)'thermalization increment, thermalization maximum'
      write(*,*) 'display increment:'
      read(*,*)g,tmp,rowst,nside,nmol,max,rcut,nrand,maxran,idisp
      write(*,*)g,tmp,rowst,nside,nmol,max,rcut,nrand,maxran,idisp

c     open (21,file='md2d.ou1',status='new')
c     open (22,file='md2d.ou2',status='new')
c     open (23,file='md2d.ou3',status='new')
      open (21,file='md2d.ou1')
      open (22,file='md2d.ou2')
      open (23,file='md2d.ou3')

      open (2, file='md2d.ran')
      read (2,*) rseed
      close(2)
      call srand(rseed)


c     initialize constants
c     units are:
c     time : tau0
c     length : sigma
c     energy : epsilon

      ndim=2
      tau0=1.0/48.0**0.5
      sigma=1.0
      epsilon=4.0
      dt=tau0*g
      boxx=(float(nmol)/rowst)**(1.0/float(ndim))
c     boxy=boxx*sqrt(3.0/4.0)
      boxy=boxx
      write(*,*) 'boxx=',boxx,'   boxy=',boxy
      twopi=8.0*atan(1.0)
      i1=1234
      i2=9643
      vbar=sqrt(tmp)
      write(*,*)'vbar=',vbar
      rnorm=2.0/float(ndim*nmol)
      sigma2=sigma*sigma
      sigma6=sigma2*sigma2*sigma2
      rcut2=rcut*rcut
      rcut6=rcut2*rcut2*rcut2
      e1=(sigma6/rcut6)
      e2=e1*e1
      delta=-epsilon*(e2-e1)

c     initialize positions

      call dqint(q,boxx,boxy,nside,nmol,v,f)

c     initialize velocities

      call vint(v,vbar,twopi,i1,i2,nmol,ndim,flg,nside)

c     start main loop for integration

c     call energy(v,pe,e,nmol,ndim)
c     write(*,*)'e1=',e,'pe1=',pe

      do i=1,max

c           plot all particles

            if (mod(i,idisp).eq.0) then
c                 write(*,*)'i=',i,'tmp=',tmpnew
c                 do j=1,ndim
c                       write(*,*)'vave(',j,')=',va(j)
c                 enddo

                  write(21,1001)real(i),tmpnew
                  write(22,1002)real(i),(q(j,1),q(j,2),j=1,nmol)
                  write(23,1002)real(i),(v(j,1),v(j,2),j=1,nmol)
1001              format(5e15.6)
1002              format(132e15.6)
            endif

c           determine force on each atom

            call dforce(f,q,boxx,boxy,nmol,ndim,pe,delta,sigma6
     1      ,epsilon,rcut2)

c           take one time step

            call advance(q,v,oq,f,dt,flg,nmol,ndim)

c           check energy

            call energy(v,pe,e,nmol,ndim)

c           write(*,*)'e=',e,'pe=',pe

c           compute temperature

            tmpnew=(e-pe)*rnorm

c           check momentum

            call vave(v,va,nmol,ndim)

c           randomize velocities every nrand steps up to maxran
            if (mod(i,nrand).eq.0.and.i.le.maxran)
     1      call vint(v,vbar,twopi,i1,i2,nmol,ndim,flg,nside)

      enddo

      iseed = (float(rand())/32767.0)*(2.D0**31-1.D0)
c     open (2, file='md2d.ran',status='new')
      open (2, file='md2d.ran')
      write(2,*) iseed
      close(2)

c     open (1, file='md2d.cfg',status='new')
      open (1, file='md2d.cfg')
      write(1,1003) (q(i,1),q(i,2),v(i,1),v(i,2),f(i,1),f(i,2),i=1,nmol)
1003  format(6e20.8)
      close(1)
      end

      subroutine energy(v,pe,e,nmol,ndim)
c     precedure to compute total energy

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension v(64,3)
      integer rand
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

      subroutine vint(v,vbar,twopi,i1,i2,nmol,ndim,flg,nside)
c     precedure to initialize velocities

c     declaration of variables
      
      implicit real*4(a-h,o-z)
      dimension v(64,3),va(3)
      integer rand
      
c     set flag to initial velocities

      flg=1.0

c     main loop for randomizing velocities

      if (nside.gt.0) then
      do j=1,ndim
            do i=1,nmol

c                 pick random velocity from gaussian distribution

                  x=float(rand())/32767.0
                  y=float(rand())/32767.0
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

      subroutine advance(q,v,oq,f,dt,flg,nmol,ndim)
c     procedure to advance by one time step

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension q(64,3),oq(64,3),v(64,3),f(64,3)
      integer rand
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

      subroutine dforce(f,q,boxx,boxy,nmol,ndim,pe,delta,sigma6
     1,epsilon,rcut2)

c     declaraction of variables

      implicit real*4(a-h,o-z)
      dimension f(64,3),q(64,3),efd2(64),efr(64,3)
      integer rand
c     initialization of constants

      boxxi=1.0/boxx
      boxyi=1.0/boxy
      sigma12=sigma6*sigma6
      a=-sigma6*epsilon
      b=sigma12*epsilon
c     write(*,*)'force'
c     write(*,*)'a=',a,'b=',b
c     write(*,*)'force','ndim=',ndim,'nmol=',nmol,'sigma6=',sigma6
c     write(*,*)'epsilon=',epsilon

c     set initial force & pe to zero

      pe=0.0
      do j=1,ndim
            do i=1,nmol
                  f(i,j)=0.0
            enddo
      enddo

c     main loop for force calculation

      do i=1,nmol-1

c           find distance between particles

            do j=i+1,nmol

                  rtmpx=q(j,1)-q(i,1)
                  rtmpy=q(j,2)-q(i,2)

c                 enforce periodic boundaries

                  efr(j,1)=rtmpx-boxx*sign(1.0,rtmpx)*
     1            int(abs(rtmpx*boxxi)+.50)
                  efr(j,2)=rtmpy-boxy*sign(1.0,rtmpy)*
     1            int(abs(rtmpy*boxyi)+.50)

                  efd2(j)=efr(j,1)**2+efr(j,2)**2
c                 write(*,*)'efr',efr(j,1),efr(j,2)
c                 write(*,*)'efd2(',j,')=',efd2(j)
            enddo

c           calculate force on each particle

            do j=i+1,nmol

c                 check for distances above cut off potential

                  if (efd2(j).le.rcut2) then

c                       distance within interaction calculate dU/dr

                        r2i=1.0/efd2(j)
                        r6i=r2i*r2i*r2i
                        r12i=r6i*r6i
                        e1=a*r6i
                        e2=b*r12i
                        pe=pe+e1+e2+delta
c                       write(*,*)'fpe=',pe
                        dudr=6.0*(e1+2.0*e2)*r2i

c                       calculate forces

                        do k=1,ndim
                              ftmp=dudr*efr(j,k)
                              f(j,k)=f(j,k)+ftmp
                              f(i,k)=f(i,k)-ftmp
                        enddo
                  endif
            enddo
      enddo
      return
      end

      subroutine dqint(q,boxx,boxy,nside,nmol,v,f)
c     procedure to initialize positions
c     in a 2-d hexagonal fashion

c     declaration of variables

      implicit real*4(a-h,o-z)
      dimension q(64,3),v(64,3),f(64,3)
      integer rand
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
