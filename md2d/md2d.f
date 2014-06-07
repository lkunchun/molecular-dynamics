c     program to simulation a system of 2-d Lennard-Jones spheres
c     Ref: David Chandler, Lecture Notes, Stat. Mech. of Liquids.

c     declare variables and types

      implicit real*4(a-h,o-z)
      integer wtype,errind,dcunit,lx,ly
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
      read (2,*) iseed
      close(2)
      jseed=rand(iseed)


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
      pe=0.d0

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

      iseed = rand(0)*(2.D0**31-1.D0)
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
