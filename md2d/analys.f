c     Analyse results from md2d

      implicit real*4(a-h,o-z)
      dimension t(5000),q(64,3,5000),v(64,3,5000)
      character*20 irsp
      dimension nwhich(64)
      dimension gofr(500)

      write(*,*) 'Total number of particles?'
      read (*,*) nmol
      open (11,file='md2d.ou2')
      open (12,file='md2d.ou3')
      do i=1,1000
         read (11,*,end=88) t(i),(q(j,1,i),q(j,2,i),j=1,nmol)
         read (12,*,end=88) t(i),(v(j,1,i),v(j,2,i),j=1,nmol)
      enddo
88    continue
      npts = i-1
c     open (21,file='analys.ou1',status='new')
      open (21,file='analys.ou1')

99    continue
      write(*,*) ''
      write(*,*) 'Choose one:'
      write(*,*) '   [1] Plot position of particle(s)'
      write(*,*) '   [2] Mean squared displacement'
      write(*,*) '   [3] Autocorrelate velocity'
      write(*,*) '   [4] g(r)'
      read (*,*) iopt
      write(*,*) ''
      write(*,*) 'You chose:',iopt,'.  Is this right?(y/n)'
      read (*,1001) irsp
1001  format(a1)
      if (irsp.eq.'n'.or.irsp.eq.'N') goto 99

      if (iopt.eq.1) then
         write(*,*) 'How many?(0 for all)'
         read (*,*) nmany
         if (nmany.ne.0) then
            write(*,*) 'Which one(s)?'
            read (*,*) (nwhich(i),i=1,nmany)
         else
            nmany=nmol
            do i=1,nmany
               nwhich(i) = i
            enddo
         endif
         do i=1,nmany
            write(21,*) ' '
            k = nwhich(i)
            write(21,1002) q(k,1,1),q(k,2,1)
1002        format(2e20.8,'   mov')
            do j=2,npts
               write(21,1003) q(k,1,j),q(k,2,j)
1003           format(2e20.8,'   lin')
            enddo
         enddo

      elseif (iopt.eq.2) then
         write(*,*) 'Upper bound?(0 for max)'
         read (*,*) nub
         if (nub.eq.0) nub=npts-1
         do k=0,nub
            tmp=0.
            do i=1,nmol
               do j=1,npts-k
                  tmp=tmp+(q(i,1,j)-q(i,1,j+k))**2
     1            +(q(i,2,j)-q(i,2,j+k))**2
               enddo
            enddo
            tmp=tmp/real(nmol)/real(npts-k)
            write(21,1004) t(1+k)-t(1),tmp
1004        format(2e20.8)
         enddo

      elseif (iopt.eq.3) then
         write(*,*) 'Upper bound?(0 for max)'
         read (*,*) nub
         if (nub.eq.0) nub=npts-1
         do k=0,nub
            tmp=0.
            do i=1,nmol
               do j=1,npts-k
                  tmp=tmp+(v(i,1,j)*v(i,1,j+k))
     1            +(v(i,2,j)*v(i,2,j+k))
               enddo
            enddo
            tmp=tmp/real(nmol)/real(npts-k)
            write(21,1004) t(1+k)-t(1),tmp
         enddo

      elseif (iopt.eq.4) then
         write(*,*) 'number of binx, bin size?'
         read (*,*) nbin,dr
         write(*,*) 'boxx,boxy?'
         read (*,*) boxx,boxy
         boxxi=1./boxx
         boxyi=1./boxy
         do k=1,npts
            do i=1,nmol
               do j=i+1,nmol
                  rtmpx=q(j,1,k)-q(i,1,k)
                  rtmpy=q(j,2,k)-q(i,2,k)

                  delx0=rtmpx-boxx*sign(1.0,rtmpx)*
     1            int(abs(rtmpx*boxxi)+.50)
                  dely0=rtmpy-boxy*sign(1.0,rtmpy)*
     1            int(abs(rtmpy*boxyi)+.50)
                  delr0=sqrt(delx0**2+dely0**2)
                  ibin=delr0/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx1=delx0
                  dely1=dely0+boxy
                  delr1=sqrt(delx1**2+dely1**2)
                  ibin=delr1/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx2=delx0+boxx
                  dely2=dely0+boxy
                  delr2=sqrt(delx2**2+dely2**2)
                  ibin=delr2/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx3=delx0+boxx
                  dely3=dely0
                  delr3=sqrt(delx3**2+dely3**2)
                  ibin=delr3/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx4=delx0+boxx
                  dely4=dely0-boxy
                  delr4=sqrt(delx4**2+dely4**2)
                  ibin=delr4/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx5=delx0
                  dely5=dely0-boxy
                  delr5=sqrt(delx5**2+dely5**2)
                  ibin=delr5/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx6=delx0-boxx
                  dely6=dely0-boxy
                  delr6=sqrt(delx6**2+dely6**2)
                  ibin=delr6/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx7=delx0-boxx
                  dely7=dely0
                  delr7=sqrt(delx7**2+dely7**2)
                  ibin=delr7/dr+.50
                  gofr(ibin)=gofr(ibin)+1

                  delx8=delx0-boxx
                  dely8=dely0+boxy
                  delr8=sqrt(delx8**2+dely8**2)
                  ibin=delr8/dr+.50
                  gofr(ibin)=gofr(ibin)+1

               enddo
            enddo
         enddo
         rho=real(nmol)/boxx/boxy
         do k=1,nbin
            r=(real(k)-.5)*dr
            rnmlz=rho/(2.*3.141592654*r*dr)
            gofr(k)=gofr(k)*rnmlz/real(npts)/9.
         enddo
         write(21,1004) ((real(k)-.5)*dr,gofr(k),k=1,nbin)

      endif

      end
