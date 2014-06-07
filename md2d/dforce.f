      subroutine dforce(f,q,boxx,boxy,nmol,ndim,pe,delta,sigma6
     1,epsilon,rcut2)

c     declaraction of variables

      implicit real*4(a-h,o-z)
      dimension f(64,3),q(64,3),efd2(64),efr(64,3)

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
