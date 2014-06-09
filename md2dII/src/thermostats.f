

      subroutine scaling_thermostat() 
      include 'md2d.h'

      ene_kinetic=ZERO     

      if(iscale_max.gt.0) then

         if(iscale_freq.eq.0) then
         iscale_max = iscale_max - 1

         do i = 1,natoms
            vbar = sqrttmp*sqrtmass(i)
            do j =1, DIM
               vel(j,i) = vbar*gausspark(myseed)
c              print*, "problem", vel(j,i)
c              stop
            enddo
         enddo

c        call velocity_scaling()

         else

            do i=1, natoms
               temp_mi = sqrtmass(i)*sqrtmass(i)
               vel(j,i) = vel(j,i) + 
     1                 0.5*force(j,i)*time_step*temp_mi
            enddo

         endif

         iscale_freq=mod(iscale_freq+1,iscale_int)

c        print*, "scaled"
c        stop
      else
         iscale_freq=1
      endif
c     print*, iscale_freq

      do i = 1, natoms

         do j=1,DIM
            temp = vel(j,i)/sqrtmass(i)
            ene_kinetic = ene_kinetic + temp*temp        
         enddo

      enddo

      end

      subroutine andersen_thermostat() 

      include 'md2d.h'


      v = ncollision*time_step

      ene_kinetic=ZERO     

      if(ncollision.gt.0) then
         do i = 1, natoms
            temp_mi = sqrtmass(i)*sqrtmass(i)
            do j = 1,DIM
               vel(j,i) = vel(j,i) + 
     1                 0.5*force(j,i)*time_step*temp_mi
            enddo
         enddo
      endif


      do i = 1, natoms
         vbar = sqrttmp*sqrtmass(i)
         if (rand().lt.v) then
c          print*, ncollision, teme_step, v
c           stop
            do j = 1, DIM

C gaussian distribution not normalized

               vel(j,i) = vbar*gausspark(myseed)

            enddo
         endif
         do j=1,DIM
            temp = vel(j,i)/sqrtmass(i)
            ene_kinetic = ene_kinetic + temp*temp        
         enddo
      enddo


      end

C--------------------------------------------------------------------

      subroutine velocity_scaling()

      include 'md2d.h'
c     declaration of variables


      dimension vc(DIM)

c     main loop for randomizing velocities

c     if (nside.gt.0) then
      do i=1,natoms
         do j=1, DIM

c           pick random velocity from gaussian distribution


C gaussian distribution not normalized

            vel(j,i) = sqrttmp*sqrtmass(i)*gausspark(myseed)
c           print*, "vel", i, vel(j,i), sqrttmp

c           x=rand(0)
c           y=rand(0)
c           vel(j,i) = sqrttmp*sqrtmass(i)*sqrt(-2.0*log(x))*cos(twopi*y)  

         enddo
      enddo

c     endif

c     determine center of mass drift velocities
c     and correct for the drift called by vint

      call vave(vc)

c     loop to set CM velocity of system to zero

      do i=1,natoms

            do j=1,DIM
                  vel(j,i)=vel(j,i)-vc(j)
c                 print*, "problem v=", vel(j,i)
            enddo

      enddo

      return
      end
c--------------------------------------------------------------


      subroutine vave(vc)
c     precedure to compute CM drift velocity

c     declaraction of variables
      include 'md2d.h'
      dimension vc(DIM)

      do i=1,DIM
         vc(i)=0.0
      enddo

      do i=1,natoms

         do j=1,DIM
            vc(j)=vc(j)+vel(j,i)
         enddo

      enddo

      do i=1,DIM
         vc(i)=vc(i)/float(natoms)
      enddo

      return

      end
