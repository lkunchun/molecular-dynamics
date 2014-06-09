      program md2d
      include 'md2d.h'
      integer parser_flag
      j=0

C     ===============================
C     INITIALIZATION ROUTINE 
C     READIN DATA AND JUMPSTART THE VELOCITY 
      call init(parser_flag)

C     -------------------------------

C     INITIALIZE THE NEIGHBOR ROUTINE
      call build_neighbor(0)


C     -------------------------------

C     INITIALIZE THE FORCES FOR THE FIRST CYCLE
      call dforce()


C     ===============================



C     LOOP STARTS HERE
C     ----------------

      do i = 0, MAXCYCLE

C     ===============================
C     ADD OUTPUT ROUTINES HERE

c        if(time.gt.22.79) then
c           print*, pos(1,45), pos(2,45), vel(1,45),vel(2,45)
c           print*, pos(1,53), pos(2,53), vel(1,53),vel(2,53)
c        endif

         if(mod(i,interval_output).eq.0) then
            j=mod(j,ENE_OUTPUT)
            if(mod(j,ENE_OUTPUT).eq.0) then   
               call print_atoms()
            else
               call print_ene()
            endif
            j=j+1
         endif
C     ===============================

         if(time.gt.stop_time) goto 191
         call ballview(natoms,pos,size,boxlength,vel)
         call move()
         call build_neighbor(1)
         call dforce()
         call andersen_thermostat()
c        call scaling_thermostat()
         time = time + time_step

C     ===============================
C     ADD ANALYSIS ROUTINES HERE


         ave_ene=ave_ene+ene_pot/natoms
         ave_ene2=ave_ene2+ene_pot*ene_pot/natoms/natoms
         pressure=pressure/2/temperature_in/
     1            boxlength(XX)/boxlength(YY)
         ave_pres=ave_pres+pressure
         ave_pres2=ave_pres2+pressure*pressure
         

C     ===============================
      enddo
 191  continue


      print*, '#@md2d finished', natoms
      print*, '@run'

      end
