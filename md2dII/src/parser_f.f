
C-----------------------------------------------------------
C     THIS PART CONTAINS THE CODE FOR THE PARSER.  IT IS NOT
C     IMPORTANT TO KNOW HOW THIS PART WORKS.  YOU SHOULD 
C     ALSO IGNORE THE C ROOUTINES.
C-----------------------------------------------------------
      subroutine parser(parser_flag)

      include 'md2d.h'
      integer inst
      integer parser_flag
      inst = PARSER_ERROR

 191  read(*,31,end=192) string_buffer

         call copyparse(inst, string_buffer,double_buffer,int_buffer)

         if(inst.eq.PARSER_ERROR) then
c           write(*,*) 'parser error', string_buffer
c           stop
         else if(inst.eq.QUITTING) then
            write(*,*) 'quitting program'
            stop
         else if (inst.eq.RUN) then
            if(int_buffer.eq.1) then
               parser_flag=int(double_buffer(1))
               goto 192
            else
               goto 192
            endif
         else if (inst.eq.SET_SEED) then
            if(int_buffer.eq.1) then
               myseed=int(double_buffer(1))
               call srand(myseed)
            else
               write(*,*) 'parse error: #arg in time.'
            endif
         else if (inst.eq.COMMENT) then
c           write(*,*) string_buffer  
         else if (inst.eq.TEMPERATURE_COEFF) then
            if(int_buffer.eq.1) then
               temperature_in=double_buffer(1)  
               sqrttmp = sqrt(temperature_in)
c             print*, 'temperature', temperature_in
c             stop
               temperature=temperature_in
            else
               write(*,*) 'parse error: #arg in time.'
            endif
         else if (inst.eq.BOX_COEFF) then
            if(int_buffer.eq.2) then
               boxlength(XX) = double_buffer(XX)    
               boxlength(YY) = double_buffer(YY)    
            else
               write(*,*) 'parse error: #arg in box.'
            endif
         else if (inst.eq.TIME_COEFF) then
            if(int_buffer.eq.3) then
               time = double_buffer(1)
               time_step = double_buffer(2)
               stop_time = double_buffer(3)
            else if(int_buffer.eq.4) then
               time = double_buffer(1)
               time_step = double_buffer(2)
               stop_time = double_buffer(3)
               interval_output=int(double_buffer(4))
               itemp=int((stop_time-time)/time_step)
C              print*, time_step, itemp/interval_output
C              stop
               if((itemp/interval_output).gt.MAXOUTPUT) then
                  write(*,*) "number of output cycle exceed max"
                  stop
               endif
            else
               write(*,*) 'parse error: #arg in box.'
            endif
         else if (inst.eq.ANDERSEN_COEFF) then

            if(int_buffer.eq.1) then
               ncollision = int(double_buffer(1))    
               if(ncollision.eq.0) then
                  iscale_freq=1
               endif
            else
               write(*,*) 'parse error: #arg in box.'
            endif

         else if (inst.eq.VELOCITY_SCALE_COEFF) then
c           write(*,*) 'this feature needs to be turn on 
c    1                  at compile time'
            if(int_buffer.eq.2) then
               iscale_int= int(double_buffer(1))    
               iscale_max= int(double_buffer(2))    
c              print*,"iscale_int=",iscale_int
c              stop
            else
               write(*,*) 'parse error: #arg in box.'
            endif
           
         else if (inst.eq.ATOM_INPUT) then
            if(int_buffer.eq.4) then
               natoms = natoms + 1
               if (natoms.gt.MAXATOM) then
                  write(*,*) 'exceed max atom. quitting ...'
                  stop
               endif
               size(natoms) = double_buffer(1)
               sqrtmass(natoms) = double_buffer(2)
               pos(XX,natoms) = double_buffer(3)
               pos(YY,natoms) = double_buffer(4)
            else if(int_buffer.eq.7.or.int_buffer.eq.9) then
               itemp = int(double_buffer(1))
               if(itemp.gt.natoms) then
                  write(*,*) '#atom ', itemp, 'is not defined'
               else
                  size(itemp) = double_buffer(2)
                  sqrtmass(itemp) = double_buffer(3)
                  pos(XX,itemp) = double_buffer(4)
                  pos(YY,itemp) = double_buffer(5)
                  vel(XX,itemp) = double_buffer(6)
                  vel(YY,itemp) = double_buffer(7)
                  irestartQ = NOTNULL
               endif
               if(itemp.eq.natoms.and.parser_flag.gt.0) then
                  goto 193
               endif
            else
               write(*,*) 'parse error ATOM: incorrect number arg.'
            endif
         else if (inst.eq.CONFIG_OUTFILE) then
            config=string_buffer
         else if (inst.eq.RESTART_OUTFILE) then
            restart=string_buffer
         else if (inst.eq.THERMO_OUTFILE) then
            thermo=string_buffer
         endif

         goto 191      

 192  continue

      parser_flag = 0

 193  continue

      end
