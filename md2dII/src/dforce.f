
C     sigma(i,j) calculates the average size of i and j
C

      double precision function sigma(i,j)

      include 'md2d.h'

         sigma = (size(i)+size(j))*0.5

      end

C-------------------------------------------------------------

      subroutine dforce()

      include 'md2d.h'

C     DEFINE LENNARD JONE PARAMETERS

      parameter(ljeps=4.0d0)
c     parameter(ljcut=1.122462048309372981433533049677)
      parameter(ljcut=3.6)
      dimension boxi(DIM)
      dimension dist(MAX_NEIGHBOR,DIM)
      dimension r2(MAX_NEIGHBOR)

c     initialization of constants

      do i = 1, DIM
         boxi(i)=1.0/boxlength(i)
      enddo

      ene_pot = ZERO
      pressure = ZERO

      do i=1,natoms
         do j=1,DIM
            force(j,i)=ZERO
         enddo
      enddo

c     main loop for force calculation

      do i=1, natoms-1
         ncounter = 0
         do k=1, MAX_NEIGHBOR

            j = neighbor_list(k,i)

            if(j.eq.NULL) then
               goto 91
            endif

c           if(j.le.i) then
c              print*, 'problem', j, i
c              stop
c           endif

            ncounter = ncounter + 1
c           find distance between particles

            r2(k) = ZERO
            do m=1, DIM

               dist(k,m)=pos(m,i)-pos(m,j)
               dist(k,m) = dist(k,m) - nint(dist(k,m)
     1                     *boxi(m))*boxlength(m)
               r2(k) = r2(k) + dist(k,m)*dist(k,m)

            enddo

         enddo

 91      continue

         do k = 1, ncounter

            j = neighbor_list(k,i)

            sig = sigma(i,j)

c           print*, 'sig=', sig 
c           stop

            temp_ljcut = sig*ljcut
            temp_ljcut = temp_ljcut*temp_ljcut
            if(r2(k).ge.temp_ljcut) goto 92 

c          if(j.eq.53.and.i.eq.45.and.time.gt.22.79) then
c             print*, time, neighbor_list(i,k),i,r2(k),
c    1        dist(k,1), dist(k,2)
c             print*, pos(1,i), pos(2,i)
c             itemp=neighbor_list(i,k)
c             print*, pos(1,itemp), pos(2,itemp)
c          endif

           if(r2(k).le.(0.49)) then
              call print_atoms()
              write(*,*) 'atoms too close', time, neighbor_list(k,i),  
     1         i, r2(k)
              stop
           endif

            sig2 = sig*sig
            sig6 = sig2*sig2*sig2
            attract = -sig6*ljeps
            sig12 = sig6*sig6
            hard_core = sig12*ljeps


            r2i = 1.0/r2(k)
            r6i=r2i*r2i*r2i
            r12i=r6i*r6i
            e1=attract*r6i
            e2=hard_core*r12i

            ene_pot=ene_pot+e1+e2

c           if((e1+e2).lt.(0.0d0)) then
c           if(i>j) then
c              print*, 'problem', i,j,(e1+e2)  
c              stop
c           endif
 

c           dudr=6.0*(e1+2.0*e2)*r2i
            dudr=6.0*(e1+2.0*e2)

            pressure = pressure + dudr           

            dudr=dudr*r2i
        
c         print*, k, j, sqrt(r2(k)), dudr, e1+e2, sig

            do m = 1, DIM

               temp=dudr*dist(k,m)
               force(m,i)=force(m,i)+temp
               force(m,j)=force(m,j)-temp

C     ===============================================
C     ADD PRESSURE CALCULATION HERE IF NEEDED
C              prssure = temp*dist
C     ===============================================

            enddo

 92         continue

         enddo

      enddo

      return

      end
