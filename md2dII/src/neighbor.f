     
 
      integer function iupdateQ()

      include 'md2d.h'
      data neighbor_update /0/
      save 
      neighbor_update = mod(neighbor_update,NEUPDATE)

c     if((ave_disp2/float(natoms)).ge.(DISPMAX)) then
         ave_disp2 = ZERO

      if(neighbor_update.eq.0) then
         iupdateQ=0
      else
         iupdateQ=1
      endif

      neighbor_update = neighbor_update+1
         
      end

C-------------------------------------------------------------

      subroutine build_neighbor(iflag)
      include 'md2d.h'
      dimension boxi(DIM)

      if(iupdateQ().eq.0.or.iflag.eq.0) then

c       write(*,*) 'nb', time, iflag, 
c    1     pos(1,1), pos(2,1)
        do i = 1, DIM
           boxi(i) = 1.0/boxlength(i)
        enddo

        do i=1, natoms

           ncount = 1

           do j=i+1, natoms

              distance = ZERO

              do k=1, DIM
                 temp = (pos(k,i)-pos(k,j))
                 temp = temp - nint(temp*boxi(k))*boxlength(k)
                 distance=distance + temp*temp
              enddo

              distance = sqrt(distance)
c             print*, distance

              if(distance.le.NEIGHBOR_CUT) then
                 neighbor_list(ncount,i) = j
                 ncount = ncount + 1

c                if(i.eq.45.and.j.eq.53) then
c                   print*, 'nb', time
c                endif

                 if(ncount.gt.MAX_NEIGHBOR) then
                    write(*,*) "density too high, try increasing
     1                          MAX_NEIGHBOR"
                    write(*,*) i, j
                    stop
                 endif
              endif
             
           enddo

           neighbor_list(ncount,i) = NULL

        enddo


      endif

      return
      
      end



C-------------------------------------------------------------
