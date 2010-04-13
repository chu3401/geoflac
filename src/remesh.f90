 
subroutine re_mesh

include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

common /remesh/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt
dimension cordo(mnz+1,mnx+1,2)
allocatable :: dummy(:,:)



! Save old mesh for interpolations
cordo = cord

! Create The New grid (cord) using cordo(nz,i,2) for the bottom and cordo(1,i,2) for the surface
call rem_cord(cordo)

! REMESHING FOR ACCRETION
if(iac_rem.eq.1) then
!  NODE WISE
do i= 1,iinj-1
   do j= 1,nz
      vel(j,i,1) = vel(j,i+1,1)
      vel(j,i,2) = vel(j,i+1,2)
enddo
enddo
do i = nx,iinj+2,-1
   do j =1,nz
      vel(j,i,1) = vel(j,i-1,1)
      vel(j,i,2) = vel(j,i-1,2)
enddo
enddo
!do j= 1,nz
!      vel(j,iinj+1,1) = vel(j,iinj,1)
!      vel(j,iinj+1,2) = vel(j,iinj,2)
!enddo
! ELEMENT WISE
do k = 1,4
do l = 1,4
do i= 1,iinj-2
   do j= 1,nz
      stress0(k,l,j,i) = stress0(k,l,j,i+1)
enddo
enddo
enddo
enddo
do k = 1,3
do i= 1,iinj-2
   do j= 1,nz
      strain(k,j,i) = strain(k,j,i+1)
enddo
enddo
enddo
   do i= 1,iinj-2
   do j= 1,nz
      visn(j,i)= visn(j,i+1)
      aps(j,i) = aps(j,i+1)
      source(j,i) = source(j,i+1)
enddo
enddo

do k = 1,4
do l = 1,4
do i= nx-1,iinj+2,-1
   do j= 1,nz
      stress0(k,l,j,i) = stress0(k,l,j,i-1)
enddo
enddo
enddo
enddo
do k = 1,3
do i= nx-1,iinj+2,-1
   do j= 1,nz
      strain(k,j,i) = strain(k,j,i-1)
enddo
enddo
enddo
   do i= nx-1,iinj+2,-1
   do j= 1,nz
      visn(j,i)= visn(j,i-1)
      aps(j,i) = aps(j,i-1)
      source(j,i) = source(j,i-1)
enddo
enddo


goto 1017
endif

! REMESHING FOR ELEMENT-WISE PROPERTIES
! Linear interpolation in baricentric coordinates defined as centers of old mesh
nxt = nx-1
nzt = nz-1
allocate( dummy(nzt,nxt) )

! Old mesh - old-element centers
! New mesh - new-element centers
do i = 1, nx-1
    do j = 1, nz-1
        cold(j,i,1) = 0.25*( cordo(j,i,1)+cordo(j+1,i,1)+cordo(j,i+1,1)+cordo(j+1,i+1,1) )
        cold(j,i,2) = 0.25*( cordo(j,i,2)+cordo(j+1,i,2)+cordo(j,i+1,2)+cordo(j+1,i+1,2) )
        cnew(j,i,1) = 0.25*( cord(j,i,1)+cord(j+1,i,1)+cord(j,i+1,1)+cord(j+1,i+1,1) )
        cnew(j,i,2) = 0.25*( cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2) )
    enddo
enddo

! Calculate parameters of old-mesh triangles
call rem_trpars

! Baricentric coordinates of new-elements centers
call rem_barcord


! Do interpolations

! Interpolate Stress (in quadralaterals) 
do k = 1,4
    do l = 1,4
        dummy(1:nzt,1:nxt) = stress0(k,l,1:nzt,1:nxt)
        call rem_interpolate( dummy )
        stress0(k,l,1:nzt,1:nxt) = dummy(1:nzt,1:nxt)
    end do
end do
!do k = 1,4
!do l = 1,4
!do j = 1, nz
!   do i = 1,  3
!	stress0(k,l,j,i) = stress0(k,l,j,4)
!enddo
!enddo
!enddo
!enddo

! HOOK
! Remeshing mode 3 - see user_luc.f90
if( mode_rem .eq. 5 ) call rem_stress_alt()


! Interpolate strains
do k = 1, 3
    dummy(1:nzt,1:nxt) = strain(k,1:nzt,1:nxt)
    call rem_interpolate( dummy )
    strain(k,1:nzt,1:nxt) = dummy(1:nzt,1:nxt)
end do


! plastic strain
dummy(1:nzt,1:nxt) = aps(1:nzt,1:nxt)
call rem_interpolate( dummy )
do i = 1, nxt
    do j = 1, nzt
        if( dummy(j,i) .ge. 0. ) then
            aps(j,i) = dummy(j,i)
        else
            aps(j,i) = 0.
        endif
!       write(*,*) i,j,aps(j,i)
    end do
end do
        
! viscosity
dummy(1:nzt,1:nxt) = visn(1:nzt,1:nxt)
call rem_interpolate( dummy )
visn(1:nzt,1:nxt) = dummy(1:nzt,1:nxt)
! phases
!iph_int = 1 
!dummy(1:nzt,1:nxt) = phasez(1:nzt,1:nxt)
!call rem_interpolate( dummy )
!phasez(1:nzt,4:nxt) = dummy(1:nzt,4:nxt)
do jj = 1,4
   do ii = 1,3
        phasez(jj,ii) = 11.  !7 is ocean coming from the left!
        aps(jj,ii) = 0.0
        do kk = 1,2
        ! Calculate triangle number in which the markers belong
         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk

                                                                                
         call newphase2marker(ii,jj,ntriang)
         enddo
         do kll = 1, nphasl
              if (lphase(kll).eq.11) then
                  phase_ratio(jj,ii,kll) = 1.0
              else
                  phase_ratio(jj,ii,kll) = 0.0
              endif
         enddo
   enddo
enddo
do jj = 5,6
   do ii = 1,3
        phasez(jj,ii) = 12.  !7 is ocean coming from the left!
        aps(jj,ii) = 0.0
        do kk = 1,2
        ! Calculate triangle number in which the markers belong
         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk

                                                                                
         call newphase2marker(ii,jj,ntriang)
         enddo
         do kll = 1, nphasl
              if (lphase(kll).eq.12) then
                  phase_ratio(jj,ii,kll) = 1.0
              else
                  phase_ratio(jj,ii,kll) = 0.0
              endif
         enddo
   enddo
enddo
do jj = 7,10
   do ii = 1,3
        phasez(jj,ii) = 6.  !7 is ocean coming from the left!
        aps(jj,ii) = 0.0
        do kk = 1,2
        ! Calculate triangle number in which the markers belong
         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk

                                                                                
         call newphase2marker(ii,jj,ntriang)
         enddo
         do kll = 1, nphasl
              if (lphase(kll).eq.6) then
                  phase_ratio(jj,ii,kll) = 1.0
              else
                  phase_ratio(jj,ii,kll) = 0.0
              endif
         enddo
   enddo
enddo
do jj = 11,12 
   do ii = 1,3
        phasez(jj,ii) = 9.  !7 is ocean coming from the left!
        aps(jj,ii) = 0.0
        do kk = 1,2
        ! Calculate triangle number in which the markers belong
         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk

                                                                                
         call newphase2marker(ii,jj,ntriang)
         enddo
         do kll = 1, nphasl
              if (lphase(kll).eq.9) then
                  phase_ratio(jj,ii,kll) = 1.0
              else
                  phase_ratio(jj,ii,kll) = 0.0
              endif
         enddo
   enddo
enddo
do jj = 13,nz-3 
   do ii = 1,3
        phasez(jj,ii) = 8.  !7 is ocean coming from the left!
        aps(jj,ii) = 0.0
        do kk = 1,2
        ! Calculate triangle number in which the markers belong
         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk

                                                                                
         call newphase2marker(ii,jj,ntriang)
         enddo
         do kll = 1, nphasl
              if (lphase(kll).eq.8) then
                  phase_ratio(jj,ii,kll) = 1.0
              else
                  phase_ratio(jj,ii,kll) = 0.0
              endif
         enddo
   enddo
enddo

!do jj = 4,10
!   do ii = 1,3
!      if (time/sec_year/1.e6.gt.1000..and.time/sec_year/1.e6.lt.1030.) then
!         icon = 1.
!      else
!         icon = 8.
!      endif
!      phasez(jj,ii) = 8.
!      aps(jj,ii) = 0.0
!      do kk = 1,2
         ! Calculate triangle number in which the markers belong
!         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk
!                                                                                 
!         call newphase2marker(i,j,ntriang)
!      enddo
!      do kl = 1, nphasl
!         if (lphase(kl).eq.8) then
!            phase_ratio(jj,ii,kl) = 1.0
!         else
!            phase_ratio(jj,ii,kl) = 0.0
!         endif
!      enddo
!1   enddo
!enddo
!do jj = 11,87 
!   do ii = 1,3
!        phasez(jj,ii) = 8.
!        aps(jj,ii) = 0.0
!        do kk = 1,2
!        ! Calculate triangle number in which the markers belong
!         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk
!
!                                                                                
!         call newphase2marker(ii,jj,ntriang)
!         enddo
!         do kll = 1, nphasl
!              if (lphase(kll).eq.8) then
!                  phase_ratio(jj,ii,kll) = 1.0
!              else
!                  phase_ratio(jj,ii,kll) = 0.0
!              endif
!         enddo
!enddo
!enddo
!

do jj = nz-2,nz-1 
   do ii = 1,nx-1
        phasez(jj,ii) = 8.
        aps(jj,ii) = 0.0
        do kk = 1,2
        ! Calculate triangle number in which the markers belong
         ntriang = 2 * ( (nz-1)*(ii-1)+jj-1) + kk

                                                                                
         call newphase2marker(ii,jj,ntriang)
         enddo
         do kll = 1, nphasl
              if (lphase(kll).eq.8) then
                  phase_ratio(jj,ii,kll) = 1.0
              else
                  phase_ratio(jj,ii,kll) = 0.0
              endif
         enddo
enddo
enddo

!do i = 1, nxt
!    do j = 1, nzt
!        iph = iphase(phasez(j,i) )
!        phasez(j,i) = float(iph)
!    end do
!end do
!iph_int = 0
! sources
dummy(1:nzt,1:nxt) = source(1:nzt,1:nxt)
call rem_interpolate( dummy )
source(1:nzt,1:nxt) = dummy(1:nzt,1:nxt)


deallocate( dummy )


! REMESHING FOR NODE-WISE PROPERTIES
! Linear interpolation in baricentric coordinates of old mesh
nxt = nx
nzt = nz
allocate( dummy(nzt,nxt) )

! Old mesh - old coordinates points
cold = cordo
temp0 = temp
! New mesh - new coordinates points
cnew = cord
if (iac_rem.eq.1) cold = cnew
! Calculate parameters of triangles of this mesh
call rem_trpars

! Baricentric coordinates of new-elements centers
call rem_barcord

! Do node-wise interpolations

! Velocities (in nodes)
do k = 1, 2
    dummy(1:nzt,1:nxt) = vel(1:nzt,1:nxt,k)
    call rem_interpolate( dummy )
    vel(1:nzt,1:nxt,k) = dummy(1:nzt,1:nxt)
end do

! Temperatures (in nodes) 
dummy(1:nzt,1:nxt) = temp(1:nzt,1:nxt)
call rem_interpolate( dummy )
temp(1:nzt,1:nxt) = dummy(1:nzt,1:nxt)
deallocate( dummy )
iynts = 20
call init_temp
!do i = 1,4 
 !  do j = 1, 8 
!	phasez(j,i) = 2. 
! enddo
!   do j = 9, nz-1 
!	phasez(j,i) = 4. 
! enddo
!enddo
! AFTER INTERPOLATIONS - RECALCULATE SOME DEPENDENT VARIABLES
1017 continue
! Calculation of areas of triangle
call init_areas
if (iac_rem.eq.1) call init_phase
if (iac_rem.eq.1) call init_temp
! reinitialize the stress in the 3 middle element if iac_rem 1
if (iac_rem.eq.1) then
do 522 i = iinj-1,iinj+1
     rogh = 0.
   do 522 j = 1,nz-1
     iph = iphase(phasez (j,i))
     if (iph.eq.0) goto 522
        tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
        densT = den(iph) * ( 1 - alfa(iph)*tmpr )
        dh1 = cord (j,i  ,2) - cord (j+1,i  ,2)
        dh2 = cord (j,i+1,2) - cord (j+1,i+1,2)
        dh  = 0.5 * (dh1+dh2)
        dPT = densT * g * dh

        dP = dPT * ( 1 - beta(iph)*rogh ) / ( 1 + beta(iph)/2*dPT )

        press = rogh + 0.5*dP
        do ii = 1,4
            stress0 (1,ii,j,i) = -press
            stress0 (2,ii,j,i) = -press
            stress0 (3,ii,j,i) = 0.
            stress0 (4,ii,j,i) = -press
        end do
        rogh = rogh + dP
        aps(j,i) = 0.
        do k = 1,3
            strain(k,j,i) = 0.
        enddo
522  continue
endif

! Distribution of masses in nodes
call rmasses

! 1) Determine Inertial Masses with a given dt_scale (idt_scale=1)
! 2) dt_mech with a given  Real Masses (idt_scale = 0)
call dt_mass

! drop the time step to the smallest one
dt = min(dt_elastic, dt_maxwell)


return 
end




!===============================================
! parameters of triangles of a grid
!===============================================
subroutine rem_trpars
include 'precision.inc'
include 'arrays.inc'

common /remesh/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt


do i = 1,nxt-1
    do j = 1,nzt-1
        do k = 1,2

            !  diagonal / :
            !   ii=1     ii=2 

            !  1---3         1
            !  | /         / |
            !  2         2---3

            if (k.eq.1) then
                x1 = cold(j  ,i  ,1)
                x2 = cold(j+1,i  ,1)
                x3 = cold(j  ,i+1,1)
                y1 = cold(j  ,i  ,2)
                y2 = cold(j+1,i  ,2)
                y3 = cold(j  ,i+1,2)
            else  !if (k.eq.2) then 
                x1 = cold(j  ,i+1,1) 
                x2 = cold(j+1,i  ,1) 
                x3 = cold(j+1,i+1,1)
                y1 = cold(j  ,i+1,2) 
                y2 = cold(j+1,i  ,2) 
                y3 = cold(j+1,i+1,2) 
            endif 
 
            det=( (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2) )

            n = 2*( (nzt-1)*(i-1)+j-1 ) + k

            !Find the parameters ONLY for 2 vertices
            pt(n,1,1)=(x2*y3-y2*x3)/det
            pt(n,1,2)=(y2-y3)/det
            pt(n,1,3)=(x3-x2)/det
            pt(n,2,1)=(x3*y1-y3*x1)/det
            pt(n,2,2)=(y3-y1)/det
            pt(n,2,3)=(x1-x3)/det
        end do
    end do
end do     

return
end



!===============================================
! baricentric coordinates of new mesh in old triangles
!===============================================
subroutine rem_barcord
include 'precision.inc'
include 'arrays.inc'

common /remesh/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt


perr = 1.e-4

do i = 1, nxt
    do j = 1, nzt
        xx = cnew(j,i,1)
        yy = cnew(j,i,2)

        amodmin = 1.e+10

        numtr(j,i) = 0
        do l = 0, max( nxt-1, nzt-1 )
            do lt = 1, max(8*l,1)
                if( lt .le. 2*l ) then
                    jo = j - l
                    io = i + lt-1 - l
                elseif( lt .le. 4*l ) then
                    jo = j + lt-2*l - l - 1
                    io = i + l
                elseif( lt .le. 6*l ) then
                    jo = j + l
                    io = i + lt-4*l - l
                elseif( lt .le. 8*l ) then
                    jo = j + lt-6*l - l
                    io = i - l
                else ! only at l=0
                    jo = j
                    io = i
                endif

                if( io.lt.1 .or. io.gt.nxt-1 ) cycle
                if( jo.lt.1 .or. jo.gt.nzt-1 ) cycle

                do k = 1, 2
                    n = 2*( (nzt-1)*(io-1)+jo-1 ) + k
                    a1 = pt(n,1,1) + xx*pt(n,1,2) + yy*pt(n,1,3)
                    a2 = pt(n,2,1) + xx*pt(n,2,2) + yy*pt(n,2,3)
                    a3 = 1. - a1 - a2
                    amod = a1*a1 + a2*a2 + a3*a3
                    if( amod .lt. amodmin ) then
                        amodmin = amod
                        nmin = n
                    endif
                    if ( a1.le.-perr .or. a1.ge.(1.0+perr) ) cycle
                    if ( a2.le.-perr .or. a2.ge.(1.0+perr) ) cycle
                    if ( a3.le.-perr .or. a3.ge.(1.0+perr) ) cycle
                    numtr(j,i) = n
                    barcord(j,i,1) = a1
                    barcord(j,i,2) = a2
                    barcord(j,i,3) = a3
                    goto 10
                end do

            end do
        end do

10      continue
    
    ! if a point is outside of the mesh
    if( numtr(j,i) .eq. 0 ) then
        numtr(j,i) = nmin
        numqu = (numtr(j,i) + 1)/2          
        io = (numqu-1)/(nzt-1) + 1
        jo = numqu - (io-1)*(nzt-1)

        if( mod(numtr(j,i),2).eq.0 ) then
            dist1 = sqrt((xx-cold(jo,io+1,1))**2 + (yy-cold(jo,io+1,2))**2)
            dist2 = sqrt((xx-cold(jo+1,io,1))**2 + (yy-cold(jo+1,io,2))**2)
            dist3 = sqrt((xx-cold(jo+1,io+1,1))**2 + (yy-cold(jo+1,io+1,2))**2)
        else
            dist1 = sqrt((xx-cold(jo,io,1))**2 + (yy-cold(jo,io,2))**2)
            dist2 = sqrt((xx-cold(jo+1,io,1))**2 + (yy-cold(jo+1,io,2))**2)
            dist3 = sqrt((xx-cold(jo,io+1,1))**2 + (yy-cold(jo,io+1,2))**2)
        endif
        
        ! select two nearest points for interpolation
        if( dist1 .gt. dist2 ) then
            if( dist2 .gt. dist3 ) then
                ! 1,2,3
                barcord(j,i,1) = 0.
                barcord(j,i,2) = 1./dist2/(1./dist2+1./dist3)
                barcord(j,i,3) = 1./dist3/(1./dist2+1./dist3)
            else
                if( dist1 .gt. dist3 ) then
                    ! 1,3,2
                    barcord(j,i,1) = 0.
                    barcord(j,i,2) = 1./dist2/(1./dist2+1./dist3)
                    barcord(j,i,3) = 1./dist3/(1./dist2+1./dist3)
                else
                    ! 3,1,2
                    barcord(j,i,1) = 1./dist1/(1./dist1+1./dist2)
                    barcord(j,i,2) = 1./dist2/(1./dist1+1./dist2)
                    barcord(j,i,3) = 0.
                endif
            endif
        else
            if( dist3 .gt. dist1 ) then
                if( dist2 .gt. dist3 ) then
                    ! 2,3,1
                    barcord(j,i,1) = 1./dist1/(1./dist1+1./dist3)
                    barcord(j,i,2) = 0.
                    barcord(j,i,3) = 1./dist3/(1./dist1+1./dist3)
                else
                    ! 3,2,1
                    barcord(j,i,1) = 1./dist1/(1./dist1+1./dist2)
                    barcord(j,i,2) = 1./dist2/(1./dist1+1./dist2)
                    barcord(j,i,3) = 0.
                endif
            else
                ! 2,1,3
                barcord(j,i,1) = 1./dist1/(1./dist1+1./dist3)
                barcord(j,i,2) = 0.
                barcord(j,i,3) = 1./dist3/(1./dist1+1./dist3)
            endif
        endif

    endif

    end do
end do

return
end                    


!===============================================
! interpolation
!===============================================
subroutine rem_interpolate( arr )
include 'precision.inc'
include 'arrays.inc'
include 'params.inc'
common /remesh/ pt(mnz*mnx*2,2,3),barcord(mnz+1,mnx+1,3), &
cold(mnz+1,mnx+1,2),cnew(mnz+1,mnx+1,2),numtr(mnz+1,mnx+1),nzt,nxt
dimension dummy(nzt,nxt),arr(nzt,nxt)


dummy = arr

do i = 1, nxt
    do j = 1, nzt

        numq = (numtr(j,i)+1) / 2
        io = (numq-1)/(nzt-1) + 1
        jo = numq - (io-1)*(nzt-1)

        !  diagonal / :
        !   kk=1       kk=2
        !
        !  1---3         1
        !  | /         / |
        !  2         2---3

        if( mod(numtr(j,i),2).eq.0 ) then
            f1 = dummy(jo  ,io+1)
            f2 = dummy(jo+1,io  )
            f3 = dummy(jo+1,io+1)
        else
            f1 = dummy(jo  ,io  )
            f2 = dummy(jo+1,io  )
            f3 = dummy(jo  ,io+1)
        endif

        arr(j,i) = barcord(j,i,1)*f1 + barcord(j,i,2)*f2 + barcord(j,i,3)*f3
!  For the phases make sure that you do not create new phases
!        iph_int = 0
!       if(iph_int.eq.1) then
!     xmax1 = max(int(f1+0.5),int(f2+0.5))
!        xmax2 = int(f3+0.5)
!        xmax = max(xmax1,xmax2)
!        xmin1 = min(int(f1+0.5),int(f2+0.5))
!        xmin2 = int(f3+0.5)
!        xmin = min(xmin1,xmin2)
!        xnorm = abs(xmax -xmin)
!        xavr = abs(xmax+xmin)*0.5
!        if (xnorm.eq.0.) goto 132
!        if(arr(j,i).lt.xavr) then 
!        xvalnorm = abs((xmin-arr(j,i))/xnorm)
!        else
!        xvalnorm = abs((xmax-arr(j,i))/xnorm)
!        endif
!        if (xvalnorm.ge.0.5) arr(j,i) = xmax-xvalnorm
!        if (xvalnorm.lt.0.5) arr(j,i) = xmin+xvalnorm
!132     arr(j,i) = xmin
!        if(arr(j,i).eq.0.) arr(j,i) = 1.
!        endif
    end do
end do

return

end
