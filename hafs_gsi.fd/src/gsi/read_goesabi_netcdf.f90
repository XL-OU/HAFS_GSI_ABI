subroutine read_goesabi_netcdf(val_img,ithin,rmesh,jsatid,gstime,&
     infile,lunout,obstype,nread,ndata,nodata,twind,sis,nobs)


!$$$  subprogram documentation block
!                .      .    .                                       .
!   2021-10-25  X. Lu and X. Wang - modified from read_goesimg.f90 to read in
!                   NCEI goes abi, POC: xuguang.wang@ou.edu
!
!
!   input argument list:
!     mype     - mpi task id
!     val_img  - weighting factor applied to super obs
!     ithin    - flag to thin data
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!
!   output argument list:
!     nread    - number of BUFR GOES imager observations read
!     ndata    - number of BUFR GOES imager profiles retained for further processing
!     nodata   - number of BUFR GOES imager observations retained for further processing
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
      checkob,finalcheck,score_crit
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,txy2ll,tll2xy,rlats,rlons
  use constants, only: deg2rad,zero,one,rad2deg,r60inv,r60
  use radinfo, only: iuse_rad,jpch_rad,nusis !,nst_gsi,nstinfo
  use gsi_nstcouplermod, only: nstinfo, nst_gsi
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen
  use deter_sfc_mod, only: deter_sfc
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth, gsi_nstcoupler_deter
  use netcdf
  use obsmod, only: iadate,time_window_max
  use satthin, only: radthin_time_info,tdiff2crit
  implicit none

  include 'netcdf.inc'
  
! Declare passed variables
  character(len=*),intent(in   ) :: infile,obstype,jsatid
  character(len=20),intent(in  ) :: sis
  integer(i_kind) ,intent(in   ) :: lunout,ithin
  integer(i_kind) ,intent(inout) :: ndata,nodata
  integer(i_kind) ,intent(inout) :: nread
  real(r_kind)    ,intent(in   ) :: rmesh,gstime,twind
  real(r_kind)    ,intent(inout) :: val_img
! Declare local parameters
  integer(i_kind),parameter:: nimghdr=13
  integer(i_kind),parameter:: maxinfo=43 !ajohnson edit
!  integer(i_kind),parameter:: maxinfo=37 !ajohnson edit
  real(r_kind),parameter:: r360=360.0_r_kind
  real(r_kind),parameter:: tbmin=50.0_r_kind
  real(r_kind),parameter:: tbmax=550.0_r_kind
  character(80),parameter:: hdrgoes  = &            ! goes imager header
        'SAID YEAR MNTH DAYS HOUR MINU SECO CLAT CLON SAZA SOZA BEARAZ SOLAZI'

! Declare local variables
  logical outside,iuse,assim
  integer(i_kind) :: nobs !ajohnson edit
  character(8) subset

!debug
integer :: ndims,dimids(20),testdim1,testdim2


  integer(i_kind) nchanl,ilath,ilonh,ilzah,iszah,irec,next
  integer(i_kind) nmind,lnbufr,idate,ilat,ilon, sthin
  integer(i_kind) ireadmg,ireadsb,iret,nreal,nele,itt
  integer(i_kind) itx,i,k,isflg,kidsat,n,iscan,idomsfc, ii
  integer(i_kind) idate5(5)
  integer(i_kind),allocatable,dimension(:)::nrec

  real(r_kind) dg2ew,sstime,tdiff,t4dv,sfcr
  real(r_kind) dlon,dlat,timedif,crit1,dist1,crit0,timeinflat
  real(r_kind) dlon_earth,dlat_earth, thislon, thislat
  real(r_kind) pred, boxcwp3,  boxcwp5
  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr
  real(r_kind),allocatable,dimension(:,:):: data_all

  real(r_double),dimension(nimghdr) :: hdrgoesarr       !  goes imager header
  real(r_double),dimension(3,6) :: dataimg              !  goes imager data

  real(r_kind) cdist,disterr,disterrmax,dlon00,dlat00
  integer(i_kind) ntest
  integer(i_kind) mins_ob,mins_an,ithin_time,n_tbin,it_mesh
  real(r_kind) ptime

!------------------
!  NETCDF-RELATED
!------------------
   INTEGER(i_kind)   :: ncdfID
   INTEGER(i_kind)   :: status
   INTEGER(i_kind)   :: datestlen, varID, DimID
   INTEGER(i_kind)   :: vardim, natts
   INTEGER(i_kind)   :: vartype, id_time
   INTEGER(i_kind)   :: numdim, numvars, numatt
   INTEGER(i_kind)   :: xx,yy,mx,my,nobsin

!------------------
!  SATELLITE DATA
!------------------
   REAL(r_kind), ALLOCATABLE, DIMENSION( :,: )       :: tb_all
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch8
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch9
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch10
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch11
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch12
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch13
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch14
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch15
!   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: ch16
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: cldcod
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: cldheight
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: cldpres
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: cldfrac
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: cldprob
   INTEGER(i_kind), ALLOCATABLE, DIMENSION( : )    :: phase 
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: vza
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: sza
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: lon
   REAL(r_kind), ALLOCATABLE, DIMENSION( : )       :: lat
   Integer(i_kind), ALLOCATABLE, DIMENSION( :,: )       :: obtime
   Integer(i_kind), DIMENSION(5)       :: obdate
!if (mype .eq. mype_root) then !@#$%

!**************************************************************************
! Initialize variables
  sthin = 2 !THIN OBSERVATIONS:  1 = no thin, 2 = every other grid point, etc

  lnbufr = 10
  disterrmax=zero
  ntest=0
  dg2ew = r360*deg2rad

  ilon=3
  ilat=4

  if (nst_gsi > 0 ) then
     call gsi_nstcoupler_skindepth(obstype, zob)         ! get penetration depth (zob) for the obstype
  endif

  rlndsea(0) = zero
  rlndsea(1) = 15._r_kind
  rlndsea(2) = 10._r_kind
  rlndsea(3) = 15._r_kind
  rlndsea(4) = 30._r_kind

  ndata=0
  nodata=0
  nchanl=10       ! the channel number
  ilath=8        ! the position of latitude in the header
  ilonh=9        ! the position of longitude in the header
  ilzah=10       ! satellite zenith angle
  iszah=11       ! solar zenith angle

! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_img=zero


! Make thinning grids

  call radthin_time_info(obstype, jsatid, sis, ptime, ithin_time)
  if( ptime > 0.0_r_kind) then
!    ptime=0.5 !every half hour?
!    ithin_time=3
    n_tbin=nint(2*time_window_max/ptime)
  else
    n_tbin=1
  endif
write(6,*) 'MAKEGRIDS: ',rmesh,ithin,n_tbin,ptime,ithin_time,time_window_max
  call makegrids(rmesh,ithin,n_tbin=n_tbin)


! OPEN NETCDF FILE
write(6,*)"infile: ",trim(infile)
write(6,*)"itxmax/nstinfo: ",itxmax,nstinfo
flush(6)
status = nf90_open(TRIM(infile), NF90_NOWRITE, ncdfID)
print*, '*** OPENING GOES RADIANCE NETCDF FILE', status
!------------------------
! Get global attributes
!-------------------------

!status = nf90_get_att( ncdfID, nf90_global, 'year', idate5(1) )
!status = nf90_get_att( ncdfID, nf90_global, 'month', idate5(2) )
!status = nf90_get_att( ncdfID, nf90_global, 'day', idate5(3) )
!status = nf90_get_att( ncdfID, nf90_global, 'hour', idate5(4) )
!status = nf90_get_att( ncdfID, nf90_global, 'minute', idate5(5) )

status = nf90_inq_varid( ncdfID, 'date', varID )
status = nf90_get_var( ncdfID, varID, idate5 )

print*, idate5

!------------------------
! Get Dimension Info (2-D)
!-------------------------
status = NF_INQ_DIMID(ncdfID, 'nobs', varID)
print *, "nobs status/varID: ",status,varID
status = NF_INQ_DIMLEN(ncdfID, varID, nobsin)

print *, "nobsin, status: ",nobsin,status

!------------------------
! Allocate data arrays
!-------------------------
ALLOCATE( obtime( nobsin,5) )
ALLOCATE( lat( nobsin ) )
ALLOCATE( lon( nobsin ) )
ALLOCATE( sza( nobsin ) )
ALLOCATE( vza( nobsin ) )
ALLOCATE( cldcod( nobsin ) )
ALLOCATE( cldheight( nobsin ) )
ALLOCATE( cldpres( nobsin ) )
ALLOCATE( cldfrac( nobsin ) )
ALLOCATE( cldprob( nobsin ) )
ALLOCATE( phase( nobsin ) )
ALLOCATE( tb_all( nobsin,10 ) )
!------------------------
! Get useful data arrays
!-------------------------
! Time
status = nf90_inq_varid( ncdfID, 'time', varID )
print *, "time status: ",status
status = nf90_get_var( ncdfID, varID, obtime )
print *, "time status: ",status
! LAT
status = nf90_inq_varid( ncdfID, 'lat', varID )
print *, "lat status: ",status
status = nf90_get_var( ncdfID, varID, lat )
print *, "lat status: ",status
! LON
status = nf90_inq_varid( ncdfID, 'lon', varID )
print *, "lon status: ",status
status = nf90_get_var( ncdfID, varID, lon )
print *, "lon status: ",status
! VZA
status = nf90_inq_varid( ncdfID, 'vza', varID )
print *, "vza status: ",status
status = nf90_get_var( ncdfID, varID, vza )
print *, "vza status: ",status

! SZA
status = nf90_inq_varid( ncdfID, 'sza', varID )
print *, "sza status: ",status
status = nf90_get_var( ncdfID, varID, sza )
print *, "sza status: ",status

! Cloud optical depth
status = nf90_inq_varid( ncdfID, 'cldcod', varID )
print *, "cldcod status: ",status
status = nf90_get_var( ncdfID, varID, cldcod )
print *, "cldcod status: ",status

! Cloud top height
status = nf90_inq_varid( ncdfID, 'cldheight', varID )
print *, "cldheight status: ",status
status = nf90_get_var( ncdfID, varID, cldheight )
print *, "cldheight status: ",status

! Cloud top pressure
status = nf90_inq_varid( ncdfID, 'cldpres', varID )
print *, "cldpres status: ",status
status = nf90_get_var( ncdfID, varID, cldpres )
print *, "cldpres status: ",status

! Cloud fraction
status = nf90_inq_varid( ncdfID, 'cldfrac', varID )
print *, "cldfrac status: ",status
status = nf90_get_var( ncdfID, varID, cldfrac )
print *, "cldfrac status: ",status

! Cloud probability
status = nf90_inq_varid( ncdfID, 'cldprob', varID )
print *, "cldprob status: ",status
status = nf90_get_var( ncdfID, varID, cldprob )
print *, "cldprob status: ",status

! Cloud Phase (use to determine clear-sky pixels)
!status = nf90_inq_varid( ncdfID, 'cloud_phase', varID )
!status = nf90_get_var( ncdfID, varID, phase )


! ABI brightness temperatures
status = nf90_inq_varid( ncdfID, 'value', varID )
print *, "tb_all status: ",status
status = nf90_inquire_variable( ncdfID, varID, ndims=ndims,dimids=dimids )
print *, "status,ndims,dimids(1:ndims): ",status,ndims,dimids(1:ndims)
status = nf90_get_var( ncdfID, varID, tb_all )
print *, "tb_all status: ",status


status = NF_INQ_DIMLEN(ncdfID, dimids(1), testdim1)
status = NF_INQ_DIMLEN(ncdfID, dimids(2), testdim2)
print *, "testdim1,testdim2: ",testdim1,testdim2

! CLOSE NETCDF FILE
status = nf90_close( ncdfID )

! SAT ID
  if(jsatid == 'g08') kidsat = 252
  if(jsatid == 'g09') kidsat = 253
  if(jsatid == 'g10') kidsat = 254
  if(jsatid == 'g11') kidsat = 255
  if(jsatid == 'g12') kidsat = 256
  if(jsatid == 'g13') kidsat = 257
  if(jsatid == 'g14') kidsat = 258
  if(jsatid == 'g15') kidsat = 259
  if(jsatid == 'g16') kidsat = 257

! Allocate arrays to hold all data for given satellite
  nreal = maxinfo + nstinfo + 5 !ajohnson: added 5 for all sky research

  nele  = nreal   + nchanl
  allocate(data_all(nele,nobsin),nrec(nobsin))

  next=0
  nrec=999999
  irec=0
  ii=1
!print *, "mxmy: ",mx,my
! Big loop over bufr file
read_loop:  do xx=1, nobsin

     irec=irec+1
     next=next+1
     
        nread=nread+nchanl
        
       ! print*, xx,yy,irec, wv68(xx,yy), cwp(xx,yy), phase(xx,yy)
        
        
!      first step QC filter out data with less clear sky fraction
        !if (hdrgoesarr(1) == 256_r_double .and. dataimg(2,3) < 70.0_r_kind) cycle read_loop
        !if (hdrgoesarr(1) == 254_r_double .and. dataimg(2,3) < 40.0_r_kind) cycle read_loop
        !if (hdrgoesarr(ilzah) >r60) cycle read_loop
     

!       Convert obs location from degrees to radians
        thislon=lon(xx)
        thislat=lat(xx)

        if (thislon >= r360) thislon=thislon-r360
        if (thislon < zero)  thislon=thislon+r360

!if ((thislat.gt.0.1).and.(thislon.gt.0.1).and.(thislat .le. 10e9).and.(thislon.le.10e9)) print *, "latlon: ",thislat,thislon

        dlon_earth=thislon*deg2rad
        dlat_earth=thislat*deg2rad

!       If regional, map obs lat,lon to rotated grid.
        if(regional)then

!          Convert to rotated coordinate.  dlon centered on 180 (pi), 
!          so always positive for limited area
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
!if((.not. outside).and.(dlat .gt. 400)) print *, "location: ",dlon_earth,dlat_earth,dlon,dlat,outside
           if(diagnostic_reg) then
              call txy2ll(dlon,dlat,dlon00,dlat00)
              ntest=ntest+1
              cdist=sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                   (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00))
              cdist=max(-one,min(cdist,one))
              disterr=acos(cdist)*rad2deg
              disterrmax=max(disterrmax,disterr)
           end if

!          Check to see if in domain.  outside=.true. if dlon_earth,
!          dlat_earth outside domain, =.false. if inside
           if(outside) cycle read_loop

!       Global case
        else
           dlon=dlon_earth
           dlat=dlat_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

!        obdate(1)=obtime(xx)/100000000
!        obdate(2)=(obtime(xx)-obdate(1)*100000000)/1000000
!        obdate(3)=(obtime(xx)-obdate(1)*100000000-obdate(2)*1000000)/10000
!        obdate(4)=(obtime(xx)-obdate(1)*100000000-obdate(2)*1000000-obdate(3)*10000)/100
!        obdate(5)=obtime(xx)-obdate(1)*100000000-obdate(2)*1000000-obdate(3)*10000-obdate(4)*100
        obdate(1)=obtime(xx,1)
        obdate(2)=obtime(xx,2)
        obdate(3)=obtime(xx,3)
        obdate(4)=obtime(xx,4)
        obdate(5)=obtime(xx,5)
        call w3fs21(obdate,mins_ob)
        call w3fs21(iadate,mins_an)
        t4dv = real((mins_ob-iwinbgn),r_kind) * r60inv
        sstime = real(mins_ob,r_kind)
        tdiff=(sstime-gstime)*r60inv
!        t4dv = 5
!        tdiff= 1
!        print *,"Xu tdiff",t4dv,sstime,tdiff
!        print *,"Xu time:",obdate,iadate
        if (l4dvar.or.l4densvar) then
           if (t4dv<zero .OR. t4dv>winlen) cycle read_loop
        else
           if (abs(tdiff)>twind) cycle read_loop
        endif
        crit0=0.01_r_kind
        timeinflat=6.0_r_kind
        call tdiff2crit(tdiff,ptime,ithin_time,timeinflat,crit0,crit1,it_mesh)
!        print *,"Xu it_mesh:",it_mesh
        call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis,it_mesh=it_mesh)
        if(.not. iuse)cycle read_loop

!        if (l4dvar) then
!           crit1=0.01_r_kind
!        else
!           timedif = mins_ob-mins_an !6.0_r_kind*abs(tdiff)        ! range:  0 to 18
!           crit1=0.01_r_kind+timedif
!        endif
!        t4dv=timedif*r60inv
!        print *,"Xu t4dv",t4dv,xx
!        call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis)
!       if(.not. iuse)cycle read_loop


!       Locate the observation on the analysis grid.  Get sst and land/sea/ice
!       mask.  
!     isflg    - surface flag
!                0 sea      LARC: 1 or 6
!                1 land     LARC: 4
!                2 sea ice  LARC: 2 or 7
!                3 snow     LARC: 0
!                4 mixed                 
       ! SET DEFAULT TO LAND
       !!isflg = 1        
       !!if ( phase(xx,yy).eq.1 .or. phase(xx,yy).eq.6 ) isflg = 0
       !!if ( phase(xx,yy).eq.4 ) isflg = 1 
       !!if ( phase(xx,yy).eq.2 .or. phase(xx,yy).eq.7 ) isflg = 2
       !!if ( phase(xx,yy).eq.0 ) isflg = 3


      ! TEMPORARY...set isflg to 1 (all land)
      isflg = 1
        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
            ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)

!if (sfcpct(0) .lt. 0) then
!  sfcpct(0) = 0.0
!else if (sfcpct(0) .gt. 1) then
!  sfcpct(0) = 1.0
!endif


!if (sfcpct(1) .lt. 0) then
!  sfcpct(1) = 0.0
!else if (sfcpct(1) .gt. 1) then
!  sfcpct(1) = 1.0
!endif


!if (sfcpct(2) .lt. 0) then
!  sfcpct(2) = 0.0
!else if (sfcpct(2) .gt. 1) then
!  sfcpct(2) = 1.0
!endif


!if (sfcpct(3) .lt. 0) then
!  sfcpct(3) = 0.0
!else if (sfcpct(3) .gt. 1) then
!  sfcpct(3) = 1.0
!endif

!if ((sfcpct(0)+sfcpct(1)+sfcpct(2)+sfcpct(3) < 0.99999999) .or. (sfcpct(0)+sfcpct(1)+sfcpct(2)+sfcpct(3) > 1.00000001))then
! write(6,*) 'sfcpct: ',sfcpct
! write(6,*) 'dlat,dlon,isflg: ',dlat,dlon,isflg
! write(6,*) 'dlat_earth,dlon_earth: ',dlat_earth,dlon_earth
! write(6,*) 'thislon,thislat: ',thislon,thislat
! write(6,*) 'idomsfc:',idomsfc
!endif 
!write(6,*) 'dlat,dlon,isflg,sfcpct: ',dlat,dlon,isflg,sfcpct

!       Set common predictor parameters
        crit1=crit1+rlndsea(isflg)
        call checkob(dist1,crit1,itx,iuse)
!        if(.not. iuse)cycle read_loop

!       Set data quality predictor 
!        pred =(10.0_r_kind-dataimg(2,1)/10.0_r_kind)+dataimg(3,3)*10.0_r_kind  ! clear sky and
                                                                 ! bt std as quality indicater
                                                                                                                    
!       Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"
!        crit1 = crit1+pred 
!        call finalcheck(dist1,crit1,itx,iuse)
!        if(.not. iuse) cycle read_loop

!       Map obs to grids
!        iscan = nint(hdrgoesarr(ilzah))+1.001_r_kind ! integer scan position
        
!
!       interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
        if ( nst_gsi > 0 ) then
           tref  = ts(0)
           dtw   = zero
           dtc   = zero
           tz_tr = one
           if ( sfcpct(0) > zero ) then
              call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
           endif
        endif
!print *,"Xu use?",iuse,outside,ii,xx

 if ((.not.outside).and. iuse) then  
!          print*, "ob",xx,yy, dlat, dlon, irec, wv68(xx,yy), cwp(xx,yy), sza(xx,yy), boxcwp5
!print *, "t4dv: ",t4dv 
           !if (wv68(xx,yy).gt.0.1) print*, xx,yy, dlat, dlon, irec, wv68(xx,yy), cwp(xx,yy), phase(xx,yy), sza(xx,yy), outside, iuse, sfcpct(:)
         
!       Transfer information to work array
        data_all( 1,ii) = kidsat                    ! satellite id
        data_all( 2,ii) = t4dv                       ! analysis relative time
        data_all( 3,ii) = dlon                       ! grid relative longitude
        data_all( 4,ii) = dlat                       ! grid relative latitude
        data_all( 5,ii) = vza(xx)*deg2rad         ! satellite zenith angle (radians)
        data_all( 6,ii) = zero             ! satellite azimuth angle (radians)
        data_all( 7,ii) = 999 !phase(xx)               ! clear sky amount
        data_all( 8,ii) = iscan                      ! integer scan position
        data_all( 9,ii) = sza(xx)                 ! solar zenith angle
        data_all(10,ii) = zero             ! solar azimuth angle
        data_all(11,ii) = sfcpct(0)                  ! sea percentage of
        data_all(12,ii) = sfcpct(1)                  ! land percentage
        data_all(13,ii) = sfcpct(2)                  ! sea ice percentage
        data_all(14,ii) = sfcpct(3)                  ! snow percentage
        data_all(15,ii)= ts(0)                       ! ocean skin temperature
        data_all(16,ii)= ts(1)                       ! land skin temperature
        data_all(17,ii)= ts(2)                       ! ice skin temperature
        data_all(18,ii)= ts(3)                       ! snow skin temperature
        data_all(19,ii)= tsavg                       ! average skin temperature
        data_all(20,ii)= vty                         ! vegetation type
        data_all(21,ii)= vfr                         ! vegetation fraction
        data_all(22,ii)= sty                         ! soil type
        data_all(23,ii)= stp                         ! soil temperature
        data_all(24,ii)= sm                          ! soil moisture
        data_all(25,ii)= sn                          ! snow depth
        data_all(26,ii)= zz                          ! surface height
        data_all(27,ii)= idomsfc + 0.001_r_kind      ! dominate surface type
        data_all(28,ii)= sfcr                        ! surface roughness
        data_all(29,ii)= ff10                        ! ten meter wind factor
        data_all(30,ii)= dlon_earth*rad2deg          ! earth relative longitude (degrees)
        data_all(31,ii)= dlat_earth*rad2deg          ! earth relative latitude (degrees)

!        data_all(36,ii) = val_img !ajohnson edit
!        data_all(37,ii) = itt !ajohnson edit

        data_all(42,ii) = val_img !ajohnson edit
        data_all(43,ii) = itt !ajohnson edit

        if ( nst_gsi > 0 ) then
           data_all(maxinfo+1,ii) = tref         ! foundation temperature
           data_all(maxinfo+2,ii) = dtw          ! dt_warm at zob
           data_all(maxinfo+3,ii) = dtc          ! dt_cool at zob
           data_all(maxinfo+4,ii) = tz_tr        ! d(Tz)/d(Tr)
        endif

!       Transfer observation location and other data to local arrays
!*        do k=1,nchanl
!*           data_all(k+31,ii) = 2.0
!*           !data_all(k+31,ii)=dataimg(3,k+1)		!STD DEVS
!*           !data_all(k+nreal,itx)=dataimg(1,k+1)     !TBs
 !*       end do
        
        !Transfer Observations to local array
        do k=1,nchanl
          data_all(k+nreal,ii) = tb_all(xx,k)
        enddo
!ajohnson: I added space for these 5 variables for all sky research:
        data_all(nreal-4,ii)=cldcod(xx)
        data_all(nreal-3,ii)=cldheight(xx)
        data_all(nreal-2,ii)=cldpres(xx)
        data_all(nreal-1,ii)=cldfrac(xx)
        data_all(nreal,ii)=cldprob(xx)

        score_crit(itx) = 1.0
        
        nrec(itx)=irec
        ndata = ndata+1
        ii = ii+1 

      endif

  enddo read_loop  !xx loop 


!  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
!     nele,itxmax,nread,ndata,data_all,score_crit,nrec)

! Write final set of "best" observations to output file
print*, obstype,sis,nreal,nchanl,ilat,ilon, k, nele, ndata
!print*,   data_all(:,:)
  
  call count_obs(ndata,nele,ilat,ilon,data_all,nobs)
  write(lunout) obstype,sis,nreal,nchanl,ilat,ilon 
  write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
!  write(lunout) data_all(1:nele,1:ndata)
 
write(6,*) 'nele,ndata:',nele,ndata
!do xx=1,ndata
!write(6,*) ndata,data_all(nreal+3,xx),data_all(nreal+5,xx),data_all(nreal+7,xx)
!enddo


! Deallocate local arrays
  deallocate(data_all,nrec)

  DEALLOCATE(lat)
  DEALLOCATE(lon)
  DEALLOCATE(vza)
  DEALLOCATE(sza)
  DEALLOCATE(cldcod)
  DEALLOCATE(cldheight)
  DEALLOCATE(cldpres)
  DEALLOCATE(cldfrac)
  DEALLOCATE(cldprob)
  DEALLOCATE(phase)
  !DEALLOCATE(nir16)
  DEALLOCATE(tb_all)
  
   print*, '*** FINISHED READING RADIANCE NETCDF FILE'
   
! Deallocate satthin arrays
900 continue
  call destroygrids

!endif !@#$%

!  if(diagnostic_reg.and.ntest>0) write(6,*)'READ_GOESIMG_NETCDF:  ',&
!     'mype,ntest,disterrmax=',mype,ntest,disterrmax

  return
end subroutine read_goesabi_netcdf
