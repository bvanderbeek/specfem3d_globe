!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
!
!
!  input is (r, theta, phi), output is the matrix cij(6x6)
!  0 <= r <= 1, 0 <= phi <= 80 , 70 <= theta <= 110
!
!  read cij from text file and interpolate them to the grid
!
!
!--------------------------------------------------------------------------------------------------

  module model_aniso_mantle_cij_par

  ! model_aniso_mantle_cij_variables
  double precision,dimension(:,:,:,:),allocatable :: AMM_V_Cij
  double precision,dimension(:,:),allocatable :: AMM_V_Cije 
  double precision,dimension(:),allocatable :: AMM_V_pro
  double precision, dimension(:), allocatable :: AMM_V_lon
  double precision, dimension(:), allocatable :: AMM_V_colat
  integer :: nx,ny,nz
 end module model_aniso_mantle_cij_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_aniso_mantle_cij_broadcast()

! standard routine to setup model

  use constants, only: myrank,IMAIN,IIN
  use model_aniso_mantle_cij_par

  implicit none

  ! local parameters
  integer :: ier
  character(len=*), parameter :: cij_model = 'DATA/CIJ/cij_model.xyz'
  
  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: cij_model (aniso_mantle from CIJ model)'
    call flush_IMAIN()
   endif


  !Read CIJ model
  open(IIN,file=cij_model,status='old',action='read',iostat=ier)

  if (ier /= 0 ) stop 'Error opening file cij_model.xyz'

! read the number of nodes in the three dimensions
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) nx, ny, nz

888 close(IIN)

  ! allocates model arrays
  ! modify these values according to your needs
  ! Uses custom iso_prem model defined by 20 elements
  allocate(AMM_V_Cij(22,nx,ny,nz), & 
           AMM_V_pro(nz),&
           AMM_V_Cije(3,nz),&
           AMM_V_lon(nx),&
           AMM_V_colat(ny), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating AMM_V arrays')

  ! the variables read are declared and stored in structure AMM_V
  if (myrank == 0) call read_aniso_mantle_model_cij()


  ! broadcast the information read on the master to the nodes
  call bcast_all_singlei(nx)
  call bcast_all_singlei(ny)
  call bcast_all_singlei(nz)
  call bcast_all_dp(AMM_V_Cij,22*nx*ny*nz) !22*nx*ny*nz
  call bcast_all_dp(AMM_V_pro,nz)  !nz
  call bcast_all_dp(AMM_V_lon,nx)  !nx
  call bcast_all_dp(AMM_V_colat,ny)  !ny
  call bcast_all_dp(AMM_V_Cije,3*nz)  !3*nz


   end subroutine model_aniso_mantle_cij_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_aniso_mantle_cij(r,theta,phi,&
                                rho,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use constants, only: PI,GRAV,EARTH_RHOAV,DEGREES_TO_RADIANS,EARTH_R,EARTH_R_KM, &
        IREGION_CRUST_MANTLE,R_UNIT_SPHERE,ZERO
  !use shared_parameters, only: R670

  use model_aniso_mantle_cij_par

  implicit none

  double precision,intent(in) :: r,theta,phi
  double precision,intent(out) :: rho  !(inout)
  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  ! local parameters
  double precision :: lon,colat,A,C,F,L,N,colat_const,lon_const
  double precision :: vp,vs,Qkappa,Qmu
  double precision :: scale_Pa,scaleval


  !Define 1D model below the CIJ model domain
  if (r <= 1.d0 - AMM_V_pro(1)/EARTH_R_KM) then

    scaleval = dsqrt(PI*GRAV*EARTH_RHOAV)
    scale_Pa =(EARTH_RHOAV)*((EARTH_R*scaleval)**2)

    call model_ak135(r,rho,vp,vs,Qkappa,Qmu,IREGION_CRUST_MANTLE)
    !Dimensionalize
    rho = rho*EARTH_RHOAV
    vp = vp*(scaleval * EARTH_R)
    vs = vs*(scaleval * EARTH_R)
    !print *,'rho kg/m3',rho,'vp m/s',vp,'vs m/s',vs

    ! c11 from vp,vs and rho in m/s and kg/m3 is in Pa 
    c11 = (rho*vp*vp)/scale_Pa
    c12 = (rho*(vp*vp-2.d0*vs*vs))/scale_Pa
    c13 = c12
    c14 = 0.d0
    c15 = 0.d0
    c16 = 0.d0
    c22 = c11
    c23 = c12
    c24 = 0.d0
    c25 = 0.d0
    c26 = 0.d0
    c33 = c11
    c34 = 0.d0
    c35 = 0.d0
    c36 = 0.d0
    c44 = (rho*vs*vs)/scale_Pa
    c45 = 0.d0
    c46 = 0.d0
    c55 = c44
    c56 = 0.d0
    c66 = c44

    ! non dimensionalized parameters
    rho = rho/EARTH_RHOAV  ! rho from text file is kg/m3, we need it non-dimens
    ! vp from text file is in m/s, we need it non-dimens 
    vp = vp/(scaleval * EARTH_R) 
    vs = vs/(scaleval * EARTH_R)

  else

    ! above 670 discontinuity
    ! from ~24.4 to 670 km, 0.8<r<1
    ! converts colat/lon to degrees if necessary
    lon = phi / DEGREES_TO_RADIANS
    colat = theta / DEGREES_TO_RADIANS

    ! it reads the model parameters from text file and interpolate them to the grid
    call build_cij_cij(AMM_V_pro,rho,r,colat,lon, &
                   c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36, &
                   c44,c45,c46,c55,c56,c66)

        ! DONT USE IT IF YOU USE ROTATION (RAD to GLOB) IN MESHFEM3D_MODELS.f90
        ! 24/08/2020  in the case of 1D model from matlab or if cij produces tensor in
        ! local (radial) coordinates system use rad_to_glob
        ! call rotate_tensor_radial_to_global(theta,phi,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                        ! c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                        ! c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                        ! c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)


    ! ----------------------------------------------------------
    ! 02/09/2020 modification to create a 1D or 3D radially anisotropic model
    !activate it only if you want to calculate approximated tensor
    IF(1==0) THEN  !calculation radial anisotropic model

        ! select constant values only for 1D cases and put them into build_cij_cij
        colat_const = 90.d0
        lon_const = 40.d0

        call build_cij_cij(AMM_V_pro,rho,r,colat,lon, &
                   c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36, &
                   c44,c45,c46,c55,c56,c66) 
        

        !inizialize values
        A = ZERO
        F = ZERO
        C = ZERO
        L = ZERO
        N = ZERO

        ! calculate Love Parameters (Love 1927)
        A = (3.d0/8.d0)*(c11+c22)+0.25d0*c12+0.5d0*c66
        C = c33
        F = 0.5d0*(c13+c23)
        L = 0.5d0*(c44+c55)
        N = (1.d0/8.d0)*(c11+c22)-0.25d0*c12+0.5d0*c66

        ! radially anisotropic elastic tensor:
        !      A       A-2N   F    0        0        0
        !      A-2N     A    F     0        0        0
        !      F        F    C     0        0        0
        !      0        0    0     L        0        0
        !      0        0    0     0        L        0
        !      0        0    0     0        0        N

        ! calculate back the 21 components of which 5 are independent
        c11=A
        c12=A-(2.d0*N)
        c13=F
        c14=0.d0
        c15=0.d0
        c16=0.d0
        c22=A
        c23=F
        c24=0.d0
        c25=0.d0
        c26=0.d0
        c33=C
        c34=0.d0
        c35=0.d0
        c36=0.d0
        c44=L
        c45=0.d0
        c46=0.d0
        c55=L
        c56=0.d0
        c66=N

        
    END IF  !end of calculation radial anisotropic model
    ! ---------------------------------------------------------

  endif

  end subroutine model_aniso_mantle_cij

!
!-------------------------------------------------------------------------------------------------
!

  subroutine build_cij_cij(pro,rho,r,colat,lon, &
                       d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,d33,d34,d35,d36, &
                       d44,d45,d46,d55,d56,d66)

  use constants, only: ZERO,EARTH_R,EARTH_R_KM,R_UNIT_SPHERE,DEGREES_TO_RADIANS,PI, &
        IREGION_CRUST_MANTLE,GRAV,EARTH_RHOAV
  use model_aniso_mantle_cij_par

  implicit none

  double precision,intent(in) :: pro(nz)  

  double precision,intent(out) :: rho

  double precision,intent(in) :: r,colat,lon
  double precision,intent(out) :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66  

  ! local parameters
 
  double precision :: d11r,d12r,d13r,d14r,d15r,d16r,d22r,d23r,d24r,d25r,d26r, &
                      d33r,d34r,d35r,d36r,d44r,d45r,d46r,d55r,d56r,d66r  
  double precision :: tet,ph,dtet,dph,depth,dout1,dout2,thick_lbz,thick_vbz
  double precision :: rhor,vpr,vsr,PwaveMod,SwaveMod
  double precision :: d1,d2,d3,d4,pc1,pc2,pc3,pc4, &
                   dpr1,dpr2,scale_GPa,scaleval
  double precision :: anispara(22,2,4),elpar(22)
  integer :: idep,ipar,icolat,ilon,out_flag,oneD_flag, &
             ict0,ict1,icp0,icp1,icz0,icz1 

  !3 options for setting the parameters in the empty grid nodes 
  !at the sides of the geodynamic model domain
  !
  !oneD_flag = 0 --> smooth transition to the 1D profile averaged from the geodyanmic model
  !oneD_flag = 1 --> smooth transition to the 1D profile present in specfem3d_globe (ak135)
  !oneD_flag = 2 --> interpolate parameters from lateral boundaries
  oneD_flag = 0

  thick_lbz = 5.d0 ! thickness in degrees of lateral buffer zone
  thick_vbz = 5.d1 ! thickness in km of vertical buffer zone

  !Set very wide buffer zone to force interpolation to stay close to the lateral boundaries
  if (oneD_flag == 2) thick_lbz = 1.d10
 
  !Scaling factors
  scaleval = dsqrt(PI*GRAV*EARTH_RHOAV)
  scale_GPa =(EARTH_RHOAV/1000.d0)*((EARTH_R*scaleval/1000.d0)**2)

  tet = colat !in degrees 
  ph  = lon   !in degrees

  ! avoid edge effects
  ! the grid has to be inside the model chunck
  dout1 = 0; dout2 = 0
  out_flag = 0 !inside the CIJ model domain
  ! West boundary
  if (ph < AMM_V_lon(1)) then 
     !outside the buffer zone
     if(ph <= AMM_V_lon(1) - thick_lbz) then
        out_flag = 2
     !inside the buffer zone
     else
        dout1 = ABS( (AMM_V_lon(1) - ph) / thick_lbz ) 
        out_flag = 1
     end if
     ph = AMM_V_lon(1)+0.00001
  end if

  ! East boundary
  if (ph > AMM_V_lon(nx)) then
     !outside the buffer zone
     if(ph >= AMM_V_lon(nx) + thick_lbz) then
        out_flag = 2
     !inside the buffer zone
     else
        !Update dout if more distant from this boundary than others
        if (dout1 < ABS( (ph - AMM_V_lon(nx)) / thick_lbz ) ) dout1 = ABS( (ph - AMM_V_lon(nx)) / thick_lbz ) 
        out_flag = 1
     end if
     ph = AMM_V_lon(nx)-0.00001
  end if

  ! South boundary
  if (tet < AMM_V_colat(1)) THEN
     !outside the buffer zone
     if(tet <= AMM_V_colat(1) - thick_lbz) then
        out_flag = 2
     !inside the buffer zone
     else
        !Update dout if more distant boundary
        if (dout1 < ABS( (AMM_V_colat(1) - tet) / thick_lbz ) ) dout1 = ABS( (AMM_V_colat(1) - tet) / thick_lbz ) 
        out_flag = 1
     end if
     tet = AMM_V_colat(1)+0.00001
  end if

  ! North boundary
  if (tet > AMM_V_colat(ny)) then
     !outside the buffer zone
     if(tet >= AMM_V_colat(ny) + thick_lbz) then
        out_flag = 2 
     !inside the buffer zone
     else
        !Update dout if more distant boundary
        if (dout1 < ABS( (tet - AMM_V_colat(ny)) / thick_lbz ) ) dout1 = ABS( (tet - AMM_V_colat(ny)) / thick_lbz ) 
        out_flag = 1
     end if
     tet = AMM_V_colat(ny)-0.00001
  end if

  if(out_flag == 1) then
     dout1 = 1.d0 - dout1
     if (dout1 < 0.d0)  call exit_MPI_without_rank('dout1 < 0')
     if (dout1 > 1.d0)  call exit_MPI_without_rank('dout1 > 1')
  end if

! dimensionalize
  depth = EARTH_R_KM*(R_UNIT_SPHERE - r)

  if (depth <= pro(nz) .or. depth >= pro(1)) print *, 'depth',depth,'pro(nz)',pro(nz),'pro(1)',pro(1)
  if (depth <= pro(nz) .or. depth >= pro(1)) call exit_MPI_without_rank('r out of range in build_cij')


  ict0 = 0
  do icolat = 1,ny-1  ! I add -1 to avoid icolat=ny and ict0>ny
    if (AMM_V_colat(icolat) < tet) ict0 = ict0 + 1 
  enddo

  icp0 = 0
  do ilon = 1,nx-1
    if (AMM_V_lon(ilon) < ph) icp0 = icp0 + 1  
  enddo

  icz0 = 0
  do idep = 1,nz     
    if (pro(idep) > depth) icz0 = icz0 + 1
  enddo

  ict1 = ict0 + 1 

  icp1 = icp0 + 1

  icz1 = icz0 + 1


  if (icp0 < 1 .or. icp0 > nx) print *,ph,icp0,nx,AMM_V_lon(1),AMM_V_lon(nx)
  if (icp1 < 1 .or. icp1 > nx) print *,ph,icp0,nx,AMM_V_lon(1),AMM_V_lon(nx)
  if (ict0 < 1 .or. ict0 > ny) print *,tet,ict0,ny,AMM_V_colat(1),AMM_V_colat(ny)
  if (ict1 < 1 .or. ict1 > ny) print *,tet,ict1,ny,AMM_V_colat(1),AMM_V_colat(ny)
  if (icz0 < 1 .or. icz0 > nz) print *,depth,icz0,nz,pro(1),pro(nz)
  if (icz1 < 1 .or. icz1 > nz) print *,depth,icz1,nz,pro(1),pro(nz)
! check that parameters make sense
  if (ict0 < 1 .or. ict0 > ny) call exit_MPI_without_rank('ict0 out of range')
  if (ict1 < 1 .or. ict1 > ny) call exit_MPI_without_rank('ict1 out of range')
  if (icp0 < 1 .or. icp0 > nx) call exit_MPI_without_rank('icp0 out of range')
  if (icp1 < 1 .or. icp1 > nx) call exit_MPI_without_rank('icp1 out of range')
  if (icz0 < 1 .or. icz0 > nz) call exit_MPI_without_rank('icz0 out of range')
  if (icz1 < 1 .or. icz1 > nz) call exit_MPI_without_rank('icz1 out of range')


  !1d density and isotropic velocity profile
  IF(out_flag > 0) THEN
     
    !Example of vertical grid
    !600 km + icz1 (model grid node)
    !       |
    !       |  drp1
    !       - specfem grid node
    !       |
    !       | dpr2
    !       |
    !670 km + icz0 (model grid node)
    !
    dpr1 = (depth - pro(icz1))/(pro(icz0) - pro(icz1))
    dpr2 = 1.0 - dpr1
    if (dpr1<0.d0)  call exit_MPI_without_rank('dpr1 < 0 for 1d profile')
    if (dpr2<0.d0)  call exit_MPI_without_rank('dpr2 < 0 for 1d profile')
    if (dpr1>1.d0)  call exit_MPI_without_rank('dpr1 > 1 for 1d profile')
    if (dpr2>1.d0)  call exit_MPI_without_rank('dpr2 > 1 for 1d profile')

    rhor= AMM_V_Cije(1,icz0)*dpr1 + AMM_V_Cije(1,icz1)*dpr2
    vpr = AMM_V_Cije(2,icz0)*dpr1 + AMM_V_Cije(2,icz1)*dpr2
    vsr = AMM_V_Cije(3,icz0)*dpr1 + AMM_V_Cije(3,icz1)*dpr2

    if(oneD_flag == 1) then
       !Get density and isotorpic Vp, Vs from 1D reference model
       !dph = Qkappa
       !dtet = Qmu 
       call model_ak135(r,rhor,vpr,vsr,dph,dtet,IREGION_CRUST_MANTLE)
       !Dimensionalize
       rhor = rhor*EARTH_RHOAV
       vpr = vpr*(scaleval * EARTH_R)
       vsr = vsr*(scaleval * EARTH_R)
       !print *,'rho kg/m3',rho,'vp m/s',vp,'vs m/s',vs
    end if

    !In GPa. Scaling to specfem GPa units is done at the end of the subroutine
    d11r = rhor*vpr*vpr/1.d9
    d12r = rhor*(vpr*vpr-2.d0*vsr*vsr)/1.d9
    d13r = d12r              
    d14r = ZERO
    d15r = ZERO
    d16r = ZERO
    d22r = d11r    
    d23r = d12r              
    d24r = ZERO
    d25r = ZERO
    d26r = ZERO
    d33r = d11r    
    d34r = ZERO
    d35r = ZERO
    d36r = ZERO
    d44r = rhor*vsr*vsr/1.d9
    d45r = ZERO
    d46r = ZERO
    d55r = d44r    
    d56r = ZERO
    d66r = d44r    

  END IF


  !Nodes in the domain or lateral buffer zone
  IF(out_flag < 2) THEN

  ! intepolate the 22 parameters from AMM_V_Cij from 8 nodes of the integration cell
  do ipar = 1,22
    anispara(ipar,1,1) = AMM_V_Cij(ipar,icp0,ict0,icz0)
    anispara(ipar,2,1) = AMM_V_Cij(ipar,icp1,ict0,icz0)
    anispara(ipar,1,2) = AMM_V_Cij(ipar,icp0,ict0,icz1)
    anispara(ipar,2,2) = AMM_V_Cij(ipar,icp1,ict0,icz1)
    anispara(ipar,1,3) = AMM_V_Cij(ipar,icp0,ict1,icz0)
    anispara(ipar,2,3) = AMM_V_Cij(ipar,icp1,ict1,icz0)
    anispara(ipar,1,4) = AMM_V_Cij(ipar,icp0,ict1,icz1)
    anispara(ipar,2,4) = AMM_V_Cij(ipar,icp1,ict1,icz1)
  enddo

  !
  ! calculation of distances between the selected point and grid points
  !
  dtet = (tet - AMM_V_colat(ict0))/(AMM_V_colat(ict1) - AMM_V_colat(ict0))
  dph  = (ph  - AMM_V_lon(icp0))/(AMM_V_lon(icp1) - AMM_V_lon(icp0)) 
  ! check that parameters make sense
  if (dtet<0.d0)  call exit_MPI_without_rank('dtet < 0')
  if (dph<0.d0)  call exit_MPI_without_rank('dph < 0')
  if (dtet>1.d0)  call exit_MPI_without_rank('dtet > 1')
  if (dph>1.d0)  call exit_MPI_without_rank('dph > 1')

  d1 = (1.0 - dtet)*(1.0 - dph)

  d2 = dtet*(1.0 - dph)

  d3 = (1.0 - dtet)*dph

  d4 = dtet*dph

  !Example of vertical grid
  !600 km + icz1 (model grid node)
  !       |
  !       |  drp1
  !       - specfem grid node
  !       |
  !       | dpr2
  !       |
  !670 km + icz0 (model grid node)
  !
  dpr1 = (depth - pro(icz1))/(pro(icz0) - pro(icz1))
  dpr2 = 1.0 - dpr1
  if (dpr1<0.d0)  call exit_MPI_without_rank('dpr1 < 0 for 3d model')
  if (dpr2<0.d0)  call exit_MPI_without_rank('dpr2 < 0 for 3d model')
  if (dpr1>1.d0)  call exit_MPI_without_rank('dpr1 > 1 for 3d model')
  if (dpr2>1.d0)  call exit_MPI_without_rank('dpr2 > 1 for 3d model')

  do ipar = 1,22
     pc1 = anispara(ipar,1,1)*dpr1+anispara(ipar,1,2)*dpr2
     pc2 = anispara(ipar,1,3)*dpr1+anispara(ipar,1,4)*dpr2
     pc3 = anispara(ipar,2,1)*dpr1+anispara(ipar,2,2)*dpr2
     pc4 = anispara(ipar,2,3)*dpr1+anispara(ipar,2,4)*dpr2
     elpar(ipar) = pc1*d1 + pc2*d2 + pc3*d3 + pc4*d4
  enddo

  !
  !   create dij
  !
  rho = elpar(1)
  d11 = elpar(2)
  d12 = elpar(3)
  d13 = elpar(4)
  d14 = elpar(5)
  d15 = elpar(6)
  d16 = elpar(7)
  d22 = elpar(8)
  d23 = elpar(9)
  d24 = elpar(10)
  d25 = elpar(11)
  d26 = elpar(12)
  d33 = elpar(13)
  d34 = elpar(14)
  d35 = elpar(15)
  d36 = elpar(16)
  d44 = elpar(17)
  d45 = elpar(18)
  d46 = elpar(19)
  d55 = elpar(20)
  d56 = elpar(21)
  d66 = elpar(22)

  !Nodes in the buffer zone
  IF(out_flag == 1) THEN 
 
    dout2 = 1.d0 - dout1

    if (dout2 < 0.d0)  call exit_MPI_without_rank('dout2 < 0')
    if (dout2 > 1.d0)  call exit_MPI_without_rank('dout2 > 1')

    rho = rho*dout1 + rhor*dout2
    d11 = d11*dout1 + d11r*dout2
    d12 = d12*dout1 + d12r*dout2
    d13 = d13*dout1 + d13r*dout2
    d14 = d14*dout1 + d14r*dout2
    d15 = d15*dout1 + d15r*dout2
    d16 = d16*dout1 + d16r*dout2
    d22 = d22*dout1 + d22r*dout2
    d23 = d23*dout1 + d23r*dout2
    d24 = d24*dout1 + d24r*dout2
    d25 = d25*dout1 + d25r*dout2
    d26 = d26*dout1 + d26r*dout2
    d33 = d33*dout1 + d33r*dout2
    d34 = d34*dout1 + d34r*dout2
    d35 = d35*dout1 + d35r*dout2
    d36 = d36*dout1 + d36r*dout2
    d44 = d44*dout1 + d44r*dout2
    d45 = d45*dout1 + d45r*dout2
    d46 = d46*dout1 + d46r*dout2
    d55 = d55*dout1 + d55r*dout2
    d56 = d56*dout1 + d56r*dout2
    d66 = d66*dout1 + d66r*dout2

  END IF

  ! -----------------------------------------------------------------------
  ! 04/06/2020 calculate isotropic elastic tensor if necessary
  IF(1==1) THEN  
        ! elastic tensor for hexagonal symmetry in reduced notation:
        !      c11 c12 c13  0   0        0
        !      c12 c11 c13  0   0        0
        !      c13 c13 c33  0   0        0
        !       0   0   0  c44  0        0
        !       0   0   0   0  c44       0
        !       0   0   0   0   0  c66=(c11-c12)/2
        ! where 
        ! c11=k+4/3G
        ! c22=k+4/3G
        ! c33=k+4/3G
        ! c44=G
        ! c55=G
        ! c66=G
        ! c12=k-2/3G
        ! c13=k-2/3G
        ! c23=k-2/3G

        PwaveMod = (3.d0/15.d0)*(d11+d22+d33)+(2.d0/15.d0)*(d23+d13+d12)+(4.d0/15.d0)*(d44+d55+d66)
        SwaveMod = (1.d0/15.d0)*(d11+d22+d33)-(1.d0/15.d0)*(d23+d13+d12)+(3.d0/15.d0)*(d44+d55+d66)

        d11 = PwaveMod
        d12 = PwaveMod-2.d0*SwaveMod
        d13 = PwaveMod-2.d0*SwaveMod
        d14 = 0.d0
        d15 = 0.d0
        d16 = 0.d0
        d22 = PwaveMod
        d23 = PwaveMod-2.d0*SwaveMod
        d24 = 0.d0
        d25 = 0.d0
        d26 = 0.d0
        d33 = PwaveMod
        d34 = 0.d0
        d35 = 0.d0
        d36 = 0.d0
        d44 = SwaveMod
        d45 = 0.d0
        d46 = 0.d0
        d55 = SwaveMod
        d56 = 0.d0
        d66 = SwaveMod

  END IF
  ! -------------------------------------------------------------

  END IF

  !Nodes outside the buffer zone
  IF(out_flag == 2) THEN 
 
    rho = rhor
    d11 = d11r
    d12 = d12r
    d13 = d13r
    d14 = d14r
    d15 = d15r
    d16 = d16r
    d22 = d22r
    d23 = d23r
    d24 = d24r
    d25 = d25r
    d26 = d26r
    d33 = d33r
    d34 = d34r
    d35 = d35r
    d36 = d36r
    d44 = d44r
    d45 = d45r
    d46 = d46r
    d55 = d55r
    d56 = d56r
    d66 = d66r

  END IF

  ! Vertical buffer zone
  if(depth > pro(1) - thick_vbz) then

    !Get density and isotorpic Vp, Vs from 1D reference model
    !dph = Qkappa
    !dtet = Qmu 
    call model_ak135(r,rhor,vpr,vsr,dph,dtet,IREGION_CRUST_MANTLE)
    !Dimensionalize
    rhor = rhor*EARTH_RHOAV
    vpr = vpr*(scaleval * EARTH_R)
    vsr = vsr*(scaleval * EARTH_R)
    !print *,'rho kg/m3',rho,'vp m/s',vp,'vs m/s',vs

    !Distance from top of vertical buffer zone
    !Example of vertical grid
    !       + pro(1) - thick_vbz (top of vertical buffer zone)    
    !       |
    !       | drp2
    !       - specfem grid node
    !       |
    !       | dpr1
    !       |
    !pro(1) + bottom of vertical buffer zone
    !
    dpr2 = (depth - (pro(1)-thick_vbz))/thick_vbz !weight for CIJ paramters
    dpr1= 1.d0 - dpr2 !weight for 1D reference model parameters
    if (dpr1<0.d0)  call exit_MPI_without_rank('dpr1 < 0 in vbf2')
    if (dpr2<0.d0)  call exit_MPI_without_rank('dpr2 < 0 in vbf2')
    if (dpr1>1.d0)  call exit_MPI_without_rank('dpr1 > 1 in vbf2')
    if (dpr2>1.d0)  call exit_MPI_without_rank('dpr2 > 1 in vbf2')
    !Inteprolate density and elastic moduli
    rho = rho*dpr1 + rhor*dpr2
    d11 = d11*dpr1 + (rhor*vpr*vpr)/1.d9*dpr2
    d12 = d12*dpr1 + (rhor*(vpr*vpr-2.d0*vsr*vsr))/1.d9*dpr2
    d13 = d13*dpr1 + (rhor*(vpr*vpr-2.d0*vsr*vsr))/1.d9*dpr2
    !d14 = 0.d0
    !d15 = 0.d0
    !d16 = 0.d0
    d22 = d22*dpr1 + (rhor*vpr*vpr)/1.d9*dpr2
    d23 = d23*dpr1 + (rhor*(vpr*vpr-2.d0*vsr*vsr))/1.d9*dpr2
    !d24 = 0.d0
    !d25 = 0.d0
    !d26 = 0.d0
    d33 = d33*dpr1 + (rhor*vpr*vpr)/1.d9*dpr2
    !d34 = 0.d0
    !d35 = 0.d0
    !d36 = 0.d0
    d44 = d44*dpr1 + (rhor*vsr*vsr)/1.d9*dpr2
    !d45 = 0.d0
    !d46 = 0.d0
    d55 = d55*dpr1 + (rhor*vsr*vsr)/1.d9*dpr2
    !d56 = 0.d0
    d66 = d66*dpr1 + (rhor*vsr*vsr)/1.d9*dpr2

  end if

  ! non-dimensionalize the elastic coefficients using
  ! the scale of GPa--[g/cm^3][(km/s)^2]

  d11 = d11/scale_GPa
  d12 = d12/scale_GPa
  d13 = d13/scale_GPa
  d14 = d14/scale_GPa
  d15 = d15/scale_GPa
  d16 = d16/scale_GPa
  d22 = d22/scale_GPa
  d23 = d23/scale_GPa
  d24 = d24/scale_GPa
  d25 = d25/scale_GPa
  d26 = d26/scale_GPa
  d33 = d33/scale_GPa
  d34 = d34/scale_GPa
  d35 = d35/scale_GPa
  d36 = d36/scale_GPa
  d44 = d44/scale_GPa
  d45 = d45/scale_GPa
  d46 = d46/scale_GPa
  d55 = d55/scale_GPa
  d56 = d56/scale_GPa
  d66 = d66/scale_GPa
  
  ! ---------------------------------------------------------------------------------
  ! non-dimensionalize
  rho = rho/EARTH_RHOAV ! rho*1000.d0/EARTH_RHOAV  rho [Kg/m3] EARTH_RHOAV 5514.3 [kg/m3]

  end subroutine build_cij_cij

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_aniso_mantle_model_cij()

  use constants, only: IIN,DEGREES_TO_RADIANS,ZERO 
  !use shared_parameters, only:RMOHO
  use model_aniso_mantle_cij_par

  implicit none

  ! local parameters
  integer :: i,k,l,param
  integer :: ier
  double precision :: PwaveMod,SwaveMod
  double precision, DIMENSION(21) :: XE 
  double precision, DIMENSION(22) :: AMM_V_Cijd
  character(len=*), parameter :: cij_model = 'DATA/CIJ/cij_model.xyz'
  
  !Read CIJ model
  open(IIN,file=cij_model,status='old',action='read',iostat=ier)

  if (ier /= 0 ) stop 'Error opening file cij_model.xyz'
 
  ! read the number of nodes in the three dimensions
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) nx, ny, nz

  !nz  Number of layers in the CIJ model
 

  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) AMM_V_lon
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) AMM_V_colat
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) AMM_V_pro

  ! REMEMBER to modify this part according to your needs
  !AMM_V_pro = (AMM_V_pro + (6371000.00 - RMOHO) - 500)/1000.00  ! convert meters to Km and add 24.4km of moho
  AMM_V_pro = AMM_V_pro / 1000.00  ! convert meters to Km and replace the crust 
  AMM_V_lon = (AMM_V_lon / DEGREES_TO_RADIANS) !+ 80.00  ! +40.00 convert radians to degree 
  AMM_V_colat = (AMM_V_colat / DEGREES_TO_RADIANS)  


  read(IIN, '(a)',end = 888)

  do l = 1,nz
    do i = 1,ny
      do k = 1,nx
           read(IIN, *,end = 888) (AMM_V_Cij(param,k,i,l), param=1,22)
           !Remove anisotropic component (however this can be done in
           !build_cij_cij)
           !XE(:) = AMM_V_Cij(2:22,k,i,l)
           !Mcur = 3D0/15D0*(XE(1)+XE(7)+XE(12))+2D0/15D0*(XE(2)+XE(3)+XE(8))+4D0/15D0*(XE(16)+XE(19)+XE(21))
           !Gcur = 1D0/15D0*(XE(1)+XE(7)+XE(12))-1D0/15D0*(XE(2)+XE(3)+XE(8))+1D0/5D0*(XE(16)+XE(19)+XE(21))
           !Lambda
           !Lcur = Mcur - 2*Gcur
           !AMM_V_Cij(2:22,k,i,l) = 0d0
           !AMM_V_Cij( 2,k,i,l) = Mcur !C11
           !AMM_V_Cij( 8,k,i,l) = Mcur !C22
           !AMM_V_Cij(13,k,i,l) = Mcur !C33
           !AMM_V_Cij( 3,k,i,l) = Lcur !C12
           !AMM_V_Cij( 4,k,i,l) = Lcur !C13
           !AMM_V_Cij( 9,k,i,l) = Lcur !C23
           !AMM_V_Cij(17,k,i,l) = Gcur !C44
           !AMM_V_Cij(20,k,i,l) = Gcur !C55
           !AMM_V_Cij(22,k,i,l) = Gcur !C66
       enddo
    enddo
  enddo
  ! Check boundaries
  do l = 1,nz
    do i = 1,ny
      do k = 1,nx
           if (AMM_V_Cij(1,k,i,l) == 0 .and. AMM_V_Cij(2,k,i,l) == 0 .and. AMM_V_Cij(3,k,i,l) == 0) then 
                  do param =1,22 
                     IF(k>1) AMM_V_Cij(param,k,i,l)=AMM_V_Cij(param,k-1,i,l)
                     IF(k==1 .AND. i> 1) AMM_V_Cij(param,k,i,l)=AMM_V_Cij(param,k,i-1,l)
                     IF(k==1 .AND. i==1 .AND. l > 1) AMM_V_Cij(param,k,i,l)=AMM_V_Cij(param,k,i,l-1)
                  enddo
           endif
           ! if an entire row of parameters is zero make it equal to the
           ! previous row
       enddo
    enddo
  enddo

888 close(IIN)

  !Find 1D model
  do l = 1,nz

    !Average of parameters for each layer
    AMM_V_Cijd = ZERO
    do i = 1,ny
      do k = 1,nx
         AMM_V_Cijd(:) = AMM_V_Cijd(:) + AMM_V_Cij(:,k,i,l)
      enddo
    enddo
    AMM_V_Cijd(:) = AMM_V_Cijd(:)/dble(nx*ny)

    !Compute isotropic velocities
    XE(:) = AMM_V_Cijd(2:22)
    PwaveMod = 3D0/15D0*(XE(1)+XE(7)+XE(12))+2D0/15D0*(XE(2)+XE(3)+XE(8))+4D0/15D0*(XE(16)+XE(19)+XE(21))
    SwaveMod = 1D0/15D0*(XE(1)+XE(7)+XE(12))-1D0/15D0*(XE(2)+XE(3)+XE(8))+1D0/5D0*(XE(16)+XE(19)+XE(21))
    AMM_V_Cije(1,l) = AMM_V_Cijd(1)
    AMM_V_Cije(2,l) = dsqrt(PwaveMod*1.d9/AMM_V_Cije(1,l))
    AMM_V_Cije(3,l) = dsqrt(SwaveMod*1.d9/AMM_V_Cije(1,l))

    !print *,AMM_V_pro(l),AMM_V_Cije(:,l)

  enddo

  end subroutine read_aniso_mantle_model_cij
