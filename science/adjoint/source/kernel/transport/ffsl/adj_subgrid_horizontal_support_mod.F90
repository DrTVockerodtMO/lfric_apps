!------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing adjoint horizontal Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module adj_subgrid_horizontal_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN

implicit none

private

! Edge interpolation routines
public :: adj_fourth_order_horizontal_edge
public :: adj_nirvana_horizontal_edge

contains

  !----------------------------------------------------------------------------
  !> @brief  Calculates the adjoint of the interpolated field at the edges of a cell
  !!         required for using horizontal PPM to estimate the quadratic subgrid
  !!         representation of a field.
  !!
  !> @param[in,out] field         Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !> @param[in,out] edge_left     Field value at left edge of cell
  !> @param[in,out] edge_right    Field value at right edge of cell
  !> @param[in]     nlayers       Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine adj_fourth_order_horizontal_edge(field, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(inout) :: field(nlayers,5)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers), edge_right(nlayers)

    real(kind=r_tran), parameter :: twelfth = 1.0_r_tran/12.0_r_tran

    field(:,1) = field(:,1) &
               - twelfth * edge_left(:)
    field(:,2) = field(:,2) &
               + twelfth * (7.0_r_tran * edge_left(:) - edge_right(:))
    field(:,3) = field(:,3) &
               + 7.0_r_tran * twelfth * (edge_left(:) + edge_right(:))
    field(:,4) = field(:,4) &
               + twelfth * (7.0_r_tran * edge_right(:) - edge_left(:))
    field(:,5) = field(:,5) &
               - twelfth * edge_right(:)
    edge_left(:) = 0.0_r_tran
    edge_right(:) = 0.0_r_tran

  end subroutine adj_fourth_order_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the adjoint of interpolated field at the edges of a cell
  !!         required for using horizontal Nirvana to estimate the quadratic subgrid
  !!         representation of a field.
  !!
  !> @param[in,out] field            Has dof map of the form | 1 | 2 | 3 |
  !> @param[in,out] edge_left        Field value at left edge of cell
  !> @param[in,out] edge_right       Field value at right edge of cell
  !> @param[in]     nlayers          Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine adj_nirvana_horizontal_edge(field, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(inout) :: field(nlayers,3)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers), edge_right(nlayers)

    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    field(:,1) = field(:,1) &
               + sixth * (2.0_r_tran * edge_left(:) - edge_right(:))
    field(:,2) = field(:,2) &
               + 5.0_r_tran * sixth * (edge_left(:) + edge_right(:))
    field(:,3) = field(:,3) &
               + sixth * (2.0_r_tran * edge_right(:) - edge_left(:))
    edge_left(:) = 0.0_r_tran
    edge_right(:) = 0.0_r_tran

  end subroutine adj_nirvana_horizontal_edge

end module adj_subgrid_horizontal_support_mod
