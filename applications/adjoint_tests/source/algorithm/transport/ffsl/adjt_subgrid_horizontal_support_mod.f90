!------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief Module containing adjoint tests for adj_subgrid_horizontal_support_mod routines
module adjt_subgrid_horizontal_support_mod

  use constants_mod,                  only : i_def, r_tran, EPS
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG

  implicit none

  public :: adjt_fourth_order_horizontal_edge
  public :: adjt_nirvana_horizontal_edge

  contains
  !=============================================================================
  !> @brief   Adjoint test for adj_fourth_order_horizontal_edge.
  !> @details Passes if adjoint is transpose of tangent linear.
  !>          Determined by testing the equality of inner products <Mx, Mx> and <AMx, x>,
  !>          where M is the tangent linear and A is the adjoint.
  !> @param[in] nlayers     Number of layers
  subroutine adjt_fourth_order_horizontal_edge(nlayers)

    use subgrid_horizontal_support_mod,     only : fourth_order_horizontal_edge
    use adj_subgrid_horizontal_support_mod, only : adj_fourth_order_horizontal_edge

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    ! Arguments for tl and adj calls
    real(kind=r_tran), dimension(:,:), allocatable :: field
    real(kind=r_tran), dimension(:),   allocatable :: edge_left
    real(kind=r_tran), dimension(:),   allocatable :: edge_right

    ! Copies of input fields used in inner products
    real(kind=r_tran), dimension(:,:), allocatable :: fld_inp
    real(kind=r_tran), dimension(:),   allocatable :: edge_l_inp
    real(kind=r_tran), dimension(:),   allocatable :: edge_r_inp

    ! Variables for allocating arrays
    integer(kind=i_def) :: order

    ! Inner products
    integer(kind=i_def) :: i_order
    real(kind=r_tran)   :: field_inner_prod
    real(kind=r_tran)   :: edge_l_inner_prod
    real(kind=r_tran)   :: edge_r_inner_prod
    real(kind=r_tran)   :: field_sf
    real(kind=r_tran)   :: edge_l_sf
    real(kind=r_tran)   :: edge_r_sf
    real(kind=r_tran)   :: inner1
    real(kind=r_tran)   :: fld_fld_inp_inner_prod
    real(kind=r_tran)   :: edge_l_edge_l_inp_inner_prod
    real(kind=r_tran)   :: edge_r_edge_r_inp_inner_prod
    real(kind=r_tran)   :: inner2

    ! Test parameters and variables
    real(kind=r_tran), parameter :: overall_tolerance = 1500.0_r_tran
    real(kind=r_tran)            :: machine_tol
    real(kind=r_tran)            :: relative_diff

    order = 2

    allocate(field(1:nlayers, 1:(2*order + 1)),   &
             edge_left(1:nlayers),                &
             edge_right(1:nlayers),               &
             fld_inp(1:nlayers, 1:(2*order + 1)), &
             edge_l_inp(1:nlayers),               &
             edge_r_inp(1:nlayers))

    ! Initialise arguments and call the forward routine.
    call RANDOM_NUMBER(field)
    fld_inp = field
    call RANDOM_NUMBER(edge_left)
    edge_l_inp = edge_left
    call RANDOM_NUMBER(edge_right)
    edge_r_inp = edge_right

    call fourth_order_horizontal_edge(field,      &
                                      edge_left,  &
                                      edge_right, &
                                      nlayers)

    field_inner_prod = 0.0_r_tran
    do i_order = 1, 2*order + 1
      field_inner_prod = field_inner_prod &
                       + dot_product(field(:, i_order), field(:, i_order))
    end do
    edge_l_inner_prod = dot_product(edge_left, edge_left)
    edge_r_inner_prod = dot_product(edge_right, edge_right)

    write( log_scratch_space, * ) "adjt_fourth_order_horizontal_edge inner products:"
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "field inner product = ", field_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_left inner product = ", edge_l_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_right inner product = ", edge_r_inner_prod
    call log_event( log_scratch_space, log_level_debug )

    field_sf = 1.0_r_tran / real(field_inner_prod + EPS, kind = r_tran)
    edge_l_sf = 1.0_r_tran / real(edge_l_inner_prod + EPS, kind = r_tran)
    edge_r_sf = 1.0_r_tran / real(edge_r_inner_prod + EPS, kind = r_tran)

    inner1 = 0.0_r_tran
    inner1 = inner1 + field_inner_prod * field_sf
    inner1 = inner1 + edge_l_inner_prod * edge_l_sf
    inner1 = inner1 + edge_r_inner_prod * edge_r_sf

    field(:,:) = field(:,:) * field_sf
    edge_left(:) = edge_left(:) * edge_l_sf
    edge_right(:) = edge_right(:) * edge_r_sf

    call adj_fourth_order_horizontal_edge(field,      &
                                          edge_left,  &
                                          edge_right, &
                                          nlayers)

    fld_fld_inp_inner_prod = 0.0_r_tran
    do i_order = 1, 2*order + 1
      fld_fld_inp_inner_prod = fld_fld_inp_inner_prod &
                             + dot_product(field(:, i_order), fld_inp(:, i_order))
    end do
    edge_l_edge_l_inp_inner_prod = dot_product(edge_left, edge_l_inp)
    edge_r_edge_r_inp_inner_prod = dot_product(edge_right, edge_r_inp)

    inner2 = 0.0_r_tran
    inner2 = inner2 + fld_fld_inp_inner_prod
    inner2 = inner2 + edge_l_edge_l_inp_inner_prod
    inner2 = inner2 + edge_r_edge_r_inp_inner_prod

    deallocate(field,          &
               edge_left,      &
               edge_right,     &
               fld_inp,        &
               edge_l_inp,     &
               edge_r_inp)

    ! Test the inner-product values for equality, allowing for the precision of the active variables
    machine_tol = spacing(max(abs(inner1), abs(inner2)))
    relative_diff = abs(inner1 - inner2) / machine_tol
    if (relative_diff < overall_tolerance) then
      write(log_scratch_space, *) "PASSED fourth_order_horizontal_edge:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
    else
      write(log_scratch_space, *) "FAILED fourth_order_horizontal_edge:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine adjt_fourth_order_horizontal_edge

  !=============================================================================
  !> @brief   Adjoint test for adj_nirvana_horizontal_edge.
  !> @details Passes if adjoint is transpose of tangent linear.
  !>          Determined by testing the equality of inner products <Mx, Mx> and <AMx, x>,
  !>          where M is the tangent linear and A is the adjoint.
  !> @param[in] nlayers     Number of layers
  subroutine adjt_nirvana_horizontal_edge(nlayers)

    use subgrid_horizontal_support_mod,     only : nirvana_horizontal_edge
    use adj_subgrid_horizontal_support_mod, only : adj_nirvana_horizontal_edge

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    ! Arguments for tl and adj calls
    real(kind=r_tran), dimension(:,:), allocatable :: field
    real(kind=r_tran), dimension(:),   allocatable :: edge_left
    real(kind=r_tran), dimension(:),   allocatable :: edge_right

    ! Copies of input fields used in inner products
    real(kind=r_tran), dimension(:,:), allocatable :: fld_inp
    real(kind=r_tran), dimension(:),   allocatable :: edge_l_inp
    real(kind=r_tran), dimension(:),   allocatable :: edge_r_inp

    ! Variables for allocating arrays
    integer(kind=i_def) :: order

    ! Inner products
    integer(kind=i_def) :: i_order
    real(kind=r_tran)   :: field_inner_prod
    real(kind=r_tran)   :: edge_l_inner_prod
    real(kind=r_tran)   :: edge_r_inner_prod
    real(kind=r_tran)   :: field_sf
    real(kind=r_tran)   :: edge_l_sf
    real(kind=r_tran)   :: edge_r_sf
    real(kind=r_tran)   :: inner1
    real(kind=r_tran)   :: fld_fld_inp_inner_prod
    real(kind=r_tran)   :: edge_l_edge_l_inp_inner_prod
    real(kind=r_tran)   :: edge_r_edge_r_inp_inner_prod
    real(kind=r_tran)   :: inner2

    ! Test parameters and variables
    real(kind=r_tran), parameter :: overall_tolerance = 1500.0_r_tran
    real(kind=r_tran)            :: machine_tol
    real(kind=r_tran)            :: relative_diff

    order = 1

    allocate(field(1:nlayers, 1:(2*order + 1)),   &
             edge_left(1:nlayers),                &
             edge_right(1:nlayers),               &
             fld_inp(1:nlayers, 1:(2*order + 1)), &
             edge_l_inp(1:nlayers),               &
             edge_r_inp(1:nlayers))

    ! Initialise arguments and call the forward routine.
    call RANDOM_NUMBER(field)
    fld_inp = field
    call RANDOM_NUMBER(edge_left)
    edge_l_inp = edge_left
    call RANDOM_NUMBER(edge_right)
    edge_r_inp = edge_right

    call nirvana_horizontal_edge(field,      &
                                 edge_left,  &
                                 edge_right, &
                                 nlayers)

    field_inner_prod = 0.0_r_tran
    do i_order = 1, 2*order + 1
      field_inner_prod = field_inner_prod &
                       + dot_product(field(:, i_order), field(:, i_order))
    end do
    edge_l_inner_prod = dot_product(edge_left, edge_left)
    edge_r_inner_prod = dot_product(edge_right, edge_right)

    write( log_scratch_space, * ) "adjt_nirvana_horizontal_edge inner products:"
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "field inner product = ", field_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_left inner product = ", edge_l_inner_prod
    call log_event( log_scratch_space, log_level_debug )
    write( log_scratch_space, * ) "edge_right inner product = ", edge_r_inner_prod
    call log_event( log_scratch_space, log_level_debug )

    field_sf = 1.0_r_tran / real(field_inner_prod + EPS, kind = r_tran)
    edge_l_sf = 1.0_r_tran / real(edge_l_inner_prod + EPS, kind = r_tran)
    edge_r_sf = 1.0_r_tran / real(edge_r_inner_prod + EPS, kind = r_tran)

    inner1 = 0.0_r_tran
    inner1 = inner1 + field_inner_prod * field_sf
    inner1 = inner1 + edge_l_inner_prod * edge_l_sf
    inner1 = inner1 + edge_r_inner_prod * edge_r_sf

    field(:,:) = field(:,:) * field_sf
    edge_left(:) = edge_left(:) * edge_l_sf
    edge_right(:) = edge_right(:) * edge_r_sf

    call adj_nirvana_horizontal_edge(field,      &
                                     edge_left,  &
                                     edge_right, &
                                     nlayers)

    fld_fld_inp_inner_prod = 0.0_r_tran
    do i_order = 1, 2*order + 1
      fld_fld_inp_inner_prod = fld_fld_inp_inner_prod &
                             + dot_product(field(:, i_order), fld_inp(:, i_order))
    end do
    edge_l_edge_l_inp_inner_prod = dot_product(edge_left, edge_l_inp)
    edge_r_edge_r_inp_inner_prod = dot_product(edge_right, edge_r_inp)

    inner2 = 0.0_r_tran
    inner2 = inner2 + fld_fld_inp_inner_prod
    inner2 = inner2 + edge_l_edge_l_inp_inner_prod
    inner2 = inner2 + edge_r_edge_r_inp_inner_prod

    deallocate(field,          &
               edge_left,      &
               edge_right,     &
               fld_inp,        &
               edge_l_inp,     &
               edge_r_inp)

    ! Test the inner-product values for equality, allowing for the precision of the active variables
    machine_tol = spacing(max(abs(inner1), abs(inner2)))
    relative_diff = abs(inner1 - inner2) / machine_tol
    if (relative_diff < overall_tolerance) then
      write(log_scratch_space, *) "PASSED nirvana_horizontal_edge:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
    else
      write(log_scratch_space, *) "FAILED nirvana_horizontal_edge:", inner1, inner2, relative_diff
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine adjt_nirvana_horizontal_edge

end module adjt_subgrid_horizontal_support_mod
