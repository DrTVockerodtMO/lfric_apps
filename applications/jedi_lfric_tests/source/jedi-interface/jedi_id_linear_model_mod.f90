!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class than handles the identity linear model
!>
!> @details This module includes a class that handles the identity linear model
!>          initialisation, time stepping and finalisation for both the Tangent
!>          Linear (TL) and Adjoint (AD). These are the required by the base
!>          interface class that uses them to define the forecastTL and
!>          forecastAD. The linear model forecasts (TL and AD) require a
!>          lineariastion-state (LS) trajectory. The LS trajectory is created
!>          by running the non-linear model and storing the result in an object
!>          of type: linear_state_trajectory_type. The set_trajectory method is
!>          included to provide the means to create and populate the LS fields.
!>
!>          An included test application (jedi_tlm_tests) performs the
!>          standard adjoint test to ensure the correctness of the code.
!>
module jedi_id_linear_model_mod

  use constants_mod,                 only : str_def
  use driver_time_mod,               only : init_time, final_time
  use field_collection_mod,          only : field_collection_type
  use driver_modeldb_mod,            only : modeldb_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_increment_mod,            only : jedi_increment_type
  use jedi_base_linear_model_mod,    only : jedi_base_linear_model_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_lfric_linear_fields_mod,  only : variable_names, &
                                            create_linear_fields
  use jedi_lfric_wind_fields_mod,    only : create_scalar_winds, &
                                            setup_vector_wind
  use jedi_state_mod,                only : jedi_state_type
  use linear_state_trajectory_mod,   only : linear_state_trajectory_type
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_ERROR
  use namelist_mod,                  only : namelist_type
  use transform_winds_mod,           only : wind_scalar_to_vector, wind_vector_to_scalar, &
                                            adj_wind_scalar_to_vector, adj_wind_vector_to_scalar
  use zero_field_collection_mod,     only : zero_field_collection

  implicit none

  private

type, public, extends(jedi_base_linear_model_type) :: jedi_id_linear_model_type
  private

  !> The model time step duration
  type ( jedi_duration_type )          :: time_step

  !> Trajectory of linear states obtained by running the non-linear model
  type( linear_state_trajectory_type ) :: linear_state_trajectory

  !> Modeldb that stores the model fields to propagate
  !> @todo: Required public for checksum but need to move to atlas checksum
  !>        so make it private when that work is done.
  type(modeldb_type), public           :: modeldb

contains

  !> Model initialiser.
  procedure, public  :: initialise

  !> Methods
  procedure, public  :: set_trajectory

  procedure, public :: model_initTL
  procedure, public :: model_stepTL
  procedure, public :: model_finalTL

  procedure, public :: model_initAD
  procedure, public :: model_stepAD
  procedure, public :: model_finalAD

  !> Finalizer
  final             :: jedi_linear_model_destructor

end type jedi_id_linear_model_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_id_linear_model_type
!>
!> @param [in] jedi_geometry A JEDI geometry object
!> @param [in] config        The linear model configuration
!> @param [in] config_filename The name of the configuration file
subroutine initialise( self, jedi_geometry, config_filename )

  use jedi_lfric_linear_modeldb_driver_mod, only : initialise_modeldb
  use jedi_lfric_timestep_mod,              only : get_configuration_timestep

  implicit none

  class( jedi_id_linear_model_type ), target, intent(inout) :: self
  type( jedi_geometry_type ),                    intent(in) :: jedi_geometry
  character(len=*),                              intent(in) :: config_filename

  ! Local
  type( namelist_type ), pointer :: jedi_lfric_settings_config
  type( namelist_type ), pointer :: jedi_linear_model_config
  character( str_def )           :: forecast_length_str
  character( str_def )           :: nl_time_step_str
  type( jedi_duration_type )     :: forecast_length
  type( jedi_duration_type )     :: nl_time_step

  ! 1. Setup modeldb

  ! 1.1 Initialise the modeldb
  call initialise_modeldb( "linear modeldb", config_filename, &
                            jedi_geometry%get_mpi_comm(), self%modeldb )

  ! 1.2 Add scalar winds that link the Atlas fields. These are used to
  ! perform interpolation between horizontally cell-centred and
  ! edge based W2 winds.
  call create_scalar_winds( self%modeldb, jedi_geometry%get_mesh() )

  ! 2. Setup time
  self%time_step = get_configuration_timestep( self%modeldb%configuration )

  jedi_lfric_settings_config => self%modeldb%configuration%get_namelist('jedi_lfric_settings')
  call jedi_lfric_settings_config%get_value( 'forecast_length', forecast_length_str )
  call forecast_length%init(forecast_length_str)

  ! 3. Setup trajactory
  jedi_linear_model_config => self%modeldb%configuration%get_namelist('jedi_linear_model')
  call jedi_linear_model_config%get_value( 'nl_time_step', nl_time_step_str )
  call nl_time_step%init(nl_time_step_str)

  call self%linear_state_trajectory%initialise( forecast_length, &
                                                nl_time_step )

end subroutine initialise

!> @brief    Set an instance of the trajectory
!>
!> @param [in] jedi_state The state to add to the trajectory
subroutine set_trajectory( self, jedi_state )

  implicit none

  class( jedi_id_linear_model_type ),   intent(inout) :: self
  type( jedi_state_type ),              intent(inout) :: jedi_state

  ! Local
  type( field_collection_type ) :: next_linear_state

  ! Create field collection that contains the linear state fields
  ! without "ls_" prepended.
  call create_linear_fields(jedi_state%geometry%get_mesh(), next_linear_state)

  ! Copy data from the input state into next_linear_state
  call jedi_state%get_to_field_collection( variable_names, &
                                           next_linear_state )

  ! Create W2 wind, interpolate from scaler winds (W3/Wtheta) then
  ! remove the scaler winds
  call setup_vector_wind(next_linear_state)

  ! Add the new linear state to the trajectory (prepends fieldnames with "ls_")
  call self%linear_state_trajectory%add_linear_state( &
                                              jedi_state%valid_time(), &
                                              next_linear_state )

end subroutine set_trajectory

!> @brief    Initialise the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model initialise
subroutine model_initTL(self, increment)

  implicit none

  class( jedi_id_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),        intent(inout) :: increment

  ! Local
  type( field_collection_type),  pointer :: prognostic_fields

  nullify(prognostic_fields)

  ! 1. Update the the modeldb prognostic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")

  ! Get Atlas field emulators to the model_prognostics
  call increment%get_to_field_collection( variable_names, prognostic_fields )

  ! Cell-centred winds to Edge based winds
  call wind_scalar_to_vector( prognostic_fields )

  ! Initialise clock and calendar
  call init_time( self%modeldb )

end subroutine model_initTL

!> @brief    Step the Tangent Linear model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepTL(self, increment)

  use jedi_lfric_linear_modeldb_driver_mod, only : identity_step_tl

  implicit none

  class( jedi_id_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),                intent(inout) :: increment

  ! Local
  type( field_collection_type ), pointer :: depository
  type( field_collection_type),  pointer :: prognostic_fields

  nullify(depository, prognostic_fields)

  ! 1. Update the LFRic modeldb linear state fields
  depository => self%modeldb%fields%get_field_collection("depository")
  call self%linear_state_trajectory%get_linear_state( increment%valid_time(), &
                                                      depository )

  ! 2. Step the linear model
  call identity_step_tl( self%modeldb )

  ! 3. Update the Atlas fields from the LFRic prognostic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")

  ! Interpolate W2 vector to W3/Wtheta scalar winds
  call wind_vector_to_scalar( prognostic_fields )

  ! Set model_prognostics to the Atlas field emulators
  call increment%set_from_field_collection( variable_names, prognostic_fields )

  ! 4. Update the increment time
  call increment%update_time( self%time_step )

end subroutine model_stepTL

!> @brief    Finalise the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalTL(self, increment)

  implicit none

  class( jedi_id_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),        intent(inout) :: increment

  ! Finalise clock and calendar
  call final_time( self%modeldb )

end subroutine model_finalTL

!> @brief    Initialise the Adjoint of the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model initialise
subroutine model_initAD(self, increment)

  implicit none

  class( jedi_id_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),        intent(inout) :: increment

  ! Local
  type( field_collection_type), pointer :: prognostic_fields

  nullify(prognostic_fields)

  ! Update the prognostic fields: zero LFRic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")
  call zero_field_collection(prognostic_fields)

  !>@todo: in ticket #267
  !>       Add call similar to init_time to setup reversable clock

end subroutine model_initAD

!> @brief    Step the Adjoint of the Tangent Linear model
!>
!> @param [inout] increment Increment object to be propagated
subroutine model_stepAD(self, increment)

  implicit none

  class( jedi_id_linear_model_type ), target, intent(inout) :: self
  type( jedi_increment_type ),                intent(inout) :: increment

  ! Local
  type( field_collection_type ), pointer :: depository
  type( field_collection_type),  pointer :: prognostic_fields

  ! 1. Update the LFRic modeldb linear state fields

  ! 1.1 Copy from the trajectory into the model_data
  depository => self%modeldb%fields%get_field_collection("depository")
  call self%linear_state_trajectory%get_linear_state( increment%valid_time(), &
                                                      depository )

  ! 2. Step the linear model
  !>@todo: in ticket #267
  !> add identity_step_ad when reversable clock is available

  ! 3. Update the Atlas fields from the LFRic prognostic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")

  ! Set model_prognostics to the Atlas field emulators
  call increment%set_from_field_collection_ad( variable_names, prognostic_fields )

  ! Adjoint of ... Interpolate W2 vector to W3/Wtheta scalar winds
  call adj_wind_vector_to_scalar( prognostic_fields )

  ! 4. Update the increment time
  call increment%update_time( self%time_step*(-1) )

end subroutine model_stepAD

!> @brief    Finalise the Adjoint of the Tangent Linear model
!>
!> @param [inout] increment Increment object to be used in the model finalise
subroutine model_finalAD(self, increment)

  implicit none

  class( jedi_id_linear_model_type ), intent(inout) :: self
  type( jedi_increment_type ),        intent(inout) :: increment

  ! Local
  type( field_collection_type),  pointer :: prognostic_fields

  nullify(prognostic_fields)

  ! 1. Update the the modeldb prognostic fields
  prognostic_fields => self%modeldb%fields%get_field_collection("prognostic_fields")

  ! Cell-centred winds to Edge based winds
  call adj_wind_scalar_to_vector( prognostic_fields )

  ! Get Atlas field emulators to the model_prognostics
  call increment%get_to_field_collection_ad( variable_names, prognostic_fields )

  !>@todo: in ticket #267
  !> add call similar to final_time to close reversable clock

end subroutine model_finalAD

!> @brief    Finalize the jedi_id_linear_model_type
!>
subroutine jedi_linear_model_destructor(self)

  use jedi_lfric_linear_modeldb_driver_mod, only : finalise_modeldb

  implicit none

  type(jedi_id_linear_model_type), intent(inout) :: self

  call finalise_modeldb( self%modeldb )

end subroutine jedi_linear_model_destructor

end module jedi_id_linear_model_mod
