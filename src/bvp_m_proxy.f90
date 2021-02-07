! # BVP_M_proxy
! (This BVP_M_proxy is part of the ODEInterface.)
!
! ## Why this BVP_M_proxy/wrapper?
!
! BVP_M-2 is written in a modern Fortran (2003) language. Fortran90 
! and above has a (completely) different conecpt of a "pointer" (or a 
! reference) than C:
!
! * In C a pointer is (more or less) a address in memory. If you pass a
!   C-pointer to a C-function, then only this address is passed.
!
! * In Fortran90 a pointer has a lot more information associated to it: e.g.
!   the type of the target, the rank (=dimension) of the target, 
!   the size, etc. If you pass a Fortran-pointer to a Fortran-function, 
!   then all this information is passed, too.
!
! At the time of writing this, there exists (e.g.) no compiler-independent
! way to create in C a Fortran90 array descriptor (a.k.a. dope vector) with
! all this information (rank, size, etc.) given above.
!
! That are the reasons why this proxy comes into play.
!
! ## What are the tasks of this proxy?
!
! 1. This proxy uses the `ISO_C_BINDING` Fortran module to export 
!    C-callable Fortran-routines which in turn transform the given input 
!    data to Fortran objects. This action is also called "wrapping".
!
! ```
! ┌──────────────────────────────┐      ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
! │Julia                         │      ┃Fortran: wrapping            ┃
! │ x = [1.0, 2.0]               │─────▶┃input: n=length(x) and ptr(x)┃
! │ proxy_func(length(x), ptr(x))│      ┃constructs Fortran Array ARR ┃
! └──────────────────────────────┘      ┃calls BVP_M-2_func(ARR)      ┃
!                                       ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
!                                           │
!                                           │
!                                           ▼
!                                       ╔══════════╗
!                                       ║BVP_M-2   ║
!                                       ║input: ARR║
!                                       ║uses ARR  ║
!                                       ╚══════════╝
! ```
!
! 2. Because BVP_M-2 needs (as input) some (boundary-value problem dependent)
!    Fortran functions/subroutines (e.g. the right-hand side of the ODEs) this
!    proxy also tries to implement such Fortran subroutines that call
!    (C-)functions supplied by julia. This can be used to even allow
!    for julia-callback functions that are closures.
!    Let's see this principle in action for a concrete example.
!    The BVP_M-2 solver has a method called `bvp_init` which can be called
!    with a Fortran-Function which is used as callback to get the
!    initial guesses at different positions with the help of this
!    Fortran-Function.
!
! ```
!   ┌─bvpm2_guess (julia) ───────────────────────────────────────────┐
!   │ inputs: ... guess_fcn ... (julia function to call)             │
!   │                                                                │
!   │ create callback-info object (cbi::Bvpm2_guess_cbi)             │
!   │   where guess_fcn is saved                                     │ 1.
!   │ call init_guess3_c( ... guess_fcn_ptr=unsafe_bvpm2_guess_cb_c, │────┐
!   │                         guess_pthrough=cbi)                    │    │
!   └────────────────────────────────────────────────────────────────┘    │
!                                                                         │
!                                                                         │
!   ┏━init_guess3_c (Fortran2003-Proxy with ISO_C_BINDING) ━━━━━━━━━━┓<───┘
!   ┃ inputs: ... guess_fcn_ptr == unsafe_bvpm2_guess_cb_c           ┃
!   ┃             guess_pthrough == cbi                              ┃
!   ┃                                                                ┃ 2.
!   ┃ call bvp_init( ... guess = guess_fcn_proxy ... )               ┃─────┐
!   ┃                                                                ┃ 3.  │
!   ┃ ┏━guess_fcn_proxy (Fortran2003-Proxy, nested!) ━━━━━━━━━━━━━━┓<╂───┐ │
!   ┃ ┃ inputs: x_point, guess_vector (needs to be calculated)     ┃ ┃   │ │
!   ┃ ┃ available because of nested Fortran func:                  ┃ ┃   │ │
!   ┃ ┃                  guess_pthrough, guess_fcn_ptr             ┃ ┃   │ │
!   ┃ ┃                                                            ┃ ┃   │ │
!   ┃ ┃ convert guess_vector to C-Pointer-Array                    ┃ ┃   │ │
!   ┃ ┃ call guess_fcn_ptr(x_point, guess_vector_c_compat,         ┃ ┃─┐ │ │
!   ┃ ┃                    guess_pthrough)                         ┃ ┃ │ │ │
!   ┃ ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛ ┃ │ │ │
!   ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛ │ │ │
!                                                                      │ │ │
!                                                                      │ │ │
!   ┌─unsafe_bvpm2_guess_cb_c (julia, generated with cfunction) ─────┐<┘ │ │
!   │ inputs: x, guess_vector_len, guess_vector_c_compat, cbi        │ 4.│ │
!   │                                                                │   │ │
!   │ use guess_vector_* to create Julia wrapper (unsafe_wrap)       │   │ │
!   │ call cbi.guess_fcn(y, y)                                       │   │ │
!   └────────────────────────────────────────────────────────────────┘   │ │
!                                                                        │ │
!                                                                        │ │
!   ╔═bvp_init =BVP_M-2.GUESS_3 (Fotran2003) ════════════════════════╗<────┘
!   ║ inputs: ... FCN==guess_fcn_proxy ...                           ║   │
!   ║                                                                ║   │
!   ║ [does a lot of things]                                         ║   │
!   ║ call FCN(X, Y)                                                 ║───┘
!   ╚════════════════════════════════════════════════════════════════╝
!
!    Legend:
!       ─────────  julia code
!       ━━━━━━━━━  Fortran2003 code in BVP_M_Proxy with ISO_C_BINDING
!       ═════════  Fortran2003 code in BVP_M-2 (without ISO_C_BINDING)
! ```
!
! 3. Because BVP_M-2 uses Fortran objects (for a user-friendly and modern
!    approach) this proxy tries to implement "handles" (see below)
!    to make it possible to have "proxy"-objects/types in Julia in such a
!    way that the objet-oriented feeling is preserved and the Fortran
!    calls are transparent/invisible to the user.
!
! ## What does "handle" mean in this context?
!
! To support the/a object-oriented look-and-feel in Julia for the
! Fortran90 objects of BVP_M-2 the following idea is used:
!
! Julia (or C) is given (via `C_LOC`) the "C-address" of an Fortran-object 
! as handle. If Julia (or C) calls some method of this proxy with
! such a handle then (via `C_F_POINTER`) the original Fortran-object
! is reobtained. This can only work, if (at the time of the call of
! `C_F_POINTER`) all additional informations (type, rank, shape, etc.) is
! known or can be determined by the "context".
! One task of this proxy (and the Julia part) is to use data structures
! in such a way, that there exists always enough data/context
! to infer the needed information:
!
!      Julia                        │      Proxy
!                                   │
!       call proxy_create_obj() ────┼─▶
!                                   │  obj=bvp_m-2_create_obj()
!                                   │  create handle from obj
!       handle                 ◀────┼─ handle
!                                   │
!       [some time later]           │
!                                   │
!       call proxy_method(handle) ──┼─▶handle
!                                   │  reobtain obj [use context]
!                                   │  result = method(obj)
!       result                 ◀────┼─ result
!
! Of course: On the Julia side, all this handles are never exposed to the
! user, but packed inside mutable structs/types and multiple-dispatch is
! used to get the object-oriented look-and-feel and to hide the
! calls to this (Fortran-)Proxy and to Fortran.
!
!
! State transition diagram for BVP_M-2 sol object, according to the 
! documentation. BVP_M-2 does *NOT* track the state. Hence this proxy takes 
! over this task.
!
!                   ┌──────────────────────────────────────┐
!                   │           STATE=0                    │
!                   │   No fields of SOL are allocated.    │
!                   └──────────────────────────────────────┘
!                    │         ▲                  │    ▲
!                    │INIT     ║                  │    │
!                    │         ║SOLVER            │    │
!                    │         ║(unsuccessful)    │    │
!                    ▼         ║                  │    │
!           ┌─────────────────────────────┐       │    │
!           │         STATE=1             │       │    │
!           │X, Y,(optionally) PARAMETERS.│       │    │
!           └─────────────────────────────┘    GET│    │TERMINATE
!                           ║   ▲                 │    │
!                     SOLVER║   ║                 │    │
!               (successful)║   ║EXTEND           │    │
!                           ║   ║                 │    │
!      ┌─────────┐          ║   ║                 │    │   ┌───────────┐
!      │   SAVE  │          ║   ║                 │    │   │   EVAL    │
!      │         │          ▼   ║                 ▼    │   │           │
!      │     ┌───────────────────────────────────────────────┐         │
!      └────▶│              STATE=2                          │◀────────┘
!            │X, Y,(optionally) PARAMETERS, IWORK, WORK.     │
!            └───────────────────────────────────────────────┘
!
!       Legend:
!     sol    ──────▶ sol      subroutine transforms BVP_M-2 sol object
!     sol_in ══════▶ sol_out  subroutine takes and *destroys* input 
!                             sol_in object and creates new sol_out object 
!                             with the shown state
!
! The Fortan routine `BVP_SOLVER` takes as input a BVP_M-2 sol object
! in STATE=1 (initial guess) and creates a *new* output object, with
!
! ```
!                                      ⎧ ┌────────────────┐
!   ┌─────────────────┐    SOLVER      ⎪ │new SOL, STATE=0│ if unsuccessful
!   │SOL_INIT, STATE=1│ ════════════▶  ⎨ └────────────────┘
!   └─────────────────┘                ⎪ ┌────────────────┐
!                                      ⎪ │new SOL, STATE=2│ if successful
!                                      ⎩ └────────────────┘
! ```
!
! The Fortran routine `BVP_EXTEND` takes as input a BVP_M2 sol object
! in STATE=2 (solution) and creates a *new* guess-object on a different
! interval:
! 
! ```
!   ┌─────────────────┐    EXTEND     ┌────────────────┐
!   │SOL_INIT, STATE=2│ ════════════▶ │new SOL, STATE=1│
!   └─────────────────┘               └────────────────┘
! ```
!

! The line above intentionally left blank because 
! the text above the blank line is interpreted as Markdown!

module bvp_m_proxy

  use, intrinsic :: iso_c_binding

  use bvp_m

  implicit none

  private

  public :: show_sol_wrapper_c, get_sol_wrapper_details_c, &
            create_sol_wrapper_c, terminate_sol_wrapper_c, &
            destroy_sol_wrapper_c, copy_sol_wrapper_c, &
            get_sol_wrapper_x_c, &
            init_guess1_c, init_guess2_c, init_guess3_c, &
            solve_c, extend_sol_s_c, extend_sol_e_c

  abstract interface
    subroutine guess_callback_c(x, y_len, y, guess_pthrough) bind(C)
      use, intrinsic :: iso_c_binding
      ! In
      real(c_double),     value, intent(in)    :: x
      integer(c_int64_t), value, intent(in)    :: y_len
      type(c_ptr),        value, intent(in)    :: guess_pthrough
      ! Out
      real(c_double), dimension(y_len), intent(out) :: y
    end subroutine guess_callback_c
  end interface

  ! RHS Callbacks    (f)
    abstract interface
      subroutine f_callback_c(x, y_len, y, f_len, f, calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        real(c_double),     value,        intent(in)  :: x
        integer(c_int64_t), value,        intent(in)  :: y_len
        real(c_double), dimension(y_len), intent(in)  :: y
        integer(c_int64_t), value,        intent(in)  :: f_len
        type(c_ptr),        value,        intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(f_len), intent(in)  :: f
      end subroutine f_callback_c
    end interface

    abstract interface
      subroutine fpar_callback_c(x, y_len, y, p_len, p, &
                            f_len, f, calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        real(c_double),     value,        intent(in)  :: x
        integer(c_int64_t), value,        intent(in)  :: y_len
        real(c_double), dimension(y_len), intent(in)  :: y
        integer(c_int64_t), value,        intent(in)  :: p_len
        real(c_double), dimension(p_len), intent(in)  :: p
        integer(c_int64_t), value,        intent(in)  :: f_len
        type(c_ptr),        value,        intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(f_len), intent(in)  :: f
      end subroutine fpar_callback_c
    end interface

  ! DRHS Callbacks   (Df)
    abstract interface
      subroutine Df_callback_c(x, y_len, y, dfdy_dim1, &
                           dfdy_dim2, dfdy, calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        real(c_double),     value,        intent(in)  :: x
        integer(c_int64_t), value,        intent(in)  :: y_len
        real(c_double), dimension(y_len), intent(in)  :: y
        integer(c_int64_t), value,        intent(in)  :: dfdy_dim1, dfdy_dim2
        type(c_ptr),        value,        intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(dfdy_dim1, dfdy_dim2), &
                                          intent(in)  :: dfdy
      end subroutine Df_callback_c
    end interface

    abstract interface
      subroutine Dfpar_callback_c(x, y_len, y, p_len, p, &
                     dfdy_dim1, dfdy_dim2, dfdy, &
                     dfdp_dim1, dfdp_dim2, dfdp, calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        real(c_double),     value,        intent(in)  :: x
        integer(c_int64_t), value,        intent(in)  :: y_len
        real(c_double), dimension(y_len), intent(in)  :: y
        integer(c_int64_t), value,        intent(in)  :: p_len
        real(c_double), dimension(p_len), intent(in)  :: p
        integer(c_int64_t), value,        intent(in)  :: dfdy_dim1, dfdy_dim2
        integer(c_int64_t), value,        intent(in)  :: dfdp_dim1, dfdp_dim2
        type(c_ptr),        value,        intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(dfdy_dim1, dfdy_dim2), &
                                          intent(in)  :: dfdy
        real(c_double), dimension(dfdp_dim1, dfdp_dim2), &
                                          intent(in)  :: dfdp
      end subroutine Dfpar_callback_c
    end interface

  ! BC Callbacks     (bc)
    abstract interface
      subroutine bc_callback_c(ya_len, ya, yb_len, yb, &
                               bca_len, bca, bcb_len, bcb, &
                               calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        integer(c_int64_t), value,         intent(in)  :: ya_len
        real(c_double), dimension(ya_len), intent(in)  :: ya
        integer(c_int64_t), value,         intent(in)  :: yb_len
        real(c_double), dimension(yb_len), intent(in)  :: yb
        integer(c_int64_t), value,         intent(in)  :: bca_len
        integer(c_int64_t), value,         intent(in)  :: bcb_len
        type(c_ptr),        value,         intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(bca_len),intent(in)  :: bca
        real(c_double), dimension(bcb_len),intent(in)  :: bcb
      end subroutine bc_callback_c
    end interface

    abstract interface
      subroutine bcpar_callback_c(ya_len, ya, yb_len, yb,  &
                               p_len, p,                   &
                               bca_len, bca, bcb_len, bcb, &
                               calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        integer(c_int64_t), value,         intent(in)  :: ya_len
        real(c_double), dimension(ya_len), intent(in)  :: ya
        integer(c_int64_t), value,         intent(in)  :: yb_len
        real(c_double), dimension(yb_len), intent(in)  :: yb
        integer(c_int64_t), value,         intent(in)  :: p_len
        real(c_double), dimension(p_len),  intent(in)  :: p
        integer(c_int64_t), value,         intent(in)  :: bca_len
        integer(c_int64_t), value,         intent(in)  :: bcb_len
        type(c_ptr),        value,         intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(bca_len),intent(in)  :: bca
        real(c_double), dimension(bcb_len),intent(in)  :: bcb
      end subroutine bcpar_callback_c
    end interface

  ! DBC Callback     (Dbc)
    abstract interface
      subroutine Dbc_callback_c(ya_len, ya, yb_len, yb, &
                                dya_dim1, dya_dim2, dya, &
                                dyb_dim1, dyb_dim2, dyb, &
                               calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        integer(c_int64_t), value,         intent(in)  :: ya_len
        real(c_double), dimension(ya_len), intent(in)  :: ya
        integer(c_int64_t), value,         intent(in)  :: yb_len
        real(c_double), dimension(yb_len), intent(in)  :: yb
        integer(c_int64_t), value,         intent(in)  :: dya_dim1, dya_dim2
        integer(c_int64_t), value,         intent(in)  :: dyb_dim1, dyb_dim2
        type(c_ptr),        value,         intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(dya_dim1, dya_dim2), &
                                           intent(out) :: dya
        real(c_double), dimension(dyb_dim1, dyb_dim2), &
                                           intent(out) :: dyb
      end subroutine Dbc_callback_c
    end interface

    abstract interface
      subroutine Dbcpar_callback_c(ya_len, ya, yb_len, yb, &
                                p_len, p,                  &
                                dya_dim1, dya_dim2, dya,   &
                                dyb_dim1, dyb_dim2, dyb,   &
                                dpa_dim1, dpa_dim2, dpa,   &
                                dpb_dim1, dpb_dim2, dpb,   &
                               calls_pthrough) bind(C)
        use, intrinsic :: iso_c_binding
        ! In
        integer(c_int64_t), value,         intent(in)  :: ya_len
        real(c_double), dimension(ya_len), intent(in)  :: ya
        integer(c_int64_t), value,         intent(in)  :: yb_len
        real(c_double), dimension(yb_len), intent(in)  :: yb
        integer(c_int64_t), value,         intent(in)  :: p_len
        real(c_double), dimension(p_len),  intent(in)  :: p
        integer(c_int64_t), value,         intent(in)  :: dya_dim1, dya_dim2
        integer(c_int64_t), value,         intent(in)  :: dyb_dim1, dyb_dim2
        integer(c_int64_t), value,         intent(in)  :: dpa_dim1, dpa_dim2
        integer(c_int64_t), value,         intent(in)  :: dpb_dim1, dpb_dim2
        type(c_ptr),        value,         intent(in)  :: calls_pthrough
        ! Out
        real(c_double), dimension(dya_dim1, dya_dim2), &
                                           intent(out) :: dya
        real(c_double), dimension(dyb_dim1, dyb_dim2), &
                                           intent(out) :: dyb
        real(c_double), dimension(dpa_dim1, dpa_dim2), &
                                           intent(out) :: dpa
        real(c_double), dimension(dpb_dim1, dpb_dim2), &
                                           intent(out) :: dpb
      end subroutine Dbcpar_callback_c
    end interface

  type, public :: bvp_sol_wrapper
    type(bvp_sol) :: sol
    integer       :: state
  end type bvp_sol_wrapper

contains

  ! take sol_wrapper-pointer and generate a "handle" which can
  ! be used/saved in C/Julia to identify the sol_wrapper.
  function sol_wrapper_to_handle(sol_wrapper) result(handle)
    ! In
    type(bvp_sol_wrapper), pointer, intent(in) :: sol_wrapper
    ! Out
    type(c_ptr) :: handle
    !
    handle = c_loc(sol_wrapper)
  end function sol_wrapper_to_handle

  ! take a "handle" generated with sol_wrapper_to_handle and
  ! do the opposite.
  function handle_to_sol_wrapper(handle) result(sol_wrapper)
    ! In
    type(c_ptr), intent(in) :: handle
    ! Out
    type(bvp_sol_wrapper), pointer  :: sol_wrapper
    !
    call c_f_pointer(handle, sol_wrapper)
  end function handle_to_sol_wrapper

  ! return (some) scalar fields of sol_wrapper (and included sol)
  ! int_slots_len >= 9 is needed.
  subroutine get_sol_wrapper_details_c(handle, int_slots_len, int_slots) &
                         bind(C, name="get_sol_wrapper_details_c")
    ! In
    type(c_ptr),        value, intent(in) :: handle
    integer(c_int64_t), value, intent(in) :: int_slots_len
    ! Out
    integer(c_int64_t), dimension(int_slots_len), intent(out) :: int_slots
    !
    type(bvp_sol_wrapper), pointer :: sol_wrapper
    type(bvp_sol) :: sol
    integer :: state

    int_slots(:) = -1

    if (int_slots_len >= 9) then
      sol_wrapper => handle_to_sol_wrapper(handle)
      sol = sol_wrapper%sol
      state = sol_wrapper%state

      int_slots(1) = sol_wrapper%state
      
      if (state > 0) then
        int_slots(2) = sol%NODE
        int_slots(3) = sol%NPAR
        int_slots(4) = sol%LEFTBC
        int_slots(5) = sol%NPTS
        int_slots(6) = sol%INFO
        int_slots(7) = sol%MXNSUB
        if (state == 2) then
          int_slots(8) = size(sol%IWORK)
          int_slots(9) = size(sol%WORK)
        end if
      end if
    end if
  end subroutine get_sol_wrapper_details_c

  ! show fields of type sol_wrapper (and the included sol type)
  subroutine show_sol_wrapper(sol_wrapper)
    ! In
    type(bvp_sol_wrapper), pointer, intent(in) :: sol_wrapper
    !
    type(bvp_sol) :: sol
    character(len=*), parameter :: fmt_i  = '(A12, 5X, I18)'
    ! character(len=*), parameter :: fmt_h  = '(A12, 5X, A, Z16.16)'
    character(len=*), parameter :: fmt_vd = '(A12, 5X, *(F5.2))'
    integer :: k, state
    !
    sol = sol_wrapper%sol
    state = sol_wrapper%state
    write(*,*) "show_sol_wrapper:"
    ! write(*, fmt_h)  "handle"    , "0x", sol_wrapper_to_handle(sol_wrapper)
    write(*, fmt_i)  "state"     , state
    if (state > 0) then
      write(*, fmt_i)  "NODE"      , sol%node
      write(*, fmt_i)  "NPAR"      , sol%npar
      write(*, fmt_i)  "LEFTBC"    , sol%leftbc
      write(*, fmt_i)  "NPTS"      , sol%npts
      write(*, fmt_i)  "INFO"      , sol%info
      write(*, fmt_i)  "MXNSUB"    , sol%mxnsub
      write(*, fmt_vd) "X"         , sol%x
      do k = 1, size(sol%y,1)
        write(*, fmt_vd) merge("Y", " ", k==1), sol%y(k,:)
      end do
      if (sol%npar > 0) then
        write(*, fmt_vd) "PARAMETERS", sol%parameters
      end if
    end if
  end subroutine show_sol_wrapper

  subroutine show_sol_wrapper_c(handle) bind(C, name="show_sol_wrapper_c")
    ! In
    type(c_ptr), value, intent(in) :: handle
    !
    type(bvp_sol_wrapper), pointer :: sol_wrapper
    sol_wrapper => handle_to_sol_wrapper(handle)
    call show_sol_wrapper(sol_wrapper)
  end subroutine show_sol_wrapper_c

  ! allocates memory for sol_wrapper.
  ! for deallocation use destroy
  function create_sol_wrapper() result (sol_wrapper)
    ! Out
    type(bvp_sol_wrapper), pointer :: sol_wrapper
    !
    allocate(sol_wrapper)
    sol_wrapper%state = 0
  end function create_sol_wrapper

  function create_sol_wrapper_c() result (handle) &
                                 bind(C, name="create_sol_wrapper_c")
    ! Out
    type(c_ptr) :: handle
    !
    type(bvp_sol_wrapper), pointer :: sol_wrapper
    !
    sol_wrapper => create_sol_wrapper()
    handle = sol_wrapper_to_handle(sol_wrapper)
  end

  function copy_sol_wrapper(sol_in) result(sol_out)
    ! In
    type(bvp_sol_wrapper), pointer, intent(in) :: sol_in
    ! Out
    type(bvp_sol_wrapper), pointer             :: sol_out
    !
    integer                                    :: no_odes, no_par, no_pts
    !
    sol_out => create_sol_wrapper()
    sol_out%state = sol_in%state
    if (sol_in%state > 0) then
      no_odes = sol_in%sol%node
      no_par = sol_in%sol%npar
      no_pts = sol_in%sol%npts
      sol_out%sol%node = no_odes
      sol_out%sol%npar = no_par
      sol_out%sol%leftbc = sol_in%sol%leftbc
      sol_out%sol%npts = no_pts
      sol_out%sol%info = sol_in%sol%info
      sol_out%sol%mxnsub = sol_in%sol%mxnsub
      allocate(sol_out%sol%x(lbound(sol_in%sol%x,1):ubound(sol_in%sol%x,1)))
      sol_out%sol%x(:) = sol_in%sol%x(:)
      allocate(sol_out%sol%y( &
        lbound(sol_in%sol%y,1):ubound(sol_in%sol%y,1),  &
        lbound(sol_in%sol%y,2):ubound(sol_in%sol%y,2) ))
      sol_out%sol%y(:,:) = sol_in%sol%y(:,:)
      if (no_par > 0) then
        allocate(sol_out%sol%parameters( &
          lbound(sol_in%sol%parameters,1):ubound(sol_in%sol%parameters,1)))
        sol_out%sol%parameters(:) = sol_in%sol%parameters(:)
      end if
      if (sol_in%state == 2) then
        allocate(sol_out%sol%iwork( &
          lbound(sol_in%sol%iwork,1):ubound(sol_in%sol%iwork,1)))
        sol_out%sol%iwork(:) = sol_in%sol%iwork(:)
        allocate(sol_out%sol%work( &
          lbound(sol_in%sol%work,1):ubound(sol_in%sol%work,1)))
        sol_out%sol%work(:) = sol_in%sol%work(:)
      end if
    end if
  end function copy_sol_wrapper

  function copy_sol_wrapper_c(handle_in) result(handle_out) &
                                     bind(C, name="copy_sol_wrapper_c")
    ! In
    type(c_ptr), value, intent(in) :: handle_in
    ! Out
    type(c_ptr)                    :: handle_out
    !
    type(bvp_sol_wrapper), pointer :: sol_in, sol_out
    !
    sol_in => handle_to_sol_wrapper(handle_in)
    sol_out => copy_sol_wrapper(sol_in)
    handle_out = sol_wrapper_to_handle(sol_out)
  end function copy_sol_wrapper_c

  ! terminates sol (in sol_wrapper). In BVP_M-2 this is only possible
  ! for STATE==2. This method extend this behaviour for every state.
  ! Unfortunately this means for STATE==1 that this method has to
  ! deallocate the (allocated) sol-fields itself ...
  subroutine terminate_sol_wrapper(sol_wrapper)
    ! In
    type(bvp_sol_wrapper), pointer, intent(in) :: sol_wrapper
    !
    select case (sol_wrapper%state)
      case (0)
        continue
      case (1)
        ! workaround: deallocate X, Y and PARAMETERS
        deallocate(sol_wrapper%sol%X, sol_wrapper%sol%Y)
        if (sol_wrapper%sol%NPAR > 0) then
          deallocate(sol_wrapper%sol%PARAMETERS)
        end if
        sol_wrapper%state = 0
      case (2)
        call bvp_terminate(sol_wrapper%sol)
    end select
    sol_wrapper%state = 0
  end subroutine terminate_sol_wrapper

  subroutine terminate_sol_wrapper_c(handle) &
                             bind(C, name="terminate_sol_wrapper_c")
    ! In
    type(c_ptr), value, intent(in) :: handle
    !
    type(bvp_sol_wrapper), pointer :: sol_wrapper
    !
    sol_wrapper => handle_to_sol_wrapper(handle)
    call terminate_sol_wrapper(sol_wrapper)
  end subroutine terminate_sol_wrapper_c

  ! calls terminate_sol_wrapper and deallocates the sol_wrapper object.
  subroutine destroy_sol_wrapper(sol_wrapper)
    ! InOut
    type(bvp_sol_wrapper), pointer, intent(inout) :: sol_wrapper
    !
    call terminate_sol_wrapper(sol_wrapper)
    deallocate(sol_wrapper)
    nullify(sol_wrapper)
  end subroutine destroy_sol_wrapper

  ! Function to get the current grid. 
  ! x_len must be sol%NPTS and x must be a vector of this length
  ! where the grid points are saved.
  !
  ! This function returns an error code:
  !   0: No error
  !  -1: no grid available because sol_wrapper has wrong state 
  !  -2: x_len has wrong size
  function get_sol_wrapper_x_c(handle, x_len, x) result(error) &
                            bind(C, name="get_sol_wrapper_x_c")
    ! In
    type(c_ptr),        value,        intent(in)  :: handle
    integer(c_int64_t), value,        intent(in)  :: x_len
    ! Out
    real(c_double), dimension(x_len), intent(out) :: x
    integer(c_int64_t)                            :: error
    !
    type(bvp_sol_wrapper), pointer                :: sol_wrapper
    !
    error = -1
    sol_wrapper => handle_to_sol_wrapper(handle)
    if (sol_wrapper%state > 0) then
      if (x_len == sol_wrapper%sol%NPTS) then
        x(:) = sol_wrapper%sol%x(:)
        error = 0
      else
        error = -2
      end if
    end if
  end function get_sol_wrapper_x_c

  ! Function to get the current parameters.
  ! p_len must be sol%NPAR and p must be a vector of this length
  ! where the parameters are saved.
  ! 
  ! This function returns an error code:
  !   0: No error
  !  -1: no parameters available
  !  -2: p_len has wrong size
  !  -3: info field in sol was /= 0
  function get_sol_wrapper_params_c(handle, p_len, p) result(error) &
                                  bind(C, name="get_sol_wrapper_params_c")
    ! In
    type(c_ptr),        value,        intent(in)  :: handle
    integer(c_int64_t), value,        intent(in)  :: p_len
    ! Out
    real(c_double), dimension(p_len), intent(out) :: p
    integer(c_int64_t)                            :: error
    !
    type(bvp_sol_wrapper), pointer                :: sol_wrapper
    !
    error = 0
    sol_wrapper => handle_to_sol_wrapper(handle)
    if (sol_wrapper%state > 0) then
      if (sol_wrapper%sol%npar > 0) then
        if (p_len == sol_wrapper%sol%npar) then
          if (sol_wrapper%sol%info == 0) then
            call bvp_eval(sol_wrapper%sol, p)
          else
            error = -3
          end if
        else
          error = -2
        endif
      end if
    else
      error = -1
    end if
  end function get_sol_wrapper_params_c

  subroutine destroy_sol_wrapper_c(handle) &
                           bind(C, name="destroy_sol_wrapper_c")
    ! In
    type(c_ptr), value, intent(in) :: handle
    !
    type(bvp_sol_wrapper), pointer :: sol_wrapper
    sol_wrapper => handle_to_sol_wrapper(handle)
    call destroy_sol_wrapper(sol_wrapper)
  end subroutine destroy_sol_wrapper_c

  subroutine init_guess1_c(handle, &
        no_odes, no_left_bc, x_len, x, y_len, y, p_len, p, &
        max_num_subintervals) bind(C, name="init_guess1_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle
    integer(c_int64_t),            value, intent(in) :: no_odes
    integer(c_int64_t),            value, intent(in) :: no_left_bc
    integer(c_int64_t),            value, intent(in) :: x_len
    real(c_double),     dimension(x_len), intent(in) :: x
    integer(c_int64_t),            value, intent(in) :: y_len
    real(c_double),     dimension(y_len), intent(in) :: y
    integer(c_int64_t),            value, intent(in) :: p_len
    real(c_double),     dimension(p_len), intent(in) :: p
    integer(c_int64_t),            value, intent(in) :: max_num_subintervals
    !
    double precision, dimension(x_len) :: x_arr
    double precision, dimension(y_len) :: y_arr
    double precision, dimension(p_len) :: p_arr
    type(bvp_sol_wrapper), pointer :: p_wrapper
    !
    p_wrapper => handle_to_sol_wrapper(handle)
    x_arr = x
    y_arr = y

    if (p_len > 0) then
      p_arr = p
      p_wrapper%sol = bvp_init(no_odes, no_left_bc, x_arr, y_arr, p_arr, &
                     max_num_subintervals=max_num_subintervals )
    else
      p_wrapper%sol = bvp_init(no_odes, no_left_bc, x_arr, y_arr, &
                     max_num_subintervals=max_num_subintervals )
    end if
    p_wrapper%state = 1
  end subroutine init_guess1_c

  subroutine init_guess2_c(handle, &
        no_odes, no_left_bc, x_len, x, y_dim1, y_dim2, y, p_len, p, &
        max_num_subintervals) bind(C, name="init_guess2_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle
    integer(c_int64_t),            value, intent(in) :: no_odes
    integer(c_int64_t),            value, intent(in) :: no_left_bc
    integer(c_int64_t),            value, intent(in) :: x_len
    real(c_double),     dimension(x_len), intent(in) :: x
    integer(c_int64_t),            value, intent(in) :: y_dim1
    integer(c_int64_t),            value, intent(in) :: y_dim2
    real(c_double),    dimension(y_dim1, &
                                 y_dim2), intent(in) :: y
    integer(c_int64_t),            value, intent(in) :: p_len
    real(c_double),     dimension(p_len), intent(in) :: p
    integer(c_int64_t),            value, intent(in) :: max_num_subintervals
    !
    double precision, dimension(x_len)         :: x_arr
    double precision, dimension(y_dim1,y_dim2) :: y_arr
    double precision, dimension(p_len)         :: p_arr
    type(bvp_sol_wrapper), pointer :: p_wrapper
    !
    p_wrapper => handle_to_sol_wrapper(handle)
    x_arr = x
    y_arr = y

    if (p_len > 0) then
      p_arr = p
      p_wrapper%sol = bvp_init(no_odes, no_left_bc, x_arr, y_arr, p_arr, &
                               max_num_subintervals=max_num_subintervals)
    else
      p_wrapper%sol = bvp_init(no_odes, no_left_bc, x_arr, y_arr, &
                               max_num_subintervals=max_num_subintervals)
    end if
    p_wrapper%state = 1
  end subroutine init_guess2_c

  subroutine init_guess3_c(handle, &
        no_odes, no_left_bc, x_len, x, guess_fcn_ptr, guess_pthrough, &
        p_len, p, max_num_subintervals) bind(C, name="init_guess3_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle
    integer(c_int64_t),            value, intent(in) :: no_odes
    integer(c_int64_t),            value, intent(in) :: no_left_bc
    integer(c_int64_t),            value, intent(in) :: x_len
    real(c_double),     dimension(x_len), intent(in) :: x
    type(c_funptr),                value, intent(in) :: guess_fcn_ptr
    type(c_ptr),                   value, intent(in) :: guess_pthrough
    integer(c_int64_t),            value, intent(in) :: p_len
    real(c_double),     dimension(p_len), intent(in) :: p
    integer(c_int64_t),            value, intent(in) :: max_num_subintervals
    !
    double precision, dimension(x_len)         :: x_arr
    double precision, dimension(p_len)         :: p_arr
    procedure(guess_callback_c), pointer       :: guess_fcn
    type(bvp_sol_wrapper), pointer :: p_wrapper 
    !
    p_wrapper => handle_to_sol_wrapper(handle)
    x_arr = x
    call c_f_procpointer(guess_fcn_ptr, guess_fcn)

    if (p_len > 0) then
      p_arr = p
      p_wrapper%sol = bvp_init(no_odes, no_left_bc, x_arr, guess_fcn_proxy, &
                               p_arr, &
                               max_num_subintervals=max_num_subintervals)
    else
      p_wrapper%sol = bvp_init(no_odes, no_left_bc, x_arr, guess_fcn_proxy, &
                               max_num_subintervals=max_num_subintervals)
    end if
    p_wrapper%state = 1

    contains
      ! nested subroutine to make it possible to pass-through the
      ! guess_pthrough pointer to the caller.
      subroutine guess_fcn_proxy(x_point, guess_vector)
        ! In
        double precision,               intent(in)        :: x_point
        ! Out
        double precision, dimension(no_odes), intent(out) :: guess_vector
        !
        call guess_fcn(x_point, no_odes, guess_vector, guess_pthrough)
      end subroutine guess_fcn_proxy
  end subroutine init_guess3_c

  ! si_dim == 0  => si_matrix not given 
  ! otherweise si_dim must be a (no_odes, no_odes) matrix
  ! handle_guess is *not* changed [because BVP_SOLVER is
  ! destructive for the SOL_INIT argument, a copy is made]
  !
  ! error
  !       0: Ok
  !      -1: handle_guess: state was neither 1 or 2
  !      -2: si_dim was neither 0 nor no_odes
  !      -3: f_fcn_ptr was NULL
  !      -4: bc_fcn_ptr was NULL
  !      -5: error_ret_len < 5
  function solve_c(handle_guess, handle_out, f_fcn_ptr, bc_fcn_ptr, &
                   si_dim, si_matrix, method, tol, &
                   dfdy_fcn_ptr, dbcdy_fcn_ptr, trace, error_control, &
                   error_ret_len, error_ret, calls_pthrough) &
                   result(error) bind(C, name="solve_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle_guess
    type(c_funptr),                value, intent(in) :: f_fcn_ptr
    type(c_funptr),                value, intent(in) :: bc_fcn_ptr
    integer(c_int64_t),            value, intent(in) :: si_dim
    real(c_double),    dimension(si_dim, &
                                 si_dim), intent(in) :: si_matrix
    integer(c_int64_t),            value, intent(in) :: method
    real(c_double),                value, intent(in) :: tol
    type(c_funptr),                value, intent(in) :: dfdy_fcn_ptr
    type(c_funptr),                value, intent(in) :: dbcdy_fcn_ptr
    integer(c_int64_t),            value, intent(in) :: trace
    integer(c_int64_t),            value, intent(in) :: error_control
    integer(c_int64_t),            value, intent(in) :: error_ret_len
    type(c_ptr),                   value, intent(in) :: calls_pthrough
    ! Out
    type(c_ptr),                          intent(out):: handle_out
    real(c_double), dimension(error_ret_len), &
                                          intent(out):: error_ret
    integer(c_int64_t)                               :: error
    !
    type(bvp_sol_wrapper), pointer         :: p_guess, p_out
    type(bvp_sol)                          :: sol_guess, sol_init
    integer                                :: fcn_flags, no_odes, no_left_bc
    integer                                :: no_par
    double precision                       :: cond, cerror, reerror
    double precision                       :: hoerror, dcerror
    procedure(f_callback_c), pointer       :: f_fcn
    procedure(fpar_callback_c), pointer    :: fpar_fcn
    procedure(Df_callback_c), pointer      :: Df_fcn
    procedure(Dfpar_callback_c), pointer   :: Dfpar_fcn
    procedure(bc_callback_c), pointer      :: bc_fcn
    procedure(bcpar_callback_c), pointer   :: bcpar_fcn
    procedure(Dbc_callback_c), pointer     :: Dbc_fcn
    procedure(Dbcpar_callback_c), pointer  :: Dbcpar_fcn
    !
    fcn_flags = 0 !   Bit 0 (1): singularity matrix given
                  !   Bit 1 (2): dfdy_fcn given
                  !   Bit 2 (4): dbcdy_fcn given
                  !   Bit 3 (8): with parameters (no_par > 0)
                  !  [pa][db][df][si]
    error = 0
    handle_out = C_NULL_PTR
    nullify(f_fcn, fpar_fcn, Df_fcn, Dfpar_fcn)
    nullify(bc_fcn, bcpar_fcn, Dbc_fcn, Dbcpar_fcn)
    p_guess => handle_to_sol_wrapper(handle_guess)
    if (p_guess%state > 0) then
      sol_guess = p_guess%sol
      no_odes = sol_guess%node
      no_par = sol_guess%npar
      no_left_bc = sol_guess%leftbc
      if (no_par > 0) then
        fcn_flags = ior(fcn_flags, 8)
        sol_init = bvp_init(sol_guess%node, sol_guess%leftbc, &
          sol_guess%x(1:sol_guess%npts),                      &
          sol_guess%y(1:sol_guess%node,1:sol_guess%npts),     &
          sol_guess%parameters, sol_guess%mxnsub)
      else
        sol_init = bvp_init(sol_guess%node, sol_guess%leftbc, &
          sol_guess%x(1:sol_guess%npts),                      &
          sol_guess%y(1:sol_guess%node,1:sol_guess%npts),     &
          max_num_subintervals=sol_guess%mxnsub)
      end if
      ! sol_init ready to use; fields will be deallocated by bvp_solver
      if (si_dim == 0) then
        continue
      else if (si_dim == no_odes) then
        fcn_flags = ior(fcn_flags, 1)
      else
        error =-2
      end if
      if (c_associated(f_fcn_ptr)) then
        if (no_par > 0) then
          call c_f_procpointer(f_fcn_ptr, fpar_fcn)
        else
          call c_f_procpointer(f_fcn_ptr, f_fcn)
        end if
      else
        error = -3
      end if
      if (c_associated(bc_fcn_ptr)) then
        if (no_par > 0) then
          call c_f_procpointer(bc_fcn_ptr, bcpar_fcn)
        else
          call c_f_procpointer(bc_fcn_ptr, bc_fcn)
        end if
      else
        error = -4
      end if
      if (c_associated(dfdy_fcn_ptr)) then
        fcn_flags = ior(fcn_flags, 2)
        if (no_par > 0) then
          call c_f_procpointer(dfdy_fcn_ptr, Dfpar_fcn)
        else
          call c_f_procpointer(dfdy_fcn_ptr, Df_fcn)
        end if
      end if 
      if (c_associated(dbcdy_fcn_ptr)) then
        fcn_flags = ior(fcn_flags, 4)
        if (no_par > 0) then
          call c_f_procpointer(dbcdy_fcn_ptr, Dbcpar_fcn)
        else
          call c_f_procpointer(dbcdy_fcn_ptr, Dbc_fcn)
        end if
      end if 
      if (error_ret_len < 5) then
        error = -5
      end if
      if (error == 0) then
        error_ret(1:5) = -1.0d0
        p_out => create_sol_wrapper()
        select case (fcn_flags)
          case (0)        ! [pa][db][df][si] = [  ][  ][  ][  ]
            p_out%sol = bvp_solver( sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              method        = method,            &
              tol           = tol,               &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (8)        ! [pa][db][df][si] = [XX][  ][  ][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              method        = method,            &
              tol           = tol,               &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (1)        ! [pa][db][df][si] = [  ][  ][  ][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (9)        ! [pa][db][df][si] = [XX][  ][  ][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (2)        ! [pa][db][df][si] = [  ][  ][XX][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Df_fcn_proxy,      &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (10)       ! [pa][db][df][si] = [XX][  ][XX][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Dfpar_fcn_proxy,   &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (3)        ! [pa][db][df][si] = [  ][  ][XX][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Df_fcn_proxy,      &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (11)       ! [pa][db][df][si] = [XX][  ][XX][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Dfpar_fcn_proxy,   &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (4)        ! [pa][db][df][si] = [  ][XX][  ][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              method        = method,            &
              tol           = tol,               &
              dbcdy         = Dbc_fcn_proxy,     &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (12)       ! [pa][db][df][si] = [XX][XX][  ][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              method        = method,            &
              tol           = tol,               &
              dbcdy         = Dbcpar_fcn_proxy,  &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (5)        ! [pa][db][df][si] = [  ][XX][  ][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              dbcdy         = Dbc_fcn_proxy,     &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (13)       ! [pa][db][df][si] = [XX][XX][  ][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              dbcdy         = Dbcpar_fcn_proxy,  &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (6)        ! [pa][db][df][si] = [  ][XX][XX][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Df_fcn_proxy,      &
              dbcdy         = Dbc_fcn_proxy,     &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (14)       ! [pa][db][df][si] = [XX][XX][XX][  ]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Dfpar_fcn_proxy,   &
              dbcdy         = Dbcpar_fcn_proxy,  &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (7)        ! [pa][db][df][si] = [  ][XX][XX][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = f_fcn_proxy,       &
              bcsub         = bc_fcn_proxy,      &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Df_fcn_proxy,      &
              dbcdy         = Dbc_fcn_proxy,     &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case (15)       ! [pa][db][df][si] = [XX][XX][XX][XX]
            p_out%sol = bvp_solver(sol_init=sol_init, &
              fsub          = fpar_fcn_proxy,    &
              bcsub         = bcpar_fcn_proxy,   &
              singularterm  = si_matrix,         &
              method        = method,            &
              tol           = tol,               &
              dfdy          = Dfpar_fcn_proxy,   &
              dbcdy         = Dbcpar_fcn_proxy,  &
              trace=trace, stop_on_fail=.false., &
              error_control=error_control, cond=cond, cerror=cerror, &
              reerror=reerror, hoerror=hoerror, dcerror= dcerror)
          case default
            error = -100
        end select
        error_ret(1) = cond
        error_ret(2) = cerror
        error_ret(3) = reerror
        error_ret(4) = hoerror
        error_ret(5) = dcerror
        p_out%state = merge(2, 0, p_out%sol%info == 0)
        handle_out = sol_wrapper_to_handle(p_out)
      end if
    else
      error = -1
    end if

    contains 
      ! nested subroutines to make it possible to pass-through the
      ! calls_pthrough pointer to the caller.

      ! RHS proxies
        subroutine f_fcn_proxy(x_point, y_state, f_values)
          ! In
          double precision,                     intent(in)  :: x_point
          double precision, dimension(no_odes), intent(in)  :: y_state
          ! Out
          double precision, dimension(no_odes), intent(out) :: f_values
          !
          call f_fcn(x_point, no_odes, y_state, &
                     no_odes, f_values, calls_pthrough)
        end subroutine f_fcn_proxy

        subroutine fpar_fcn_proxy(x_point, y_state, params, f_values)
          ! In
          double precision,                     intent(in)  :: x_point
          double precision, dimension(no_odes), intent(in)  :: y_state
          double precision, dimension(no_par),  intent(in)  :: params
          ! Out
          double precision, dimension(no_odes), intent(out) :: f_values
          !
          ! write(*,*) "fpar_fcn_proxy xpoint=", x_point, " y_state=", y_state, &
          !            " params=", params
          call fpar_fcn(x_point, no_odes, y_state, &
                     no_par, params, &
                     no_odes, f_values, calls_pthrough)
        end subroutine fpar_fcn_proxy

      ! DRHS proxies
        subroutine Df_fcn_proxy(x_point, y_state, dfdy_values)
          ! In
          double precision,                     intent(in)  :: x_point
          double precision, dimension(no_odes), intent(in)  :: y_state
          ! Out
          double precision, dimension(no_odes, no_odes), &
                                                intent(out) :: dfdy_values
          call Df_fcn(x_point, no_odes, y_state, &
                      no_odes, no_odes, dfdy_values, calls_pthrough)
        end subroutine Df_fcn_proxy

        subroutine Dfpar_fcn_proxy(x_point, y_state, params, &
                                   dfdy_values, dfdp_values)
          ! In
          double precision,                     intent(in)  :: x_point
          double precision, dimension(no_odes), intent(in)  :: y_state
          double precision, dimension(no_par),  intent(in)  :: params
          ! Out
          double precision, dimension(no_odes, no_odes), &
                                                intent(out) :: dfdy_values
          double precision, dimension(no_odes, no_par), &
                                                intent(out) :: dfdp_values
          call Dfpar_fcn(x_point, no_odes, y_state, &
                      no_par, params, no_odes, no_odes, dfdy_values, &
                      no_odes, no_par, dfdp_values, calls_pthrough)
        end subroutine Dfpar_fcn_proxy

      ! BC proxies
        subroutine bc_fcn_proxy(ya, yb, bca, bcb)
          ! In
          double precision, dimension(no_odes), intent(in)  :: ya
          double precision, dimension(no_odes), intent(in)  :: yb
          ! Out
          double precision, dimension(no_left_bc), intent(out) :: bca
          double precision, dimension(no_odes+no_par-no_left_bc), &
                                                   intent(out) :: bcb
          !
          call bc_fcn(no_odes, ya, no_odes, yb, no_left_bc, bca, &
            no_odes+no_par-no_left_bc, bcb, calls_pthrough)
        end subroutine bc_fcn_proxy

        subroutine bcpar_fcn_proxy(ya, yb, params, bca, bcb)
          ! In
          double precision, dimension(no_odes), intent(in)  :: ya
          double precision, dimension(no_odes), intent(in)  :: yb
          double precision, dimension(no_par),  intent(in)  :: params
          ! Out
          double precision, dimension(no_left_bc), intent(out) :: bca
          double precision, dimension(no_odes+no_par-no_left_bc), &
                                                   intent(out) :: bcb
          !
          call bcpar_fcn(no_odes, ya, no_odes, yb, &
            no_par, params,  no_left_bc, bca, &
            no_odes+no_par-no_left_bc, bcb, calls_pthrough)
        end subroutine bcpar_fcn_proxy

      ! DBC proxies
        subroutine Dbc_fcn_proxy(ya, yb, dya, dyb)
          ! In
          double precision, dimension(no_odes), intent(in)  :: ya
          double precision, dimension(no_odes), intent(in)  :: yb
          ! Out
          double precision, dimension(no_left_bc, no_odes), &
                                                intent(out) :: dya
          double precision, dimension(no_odes+no_par-no_left_bc, no_odes), &
                                                intent(out) :: dyb
          !
          call Dbc_fcn(no_odes, ya, no_odes, yb, &
                       no_left_bc, no_odes, dya, &
                       no_odes+no_par-no_left_bc, no_odes, dyb, &
                       calls_pthrough)
        end subroutine Dbc_fcn_proxy

        subroutine Dbcpar_fcn_proxy(ya, yb, dya, dyb, p, dpa, dpb)
          ! In
          double precision, dimension(no_odes), intent(in)  :: ya
          double precision, dimension(no_odes), intent(in)  :: yb
          double precision, dimension(no_par),  intent(in)  :: p
          ! Out
          double precision, dimension(no_left_bc, no_odes), &
                                                intent(out) :: dya
          double precision, dimension(no_odes+no_par-no_left_bc, no_odes), &
                                                intent(out) :: dyb
          double precision, dimension(no_left_bc, no_par), &
                                                intent(out) :: dpa
          double precision, dimension(no_odes+no_par-no_left_bc, no_par), &
                                                intent(out) :: dpb
          !
          call Dbcpar_fcn(no_odes, ya, no_odes, yb, &
                       no_par, p, &
                       no_left_bc, no_odes, dya, &
                       no_odes+no_par-no_left_bc, no_odes, dyb, &
                       no_left_bc, no_par, dpa, &
                       no_odes+no_par-no_left_bc, no_par, dpb, &
                       calls_pthrough)
        end subroutine Dbcpar_fcn_proxy

  end function solve_c

  ! error
  !       0: Ok
  !      -1: handle: state /= 2
  !      -2: z_len /= no_odes
  !      -3: dz_len /= 0  and dz_len /= no_odes
  function eval_s_sol_c(handle, x, z_len, z, dz_len, dz) result(error) &
          bind(C, name="eval_s_sol_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle
    real(c_double),                value, intent(in) :: x
    integer(c_int64_t),            value, intent(in) :: z_len
    integer(c_int64_t),            value, intent(in) :: dz_len
    ! Out
    real(c_double), dimension(z_len),     intent(out):: z
    real(c_double), dimension(dz_len),    intent(out):: dz
    integer(c_int64_t)                               :: error
    !
    type(bvp_sol_wrapper), pointer                   :: sol_wrapper
    integer                                          :: no_odes
    !
    error = 0
    sol_wrapper => handle_to_sol_wrapper(handle)
    if (sol_wrapper%state == 2) then
      no_odes = sol_wrapper%sol%node
      if (z_len /= no_odes) then
        error = -2
      end if
      if (dz_len /= 0) then
        if (dz_len /= no_odes) then
          error = -3
        end if
      end if
      if (error == 0) then
        if (dz_len == 0) then
          call bvp_eval(sol_wrapper%sol, x, z)
        else
          call bvp_eval(sol_wrapper%sol, x, z, dz)
        end if
      end if
    else
      error = -1
    end if
  end function eval_s_sol_c

  ! error
  !       0: Ok
  !      -1: handle: state /= 2
  !      -2: z_dim1 /= no_odes or z_dim2 /= x_len
  !      -3: (dz_dim1 /= 0 or dz_dim2 /= 0) and 
  !          (dz_dim1 /= no_odes or dz_dim2 /= x_len)
  function eval_v_sol_c(handle, x_len, x, z_dim1, z_dim2, z,  &
                        dz_dim1, dz_dim2, dz) result (error) &
                        bind(C, name="eval_v_sol_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle
    integer(c_int64_t),            value, intent(in) :: x_len
    real(c_double), dimension(x_len),     intent(in) :: x
    integer(c_int64_t),            value, intent(in) :: z_dim1
    integer(c_int64_t),            value, intent(in) :: z_dim2
    integer(c_int64_t),            value, intent(in) :: dz_dim1
    integer(c_int64_t),            value, intent(in) :: dz_dim2
    ! Out
    real(c_double), dimension(z_dim1, z_dim2), &
                                          intent(out):: z
    real(c_double), dimension(dz_dim1, dz_dim2), &
                                          intent(out):: dz
    integer(c_int64_t)                               :: error
    !
    type(bvp_sol_wrapper), pointer                   :: sol_wrapper
    integer                                          :: no_odes
    !
    error = 0
    sol_wrapper => handle_to_sol_wrapper(handle)
    if (sol_wrapper%state == 2) then
      no_odes = sol_wrapper%sol%node
      if (z_dim1 /= no_odes .or. z_dim2 /= x_len) then
        error = -2
      end if
      if (dz_dim1 /= 0 .or. dz_dim2 /= 0) then
        if (dz_dim1 /= no_odes .or. dz_dim2 /= x_len) then
          error = -3
        end if
      end if
      if (error == 0) then
        if (dz_dim1 == 0) then
          call bvp_eval(sol_wrapper%sol, x, z)
        else
          call bvp_eval(sol_wrapper%sol, x, z, dz)
        end if
      end if
    else
      error = -1
    end if
  end function eval_v_sol_c

  ! extend solution to new interval and new states.
  ! error
  !       0: Ok
  !      -1: handle_in has not state==2
  !      -2: yanew_len /= no_odes
  !      -3: ybnew_len /= no_odes
  !      -4: plen /= 0 and plen /= no_par
  function extend_sol_s_c(handle_in, handle_out, anew, bnew, &
                        yanew_len, yanew, ybnew_len, ybnew, p_len, p, &
                        max_num_subintervals) result(error) &
                        bind(C, name="extend_sol_s_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle_in
    real(c_double),                value, intent(in) :: anew, bnew
    integer(c_int64_t),            value, intent(in) :: yanew_len
    real(c_double), dimension(yanew_len), intent(in) :: yanew
    integer(c_int64_t),            value, intent(in) :: ybnew_len
    real(c_double), dimension(ybnew_len), intent(in) :: ybnew
    integer(c_int64_t),            value, intent(in) :: p_len
    real(c_double), dimension(p_len),     intent(in) :: p
    integer(c_int64_t),            value, intent(in) :: max_num_subintervals
    ! Out
    type(c_ptr),                          intent(out):: handle_out
    integer(c_int64_t)                               :: error
    !
    type(bvp_sol_wrapper), pointer                   :: p_in, p_out
    integer                                          :: no_odes, no_par
    integer                                          :: max_subint
    !
    error = 0
    handle_out = C_NULL_PTR
    nullify(p_in, p_out)
    p_in => handle_to_sol_wrapper(handle_in)
    if (p_in%state == 2) then
      no_odes = p_in%sol%node
      no_par = p_in%sol%npar
      max_subint = p_in%sol%mxnsub
      if (max_num_subintervals > 1) then
        max_subint = max_num_subintervals
      endif
      if (yanew_len /= no_odes) then
        error = -2
      end if
      if (ybnew_len /= no_odes) then
        error = -3
      end if
      if (p_len /= 0 .and. p_len /= no_par) then
        error = -4
      end if
      if (error == 0) then
        p_out => create_sol_wrapper()
        if (p_len > 0) then
          p_out%sol = bvp_extend(p_in%sol, anew, bnew, yanew, ybnew, &
                                 p, max_subint)
        else
          p_out%sol = bvp_extend(p_in%sol, anew, bnew, yanew, ybnew, &
                                 max_num_subintervals=max_subint)
        end if
        p_in%state = 0
        p_out%state = 1
        handle_out = sol_wrapper_to_handle(p_out)
      end if
    else
      error = -1
    end if
  end function extend_sol_s_c

  ! extend solution to new interval using extrapolation.
  ! error
  !       0: Ok
  !      -1: handle_in has not state==2
  !      -4: plen /= 0 and plen /= no_par
  function extend_sol_e_c(handle_in, handle_out, anew, bnew, order, &
                        p_len, p, max_num_subintervals) result(error) &
                        bind(C, name="extend_sol_e_c")
    ! In
    type(c_ptr),                   value, intent(in) :: handle_in
    real(c_double),                value, intent(in) :: anew, bnew
    integer(c_int64_t),            value, intent(in) :: order
    integer(c_int64_t),            value, intent(in) :: p_len
    real(c_double), dimension(p_len),     intent(in) :: p
    integer(c_int64_t),            value, intent(in) :: max_num_subintervals
    ! Out
    type(c_ptr),                          intent(out):: handle_out
    integer(c_int64_t)                               :: error
    !
    type(bvp_sol_wrapper), pointer                   :: p_in, p_out
    integer                                          :: no_odes, no_par
    integer                                          :: max_subint
    !
    error = 0
    handle_out = C_NULL_PTR
    nullify(p_in, p_out)
    p_in => handle_to_sol_wrapper(handle_in)
    if (p_in%state == 2) then
      no_odes = p_in%sol%node
      no_par = p_in%sol%npar
      max_subint = p_in%sol%mxnsub
      if (max_num_subintervals > 1) then
        max_subint = max_num_subintervals
      endif
      if (p_len /= 0 .and. p_len /= no_par) then
        error = -4
      end if
      if (error == 0) then
        p_out => create_sol_wrapper()
        if (p_len > 0) then
          p_out%sol = bvp_extend(p_in%sol, anew, bnew, order, p, max_subint)
        else
          p_out%sol = bvp_extend(p_in%sol, anew, bnew, order, &
                                 max_num_subintervals=max_subint)
        end if
        p_in%state = 0
        p_out%state = 1
        handle_out = sol_wrapper_to_handle(p_out)
      end if
    else 
      error = -1
    end if
  end function extend_sol_e_c

end module bvp_m_proxy

! vim:cc=79:fdm=indent:
