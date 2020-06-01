! this test cannot be successfully translated in compiler
! it checks if statements determined correctly 

module a                                        ! module stmt
integer k
end
module b                                        ! module stmt
integer j
end
module procedure                                ! module stmt
end
modulec () = something                          ! not module stmt
moduled = something                             ! not module stmt


interface
  interface name1
    interface assignment (=)
      interface operator (.OR.)
        interface operator (==)
          module procedure name2                ! not module stmt
          use a, ONLY: k                        ! correct
          use b, j => j1                        ! correct
          use procedure                         ! correct
          use c() = something                   ! not use stmt
          use d = something                     ! not use stmt
          interface sdfadf (asd) = sdfdf        ! not interface stmt
        end interface operator (==)
      end interface operator (.OR.)
    end interface assignment (=)
  end interface name
end interface

print *, "next not include stmt"; include "file17_inc.f90"
include "file17_inc.f90" ; print *, "this is also not include stmt"
include "file17_inc.f90"                            ! include stmt