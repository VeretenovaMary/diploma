! The character ";" terminates a statement, except when the ";" appears in a character context or in a comment
PROGRAM ! this one ; in a comment
PRINT *, "this ; won't terminate stmt" ; ! this will be ignored
PRINT *, "but next " ; PRINT *, "will terminate stmt"
! If a ";" separator is followed by zero or more blanks and one or more ";" separators, the sequence from the first ";" to 
!        the last, inclusive, is interpreted as a single ";" separator.
PRINT *, "just two " ; ;  ;;; ;; PRINT *, "statements"
END