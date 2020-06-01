! If character position 6 contains any character other than blank or zero, character positions 7–72 
!        of the line constitute a continuation of the preceding noncomment line
      PROGRAM
   ! there is no continuation
      PRINT *,
     a "there is a line continuation"
      PRINT *,
     ! "there is a line continuation"
      PRINT *,
     ; "there is a line continuation"
     0PRINT *, "there is no line continuation"
      PRINT *, "character context
     . line continuation"
      END