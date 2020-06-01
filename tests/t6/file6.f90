! If a character context is to be continued, the "&" shall be the last nonblank character on the line and shall not be followed 
!        by commentary. An "&" shall be the first nonblank character on the next line that is not a comment line 
!        and the statement continues with the next character following the "&".
! A label may precede any statement not forming part of another statement.
PROGRAM
PRINT *, "this string &
&won't be separated"
PRINT *, "this also &
! comment line
&won't be"
PRINT *, "this character ! &
&will be printed" ! comment
PRINT *, "this characters !&
! won't be characters !& printed 
&! will be printed" ! comment
PRINT *, "this character &
&! will be printed" ! comment
PRINT *, "this character &&
& will be printed"
PRINT *, "this characters &&
&& will be printed"
PRINT *, "this character &
&& will be printed"
END