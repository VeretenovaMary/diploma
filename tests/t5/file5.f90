! statement continuation. The character "&" is used to indicate that the current statement is continued on the next line that isnot a comment line
! statement continuation. Comments may occur within a continued statement. 
! If an "&" not in a comment is the last nonblank character on a line or the last nonblank character
!        before an "!", the statement is continued on the next line that is not a comment line. 
! If the first nonblank character on the next noncomment line is an "&", the statement continues at the next
!        character position following the "&"; otherwise, it continues with the first character position of the next noncomment line.
PROGRAM
PRINT *, & ! comment
   "statement"
PRINT *, &  ! & this character won't be printed
&   "statement"
PRINT *,
! comment
& "statement"
end
