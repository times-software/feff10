!
!
!    \ | /            /##|    @@@@  @   @@@@@   |      |             @    @
!     \|/ STAR       /###|   @      @   @     __|__    |             @    @
!  ----*----        /####|  @       @   @@@@    |      |___  __  __  @@@@@@
!     /|\          /#####|   @      @   @       |      |   \   \/         @
!    / | \         |#####|    @@@@  @   @       \___/  \___/ __/\__       @
!                  |#####|________________________________________________
!                 ||#####|                 ___________________            |
!        __/|_____||#####|________________|&&&&&&&&&&&&&&&&&&&||          |
!<\\\\\\\\_ |_____________________________|&&& 29 Nov 2009  &&||          |
!          \|     ||#####|________________|&&&&&&&&&&&&&&&&&&&||__________|
!                  |#####|
!                  |#####|                Version 4.1.0 Release
!                  |#####|
!                 /#######\
!                |#########|
!                    ====
!                     ||
!           An extended tool box of fortran routines for manipulating CIF data.
!                     ||
!                     ||  CIFtbx Version 4
!                     ||        by
!                     ||
!                     ||  Sydney R. Hall (syd@crystal.uwa.edu.au)
!                     ||  Crystallography Centre
!                     ||  University of Western Australia
!                     ||  Nedlands 6009, AUSTRALIA
!                     ||
!                     ||       and
!                     ||
!                     ||  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
!                     ||  Bernstein + Sons
!                     ||  5 Brewster Lane
!                     ||  Bellport, NY 11713, U.S.A.
!                     ||
! The latest program source and information is available from:
!                     ||
! Em: syd@crystal.uwa.edu.au       ,-_|\      Sydney R. Hall
! sendcif@crystal.uwa.edu.au      /     \     Crystallography Centre
! Fx: +61 9 380 1118  ||      --> *_,-._/     University of Western Australia
! Ph: +61 9 380 2725  ||               v      Nedlands 6009, AUSTRALIA
!                     ||
!                     ||
!_____________________||_____________________________________________________
!
! This is a version of CIFtbx which has been extended to work with CIF2, DDLm,
! DDL 2 and mmCIF as well as with DDL 1.4 and core CIF dictionaries.  CIFtbx
! version 1 was written by Sydney R. Hall (see Hall, S. R., "CIF Applications
! IV.  CIFtbx: a Tool Box for Manipulating CIFs,"  J. Appl. Cryst (1993). 26,
! 482-494.  The revisions for version 2 were done by Herbert J. Bernstein
! and Sydney R. Hall (see Hall, S. R. and Bernstein, H. J., "CIFtbx 2:
! Extended Tool Box for Manipulating CIFs," J. Appl. Cryst.(1996). 29,
! 598-603)
!
! The revisions for releases 3 and 4 were done by Herbert J. Bernstein, work
! funded in part by the International Union of Crystallography
!
!___________________________________________________________________________
!
!
!    GENERAL TOOLS
!
!
!    init_      Sets the device numbers of files.   (optional)
!               [logical function always returned .true.]
!
!               <input CIF dev number> Set input CIF device     (def=1)
!
!               <output CIF dev number>Set output CIF device    (def=2)
!
!               <diracc dev number>    Set direct access formatted
!                                      scratch device number    (def=3)
!
!               <error  dev number>    Set error message device (def=6)
!
!
!
!    dict_      Requests a CIF dictionary be used for various data checks.
!               [logical function returned as .true. if the name dictionary
!               was opened and if the check codes are recognisable.  The
!               data item names used in the first dictionary loaded are
!               considered to be preferred by the user to aliases found
!               in dictionaries loaded in later calls.  On exit from dict_
!               the variable dicname_ is either equal to the filename, or,
!               if the dictionary had a value for the tag dictionary_name
!               of dictionary.title, dicname_ is set to that value.
!               The variable dicver_ is blank or set to the value of
!               _dictionary_version or of _dictionary.version  The check codes
!               'catck' and 'catno' turn on and off checking of dictionary
!               catgeory conventions.  The default is 'catck'.  The check
!               codes 'parck' and 'parno' turn on and off checking of
!               parent-child relationships.  The default of 'parck'.  Three check
!               codes control the handling of tags from the current dictionary
!               which duplicate tags from a dictionary loaded earlier.  These
!               codes ('first', 'final' and 'nodup') have effect only for the
!               current call to dict_  The default is 'first'.]
!
!               <dictionary filename>  A CIF dictionary in DDL format
!                                      or blank if just setting flags
!                                      or resetting the dictionary
!
!               <check code string>    The codes specifying the types of
!                                      checks to be applied to the CIF.
!
!                                      'valid'  data name validation check.
!                                      'dtype'  data item data type check.
!                                      'catck'  check datanames against
!                                               categories
!                                      'catno'  don't check datanames against
!                                               categories
!                                      'parck'  check datanames against
!                                               parent-child relationships
!                                      'parno'  don't check datanames against
!                                               parent-child relationships
!                                      'first'  accept first dictionary's
!                                               definitions of duplicate tags
!                                      'final'  accept final dictionary's
!                                               definitions of duplicate tags
!                                      'nodup'  do not accept duplicate tag
!                                               definitions
!                                      'parck'  check datanames against parent-
!                                               child relationahips
!                                      'parno'  don't check datanames against
!                                               parent-child relationships
!                                      'reset'  switch off checking flags
!                                      'close'  close existing dictionaries
!
!___________________________________________________________________________
!
!
!   CIF ACCESS TOOLS  ("the get_ing commands")
!
!
!
!    ocif_      Opens the CIF containing the required data.
!               [logical function returned .true. if CIF opened]
!
!               <CIF filename>        A blank name signals that the
!                                     currently open input CIF file
!                                     will be read.
!
!
!
!    data_      Identifies the data block containing the data to be requested.
!               [logical function returned .true. if block found]
!
!               <data block name>     A blank name signals that the next
!                                     encountered block is used (the block
!                                     name is stored in the variable bloc_).
!
!
!    bkmrk_     Saves or restores the current position so that data from
!               elsewhere in the cif can be examined.
!               [logical function returned as .true. on save if there was
!               room in internal storage to hold the current position, .true.
!               on restore if the bookmark number used was valid.  If the
!               argument is zero, the call is to save the position and return
!               the bookmark number in the argument.  If the argument is
!               non-zero, the call is to restore the position saved for the
!               bookmark number given.  The bookmark and the argument are
!               cleared.  The position set on return allow reprocessing of
!               the data item or loop row last processed when the bookmark
!               was placed.
!
!               NOTE:  All bookmarks are cleared by a call to data_]
!
!               <integer variable>    Bookmark number
!
!
!    find_      Find the location of the requested item in the CIF.
!               [The argument "name" may be a data item name, blank
!               for the next such item.  The argument "type" may be
!               blank for unrestricted acceptance of any non-comment
!               string (use cmnt_ to see comments), including loop headers,
!               "name" to accept only the name itself and "valu"
!               to accept only the value, or "head" to position to the
!               head of the CIF.  Except when the "head" is requested,
!               the position is left after the data item provided.  If the
!               item found is of type "name", posnam_ is set, otherwise,
!               posval_]
!
!               <data item name>      A blank name signals that the next
!                                     item of the type specified is needed
!
!               <data item type>      blank, 'head', 'name' or 'valu'
!
!               <character variable>  Returned string is of length long_.
!
!
!
!    test_      Identify the data attributes of the named data item.
!               [logical function returned as .true. if the item is present or
!               .false. if it is not. The data attributes are stored in the
!               common variables list_, type_, dictype_, diccat_ and dicname_.
!               The list, array, tuple or table attribites are stored in
!               ttype_, depth_ index_.
!
!               The values in dictype_, diccat_ and dicname_ are valid
!               whether or not the data item is found in the input CIF, as
!               long as the named data item is found in the dictionaries
!               declared by calls to dict_.  The data item name found
!               in the input CIF is stored in tagname_.  The appropriate
!               column numbers are stored in posnam_, posval_, posend_ and (for
!               numbers) in posdec_.  The quotation mark, if any, used is
!               stored in quote_.
!
!               list_ is an integer variable containing the sequential number
!               of the loop block in the data block. If the item is not within
!               a loop structure this value will be zero.
!
!               type_ is a character*4 variable with the possible values:
!                      'numb'  for number data
!                      'char'  for character data
!                      'text'  for text data
!                      'null'  if data missing or '?' or '.'
!                              also used for blank quoted fields if
!                              nblank_ is true
!
!               ttype_ is a character*4 variable with the container type:
!                      'list'  for list or array data   [item,...]
!                      'tupl'  for tuple data           (item,...)
!                      'tabl'  for table data           {item,...}
!               The meanings change if rdbkt_, rdbrc_ or rdprn_ are
!               false.  If rdbkt_ is false, the meanings are
!                      'tupl'  for tuple data           (item,...)
!                      'list'  for list or table data   {item,...}
!               If rdprn_ is false, then 'list' is used for all
!               container types.  If depth_ is 0, then ttype_ is not
!               valid and will contain '    '
!
!               depth_ is an integer variable with the depth into a
!               list, array, tuple or table.  A depth of zero means that
!               no list, array, tuple or table is being processed.
!
!               index_ is an integer variable with the index (from 1)
!               across a list, array, tuple or table.  An index of zero
!               means that no list, array, tuple or table is being processed.
!
!               dictype_ is a character*(NUMCHAR) variable with the type code
!               given in the dictionary entry for the named data item.  If
!               no dictionary was used, or no type code was specified, this
!               field will simply agree with type_.  If a dictionary was used,
!               this type may be more specific than the one given by type_.
!
!               diccat_ is a character*(NUMCHAR) variable with the category
!               of the named data item, or '(none)'
!
!               dicname_ is a character*(NUMCHAR) variable with the name of
!               the data item which is found in the dictionary for the
!               named data item.  If alias_ is .true., this name may
!               differ from the name given in the call to test_.  If alias_
!               is .false. or no preferred alias is found, dicname_ agrees with
!               the data item name.
!
!               tagname_ is a character*(NUMCHAR) variable with the name
!               of the data item as found in the input CIF.  It will be
!               blank if the data item name requested is not found in the
!               input CIF and may differ from the data item name provided
!               by the user if the name used in the input CIF is an
!               alias of the data item name and alias_ is .true.
!
!               posnam_, posval_, posend_  and posdec_ are integer variables
!               which may be examined if information about the horizontal
!               position of the name and data read are needed.  posnam_ is the
!               starting column of the data name found (most often 1).
!               posval_ is the starting column of the data value.  If the
!               field is numeric, then posdec_ will contain the effective
!               column number of the decimal point.  For whole numbers, the
!               effective position of the decimal point is one column to the
!               right of the field.  posend_ contains the ending column of the
!               data value.
!
!               quote_ is a character*3 variable which may be examined to
!               determine if a quotation character was used on character data.]
!
!               <data name>           Name of the data item to be tested.
!
!
!    dtype_     Return the dictionary type of a data name, if any.
!               [logical function returned as .true. if the item has a type
!               in the dctionary, .false. if not.  The type returned is
!               one of the base type used by type_ (see above), if possible]
!
!               <data name>          Name of the item for which a type is needed
!               <data type>          Returned type from the dictionary
!
!
!    name_      Get the NEXT data name in the current data block.
!               [logical function returned as .true. if a new data name exists
!               in the current data block, and .false. when the end of the data
!               block is reached.]
!
!               <data name>           Returned name of next data item in block.
!
!
!
!    numb_      Extracts the number and its standard deviation (if appended).
!               [logical function returned as .true. if number present. If
!               .false. arguments 2 and 3 are unaltered. If the esd is not
!               attached to the number argument 3 is unaltered.]
!
!               <data name>           Name of the number sought.
!
!               <real variable>       Returned number.
!
!               <real variable>       Returned standard deviation.
!
!
!
!    numd_      Extracts the number and its standard deviation (if appended)
!               as double precision variables.
!               [logical function returned as .true. if number present. If
!               .false. arguments 2 and 3 are unaltered. If the esd is not
!               attached to the number argument 3 is unaltered.]
!
!               <data name>           Name of the number sought.
!
!               <double precision variable>
!                                     Returned number.
!
!               <double precision variable>
!                                     Returned standard deviation.
!
!
!
!    char_      Extracts character and text strings.
!               [logical function returned as .true. if the string is present.
!               Note that if the character string is text this function is
!               called repeatedly until the logical variable text_ is .false.
!               Non-text blank (quoted blanks) or empty ('' or "") fields
!               are converted by char to a null field, if nblank_ is true.]
!
!               <data name>           Name of the string sought.
!
!               <character variable>  Returned string is of length long_.
!
!    charnp_    Extracts character and text strings.
!               [logical function returned as .true. if the string is present.
!               Note that if the character string is text this function is
!               called repeatedly until the logical variable text_ is .false.
!               If the value is found in a container, then charnp_ should
!               be called repeatedly until both text_ is false and depth_
!               is zero. 
!
!               Non-text blank (quoted blanks) or empty ('' or "") fields
!               are converted by char to a null field, if nblank_ is true.
!               Only the number of characters returned in the third argument
!               are set.  This value is never less than 1, but may be less
!               than the allocated length of the returned string.]
!
!               <data name>           Name of the string sought.
!
!               <character variable>  Returned string is of length long_.
!
!               <integer variable>    Returned length of valid characters.
!
!
!    cmnt_      Extracts the next comment from the data block.
!               [logical function returned as .true. if a comment is present.
!               The initial comment character "#" is _not_ included in the
!               returned string.  A completely blank line is treated as
!               a comment.  A comment may be extracted while reading a list,
!               array, tuple or table]
!
!               <character variable>  Returned string is of length long_.
!
!
!    delim_     Reports the most recently seen delimiter prior to the
!               most recently extracted tag or value at the specified
!               depth.  Outside of bracketed constructs, only delimiters
!               at depth 0 (top level) can be seen.  This is not the
!               quoting character for a quoted string or text field.
!               See the variable quote_.
!               [logical function returned as .true. if the depth is
!               not negative and greater than or equal to the current
!               depth.  At depth 0, in a correctly formatted CIF, the
!               delimiter returned is always a blank,]
!
!               <integer variable>    Depth
!                    
!               <character variable>  Returned string is of length 1
!
!               <integer variable>    column position of delimiter
!
!               <integer variable>    record position of delimiter
!
!
!
!    purge_     Closes existing data files and clears tables and pointers.
!               [subroutine call]
!
!____________________________________________________________________________
!
!
!
!   CIF CREATION TOOLS ("the put_ing commands")
!
!
!
!    pfile_     Create a file with the specified file name.
!               [logical function returned as .true. if the file is opened.
!               The value will be .false. if the file already exists.]
!
!               <file name>           Blank for use of currently open file
!
!
!
!    pdata_     Put a data block command into the created CIF.
!               [logical function returned as .true. if the block is created.
!               The value will be .false. if the block name already exists.
!               Produces a save frame instead of a data block if the
!               variable saveo_ is true during the call.  No block duplicate
!               check is made for a save frame.]
!
!               <block name>
!
!    pdelim_    Emit a specific delimiter
!               [logical function returned as .true. if the delimiter is
!               appropriate to the context.  Emitting a '(', '{' or '['
!               increases the output depth by one.  Emitting a ')', '}'
!               or ']' decreases the output depth by one.  Emitting a ' ',
!               ',' or ':' does not change the depth.  Emitting a ','
!               or ':' at depth_ 0 is an error that can be overridden
!               by the second argument being .true..  Emitting a ' ' at
!               a depth_ greater than 0 is an error that can be overridden
!               by the second argument being .true.. ]
!
!               <character variable>   The one-character delimiter string
!
!               <logical variable>     .true. if an invalid delimiter is
!                                      to be forced out
!
!               <integer variable>     Column position at which to write
!                                      the delimiter or 0 if not specified
!
!
!
!    ploop_     Put a loop_ data name into the created CIF.
!               [logical function returned as .true. if the invocation
!               conforms with the CIF logical structure.  If pposval_
!               is non-zero, the "loop_" header is positioned to
!               that column.  If pposnam_ is non-zero, the data name is
!               positioned to that column.]
!
!               <data name>         If the name is blank on the first call
!                                   of a loop, only the "loop_" is placed.
!
!
!
!    pchar_     Put a character string into the created CIF.
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary,
!               AND, if the invocation conforms to the CIF logical structure.
!               The action of pchar_ is modified by the variables pquote_ and
!               nblanko_.  If pquote_ is non-blank, it is used as a quotation
!               character for the string written by pchar_.  The valid values
!               are '''', '"', ';', '(', '{', '[', '''''''', and '"""'.
!               In the last six cases a text field, bracketed construct or
!               multi-line triple-quoted string is written.  If the string
!               contains a matching character to the value of quote_, or if
!               quote_ is not one of the valid quotation characters, a valid,
!               non-conflicting quotation character is used or the line-folding
!               conventions are used to prevent the close-quote from being
!               followed by white space.  Except when writing a text field, if
!               nblanko_ is true, pchar_ converts a blank string to
!               an unquoted period.]
!
!               <data name>         If the name is blank, do not output name.
!
!               <character string>  A character string of MAXBUF chars or less.
!
!
!
!    pcmnt_     Puts a comment into the created CIF.
!               [logical function returned as .true.  The comment character
!               "#" should not be included in the string.  A blank comment
!               is presented as a blank line without the leading "#"].
!
!               <character string>  A character string of MAXBUF chars or less.
!
!
!    pnumb_     Put a single precision number and its esd into the created CIF.
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary,
!               AND, if the invocation conforms to the CIF logical structure.
!               The number of esd digits is controlled by the variable
!               esdlim_]
!
!               <data name>         If the name is blank, do not output name.
!
!               <real variable>     Number to be inserted.
!
!               <real variable>     Esd number to be appended in parentheses.
!
!
!    pnumd_     Put a double precision number and its esd into the created CIF.
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary,
!               AND, if the invocation conforms to the CIF logical structure.
!               The number of esd digits is controlled by the variable
!               esdlim_]
!
!               <data name>         If the name is blank, do not output name.
!
!               <double precision variable>
!                                   Number to be inserted.
!
!               <double precision variable>
!                                   Esd number to be appended in parentheses.
!
!
!
!    ptext_     Put a character string into the created CIF.
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary,
!               AND, if the invocation conforms to the CIF logical structure.
!               ptext_ is invoked repeatedly until the text is finished. Only
!               the first invocation will insert a data name.
!
!               If used when pclipt_ is .true. if the first character of the
!               text field is blank, it is removed.
!
!               If used when pfold_ is non-zero, the text field will be marked
!               as folded even if the first line is small enough to fit.
!               In order to produce a non-folded text field in the midst
!               of generally folded items, pfold_ should be set to 0 before
!               calling ptext_ and then restored to the previous value.]
!
!               <data name>         If the name is blank, do not output name.
!
!               <character string>  A character string of MAXBUF chars or less.
!
!
!    prefx_     Puts a prefix onto subsequent lines of the created CIF.
!               [logical function returned as .true.  The second argument
!               may be zero to suppress a previously used prefix, or
!               greater than the non-blank length of the string to force
!               a left margin.  Any change in the length of the prefix
!               string flushes pending partial output lines, but does _not_
!               force completion of pending text blocks or loops.
!               This function allows the CIF output functions to be used
!               within what appear to be text fields to support annotation
!               of a CIF. ]
!
!               <character string>  A character string of MAXBUF chars or less.
!
!               <integer variable>  The length of the prefix string to use.
!
!
!
!
!    close_     Close the creation CIF. MUST be used if pfile_ is used.
!               [subroutine call]
!
!
!____________________________________________________________________________
!
!
!
!....The CIF tool box also provides variables for data access control:
!
!
!    alias_      Logical variable: if left .true. then all calls to
!                CIFtbx functions may use aliases of data item names.
!                The preferred synonym from the dictionary will be
!                subsituted internally, provided aliased data names were
!                supplied by an input dictionary (via dict_).  The
!                default is .true., but alias_ may be set to .false.
!                in an application.
!
!    aliaso_     Logical variable: if set .true. then cif output
!                routines will convert aliases to the names to preferred
!                synonyms from the dictionary.  The default is .false., but
!                aliaso_ may be set to .true. in an application.  The
!                setting of aliaso_ is independent of the setting of
!                alias_.
!
!    align_      Logical variable signals alignment of loop_ lists during
!                the creation of a CIF. The default is .true.
!
!    append_     Logical variable:  if set .true. each call to ocif_ will
!                append the information found to the current cif.  The default
!                is .false.
!
!    bloc_       Character*(NUMCHAR) variable: the current block name.
!
!    clipt_      Logical variable: if set .true., when reading text fields,
!                an extra blank is inserted before the character string
!                returned for the first line of a text field, emulating
!                the behavior of CIFtbx versions prior to version 4.
!
!    decp_       Logical variable: set when processing numeric input, .true.
!                if there is a decimal point in the numeric value, .false.
!                otherwise
!
!    depth_      Integer variable: set to the depth within a list, array, tuple
!                or table
!
!    dictype_    Character*(NUMCHAR) variable: the precise data type code
!                (see test_)
!
!    diccat_     Character*(NUMCHAR) variable: the category (see test_)
!
!    dicname_    Character*(NUMCHAR) variable: the root alias (see test_) or
!                the name of the dictionary just loaded (see dict_)
!
!    dicpname_   Character*(NUMCHAR) variable: the parent (see test_)
!
!    dicver_     Character*(NUMCHAR) variable: the version of the dictionary
!                just loaded (see dict_)
!
!    esdlim_     Integer variable:  Specifies the upper limit of esd's
!                produced by pnumb_, and, implicitly, the lower limit.
!                The default value is 19, which limits esd's to the range
!                2-19.  Typical values of esdlim_ might be 9 (limiting
!                esd's to the range 1-9), 19, or 29 (limiting esd's
!                to the range 3-29).  If esdlim_ is given as a negative
!                value, the upper limit of esd's is the absolute value
!                of esdlim_ and the lower limit is 1.
!
!    esddig_     Integer variable:  The number of esd digits in the last
!                number read from a CIF.  Will be zero if no esd
!                was given.
!
!    file_       Character*(MAXBUF) variable: the filename of the current file.
!                Warning:  only file_(1:longf_) is valid
!
!    fold_       Logical variable signals that the current text block
!                began with the ';\' fold indicator. Only meaningful
!                when text_ is .true. and type_ is 'text'.
!                (fold_ is .true. if the indicator is present)
!
!    glob_       Logical variable signals that the current data block
!                is actually a global block (.true. for a global block).
!
!    globo_      Logical variable signals that the output data block from
!                pdata_ is actually a global block (.true. for a global block).
!
!    index_      Integer variable: Specifies the one-based index of the current
!                item in a list, array, tuple or table
!
!    line_       Integer variable: Specifies the input/output line limit
!                for processing a CIF. The default value is 80 characters.
!                This may be set by the program. The max value is MAXBUF
!                which has a default value of 2048.  In order to use
!                the CIF 1.1 line folding protocol for lines that
!                cannot be fit into line_ characters, the variable
!                pfold_ must be set to a non-zero value less than
!                or equal to line_
!
!    list_       Integer variable: the loop block number (see test_).
!
!    long_       Integer variable: the length of the data string in strg_.
!
!    longf_      Integer variable: the length of the filename in file_.
!
!    loop_       Logical variable signals if another loop packet is present.
!
!    lzero_      Logical variable: set when processing numeric input, .true.
!                if the numeric value is of the form [sign]0.nnnn rather than
!                [sign].nnnn, .false. otherwise
!
!    nblank_     Logical variable: if set .true. then all calls to
!                to char_ or test_ which encounter a non-text quoted blank
!                will return the type as 'null' rather than 'char'.
!
!    nblanko_    Logical variable: if set .true. then cif output
!                routines will convert quoted blank strings to an
!                unquoted period (i.e. to a data item of type null).
!
!    pclipt_     Logical variable: if set .true., when writing text fields,
!                if there is a blank as the first character to be written,
!                it is removed, emulating the behavior of CIFtbx versions
!                prior to version 4.
!
!    pdecp_      Logical variable: if set .true. then cif numeric output
!                routines will insert a decimal point in all numbers written by
!                pnumb_ or pnumbd_.  If set .false. then a decimal point will be
!                written only when needed.  The default is .false.
!
!    pesddig_    Integer variable: if set non-zero, and esdlim_ is negative,
!                controls the number of digits for esd's produced by
!                pnumb_ and pnumd_
!
!    pfold_      Integer variable:  If set non-zero, specifies a column
!                on which output lines are to be folded.  The default is 0.
!                If pfold_ is set to a value greater than line_ the
!                value of line_ will be used instead.  Non-zero values of
!                pfold_ less than 4 are not valid and will be reset to 4.
!                Non-zero values of pfold_ less than 80 can cause conflict
!                with the syntactic requirements of creating a valid CIF.
!
!    plzero_     Logical variable: if set .true. then cif numeric output
!                routines will insert a zero before a leading decimal point,
!                The default is .false.
!
!    pposdec_    Integer variable giving the position of the decimal point
!                for the next number to be written.  This acts very much like
!                a decimal centered tab in a word processor, to help align
!                columns of number on a decimal point, if a decimal point
!                is present.
!
!    pposend_    Integer variable giving the ending column of the next
!                number or quoted character value to be written.  Used to
!                pad with zeros or blanks.
!
!    pposnam_    Integer variable giving the starting column of the next
!                name or comment or data block to be written.
!
!    pposval_    Integer variable giving the starting column of the next
!                data value to be written by pchar_, pnumb_ or pnumd_.
!                Also used to set the position of the initial "loop_"
!                in a ploop_ call or to set the position of a terminal "save_"
!                for a save frame in a pdata_ call for which saveo_ is .true.
!
!    posdec_     Integer variable giving the position of the decimal point
!                for the last number read, if a decimal point was present.
!
!    posend_     Integer variable giving the ending column of the last
!                data value read, not including a terminal quote.
!
!    posnam_     Integer variable giving the starting column of the last
!                name or comment or data block read.
!
!    posval_     Integer variable giving the starting column of the last
!                data value read.  Also reports the column of the
!                terminal "save_" of a save frame.
!
!    pquote_     Character variable giving the quotation symbol to be
!                used for the next string written, or the comment
!                flag for the next comment written.
!
!    precn_      Integer variable:  Reports the record number of the last
!                line written to the output cif.  Set to zero by init_.  Also
!                set to zero by pfile_ and close_ if the output cif file name
!                was not blank.
!
!    ptabx_      Logical variable signals tab character expansion to blanks
!                during the creation of a CIF. The default is .true.
!
!    quote_      Character variable giving the quotation symbol found
!                delimiting the last string read or the comment flag
!                for the last comment read.  The possible valid values
!                are '''', '"', ';', '''''''', and '"""'.
!                The treble quotes are recognized only if rdtq_ is .true.
!
!    rdbrc_      Logical variable:  control recognition of { ... } constructs
!                on read.  The default is .false.
!
!    rdbkt_      Logical variable:  controls recognition of [ ... ] constructs
!                on read.  The default is .false.
!
!    rdprn_      Logical variable:  controls recognition of ( ... ) constructs
!                on read.  The default is .false.
!
!    rdtq_       Logical variable:  controls recognition of """ ... """ and
!                ''' ... ''' constructs on read.  The default is .false.
!
!    rdrcqt_     Logical variable:  controls recognition of trailing punctuation
!                after a closing quote.  If .true. a closing quotation mark is
!                recognized immediately, no matter what follows the closing
!                quoation mark (the CIF 2 convention).  If .false., a closing
!                quotation mark is only effective if followed by a blank, or,
!                in bracketed constructs by a blank, a colon, a comma or 
!                the closing bracket.
!
!    recbeg_     Integer variable:  Gives the record number of the first
!                record to be used.  May be changed by the user to restrict
!                access to a CIF.
!
!    recend_     Integer variable:  Gives the record number of the last
!                record to be used.  May be changed by the user to restrict
!                access to a CIF.
!
!    recn_       Integer variable:  Reports the record number of the last
!                line read from the direct access copy of the input cif.
!
!    save_       Logical variable signals that the current data block
!                is actually a save-frame (.true. for a save-frame).
!
!    saveo_      Logical variable signals that the output data block from
!                pdata_ is actually a save-frame (.true. for a save-frame).
!
!    strg_       Character*(MAXBUF) variable: the current data item.
!
!    tabl_       Logical variable signals tab-stop alignment of output
!                during the creation of a CIF. The default is .true.
!
!    tabx_       Logical variable signals tab character expansion to blanks
!                during the reading of a CIF. The default is .true.
!
!    tbxver_     Character*32 variable: the CIFtbx version and date
!                in the form 'CIFtbx version N.N.N, DD MMM YY '
!
!    text_       Logical variable signals if another text line or is present.
!
!    type_       Character*4 variable: the data type code (see test_).
!
!    ttype_      Character*4 variable: the list, array, tuple or table type code (see test_).
!
!    unfold_     Logical variable signals that input lines are to be
!                unfolded before presentation of data.  The default
!                is .false.
!
!    xmlout_     Logical variable:  Set by the user to change the output
!                style to XML conventions.  Note that this is not a
!                cml output, but a literal translation from the input CIF.
!                The default is .false.
!
!    xmlong_     Logical variable:  Set by the user to change the style of
!                xml output if xmlout_ is .true.  When .true. (the default)
!                xml tag names are the full CIF tag names with the leading
!                '_' removed.  When .false. an attempt is made to strip
!                the leading category name as well.
!
!
!_____________________________________________________________________________
!
!
! >>>>>> Set the device numbers.
!
         function init_(devcif,devout,devdir,deverr)
!
         logical   init_
         include   'ciftbx.sys'
         integer   devcif,devout,devdir,deverr
         integer   ii,kdig
         real      ytest
         double precision ztest
         double precision tbxxdble
         real      tbxxsngl
!
         init_=.true.
         cifdev=devcif
         outdev=devout
         dirdev=devdir
         errdev=deverr

         recn_=0
         precn_=0
         plcat = ' '
         plxcat = ' '
         plhead(1) = ' '
         plxhead(1) = ' '
         pdblok = ' '
         ploopn = 0
         nstable = 0
         nivt = 0
!
!        recompute decimal single precision precision
!        This is found by computing the smallest power of
!        10 which, when added to 1, produces a change
!        and then backing off by 1
!
         decprc = .1
         do ii = 1,8
         ytest = tbxxsngl(1.+decprc/10.)
         if (ytest.eq.1.) go to 100
         decprc = decprc/10.
         enddo
100      continue
         decprc=decprc*10.
!
!        recompute decimal double precision precision
!
         kdig = 1
         dpprc = .1D0
         do ii = 1,17
         ztest = tbxxdble(1.D0+dpprc/10.)
         if (ztest.eq.1.D0) go to 200
         dpprc = dpprc/10.D0
         kdig = kdig+1
         enddo
200      continue
         dpprc=dpprc*10.D0
         write(ndpfmt,'(5h(d30.,i2,1h))') kdig-1
!
!        recompute decimal single precision minimum power of ten
!
         decmin = .1
         do ii = 1,39
         ytest = decmin/10.
         if (ytest.eq.0.) go to 300
         decmin = decmin/10.
         enddo
300      continue
!
!        recompute decimal double precision minimum power of 10
!        and its log base 10 (minexp)
!
         dpmin = .1D0
         minexp = -1
         do ii = 1,309
         ztest = dpmin/10.
         if (ztest.eq.0.D0) go to 400
         dpmin = dpmin/10.D0
         minexp = minexp-1
         enddo
400      continue
         call clearfp
         return
         end
!
!
! >>>>>> Function to defeat the optimizer
!
!     
         function tbxxdble(x)
         double precision x
         double precision tbxxdble
         tbxxdble = x
         return
         end
!
!
! >>>>>> Function to defeat the optimizer
!
!     
         function tbxxsngl(x)
         real x
         real tbxxsngl
         tbxxsngl = x
         return
         end
!
!
!
!
!
! >>>>>> Read a CIF dictionary and prepare for checks
!
         function dict_(fname,checks)
!
         logical   dict_
         logical   ocif_
         logical   data_
         logical   charnp_
         logical   test_
         integer   lastnb
         include  'ciftbx.sys'
         logical   tbxxnewd, tbxxoldd
         logical   nresult
         character fname*(*),checks*(*)
         character temp*(MAXBUF)
         character codes(11)*5,name*(MAXBUF),bxname*(NUMCHAR)
         character bpname*(NUMCHAR)
         character bcname*(NUMCHAR),bname*(NUMCHAR)
         character baname*(NUMCHAR),ganame*(NUMCHAR),btname*(NUMCHAR)
         character batag*(NUMCHAR)
         character mcstrg*(NUMCHAR)
         character riname*(NUMCHAR),rfname*(NUMCHAR)
         character xdicnam*(NUMCHAR)
         character xdicver*(NUMCHAR)
         character xmtoken*(NUMCHAR),xmtarg*(XMLCHAR),xmtyp*(NUMCHAR)
         character xxxtemp*(NUMCHAR)
         character*3 ovchk, otchk
         integer   nrecds,recends,recbegs
         integer   lchecks,lbpname,lbcname,lbaname,lbname
         integer   kdict,ifind,jfind,iafind,ick
         integer   i,j,nmatch,mycat,ksmatch,ii,jj,idstrt,kdup
         integer   nmycat,ixmtyp,nxmc,kxmc
         integer   lstrg,lxmtoken,lxmtarg,lxmtyp,kvrtp,kstrg,sindex
         integer   lbloc,kivt

!
!        Control flags for matching categories, names and types
!
!        icloop is the loop number of the block for the
!        current category
!        ictype is the type of the current category
!          0 - none found yet
!          1 - _item.category_id             (DDL2)
!          2 - _category                     (DDL1)
!          3 - _category.id                  (DDL2)
!          4 - _name.category_id             (DDLm)
!        the last ictype entry is not a type, but a tag
!        whose value may specify that this is a category
!        with the category name given under intype
!          5 - _definition.scope             (DDLm)
!        inloop is the loop number of the block for the
!        current name
!        intype is the type of the current name
!          0 - none found yet
!          1 - _item.name                    (DDL2)
!          2 - _name                         (DDL1)
!          3 - _definition.id                (DDLm)
!        ialoop is the loop number of the block for the
!        current alias
!        iatype is the type for the current alias
!          0 - none found yet
!          1 - _item_aliases.alias_name      (DDL2)
!          2 - _aliases.definition_id        (DDL2)
!        imloop is the loop number of the block for the
!        current parent
!        imtype is the type for a mandatory item
!          0 - none found yet
!          1 - _item.mandatory_code          (DDL2)
!          2 - _category_mandatory.item_id   (DDLm)
!        iptype is the type for the current parent
!          0 - none found yet
!          1 - _item_linked.parent_name      (DDL2)
!          2 - _item_link_parent             (DDL1)
!          3 - _category.parent_id           (DDLm)
!          4 - _name.linked_item_id          (DDLm)
!        itloop is the loop number of the block for the
!        current type
!        ittype is the type of the current type
!          0 - none found yet
!          1 - _item_type.code               (DDL2)
!          2 - _type                         (DDL1)
!          3 - _type.contents                (DDLm)
!        iritype is the type of the current related item
!          0 - none found yet
!          1 - _item_related.related_name    (DDL2)
!          2 - _related_item                 (DDL1)
!          3 - _type.purpose                 (DDLm)
!        irftype is the type of the current related item function
!          0 - none found yet
!          1 - _item_related.function_code   (DDL2)
!          2 - _related_function             (DDL1)
!          3 - _type.purpose                 (DDLm)
!

         integer icloop,ictype,inloop,intype,ialoop,iatype,             &
     & imloop,imtype,iptype,itloop,ittype,                              &
     & iritype,irftype,icktype
!
         character*4 map_type(19),map_to(19),mapped
         character*(NUMCHAR) dt(2),dv(2),ct(5),nt(3),at(2),tt(3)
         character*(NUMCHAR) ri(3),rf(3),ck(4),pt(4),pc(2),mc(3)
         character*(NUMCHAR) ve(3),vr(4)
         data map_type                                                  &
     &   /'floa','int ','yyyy','symo','ucha','ucod','name','idna',      &
     &    'any ','code','line','ulin','atco','fax ','phon','emai',      &
     &    'real','inte','coun'/
         data map_to                                                    &
     &   /'numb','numb','char','char','char','char','char','char',      &
     &    'char','char','char','char','char','char','char','char',      &
     &    'numb','numb','numb'/
         data ri                                                        &
     &      /'_item_related.related_name      ',                        &
     &       '_related_item                   ',                        &
     &       '_type.purpose                   '/
         data rf                                                        &
     &      /'_item_related.function_code     ',                        &
     &       '_related_function               ',                        &
     &       '_type.purpose                   '/
         data dt                                                        &
     &      /'_dictionary.title               ',                        &
     &       '_dictionary_name                '/
         data dv                                                        &
     &      /'_dictionary.version             ',                        &
     &       '_dictionary_version             '/
         data ct                                                        &
     &      /'_item.category_id               ',                        &
     &       '_category                       ',                        &
     &       '_category.id                    ',                        &
     &       '_name.category_id               ',                        &
     &       '_definition.scope               '/
         data nt                                                        &
     &      /'_item.name                      ',                        &
     &       '_name                           ',                        &
     &       '_definition.id                  '/
         data at                                                        &
     &      /'_item_aliases.alias_name        ',                        &
     &       '_aliases.definition_id          '/
         data tt                                                        &
     &      /'_item_type.code                 ',                        &
     &       '_type                           ',                        &
     &       '_type.contents                  '/
         data ck                                                        &
     &      /'_category_key.name              ',                        &
     &       '_list_reference                 ',                        &
     &       '_category_key.generic           ',                        &
     &       '_category_key.primitive         '/
         data pt                                                        &
     &      /'_item_linked.parent_name        ',                        &
     &       '_item_link_parent               ',                        &
     &       '_category.parent_id             ',                        &
     &       '_name.linked_item_id            '/
         data pc                                                        &
     &      /'_item_linked.child_name         ',                        &
     &       '_item_link_child                '/
         data mc                                                        &
     &      /'_item.mandatory_code            ',                        &
     &       '_mandatory                      ',                        &
     &       '_category_mandatory.item_id     '/
         data ve                                                        &
     &      /'_item_enumeration.value         ',                        &
     &       '_enumeration                    ',                        &
     &       '_enumeration_set.state          '/
         data vr                                                        &
     &      /'_item_range.minimum             ',                        &
     &       '_enumeration_range              ',                        &
     &       '_item_range.maximum             ',                        &
     &       '_enumeration.range              '/

!
         data codes /'valid','dtype','reset','close',                   &
     &       'catck','catno','nodup','final','first',                   &
     &       'parck','parno'/
!
         nrecds=nrecd
         recbegs=recbeg_
         recends=recend_
         if(append_) then
           recbeg_=nrecd
         endif
!
!        Initialize kdup to 0 ('final')
!
         kdup = 0
!
!        initialize both xdicnam and xdicver to blank
!
         xdicnam = ' '
         xdicver = ' '
!
!        preserve entry values of tcheck and vcheck in case dict fails
!
         otchk = tcheck
         ovchk = vcheck
!
!....... Are the codes OK
!
         lchecks=min(len(temp),len(checks))
         call tbxxnlc(temp(1:lchecks),checks)
         i=0
120      i=i+1
         if(i.ge.lchecks)            goto 190
         if(temp(i:i).eq.' ')        goto 120
         do 150 j=1,11
         if(temp(i:i+4).eq.codes(j)) goto 170
150      continue
         dict_=.false.
         goto 500
170      i=i+4
         if(j.eq.1) then
           vcheck='yes'
           goto 120
         endif
         if(j.eq.2) then
           tcheck='yes'
           goto 120
         endif
         if(j.eq.3) then
           vcheck = 'no '
           tcheck = 'no '
           goto 120
         endif
         if(j.eq.4) then
           vcheck = 'no '
           tcheck = 'no '
           catchk = 'yes'
           ndcname = 0
           ndict = 0
           if(nname.gt.0) then
           do 180 i = 1,nname
             dtype(i)=' '
             dxtyp(i)=' '
             cindex(i)=0
             ddict(i)=0
180        continue
           endif
           dict_=.true.
           goto 500
         endif
         if (j.eq.5) then
           catchk = 'yes'
           goto 120
         endif
         if (j.eq.6) then
           catchk = 'no '
           goto 120
         endif
         if (j.eq.10) then
           parchk = 'yes'
           goto 120
         endif
         if (j.eq.11) then
           parchk = 'no '
           goto 120
         endif
         kdup=j-8
         goto 120
!
!        if no category names have been loaded, clean up
!        the hash table for dictionary category names
!
190      if(ndcname.eq.0) then
           call hash_init(dcname,dcchain,NUMDICT,ndcname,dchash,        &
     &     NUMHASH)
         endif
!
!        if no dictionary names have been loaded, clean up
!        the hash table for dictionary names
!
         if(ndict.eq.0) then
           call hash_init(dicnam,dicchain,NUMDICT,ndict,dichash,        &
     &     NUMHASH)
         endif
         idstrt=ndict
!
!....... Open and store the dictionary
!
         dict_=.true.
         if(fname.eq.' ')            goto 500
         if(nname.gt.0) call tbxxerr(' Dict_ must precede ocif_')
         dict_=ocif_(fname)
         if(.not.dict_)              goto 500
         dictfl='yes'
!
!        At this point is is proper to update xdicnam to fname
!
         xdicnam = fname
!
!....... Loop over data blocks; extract _name's, _type etc.
!
200      if(.not.data_(' '))         goto 400
         lbloc = lastnb(bloc_)
         if(bloc_(1:1).eq.'_'.or.glob_.or.bloc_.eq.' ') then
           call tbxxclc(bname,lbname,bloc_(1:lbloc),lbloc)
         else
           call tbxxclc(bname,lbname,'_'//bloc_(1:lbloc),lbloc+1)
         endif
!
!        see if this is a dictionary defining block
!
         do i = 1,2
           if(charnp_(dt(i),name,lstrg)) then
             xdicnam = name(1:lstrg)
             do j = 1,2
               if(test_(dv(j))) then
                 xdicver = strg_(1:max(1,long_))
                 goto 200
               endif
             enddo
             goto 200
           endif
         enddo
!
!dbg     WRITE(6,*) ndict,bloc_
!
!        Analyze loop structure for categories, names, types and parents
!
!
!        initalize loop info
!
         icloop = -1
         inloop = -1
         ialoop = -1
         imloop = -1
         itloop = -1
         ictype = 0
         intype = 0
         iatype = 0
         imtype = 0
         iptype = 0
         ittype = 0
         iritype = 0
         irftype = 0
         icktype = 0
         ixmtyp = 0
         bcname = ' '
         bpname = ' '
         lbcname = 1
         lbpname = 1
         baname = ' '
         batag = ' '
         lbaname = 1
         btname = ' '
         mycat=0
         loop_=.false.
         loopnl=0
         nmatch=0
         ksmatch=0
         riname = ' '
         rfname = ' '
!
!        Pick up category_keys and list_references
!
         do i = 1,4
210        if(charnp_(ck(i),name,lstrg)) then
             if (icktype.ne.0 .and. icktype.ne.i)                       &
     &         call tbxxwarn                                            &
     &         (' Multiple DDL 1, 2 or m related key definitions ')
             icktype = i
             if (tbxxnewd(name(1:lstrg),ick)) then
               catkey(ick) = .true.
             else
               if(.not.catkey(ick)) then
                 ifind = aroot(ick)
215              catkey(ifind) = .true.
                 ifind = alias(ifind)
                 if (ifind.ne.0) go to 215
               endif
             endif
             if (loop_) go to 210
           endif
         enddo

!
!        Process related items
!
         do i = 1,2
           if(charnp_(ri(i),name,lstrg)) then
             if (iritype.ne.0)                                          &
     &         call tbxxwarn                                            &
     &         (' Multiple DDL 1 and 2 related item definitions ')
             iritype = i
             call tbxxnlc(riname,name(1:lstrg))
!
!            Seek the matching function, may be in the same loop or not
!
             if(charnp_(rf(i),name,lstrg)) then
               if (irftype.ne.0)                                        &
     &           call tbxxwarn                                          &
     &           (' Multiple DDL 1 and 2 related item functions ')
               irftype = i
               call tbxxnlc(rfname,name(1:lstrg))
             endif
           endif
         enddo
         loop_ = .false.
         loopnl = 0
!
!        Process categories
!
         do i = 1,5
           if(charnp_(ct(i),name,lstrg)) then
             if(i.eq.5) then
!
!              if this is a DDLm _defintion.scope with a value of
!              category, we need to get the name from _defintion.id
!
               call tbxxnlc(bcname,name(1:lstrg))
               if(bcname.eq.'category') then
                 if(.not.charnp_(nt(3),name,lstrg)) then
                   call tbxxwarn(                                       &
     &             ' DDLm category defintion without _definition.id ')
                 else
                   go to 216
                 endif
               endif
             endif

             if(ictype.ne.0)                                            &
     &         call tbxxwarn(                                           &
     &           ' Multiple DDL 1, 2 or m category definitions ')
             ictype = i
             if(loop_) icloop = loopnl
             call tbxxnlc(bcname,name(1:lstrg))
             lbcname=long_
             nmycat = ndcname+1
             call hash_store(bcname(1:long_),                           &
     &         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
             if(mycat.eq.0) then
               call tbxxerr(' Dictionary category names > NUMDICT ')
             endif
             if (mycat.eq.nmycat) then
               ccatkey(mycat) = 0
               xmcind(mycat)=0
             endif
!
!            if this is not a loop of categories, we expect a match
!            against the block name, unless we are doing replacements
!
             if(.not.loop_) then
               if(ictype.eq.1) then
                 if(bname(1:min(lbname,lbcname+2)).ne.                  &
     &            '_'//bcname(1:lbcname)//'.'                           &
     &            .and. catchk.eq.'yes'                                 &
     &            .and. (rfname(1:7).ne.'replace')) then
                 call tbxxwarn(' Category id does not match block name')
                 endif
               else
                 if(ictype.eq.2) then
                   if(bcname.ne.'dictionary_definition' .and.           &
     &                bcname.ne.'category_overview') then
                   if(bname(1:min(lbname,lbcname+2)).ne.                &
     &               '_'//bcname(1:lbcname)//'_') then
                   if(bname(1:min(lbname,lbcname+1)).ne.                &
     &               '_'//bcname(1:lbcname)                             &
     &            .and. catchk.eq.'yes'                                 &
     &            .and. (rfname(1:7).ne.'replace')) then
                   call tbxxwarn(                                       &
     &               ' Category id does not match block name')
                   endif
                   endif
                   endif
                 endif
               endif
             endif
           endif
           loop_ = .false.
           loopnl = 0
         enddo
!
!        Process XML translations
!
216      loop_ = .false.
         loopnl = 0
         if(charnp_('_xml_mapping.token',xmtoken,lxmtoken)) then
230        if(charnp_('_xml_mapping.token_type',xmtyp,lxmtyp)) then
             if(charnp_('_xml_mapping.target',xmtarg,lxmtarg)) then
               if (xmnxlat.ge.XMLDEFS) then
                 call tbxxerr(' XML translations > XMLDEFS')
               else
                 xmnxlat=xmnxlat+1
                 xmlate(xmnxlat)=xmtarg(1:lxmtarg)
               endif
               if (xmtyp.eq.'data') then
                 ixmtyp = 1
                 if (xmdata.eq.0) then
                   xmdata = xmnxlat
                 else
                   call tbxxwarn(' XML duplicate DATA_ translation')
                 endif
               endif
               if (xmtyp(1:lxmtyp).eq.'category') then
                 ixmtyp = 2
                 nxmc = ndcname+1
                 call tbxxnlc(xxxtemp,xmtoken(1:lxmtoken))
                 call hash_store(xxxtemp,                               &
     &           dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,kxmc)
                 if( kxmc.eq.nxmc) then
                   ccatkey(kxmc) = 0
                   xmcind(kxmc) = xmnxlat
                 else
                   if (xmcind(kxmc).ne.0) then
                     call tbxxwarn(                                     &
     &                 ' XML duplicate category translation')
                   else
                     xmcind(kxmc) = xmnxlat
                   endif
                 endif
               endif
               if (xmtyp.eq.'item') then
                 ixmtyp = 3
                 if (tbxxnewd(xmtoken(1:lxmtoken),ifind)) then
                   xmindex(ifind) = xmnxlat
                 else
                   if (xmindex(ifind).ne.0) then
                     call tbxxwarn(' XML duplicate item translation')
                   else
                     ifind = aroot(ifind)
 235                 xmindex(ifind) = xmnxlat
                     ifind = alias(ifind)
                     if (ifind.ne.0) go to 235
                   endif
                 endif
               endif
               if(loop_) then
                 if(charnp_('_xml_mapping.token',xmtoken,lxmtoken)) then
                   go to 230
                 else
                   call tbxxerr(' XML dictionary logic error')
                 endif
               endif
             else
               call tbxxerr(' XML target missing')
             endif
           else
             call tbxxerr(' XML token_type missing')
           endif
         else
           xmtoken = bname(1:lbname)
           lxmtoken=lbname
           if(charnp_('_xml_mapping.token_type',xmtyp,lxmtyp)) then
             if(charnp_('_xml_mapping.target',xmtarg,lxmtarg)) then
               if (xmnxlat.ge.XMLDEFS) then
                 call tbxxerr(' XML translations > XMLDEFS')
               else
                 xmnxlat=xmnxlat+1
                 xmlate(xmnxlat)=xmtarg(1:lxmtarg)
               endif
               if (xmtyp(1:lxmtyp).eq.'data') then
                 ixmtyp = 1
                 if (xmdata.eq.0) then
                   xmdata = xmnxlat
                 else
                   call tbxxwarn(' XML duplicate DATA_ translation')
                 endif
               endif
               if (xmtyp.eq.'category') then
                 ixmtyp = 2
                 nxmc = ndcname+1
                 call tbxxnlc(xxxtemp,xmtoken(1:lxmtoken))
                 call hash_store(xxxtemp,                               &
     &           dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,kxmc)
                 if( kxmc.eq.nxmc) then
                   ccatkey(kxmc) = 0
                   xmcind(kxmc) = xmnxlat
                 else
                   if (xmcind(kxmc).ne.0) then
                     call tbxxwarn(                                     &
     &                 ' XML duplicate category translation')
                   else
                     xmcind(kxmc) = xmnxlat
                   endif
                 endif
               endif
               if (xmtyp.eq.'item') then
                 ixmtyp = 3
                 if (tbxxnewd(xmtoken(1:lxmtoken),ifind)) then
                   xmindex(ifind) = xmnxlat
                 else
                   if (xmindex(ifind).ne.0) then
                     call tbxxwarn(' XML duplicate item translation')
                   else
                     ifind = aroot(ifind)
 240                 xmindex(ifind) = xmnxlat
                     ifind = alias(ifind)
                     if (ifind.ne.0) go to 240
                     xmindex(ifind) = xmnxlat
                   endif
                 endif
               endif
               if(loop_) then
                 call tbxxerr(' XML dictionary logic error')
               endif
             else
               call tbxxerr(' XML target missing')
             endif
           endif
         endif
!
!        Process names
!
         bxname = ' '
         do i = 1,3
         if(charnp_(nt(i),name,lstrg)) then
           if(intype.ne.0)                                              &
     &       call tbxxwarn(                                             &
     &         ' Multiple DDL 1 and 2 or m name definitions ')
           intype = i
           call tbxxnlc(bxname,name(1:lstrg))
           if(loop_) inloop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         enddo
         if(intype.eq.0.and.ictype.lt.3.and.(.not.glob_)                &
     &     .and.bname(1:lbname).ne.' '.and.ixmtyp.eq.0)                 &
     &     call tbxxwarn (' No name defined in block')
         loop_ = .false.
         if(charnp_(at(1),name,lstrg)) then
           iatype=1
           call tbxxnlc(baname,name(1:lstrg))
           batag = name(1:lstrg)
           lbaname = lstrg
           if(loop_) ialoop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         mcstrg = "no"
         if(ictype.ne.3) then
           do i=1,3
             if(charnp_(tt(i),name,lstrg)) then
               if(ittype.ne.0)                                          &
     &           call tbxxwarn(                                         &
     &             ' Multiple DDL 1 and 2 type definitions ')
               ittype = i
               call tbxxnlc(btname,name(1:lstrg))
               if(loop_) itloop = loopnl
             endif
             loop_ = .false.
             loopnl=0
           enddo
           do i = 1,2
             if(charnp_(mc(i),name,lstrg)) then
               if (imtype.ne.0)                                         &
     &           call tbxxwarn(' Multiple DDL 1 and 2 mandatory codes ')
               imtype = i
               call tbxxnlc(mcstrg,name(1:lstrg))
               if (loop_) imloop = loopnl
             endif
             loop_ = .false.
             loopnl=0
           enddo
         endif
!
!        Now test for consistent combinations
!
         if(inloop.ne.-1) then
           if(icloop.ne.-1.and.icloop.ne.inloop                         &
     &            .and. catchk.eq.'yes')                                &
     &       call tbxxwarn(                                             &
     &       ' Categories and names in different loops')
           if(iatype.ne.0.and.ialoop.ne.inloop) then
             if(ialoop.eq.-1) then
               if(bxname.ne.bname(1:lbname))                            &
     &          call tbxxwarn(                                          &
     &         ' One alias, looped names, linking to first')
             else
               call tbxxwarn(                                           &
     &         ' Aliases and names in different loops '                 &
     &         //' only using first alias ')
             endif
           endif
           if(itloop.ne.-1.and.itloop.ne.inloop)                        &
     &       call tbxxwarn(                                             &
     &       ' Types and names in different loops')
           if(imloop.ne.-1.and.imloop.ne.inloop)                        &
     &       call tbxxwarn(                                             &
     &       ' Mandatory codes and names in different loops')
         else
           if(icloop.ne.-1)                                             &
     &       call tbxxwarn(                                             &
     &         ' Multiple categories for one name')
           if(itloop.ne.-1)                                             &
     &       call tbxxwarn(                                             &
     &         ' Multiple types for one name')
           if(imloop.ne.-1)                                             &
     &       call tbxxwarn(                                             &
     &         ' Multiple madatory codes for one name')
         endif
!
!        Pick up parents
!
         do i = 1,2
220        if(charnp_(pt(i),name,lstrg)) then
             if (iptype.ne.0 .and. iptype.ne.i)                         &
     &         call tbxxwarn                                            &
     &         (' Multiple DDL 1 and 2 parent definitions ')
             iptype = i
             call tbxxnlc(bpname,name(1:lstrg))
             lbpname=long_
!
!            Seek the matching child, may be in the same loop or not
!
             if (charnp_(pc(i),name,lstrg)) then
               nresult = tbxxnewd(name(1:lstrg),ifind)
               nresult = tbxxnewd(bpname(1:lbpname),dpindex(ifind))
               bpname = ' '
               lbpname = 1
             endif
             if (loop_) go to 220
           endif
         enddo

!
!        Now we need to process value enumerations and ranges
!        and load them into item value table
!
         if (tcheck .eq. 'yes' .and. bxname.ne.' ') then
         loop_ = .false.
         nresult = tbxxnewd(bxname,ifind)

         do i = 1,2
5400       if(charnp_(ve(i),name,lstrg) .and. nivt.lt.NUMIVALS) then
             call tbxxsstb(name(1:lstrg),sindex)
             if (sindex.gt.0) then
               if (deindex(ifind).eq.0) then
                 deindex(ifind)=nivt+1
               else
                 kivt = deindex(ifind)
5410             if (ivtnxt(kivt).ne.0) then
                   kivt = ivtnxt(kivt)
                   go to 5410
                 endif
                 ivtnxt(kivt)=nivt+1
               endif
               nivt = nivt+1
               ivtnxt(nivt)=0
               ivtvet(nivt)=0
               ivtsbp(nivt)=sindex
             endif
           endif
           if (loop_) go to 5400
         enddo


         do i = 1,2
         loop_ = .false.
5420     strg_=' '
         long_=1
         nresult = test_(vr(i))
         if (strg_(1:long_).ne.' '.and.type_.eq.'null')                 &
     &     nresult = .true.
         if (nresult .and. nivt.lt.NUMIVALS) then
           nresult = charnp_(vr(i),name,lstrg)
           if (type_.ne.'char'.and.type_.ne.'numb') then
             name = '.'
             lstrg = 1
           endif
           kvrtp = -1
           if(i.eq.1 .and. lstrg<len(name)-2) then
             strg_=' '
             long_=1
             nresult = test_(vr(3))
             if (strg_(1:long_).ne.' '.and.                             &
     &         type_.eq.'null') nresult = .true.
             if (nresult) then
               nresult = charnp_(vr(3),                                 &
     &           name(lstrg+2:len(name)),kstrg)
               if (type_.ne.'char'.and.type_.ne.'numb') then
                 name(lstrg+2:len(name)) = '.'
                 kstrg = 1
               endif
               if (name(1:lstrg).ne.name(lstrg+2:lstrg+1+kstrg))        &
     &           kvrtp = 1
               name(lstrg+1:lstrg+1)=':'
               lstrg = lstrg+kstrg+1
               endif
             endif
             if (name(1:lstrg).eq.'.:.') then
               loop_=.false.
             else
               call tbxxsstb(name(1:lstrg),sindex)
               if (sindex.gt.0) then
                 if (deindex(ifind).eq.0) then
                   deindex(ifind)=nivt+1
                 else
                   kivt = deindex(ifind)
5430               if (ivtnxt(kivt).ne.0) then
                     kivt = ivtnxt(kivt)
                     go to 5430
                   endif
                   ivtnxt(kivt)=nivt+1
                 endif
                 nivt = nivt+1
                 ivtnxt(nivt)=0
                 ivtvet(nivt)=kvrtp
                 ivtsbp(nivt)=sindex
               endif
               if(loop_) go to 5420
             endif
           endif
         enddo
         endif

!
!        This is the main loop
!
250      if(ictype.eq.5.or.intype.eq.0) goto 200
         if(.not.charnp_(nt(intype),name,lstrg)) goto 200
         kdict=ndict+1
251      nresult = tbxxnewd(name(1:lstrg),ifind)
         if (bpname .ne. ' ') then
           nresult=tbxxnewd(bpname(1:lbpname),dpindex(ifind))
           bpname = ' '
           lbpname = 1
         endif
         nresult = tbxxoldd(bxname,jfind)
         if (nresult.and.jfind.ne.ifind.and.deindex(ifind).eq.0)        &
     &     deindex(ifind) = deindex(jfind)
         if(ifind.le.idstrt) then
           if (kdup .lt. 0) then
             call tbxxerr(' Duplicate name in dictionary '//            &
     &       dictag(ifind)(1:lastnb(dictag(ifind))))
           endif
           if (kdup .gt.0) go to 254
           dicnam(ifind)=char(0)
           goto 251
254        continue
         endif
         if(dicnam(ifind).eq.bname(1:lbname)) nmatch=ifind
         if(dicnam(ifind)(1:lbname).eq.bname(1:lbname)) ksmatch=ifind
!dbg     if(dicnam(ifind).ne.bname(1:lbname))
!dbg *   call tbxxwarn (' Name mismatch: '//dicnam(ifind)//bname(1:lbname))
         if(inloop.ge.0)then
!
!          We are in a loop of names.  If it is the same loop as
!          for categories, we need to extract the matching category
!
           if(inloop.eq.icloop) then
             mycat=0
             if(charnp_(ct(ictype),name,lstrg)) then
               call tbxxnlc(bcname,name(1:lstrg))
               lbcname=lstrg
               nmycat=ndcname+1
               call hash_store(bcname,                                  &
     &         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call tbxxerr(' Dictionary category names > NUMDICT ')
               endif
               if(mycat.eq.nmycat) ccatkey(mycat)=0
             endif
           endif
!
!          If it is the same loop as for types, we need to extract
!          the matching type
!
           if(inloop.eq.itloop) then
             btname=' '
             if(charnp_(ct(ittype),name,lstrg)) then
               call tbxxnlc(btname,name(1:lstrg))
             endif
           endif
!
!          If it is the same loop as for mandatory codes, we need to extract
!          the matching mandatory
!
           if(inloop.eq.imloop) then
             mcstrg='no'
             if(charnp_(mc(imtype),name,lstrg)) then
               call tbxxnlc(mcstrg,name(1:lstrg))
             endif
           endif
!
!          If it is the same loop as for aliases, we need to extract
!          the matching alias
!
           if(inloop.eq.ialoop) then
             baname=' '
             batag=' '
             if(charnp_(at(1),name,lstrg)) then
               call tbxxnlc(baname,name(1:lstrg))
               batag = name(1:lstrg)
               lbaname = lstrg
             endif
           endif
         endif
!
!        now we have a name stored in dicnam at location ifind
!        the index of the category in mycat, the type in btname,
!        the alias in baname, and the mandatory code in mcstrg
!
!        First verify match between the name and category, if
!        we have one, or extract from the block name
!
         if (mycat.eq.0) then
         if (dcindex(ifind).eq.0) then
           if (dicnam(ifind).eq.bloc_) then
             call tbxxcat(dicnam(ifind),bcname,lbcname)
!dbg         call tbxxwarn(' Extracting category name from block name '
!dbg *       //bloc_(1:max(1,lastnb(bloc_))))
             if(bcname(1:1).ne.' ') then
               ictype = 1
               nmycat = ndcname+1
               call hash_store(bcname,                                  &
     &         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call tbxxerr(' Dictionary category names > NUMDICT ')
               endif
               if (mycat.eq.nmycat) then
                 ccatkey(mycat) = 0
                 xmcind(mycat) = 0
               endif
             else
               if(catchk.eq.'yes')                                      &
     &         call tbxxwarn(' No category defined in block '           &
     &       //bloc_(1:max(1,lastnb(bloc_)))//' and name '              &
     &       //dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))            &
     &       //' does not match')
             endif
           endif
         endif
         else
         if (bcname(1:lbcname).ne.'dictionary_definition' .and.         &
     &     bcname(1:lbcname).ne.'category_overview') then
           if (dicnam(ifind)(1:lbcname+1).ne.'_'//bcname(1:lbcname)     &
     &        .or.( dicnam(ifind)(lbcname+2:lbcname+2).ne.'_' .and.     &
     &          dicnam(ifind)(lbcname+2:lbcname+2).ne.'.' .and.         &
     &          dicnam(ifind)(lbcname+2:lbcname+2).ne.' ' )) then
                if (catchk.eq.'yes'.and.rfname(1:7).ne.'replace')       &
     &          call tbxxwarn(' Item name '//                           &
     &          dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))//' '//    &
     &       ' does not match category name '//bcname(1:lbcname))
           endif
         endif
         endif
!
!        We will need the type in what follows.  cif_mm.dic defines
!        some higher level types.  We map them to primitive types
!
         mapped = btname(1:4)
         do i = 1,19
           if (btname(1:4).eq.map_type(i)) mapped = map_to(i)
         enddo
         if (mapped.ne.'char' .and.                                     &
     &       mapped.ne.'text' .and.                                     &
     &       mapped.ne.'null' .and.                                     &
     &       mapped.ne.'numb' .and.                                     &
     &       mapped.ne.'    ' ) then
             if (tcheck .eq. 'yes') then
               call tbxxwarn (' Item type '//                           &
     &           btname(1:max(1,lastnb(btname)))//' not recognized')
             endif
             mapped = 'char'
         endif

!
!        There are two cases to consider, one if the name is new to
!        the dictionary, the other, if it is not
!
         if(ifind.eq.kdict) then
           aroot(ifind)=ifind
           alias(ifind)=0
           dcindex(ifind)=mycat
           dictyp(ifind)=mapped
           dicxtyp(ifind)=btname
           dmcode(ifind) = 0
           if (mcstrg .eq. 'yes') dmcode(ifind) = 1
           if (mcstrg .eq. 'implicit') dmcode(ifind) = -1
         else
           if(dcindex(ifind).ne.mycat) then
             if(dcindex(ifind).eq.0) then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=aroot(ifind)
255            continue
               dcindex(jfind)=mycat
               jfind=alias(jfind)
               if(jfind.ne.0) goto 255
             else
               if(mycat.ne.0.and.                                       &
     &           (vcheck.eq.'yes'.or.tcheck.eq.'yes')                   &
     &           .and.catchk.eq.'yes')  then
                 if(rfname(1:7).ne.'replace')                           &
     &           call tbxxwarn(' Attempt to redefine category for item')
                 endif
             endif
           endif
           if(dictyp(ifind).ne.mapped .or.                              &
     &       dicxtyp(ifind).ne.btname) then
             if(dictyp(ifind).eq.' ') then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=aroot(ifind)
256            continue
               dictyp(jfind)=mapped
               dicxtyp(jfind)=btname
               jfind=alias(jfind)
               if(jfind.ne.0) go to 256
             else
               if(mapped.ne.' '.and.tcheck.eq.'yes')                    &
     &           call tbxxwarn(' Attempt to redefine type for item')
             endif
           endif
           if(dmcode(ifind).eq.0) then
             jfind = ifind
             if (aroot(ifind).ne.0) jfind = aroot(ifind)
257          continue
             dmcode(jfind) = 0
             if (mcstrg.eq.'yes') dmcode(jfind) = 1
             if (mcstrg.eq.'implicit') dmcode(jfind) = -1
             jfind=alias(jfind)
             if(jfind.ne.0) go to 257
           else
             if((mcstrg.eq.'yes' .and. dmcode(ifind).lt.0) .or.         &
     &         (mcstrg.eq.'implicit' .and. dmcode(ifind).gt.0))         &
     &         call tbxxwarn(                                           &
     &           ' Attempt to redefine mandatory code for item')
           endif
         endif
!
!        now deal with alias, if any.
!
         if(baname.ne.' ') then
           if (tbxxnewd(baname(1:lbaname),iafind)) then
             dictag(iafind)    =batag
             aroot(iafind)     =aroot(ifind)
             if(aroot(iafind).eq.0) aroot(iafind)=ifind
             catkey(iafind)    =catkey(ifind)
             alias(ifind)      =iafind
             dcindex(iafind)   =dcindex(ifind)
             dictyp(iafind)    =dictyp(ifind)
             dicxtyp(iafind)   =dicxtyp(ifind)
             xmindex(iafind)   =xmindex(ifind)
             dmcode(iafind)    =dmcode(ifind)
             dpindex(iafind)   =dpindex(ifind)
             deindex(iafind)   =deindex(ifind)
           else
             if(aroot(iafind).ne.0 .and.                                &
     &         aroot(iafind).ne.iafind) then
               if(aroot(iafind).eq.ifind .or.                           &
     &           aroot(iafind).eq.aroot(ifind)) then
                 call tbxxwarn(' Duplicate definition of same alias')
               else
                 call tbxxwarn(' Conflicting definition of alias')
               endif
             else
               if((dcindex(iafind).eq.0.or.                             &
     &           dcindex(iafind).eq.dcindex(ifind)).and.                &
     &           (dictyp(iafind).eq.' '.or.                             &
     &           (dictyp(iafind).eq.dictyp(ifind) .and.                 &
     &            dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
               endif
               if(xmindex(iafind).eq.0)                                 &
     &           xmindex(iafind)=xmindex(ifind)
               if(xmindex(ifind).eq.0)                                  &
     &           xmindex(ifind)=xmindex(iafind)
               if (dmcode(iafind).eq.0)                                 &
     &           dmcode(iafind)=dmcode(ifind)
               if (dmcode(ifind).eq.0)                                  &
     &           dmcode(ifind)=dmcode(iafind)
               if (dpindex(iafind).eq.iafind                            &
     &           .and. dpindex(ifind).ne.ifind)                         &
     &           dpindex(iafind) = dpindex(ifind)
               if (dpindex(ifind).eq.ifind                              &
     &           .and. dpindex(iafind).ne.iafind)                       &
     &           dpindex(ifind) = dpindex(iafind)
               if (deindex(ifind).eq.0)                                 &
     &           deindex(ifind)=deindex(iafind)
               if (deindex(iafind).eq.0)                                &
     &           deindex(iafind)=deindex(ifind)
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               alias(ifind)      =iafind
               if (catkey(iafind)) catkey(ifind) = .true.
               if (catkey(ifind)) catkey(iafind) = .true.
             endif
           endif
         endif
         if(inloop.ge.0) then
           baname = ' '
           batag = ' '
         endif
!
         if(inloop.ge.0.and.loop_) go to 250
         if(nmatch.eq.0) then
         if ((ksmatch.eq.0.or.inloop.lt.0)                              &
     &     .and.(rfname(1:7).ne.'replace')) then
         call tbxxwarn(' No name in the block matches the block name')
         endif
         endif
!
!        check for aliases
!        we execute this loop only in the case of unlooped name
!        with looped alias
!
         if(inloop.lt.0.and.ialoop.ge.0) then
           loop_=.false.
           loopnl=0
           ganame=baname
260        if(.not.charnp_(at(iatype),name,lstrg)) goto 200
           call tbxxnlc(baname,name(1:lstrg))
           batag=name(1:lstrg)
           lbaname=lstrg
           if(baname.eq.ganame) then
             if(loop_) go to 260
             go to 200
           endif
           if(baname.ne.' ') then
             if (tbxxnewd(baname(1:lbaname),iafind)) then
             if(iafind.eq.0) call tbxxerr(' CIFdic names > NUMDICT')
               dictag(iafind)    =batag
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               catkey(iafind)    =catkey(ifind)
               alias(ifind)      =iafind
               dcindex(iafind)   =dcindex(ifind)
               dictyp(iafind)    =dictyp(ifind)
               dicxtyp(iafind)   =dicxtyp(ifind)
               xmindex(iafind)   =xmindex(ifind)
               dmcode(iafind)    =dmcode(ifind)
               dpindex(iafind)   =dpindex(ifind)
               deindex(iafind)   =deindex(ifind)
               ifind=iafind
             else
               if(aroot(iafind).ne.0 .and.                              &
     &           aroot(iafind).ne.iafind) then
                 if(aroot(iafind).eq.ifind .or.                         &
     &             aroot(iafind).eq.aroot(ifind)) then
                   call tbxxwarn(' Duplicate definition of same alias')
                 else
                   call tbxxwarn(' Conflicting definition of alias')
                 endif
               else
                 if((dcindex(iafind).eq.0.or.                           &
     &           dcindex(iafind).eq.dcindex(ifind)).and.                &
     &           (dictyp(iafind).eq.' '.or.                             &
     &           (dictyp(iafind).eq.dictyp(ifind) .and.                 &
     &            dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
                 ifind=iafind
                 endif
                 if(xmindex(iafind).eq.0)                               &
     &             xmindex(iafind)=xmindex(ifind)
                 if(xmindex(ifind).eq.0)                                &
     &             xmindex(ifind)=xmindex(iafind)
                 if (dmcode(iafind).eq.0)                               &
     &             dmcode(iafind)=dmcode(ifind)
                 if (dmcode(ifind).eq.0)                                &
     &             dmcode(ifind)=dmcode(iafind)
                 if (dpindex(iafind).eq.iafind                          &
     &             .and. dpindex(ifind).ne.ifind)                       &
     &             dpindex(iafind) = dpindex(ifind)
                 if (dpindex(ifind).eq.ifind                            &
     &             .and. dpindex(iafind).ne.iafind)                     &
     &             dpindex(ifind) = dpindex(iafind)
                 if (deindex(ifind).eq.0)                               &
     &             deindex(ifind) = deindex(iafind)
                 if (deindex(iafind).eq.0)                              &
     &             deindex(iafind) = deindex(ifind)
                 aroot(iafind)     =aroot(ifind)
                 if(aroot(iafind).eq.0) aroot(iafind)=ifind
                 alias(ifind)      =iafind
                 if (catkey(iafind)) catkey(ifind) = .true.
                 if (catkey(ifind)) catkey(iafind) = .true.
               endif
             endif
           endif
           if(loop_) go to 260
         endif
         go to 200
!
400      bloc_=' '
         if (ndcname.ne.0) then
         do ii = idstrt+1,ndict
         keychain(ii) = 0
         if (aroot(ii).eq.0.and.dcindex(ii).eq.0                        &
     &     .and.catchk.eq.'yes')                                        &
     &     call tbxxwarn(' No category specified for name '//           &
     &       dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
         enddo
         endif
         do ii = idstrt+1,ndict
         if (dicxtyp(ii).eq.' ') then
           if (dpindex(ii).ne.ii                                        &
     &       .and. dicxtyp(dpindex(ii)).ne.' ') then
             dicxtyp(ii) = dicxtyp(dpindex(ii))
             dictyp(ii) = dicxtyp(dpindex(ii))(1:4)
           else
             dicxtyp(ii) = 'null'
             dictyp(ii) = 'null'
             if (tcheck.eq.'yes')  then
               jj = lastnb(dicnam(ii))
               if (jj.gt.0) then
               if (dicnam(ii)(jj:jj).ne.'_')                            &
     &         call tbxxwarn(' No type specified for name '//           &
     &           dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
               endif
             endif
           endif
         endif
         if (catkey(ii) .or. dmcode(ii).gt.0) then
           ifind = aroot(ii)
           mycat = dcindex(ifind)
           if (mycat.ne.0) then
             jj = ccatkey(mycat)
             if (jj.eq.0) then
               ccatkey(mycat) = ifind
             else
410            if (keychain(jj).eq.0) then
                 keychain(jj) = ifind
                 keychain(ifind) = 0
               else
                 if(keychain(jj).ne.ifind) then
                   jj = keychain(jj)
                   goto 410
                 endif
               endif
             endif
           endif
         endif
         enddo
         if (.not.append_) then
           close(dirdev)
           nrecd=0
         endif
         dictfl='no '
500      continue
         if (append_) then
           nrecd=nrecds
           recend_=recends
           recbeg_=recbegs
         endif
         if(dict_) then
           dicname_=xdicnam
           dicver_ =xdicver
         else
           tcheck = otchk
           vcheck = ovchk
         endif
         if(tcheck.eq.'yes') vcheck='yes'
!dbg     WRITE(6,'(i5,3x,a,2x,a)') (i,dicnam(i),dictyp(i),i=1,ndict)
         return
         end
!
!
!
!
!
! >>>>>> Create a new dictionary entry, or find a matching existing one
!
         function tbxxnewd(xname,ick)
         logical   tbxxnewd
         include  'ciftbx.sys'
         character xname*(*)
         character xxxtemp*(NUMCHAR)
         integer   jck, ick, ilen
         integer   lastnb
         tbxxnewd = .true.
         ilen = lastnb(xname)
         jck = ndict
         call tbxxnlc(xxxtemp,xname(1:ilen))
         call hash_store(xxxtemp,                                       &
     &     dicnam,dicchain,                                             &
     &     NUMDICT,ndict,dichash,NUMHASH,ick)
         if(ick.eq.0) call tbxxerr(' CIFdic names > NUMDICT')
         if(ick .eq. jck+1) then
           dictag(ick) = xname(1:ilen)
           dictyp(ick) = ' '
           dicxtyp(ick) = ' '
           catkey(ick) = .false.
           dpindex(ick) = ick
           deindex(ick) = 0
           alias(ick) = 0
           aroot(ick) = ick
           keychain(ick) = 0
           dcindex(ick) = 0
           xmindex(ick) = 0
           dmcode(ick) = 0
         else
           tbxxnewd = .false.
         endif
         return
         end
!
!
!
!
!
! >>>>>> Find matching existing dictionary entry if any
!
         function tbxxoldd(xname,ick)
         logical   tbxxoldd
         include  'ciftbx.sys'
         character xname*(*)
         character xxxtemp*(NUMCHAR)
         integer   ick, ilen
         integer   lastnb
         tbxxoldd = .true.
         ilen = lastnb(xname)
         call tbxxnlc(xxxtemp,xname(1:ilen))
         call hash_find(xxxtemp,                                        &
     &     dicnam,dicchain,                                             &
     &     NUMDICT,ndict,dichash,NUMHASH,ick)
         if(ick.eq.0) tbxxoldd = .false.
         return
         end
!
!
!
!
!
! >>>>>> Find position of last non_blank in a string
!        but never less than 1
!
         function lastnb(str)
!
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) str
         integer lenn,ihi,itestl
         lenn = len(str)
!
         ihi = lenn
         if(ihi.eq.0) then
           ihi = 1
           go to 200
         endif
         itestl = ihi/4
         if (itestl.lt.4) go to 200
!
100      if (ihi.gt.itestl) then
         if (str(ihi-itestl+1:ihi-itestl+1).eq.' ') then
           if (str(ihi-itestl+1:ihi).eq.' ') then
             ihi = ihi-itestl
             go to 100
           endif
         endif
         endif
         itestl = itestl/2
         if (itestl.gt.3) go to 100
!
200      if (ihi.gt.1 .and. str(ihi:ihi).eq.' ') then
           ihi = ihi-1
           go to 200
         endif
         if (ihi.eq.0) ihi = 1
         lastnb = ihi
         return
         end
!
!
!
!
!
! >>>>>> Convert a character to a radix XXRADIX digit
!
!        given a character c, return a decimal value
!
         function tbxxc2dig(c)
         integer   tbxxc2dig
         character*(*) c
         include  'ciftbx.sys'
!
         tbxxc2dig = ichar(c)-ichar(' ')
!
!        The code above may not be portable, especially to non-ascii
!        computer systems.  In that case, comment out the line above
!        and uncomment the following lines.  Be sure to make the
!        matching change in tbxxd2chr.  Be certain to have at least
!        XXRADIX characters in the search string.
!
!         tbxxc2dig = index(
!     *   '+-01234567890'//
!     *   'abcdefghijlmnopqrstuvwxyz'//
!     *   'ABCDEFGHIJKLMNOPQRSTUVWXYZ',c)-1
         return
         end
!
!
!
!
!
! >>>>>> Convert a radix XXRADIX digit to a character
!
!        given an integer value, return a character
!
         function tbxxd2chr(d)
         character*1 tbxxd2chr
         integer   d
         include  'ciftbx.sys'
!
         tbxxd2chr = char(d+ichar(' '))
!
!        The code above may not be portable, especially to non-ascii
!        computer systems.  In that case, comment out the line above
!        and uncomment the following lines.  Be sure to make the
!        matching change in tbxxc2dig.  Be certain to have at least
!        XXRADIX characters in the search string.
!
!         character*(XXRADIX) digits
!         digits =
!     *   '+-01234567890'//
!     *   'abcdefghijlmnopqrstuvwxyz'//
!     *   'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!         tbxxd2chr = digits(d+1:d+1)
         return
         end
!
!
!
!
!
! >>>>>> Convert a string to Run Length Encoded version
!
         subroutine tbxxrle(astr,bstr,mlen)
!
!        astr is the raw input string
!        bstr is the run-length-encoded string
!          beginning with the compressed length in
!            in base-XXRADIX in the first four characters
!          followed by either individual characters or run
!          flagged by XXFLAG
!        XXFLAG//tbxxd2chr(n)//c represents n copies of c
!
         character*(*) astr, bstr
         include  'ciftbx.sys'
         character*1 c
         character*1 tbxxd2chr
         integer tbxxc2dig
         integer klen, krep, ialen, iblen, mode, ii
         integer mlen
!
         ialen = len(astr)
         iblen = len(bstr)
         mode = 0
         klen = 4
         bstr(1:4) = tbxxd2chr(0)//tbxxd2chr(0)                         &
     &     //tbxxd2chr(0)//tbxxd2chr(0)
         do ii = 1,ialen
           c = astr(ii:ii)
           if (mode .eq. -2) then
             krep = tbxxc2dig(bstr(klen-1:klen-1))
             if (c.eq.bstr(klen:klen).and.krep.lt.XXRADIX-1) then
               bstr(klen-1:klen-1) = tbxxd2chr(krep+1)
             else
               mode = 0
               if (c.eq.bstr(klen:klen)) mode=-1
             endif
           endif
           if (klen.ge.iblen) go to 100
           if (mode .ge.-1 .and. mode .le.2) then
             klen = klen+1
             bstr(klen:klen) = c
             if (klen .gt. 5) then
               if (c.eq.bstr(klen-1:klen-1)) mode=mode+1
               if (c.ne.bstr(klen-1:klen-1)) mode=0
             endif
             if (c.eq.XXFLAG .and. klen.lt.iblen-1) then
               bstr(klen+1:klen+2) = tbxxd2chr(1)//c
               mode = -2
               klen = klen+2
             endif
           endif
           if (mode.eq.2) then
             bstr(klen-2:klen-1) = XXFLAG//tbxxd2chr(3)
             mode = -2
           endif
         enddo
 100     mlen = klen
         do ii = 4,1,-1
           bstr(ii:ii) = tbxxd2chr(mod(klen,XXRADIX))
           klen = klen/XXRADIX
         enddo
         return
         end
!
!
!
!
!
! >>>>>> Decode a string from  Run Length Encoded version
!
         function tbxxrld(astr,bstr,fill)
!
!        astr is the raw output string
!        bstr is the run-length-encoded string
!          beginning with the compressed length in
!            in base-XXRADIX in the first four characters
!          followed by either individual characters or run
!          flagged by char(0)
!        char(0)//char(n)//c represents n copies of c
!        fill is a logical variable, .true. to fill astr with blanks
!        the return value is the number of valid characters in astr
!        never less than 1
!
!
         character*(*) astr, bstr
         logical fill
         integer tbxxrld
         include  'ciftbx.sys'
         character*1 c
         integer tbxxc2dig
         integer klen, krep, ialen, mode, ipos
         integer ii, jj
!
         tbxxrld = 1
         krep = 0
         ialen = len(astr)
         if (fill) then
           astr = ' '
         else
           astr(1:1) = ' '
         endif
         mode = 0
         klen = 0
         do ii = 1,4
           klen = klen*XXRADIX+tbxxc2dig(bstr(ii:ii))
         enddo
         mode = 0
         ipos = 0
         do ii = 5,klen
           c = bstr(ii:ii)
           if(mode.eq.0) then
             if(c.ne.XXFLAG) then
               if (ipos.ge.ialen) then
                 tbxxrld = ialen
                 return
               endif
               ipos = ipos+1
               astr(ipos:ipos) = c
             else
               mode = 1
             endif
           else
             if (mode.eq.1) then
               krep = tbxxc2dig(c)
               mode = -1
             else
               do jj = 1,krep
                 if (ipos.ge.ialen) return
                 ipos=ipos+1
                 astr(ipos:ipos) = c
               enddo
               mode = 0
             endif
           endif
         enddo
         if(ipos .lt. ialen) astr(ipos+1:ipos+1) = ' '
         tbxxrld = max(ipos,1)
         return
         end
!
!
!
!
!
! >>>>>> Extract the item.category_id from a save frame name
!
         subroutine tbxxcat(sfname,bcname,lbcname)
!
         character*(*) sfname,bcname
         integer lbcname,ii,ic,lastnb,lenn
!
!        Note that this logic works only for item.category_id
!        not for category.id
!
         lenn = lastnb(sfname)
         bcname = ' '
         lbcname = 1
         if (lenn.eq.0.or.sfname(1:1).ne.'_') return
         do ii = 1,lenn-2
         ic = 1+lenn-ii
         if (sfname(ic:ic).eq.'.') then
           bcname = sfname(2:ic-1)
           lbcname = ic-2
           return
         endif
         enddo
         return
         end
!
!
!
!
!
! >>>>>> Fetch line from direct access file
!
         subroutine tbxxflin(linno,lip,lipag,lipof,kip,ip,mip,mis)
!
         include   'ciftbx.sys'
         integer    linno,lip,kip,ip,mip,mis,i,mipno,miprno, kzero
         integer    lipag,lipof,kmode
!
!        linno -- the line number to locate
!        lip   -- the location of the line
!                   (page*(NUMCPP/NUMCIP)+offset)
!        lipag -- the page number (1...)
!        lipof -- the offset (1...)
!        kip   -- subindex number
!        ip    -- subindex offset
!        mip   -- master index number
!        mis   -- master index offset

         kip = (linno-1)/NUMSIP + 1
         ip = mod(linno-1,NUMSIP) + 1
         mip = (kip-1)/NUMMIP + 1
         mis = mod(kip-1,NUMMIP) + 1
!
!        test subindex page number against number in memory
!
         if (kip.ne.iabs(ipim)) then
!
!          save the current subindex page if it has been written
!
           if (ipim.lt.0) then
             do i = 1,NUMSIP
               write(scrbuf(NUMCIP*(i-1)+1:NUMCIP*i),'(i8)')            &
     &           ippoint(i)
             enddo
             write(dirdev,'(a)',rec=iabs(iprim)) scrbuf
             ipim = -ipim
           endif
!
!          find the appropriate master index page and slot
!
           if (mip.ne.iabs(mipim)) then
!
!            save the current master index page if it has been written
!
             if (mipim.lt.0) then
               write(scrbuf(1:NUMCIP),'(i8)')mipcp
               do i = 1,NUMMIP
                 write(scrbuf(NUMCIP*i+1:NUMCIP*(i+1)),'(i8)')          &
     &             mippoint(i)
               enddo
               write(dirdev,'(a)',rec=iabs(miprim))scrbuf
               mipim = -mipim
             endif
!
!            search the master index pages for a match
!
             mipno = 0
             miprno = 1
             kzero = 0
             kmode = 1
 10          read(dirdev,'(a)',rec=miprno) scrbuf
             mipno = mipno+1
             read(scrbuf(1:NUMCIP),'(i8)') mipcp
             if (mipno.ne.mip) then
               if (mipcp.eq.0) then
                 if (nfword.gt.1) then
                   nfblock = nfblock+1
                   nfword = 1
                 endif
                 mipcp = nfblock
                 nfblock = nfblock+1
                 write(scrbuf(1:NUMCIP),'(i8)') mipcp
                 write(dirdev,'(a)',rec=miprno) scrbuf
                 scrbuf = ' '
                 write(scrbuf(1:NUMCIP),'(i8)') kzero
                 write(dirdev,'(a)',rec=mipcp) scrbuf
                 kmode = -1
               endif
               miprno = mipcp
               go to 10
             endif
!
!            Have the master index in scrbuf, copy to mippoint
!
             do i = 1,NUMMIP
               if (scrbuf(NUMCIP*i+1:NUMCIP*(i+1)).eq.' ') then
                 mippoint(i) = 0
               else
                 read(scrbuf(NUMCIP*i+1:NUMCIP*(i+1)),'(i8)')           &
     &             mippoint(i)
               endif
             enddo
             mipim =kmode* mip
             miprim = miprno
           endif
!
!          See if the subindex page exists
!
           if (mippoint(mis).eq.0) then
             do i = 1,NUMSIP
               ippoint(i) = 0
             enddo
             if (nfword.gt.1) then
               nfblock=nfblock+1
               nfword = 1
             endif
             mippoint(mis) = nfblock
             mipim = -iabs(mipim)
             ipim = -kip
             iprim = -nfblock
             scrbuf = ' '
             write(dirdev,'(a)', rec=nfblock) scrbuf
             nfblock = nfblock+1
           else
             read(dirdev,'(a)', rec=mippoint(mis)) scrbuf
             do i = 1,NUMSIP
               if (scrbuf(NUMCIP*(i-1)+1:NUMCIP*i).eq.' ') then
                 ippoint(i) = 0
               else
               read(scrbuf(NUMCIP*(i-1)+1:NUMCIP*i),'(i8)')             &
     &           ippoint(i)
               endif
             enddo
             ipim = kip
             iprim = mippoint(mis)
           endif
         endif
         lip = ippoint(ip)
         lipag = (lip-1)/(NUMCPP/NUMCIP) + 1
         lipof = mod(lip-1,NUMCPP/NUMCIP) + 1
         lipof = (lipof-1)*NUMCIP + 1
         return
         end

!
!
!
!
!
! >>>>>> Store a string in the string table
!
         subroutine tbxxsstb(astrg,sindex)
!
!        store string astrg in the string table, returning the
!        index in sindex
!
         character *(*) astrg
         integer sindex
         include  'ciftbx.sys'
         character *(MAXBUF) temp
         integer mlen, ii, ibstb, icstb, ikstb, rlen
         integer iestb

         call tbxxrle(astrg,temp,mlen)
         icstb = mod(nstable,NUMCSTB)+1
         ibstb = (nstable+NUMCSTB)/NUMCSTB
         iestb = min(NUMCSTB,icstb+mlen-1)
         ikstb = iestb-icstb+1
         if (mlen+nstable .le. NUMCSTB*NUMSTB) then
           stable(ibstb)(icstb:iestb)=temp(1:ikstb)
           sindex = nstable+1
           nstable = nstable+mlen
           rlen = mlen - ikstb
           if (rlen .gt. 0) then
             do ii = ikstb+1,mlen,NUMCSTB
               ibstb = ibstb+1
               iestb = min(NUMCSTB,rlen)
               stable(ibstb)(1:iestb) = temp(ii:ii+iestb-1)
               rlen = rlen - iestb
             enddo
           endif
         else
           sindex = 0
           call tbxxwarn(                                               &
     &      ' More than NUMCSTB*NUMSTB stable characters needed')
         endif
         return
         end

!
!
!
!
!
! >>>>>> Fetch a string from the string table
!
         function tbxxfstb(astrg,sindex,fill)
!
!        fetch string astrg from the string table, starting at the
!        index in sindex, and returning the valid length.
!
!        fill is a logical variable, .true. to fill astr with blanks
!        the return value is the number of valid characters in astr
!        never less than 1, unless there is no valid string

         integer tbxxfstb
         character *(*)astrg
         integer sindex
         logical fill
         integer tbxxc2dig, tbxxrld
         integer rlen
         integer icstb, ibstb, iestb, ikstb, klen, ii

         include  'ciftbx.sys'
         character *(MAXBUF) temp

         tbxxfstb = 0

         if (sindex.le.0.or.nstable+3.gt.NUMCSTB*NUMSTB) return

         icstb = mod(sindex-1,NUMCSTB)+1
         ibstb = (sindex-1+NUMCSTB)/NUMCSTB
         iestb = min(NUMCSTB,icstb+3)
         ikstb = iestb-icstb+1
         temp(1:ikstb)=stable(ibstb)(icstb:iestb)

         rlen = 4-ikstb
         if (rlen .gt. 0) then
         temp(ikstb+1:4)=stable(ibstb+1)(1:rlen)
         endif

         klen = 0
         do ii = 1,4
           klen = klen*XXRADIX+tbxxc2dig(temp(ii:ii))
         enddo

         if (klen.gt.MAXBUF.or.klen.le.0) return
         if (sindex+klen-1.gt.NUMCSTB*NUMSTB) return

         if (klen.gt.4) then
           icstb = mod(sindex+3,NUMCSTB)+1
           ibstb = (sindex+3+NUMCSTB)/NUMCSTB
           iestb = min(NUMCSTB,icstb+klen-5)
           ikstb = iestb-icstb+1
           temp(5:ikstb+4) = stable(ibstb)(icstb:iestb)
           rlen = klen - ikstb - 4
           if (rlen .gt. 0) then
             do ii = ikstb+1,ikstb+rlen,NUMCSTB
               ibstb = ibstb+1
               iestb = min(NUMCSTB,rlen)
               temp(ii:ii+iestb-1) = stable(ibstb)(1:iestb)
               rlen = rlen - iestb
             enddo
           endif
         endif

         tbxxfstb = tbxxrld(astrg,temp(1:klen),fill)
         return
         end


!
!
!
!
!
! >>>>>> Open a CIF and copy its contents into a direct access file.
!
         function ocif_(fname)
!
         logical   ocif_
         integer   lastnb
         include  'ciftbx.sys'
         logical   test
         character fname*(*)
         integer   lfname
         integer   case,i,kp,lp,mp
         integer   klen, mlen, lip, lppag, lipof, kip, ip, mip, mis
!
         save_=.false.
         glob_=.false.
         depth_=0
         jchar=MAXBUF
         lastch=0
         if(line_.gt.MAXBUF) call tbxxerr(' Input line_ value > MAXBUF')
         if(nrecd.ne.0 .and. (.not.append_)) then
           close(dirdev)
           nrecd=0
           lrecd=0
         endif
!
!        clear the memory resident page buffer
!
         do i = 1,NUMPAGE
         mppoint(i)=0
         enddo
!
         case=ichar('a')-ichar('A')
         tab=char(05)
         if(case.lt.0) goto 100
         tab=char(09)
         bloc_=' '
!
!....... Make sure the CIF is available to open
!
100      file_(1:longf_)=' '
         lfname = len(fname)
         file_(1:lfname) = fname
         do 120 i=1,lfname
         if(file_(i:i).eq.' ' .or. file_(i:i).eq.char(0) ) goto 140
120      continue
140      longf_=i-1
         if (longf_.gt.0) then
           inquire(file=file_(1:longf_),exist=test)
           ocif_=test
           if(.not.ocif_)      goto 200
         else
           file_(1:1) = ' '
           longf_ = 1
           ocif_ = .true.
         endif
!
!....... Open up the CIF and a direct access formatted file as scratch
!
         if (file_(1:1).ne.' ')                                         &
     &   open(unit=cifdev,file=file_(1:longf_),status='OLD',            &
     &                    access='SEQUENTIAL',                          &
     &                    form='FORMATTED')
         if(nrecd.eq.0)  then
           open(unit=dirdev,file='test',status='UNKNOWN',access='DIRECT',           &
     &                    form='FORMATTED',recl=NUMCPP)
           mipim = -1
           miprim = 1
           mipcp = 0
           ipim = -1
           iprim = 2
           do i = 1,NUMPAGE
             mppoint(i) = 0
           enddo
           do i = 1,NUMMIP
             mippoint(i) = 0
           enddo
           mippoint(1)=2
           do i = 1,NUMSIP
             ippoint(i) = 0
           enddo
           nfblock = 3
           nfword = 1
         endif
         if (mppoint(1).lt.0) then
            write(dirdev,'(a)',rec=-mppoint(1)) pagebuf(1)
            mppoint(1) = 0
         endif
         if(append_ .and. nrecd.ne.0) then
           kp = 1
           lp = nfblock
           nfblock = nfblock+1
           mppoint(kp) = lp
           mp = 1
         else
           do kp = 1,NUMPAGE
             mppoint(kp)=0
           enddo
           kp = 1
           lp = 3
           nfblock = 4
           mp = 1
         endif
!
!....... Copy the CIF to the direct access file
!
160      read(cifdev,'(a)',end=180) buffer
         nrecd=nrecd+1
         irecd=nrecd
         klen = lastnb(buffer(1:MAXBUF))
         if (klen.gt.line_)                                             &
     &     call tbxxwarn(' Input line length exceeds line_')
         call tbxxrle(buffer(1:klen),scrbuf,mlen)
         if (mp+mlen-1 .gt. NUMCPP) then
           if (mp.lt.NUMCPP) pagebuf(kp)(mp:NUMCPP) = ' '
!          write(dirdev,'(a)',rec=lp) pagebuf(kp)
           mppoint(kp)=-lp
           if (nfword.gt.1) then
             nfblock = nfblock+1
             nfword = 1
           endif
           lp = nfblock
           nfblock=nfblock+1
           kp = kp+1
           if(kp.gt.NUMPAGE) kp=1
           if (mppoint(kp).lt.0) then
             write(dirdev,'(a)',rec=-mppoint(kp)) pagebuf(kp)
           endif
           mppoint(kp)=0
           mp=1
         endif
         pagebuf(kp)(mp:mp+mlen-1)=scrbuf(1:mlen)
         mppoint(kp) = -lp
         mlen = ((mlen+NUMCIP-1)/NUMCIP)
         mlen = mlen*NUMCIP
         call tbxxflin(nrecd,lip,lppag,lipof,kip,ip,mip,mis)
         ippoint(ip) = (mp-1)/NUMCIP+(lp-1)*(NUMCPP/NUMCIP)+1
         ipim = -iabs(ipim)
         mp = mp+mlen
         goto 160
!
180      if (mp.lt.NUMCPP) pagebuf(kp)(mp:NUMCPP) = ' '
         if (mp.gt.1) then
!          write(dirdev,'(a)',rec=lp) pagebuf(kp)
           mppoint(kp)=-lp
         endif
         lrecd=max(0,recbeg_-1)
         jrecd=max(0,recbeg_-1)
         jrect=-1
         irecd=max(0,recbeg_-1)
         recn_=irecd
         recend_=nrecd
         if (file_(1:1).ne.' ') close(cifdev)
200      return
         end
!
!
!
!
!
! >>>>>> Close off direct access file of the current CIF
!         and reset all data name tables and pointers
!
         subroutine purge_
!
         include   'ciftbx.sys'
!
         integer i
         if(nrecd.ne.0) close(dirdev)
         do i = 1,NUMPAGE
           mppoint(i)=0
         enddo
         do i = 1,MAXBOOK
           ibkmrk(1,i)=-1
           ibkmrk(2,i)=-1
           ibkmrk(3,i)=-1
           ibkmrk(4,i)=-1
           ibkmrk(5,i)=-1
           ibkmrk(6,i)=-1
         enddo
         recn_=0
         save_=.false.
         glob_=.false.
         jchar=MAXBUF
         depth_=0
         lastch=0
         nrecd=0
         lrecd=0
         irecd=0
         nname=0
         nhash=0
         iname=0
         loopct=0
         loopnl=0
         loop_=.false.
         text_=.false.
         textfl='no '
         append_=.false.
         recbeg_=0
         recend_=0
         nivt = 0
         nstable = 0
         return
         end
!
!
!
!
!
! >>>>>> Store the data names and pointers for the requested data block
!
         function data_(name)
!
         logical   data_
         logical   wasave
         logical   tbxxoldd
         integer   lastnb
         include  'ciftbx.sys'
         character name*(*),temp*(NUMCHAR),ltype*4
         character ctemp*(NUMCHAR)
         character xdname*(NUMCHAR)
         character ydname*(NUMCHAR)
         character isbuf*(MAXBUF),lsbuf*(MAXBUF)
         logical   ixcat(NUMDICT)
         integer   ndata,idata,nitem,npakt,i,ii,j,k,kchar,krecd
         integer   jj,icc,idd
         integer   fcatnum,lctemp,isrecd,isjchr,islast
         integer   lsrecd,lsjchr,lslast
         integer   pnname,itpos,ipp,ipj
         integer   ltemp
!DBG     if(dictfl.eq.'no ')
!DBG *     print *,' ***>>>> Entering data_ ',name
!
         jchar=MAXBUF
         depth_=0
         nname=0
         ndata=0
         nhash=0
         nitem=0
         idata=0
         iname=0
         loopct=0
         loopnl=0
         ltype=' '
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         kchar = 0
         krecd = 0
         fcatnum = 0
         data_=.false.
         wasave=.false.
         loop_=.false.
         text_=.false.
         textfl='no '
         glob_=.false.
         do ii = 1,MAXBOOK
         ibkmrk(1,ii)=-1
         enddo
         irecd=lrecd
         lrecd=min(nrecd,recend_)
         if(name(1:1).ne.' ') irecd=max(0,recbeg_-1)
         call hash_init(dname,dchain,NUMBLOCK,nname,dhash,              &
     &     NUMHASH)
         call hash_init(cname,cchain,NUMBLOCK,ncname,chash,             &
     &     NUMHASH)
         isrecd=irecd
         isjchr=jchar
         islast=lastch
         lsrecd=isrecd
         lsjchr=isjchr
         lslast=islast
         isbuf=' '
         if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         lsbuf=' '
         if(lastch.gt.0)lsbuf(1:lastch)=isbuf(1:lastch)
         call tbxxnlc(xdname,name)
!
!....... Find the requested data block in the file
!
100      lsjchr=isjchr
         call getstr
         isjchr=jchar
         if(irecd.ne.isrecd) then
           lsrecd=isrecd
           lslast=islast
           lsbuf=' '
           if(islast.gt.0)lsbuf(1:islast)=isbuf(1:islast)
           isrecd=irecd
           islast=lastch
           isbuf=' '
           if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         endif
         if(type_.eq.'fini')           goto 500
         if(text_.or.depth_.gt.0)      goto 110
         goto 120
110      call getstr
         if (type_.eq.'fini')                                           &
     &      call tbxxerr(' Unexpected termination of file')
         if (text_.or.depth_.gt.0)     goto 100
         goto 100
120      continue
         if(type_.eq.'save') then
           if(long_.lt.6) then
             if(.not.save_)                                             &
     &         call tbxxerr(                                            &
     &           ' Save frame terminator found out of context ')
             wasave=.true.
             save_=.false.
             goto 100
           else
             if(save_)                                                  &
     &         call tbxxerr(' Prior save frame not terminated ')
             save_=.true.
             if(name.eq.' ')          goto 150
             call tbxxnlc(ydname,strg_(6:long_))
             if(ydname.ne.xdname) goto 100
             goto 150
           endif
         endif
         if(type_.eq.'glob') then
           if(name.ne.' ')            goto 100
           glob_=.true.
           goto 150
         endif
         if(type_.eq.'name'.or.type_.eq.'loop') then
           if(name.ne.' ')            goto 100
           if(.not.wasave)                                              &
     &       call tbxxwarn(' Data block header missing ')
           isrecd=lsrecd
           islast=lslast
           isjchr=lsjchr
           isbuf=' '
           if(islast.gt.0)isbuf(1:islast)=lsbuf(1:islast)
           data_=.true.
           bloc_=' '
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posnam_=itpos
           goto 204
         endif
         if(type_.ne.'data')          goto 100
         if(name.eq.' ')              goto 150
         call tbxxnlc(ydname,strg_(6:long_))
         if(ydname.ne.xdname)   goto 100
150      data_=.true.
         bloc_=strg_(6:long_)
!
!DBG     if(dictfl.eq.'no ')
!DBG *     print *, 'bloc_: '//bloc_
         itpos=jchar-long_
         if(tabx_) then
         itpos=0
         do ipp=1,jchar-long_
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         endif
         posnam_=itpos
!
!....... Get the next token and identify
!        ltype is the previous type
!
200      call getstr
!DBG     WRITE(6,*) ltype,type_,loop_,nitem,ndata,idata,iname,nname
!
         if(ltype.ne.'name')                goto 201
         if(type_.eq.'numb')                goto 203
         if(type_.eq.'char')                goto 203
         if(type_.eq.'text')                goto 203
         if(type_.eq.'null')                goto 203
         if(type_.eq.'tupl'                                             &
     &     .or.type_.eq.'tabl'                                          &
     &     .or.type_.eq.'arra')             goto 203
         if(type_.eq.'name'.and.loop_)      goto 204
!DBG     WRITE(6,*) ltype,type_,loop_,nitem,ndata,idata,iname,nname
         call tbxxerr(                                                  &
     &     ' Illegal tag/value construction: tag followed by '          &
     &     //type_)    
!
!        The prior type was not a name (not a tag)
201      if(ltype.ne.'valu')                goto 204
!
!        The prior type was a data value
!
         if(type_.eq.'numb')                goto 202
         if(type_.eq.'char')                goto 202
         if(type_.eq.'text')                goto 202
         if(type_.eq.'null')                goto 202
         if(type_.eq.'tupl'                                             &
     &     .or.type_.eq.'tabl'                                          &
     &     .or.type_.eq.'arra')             goto 202
         goto 204
!
!        If we have a vaue followed by a value, we need to be
!        in a loop (item > 0)
!
202      if(nitem.gt.0)                     goto 205
         call tbxxerr(                                                  &
     &     ' Illegal tag/value construction: value followed by '        &
     &     //type_)
!
!        The prior item was a tag and this is a value
!
203      ltype='valu'
!DBG     if(dictfl.eq.'no ')
!DBG *     print *, ' ***>>>>> data_ value ',strg_(1:long_)
         goto 205
!
!        Cases that get us here
!          The prior item was a tag and this is a tag in a loop
!          The prior item was neither a tag nor a value 
204      ltype=type_
!
!        We are in a loop and have a value after a value
!        or a name after a value or come from above cases
!
205      if(type_.eq.'name')           goto 206
         if(type_.eq.'loop')           goto 210
         if(type_.eq.'data')           goto 210
         if(type_.eq.'save')           goto 210
         if(type_.eq.'glob')           goto 210
         if(type_.ne.'fini')           goto 220
206      if(loop_)                     goto 270
210      if(nitem.eq.0)                goto 215
!
!....... End of loop detected; save pointers
!        loopni(loopct) -- number of items in a row
!        loopnp(loopct) -- number of rows
!
         npakt=idata/nitem
         if(npakt*nitem.ne.idata) call tbxxerr(' Item miscount in loop')
         loopni(loopct)=nitem
         loopnp(loopct)=npakt
         nitem=0
         idata=0
215      if(type_.eq.'name')           goto 270
         if(type_.eq.'data')           goto 300
         if(type_.eq.'save')           goto 300
         if(type_.eq.'glob')           goto 300
         if(type_.eq.'fini')           goto 300
!
!....... Loop_ line detected; incr loop block counter
!        record the character position in loopos(loopct)
!        record the line number in        loorec(loopct)
!        record the detabbed char pos in  loopox(loopct)
!
         loop_=.true.
!DBG     print *,' in data_ loop_ set, type_', type_
         loopct=loopct+1
         if(loopct.gt.NUMLOOP) call tbxxerr(                            &
     &     ' Number of loop_s > NUMLOOP')
         loorec(loopct)=irecd
         loopos(loopct)=jchar-long_
         if(quote_.ne.' ') then
           if (quote_.eq.';') then
             loopos(loopct) = 1
           else
             if (quote_.eq.''''''''.or.quote_.eq.'"""') then
               loopos(loopct)=jchar-long_-3
             else
               loopos(loopct)=jchar-long_-1
             end if
           end if
         end if
         itpos=0
         do ipp=1,loopos(loopct)
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         loopox(loopct)=itpos
         goto 200
!
!....... This is the data item; store char position and length
!
220      if(loop_ .and. nitem.eq.0)                                     &
     &   call tbxxerr(' Illegal tag/value construction')
         loop_=.false.
!
         i=nname
         if(nitem.gt.0) i=i-nitem+mod(idata,nitem)+1
         if(i.lt.1) call tbxxerr(' Illegal tag/value construction')
         if(dtype(i).ne.'test')       goto 223
         if(dictfl.eq.'yes')          goto 223
         if(tcheck.eq.'no ')          goto 223
!>>>>    if(long_.eq.1.and.strg_(1:1).eq.'?') goto 223
!>>>>    if(long_.eq.1.and.strg_(1:1).eq.'.') goto 223
         if(type_.eq.'null')          goto 223
         if(type_.eq.'numb')          goto 223
         call tbxxwarn( ' Numb type violated  '//dname(i))
223      if(nitem.le.0)               goto 224
         idata=idata+1
         if(dtype(i).eq.'null') dtype(i)=type_
         if(dtype(i).eq.'numb' .and.                                    &
     &     (type_.eq.'char'.or.type_.eq.'text')) dtype(i)='char'
224      if(nname.eq.ndata)           goto 230
         ndata=ndata+1
         if(iloop(ndata).gt.1)        goto 225
         krecd=irecd
         kchar=jchar-long_-1
         if(quote_.ne.' ') then
           kchar=kchar-1
           if (quote_(2:3).ne.'  ') kchar=kchar-2
         end if
225      continue
         if(dtype(ndata).eq.'    ') dtype(ndata)=type_
         drecd(ndata)=krecd
         dchar(ndata)=kchar
         if (depth_.gt.0) then
!DBG     print *,' Setting bracket start at ',
!DBG *     'char: ', posbrkstk(1)-1, 'rec: ',srecd
         dchar(ndata) = posbrkstk(1)-1
         drecd(ndata) = srecd

         end if
         if(nloop(ndata).gt.0)        goto 230
         nloop(ndata)=0
         iloop(ndata)=long_
         if (depth_.gt.0) iloop(ndata) = 1
!
!....... Skip text lines if present
!
230      if(type_.ne.'text')           goto 250
!DBG     print *,' text field detected at 230 '
         if(nloop(ndata).eq.0.and.depth_.eq.0) dchar(ndata)=0
         if(nloop(ndata).eq.0.and.depth_.eq.0) iloop(ndata)=long_
240      call getstr
         if(type_.eq.'fini') call tbxxerr(' Unexpected end of data')
         if (type_.ne.'text'.or..not.text_) then
           if (depth_.eq.0)            goto 200
           goto 260
         endif
         goto 240
!
!....... Skip bracketed construct if present
!
250      if(depth_.eq.0)           goto 200
         
260      call getstr
         if(depth_.eq.0) goto 200
         if(type_.eq.'fini') call tbxxerr(' Unexpected end of data')
         if(type_.eq.'text') goto 240
         goto 260
!
!....... This is a data name; store name and loop parameters
!
270      call tbxxclc(temp,ltemp,strg_(1:long_),long_)
         k=0
         if(dictfl.ne.'yes' .and. ndict.gt.0) then
           tbxxrslt = tbxxoldd(temp(1:ltemp),k)
           if(k.ne.0) then
             if(alias_ .and. aroot(k).ne.0) then
               temp=dicnam(aroot(k))
               ltemp = lastnb(temp)
             endif
           endif
         endif
         pnname=nname
         call hash_store(temp(1:ltemp),                                 &
     &   dname,dchain,NUMBLOCK,nname,dhash,                             &
     &     NUMHASH,j)
!DBG     if(dictfl.eq.'no ')
!DBG *     print *,' ***>>>>> data_ name: ',temp(1:ltemp)
         if(j.eq.pnname+1) then
           dtag(j)=strg_(1:long_)
           if(k.ne.0) dtag(j)=dictag(k)
           trecd(j)=irecd
           tchar(j)=jchar-long_
           if(quote_.ne.' '.and.quote_.ne.';')                          &
     &       tchar(j)=jchar-long_-1
           itpos=0
           do ipp=1,tchar(j)
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           xchar(j)=itpos
         endif
         if(j.eq.0)                                                     &
     &     call tbxxerr(' Number of data names > NUMBLOCK')
         if(k.ne.0) then
           ltemp = lastnb(dicnam(k))
           temp(1:ltemp) = dicnam(k)(1:ltemp)
         endif
         if(j.ne.pnname+1) then
           call tbxxwarn(' Duplicate data item '//                      &
     &     temp(1:ltemp))
           goto 200
         endif
         dtype(nname)=' '
         dxtyp(nname)=' '
         cindex(nname)=0
         ddict(nname)=0
         ctemp(1:6)='(none)'
         lctemp=6
!
         if(dictfl.eq.'yes' .or. vcheck.eq.'no ') goto 290
         j=k
         if(j.ne.0) then
           ddict(nname)=j
           cindex(nname)=dcindex(j)
           dxtyp(nname)=dicxtyp(j)
           dtype(nname)=dictyp(j)
           if(vcheck.eq.'no ')          goto 280
           if(dictyp(j).eq.'numb') then
             dtype(nname)='test'
           endif
           if(cindex(nname).ne.0) then
             lctemp=lastnb(dcname(cindex(nname)))
             ctemp(1:lctemp)=dcname(cindex(nname))(1:lctemp)
             goto 290
           endif
           goto  280
         endif
         call tbxxwarn(' Data name '//                                  &
     &               temp(1:ltemp)                                      &
     &               //' not in dictionary!')
280      call tbxxcat(temp(1:ltemp),ctemp,lctemp)
         if (ctemp(1:lctemp).eq.' '.or.                                 &
     &     ('_'//ctemp(1:lctemp).eq.temp(1:ltemp))) then
           ctemp = '(none)'
           lctemp= 6
           if (ndcname.ne.0.and.vcheck.eq.'yes')                        &
     &       call tbxxwarn(' No category defined for '                  &
     &       //temp(1:ltemp))
         else
           call hash_find(ctemp(1:lctemp),                              &
     &       dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,j)
           if(j.ne.0) then
             cindex(nname) = j
           else
             ipj=ncname
             call hash_store(ctemp(1:lctemp),                           &
     &         cname,cchain,NUMBLOCK,ncname,chash,NUMHASH,j)
             if (j.eq.0)                                                &
     &         call tbxxerr(' Number of categories > NUMBLOCK ')
             cindex(nname) = -j
             if (ndcname.gt.0.and.j.eq.ipj+1.and.vcheck.eq.'yes'        &
     &         .and.catchk.eq.'yes')                                    &
     &         call tbxxwarn(' Category '//                             &
     &         ctemp(1:lctemp)//' first implicitly defined in cif ')
           endif
         endif
!
290      lloop(nname)=0
         nloop(nname)=0
         iloop(nname)=0
         if (nitem.eq.0) fcatnum=cindex(nname)
         if(.not.loop_)               goto 200
         nitem=nitem+1
         if(nitem.gt.NUMITEM)                                           &
     &     call tbxxerr(' Items per loop packet > NUMITEM')
         nloop(nname)=loopct
         iloop(nname)=nitem
         if (fcatnum.ne.cindex(nname)) then
           temp = '(none)'
           if (fcatnum.gt.0) temp=dcname(fcatnum)
           if (fcatnum.lt.0) temp=cname(-fcatnum)
           ltemp = lastnb(temp)
           if (ctemp(1:lctemp).ne.temp(1:ltemp)                         &
     &     .and.catchk.eq.'yes')                                        &
     &     call tbxxwarn (' Heterogeneous categories in loop '//        &
     &     ctemp(1:lctemp)//' vs '//                                    &
     &     temp(1:ltemp))
           fcatnum=cindex(nname)
         endif
         goto 200
300      continue
!
!....... Are names checked against dictionary?
!
         if(dictfl.eq.'yes')          goto 500
         if(vcheck.eq.'no '.or.ndict.eq.0) goto 500
         do i=1,nname
           if(dtype(i).eq.'test') dtype(i)='numb'
         enddo

!
!        prepare for category and parent checks
!
         if ((catchk.eq.'yes'.or.parchk.eq.'yes')                       &
     &      .and. ndict.gt.0) then
         do i = 1,ndict
           ixcat(i) = .false.
         enddo
!
!        make a pass marking all used tags and their aliases
!
         do i = 1,nname
           icc=cindex(i)
           idd=ddict(i)
           if(icc.ne.0.and.idd.ne.0) then
             icc = aroot(idd)
310          ixcat(icc) = .true.
             icc = alias(icc)
             if (icc.ne.0) goto 310
           endif
         enddo
         endif
!
!        check for category keys
!
!
!
!        now make a pass making certain the keys are
!        used
!
         if(catchk.eq.'yes' .and. ndict.gt.0) then
         do i = 1,nname
           idd=cindex(i)
           if (idd.gt.0) then
             icc=ccatkey(idd)
             if(icc.ne.0) then
             if(aroot(icc).ne.0) icc=aroot(icc)
320          if(icc.ne.0) then
               if(.not.ixcat(icc)) then
                 jj = irecd
                 irecd = drecd(i)
                 if (catkey(icc)) then
				 
				 !KJ 11-2011 for FEFF9 code : removing the warnings below.  They seem irrelevant to my case, and are confusing to the user.
    !               call tbxxwarn(' Category key '//                     &
    ! &               dictag(icc)(1:lastnb(dictag(icc)))//               &
    ! &               ' not given for '//                                &
    ! &               dcname(idd)(1:lastnb(dcname(idd))))
                 else
                   call tbxxwarn(' Mandatory item '//                   &
     &               dictag(icc)(1:lastnb(dictag(icc)))//               &
     &               ' not given for '//                                &
     &               dcname(idd)(1:lastnb(dcname(idd))))
                 endif
                 ixcat(icc) = .true.
                 irecd = jj
               endif
               icc = keychain(icc)
               if(icc.ne.0) go to 320
             endif
             endif
           endif
         enddo
         endif
!
!        check for parents of tags that are used
!
         if(parchk.eq.'yes' .and. ndict.gt.0) then
         do i = 1,nname
           if (ddict(i).ne.0) then
             if (dpindex(ddict(i)).ne.ddict(i)) then
               if (.not.ixcat(dpindex(ddict(i)))) then
                 call tbxxwarn(' Parent '//                             &
     &             dicnam(dpindex(ddict(i)))                            &
     &             (1:lastnb(dicnam(dpindex(ddict(i)))))//              &
     &             ' of '//                                             &
     &             dname(i)(1:lastnb(dname(i))) //                      &
     &             ' not given')
               endif
             endif
           endif
         enddo
         endif

!
!....... End of data block; tidy up loop storage
!
500      lrecd=irecd-1
         if(type_.eq.'save'.and.long_.lt.6) then
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posval_=itpos
         endif
         irecd=isrecd
         jchar=isjchr
         lastch=islast
         recn_=irecd
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=isbuf(1:lastch)
         jrecd=irecd
         loop_=.false.
         loopct=0
         if(ndata.ne.nname) call tbxxerr(' Syntax construction error')
!
!dbg     WRITE(6,'(a)')
!dbg *   ' data name                       type recd char loop leng'
!dbg     WRITE(6,'(a,1x,a,4i5)') (dname(i),dtype(i),drecd(i),dchar(i),
!dbg *              nloop(i),iloop(i),i=1,nname)
!dbg     WRITE(6,'(3i5)') (i,loopni(i),loopnp(i),i=1,loopct)
!
         return
         end
!
!
!
!
!
! >>>>>> Check dictionary for data name validation
!
         function dtype_(name,type)
!
         logical    dtype_, tbxxoldd
         include   'ciftbx.sys'
         integer    nln, ii
         character  name*(*),temp*(NUMCHAR),                            &
     &              type*4
!
         character*4 map_type(19),map_to(19),mapped
         data map_type                                                  &
     &   /'floa','int ','yyyy','symo','ucha','ucod','name','idna',      &
     &    'any ','code','line','ulin','atco','fax ','phon','emai',      &
     &    'real','inte','coun'/
         data map_to                                                    &
     &   /'numb','numb','char','char','char','char','char','char',      &
     &    'char','char','char','char','char','char','char','char',      &
     &    'numb','numb','numb'/

!
         type = ' '
         dtype_ = .false.
         nln = min(len(name),len(temp))
         call tbxxnlc(temp(1:nln),name)
         if (ndict.eq.0) go to 200
         tbxxrslt = tbxxoldd(temp(1:nln),xdchk)
         if(xdchk.eq.0) go to 200
         mapped = dictyp(xdchk)(1:4)
         do ii = 1,19
           if (dictyp(xdchk)(1:4).eq.map_type(ii)) mapped = map_to(ii)
         enddo
         if (mapped.ne.'char'.and.mapped.ne.'numb'                      &
     &     .and.mapped.ne.'null'.and.mapped.ne.'text') then
           call tbxxwarn(' Item type '                                  &
     &        //dictyp(xdchk)(1:max(1,lastnb(dictyp(xdchk))))//         &
     &       ' for item '//                                             &
     &       name(1:max(1,lastnb(name)))//' not recognized ')
           mapped = 'char'
         endif
         type = mapped
         dtype_ = .true.
200      continue
         return
         end
!
!
!
!
!
!
! >>>>>> Get the attributes of data item associated with data name
!
         function test_(temp)
!
         logical   test_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character*4 otype
         integer lname
         character  otestf*3
!
         otestf=testfl
         otype = type_
         testfl='yes'
         test_ = .false.
         call tbxxclc(name,lname,temp,len(temp))
!DBG     print *,' Entering test_ ',name(1:lname)
         if (depth_.eq.0) go to 100
         if (name(1:1).ne.' '.and.name(1:1).ne.char(0).and.             &
     &     name(1:lname).ne.nametb(1:lnametb))     goto 120
         call getstr
         test_=.true.
         if (type_.eq.'null') test_=.false.
         if (otype.eq.'text' .and. (.not. text_) .and.long_.eq.0) then
           quote_=' '
           textfl = 'no'
           type_ = 'null'
           test_ = .false.
           goto 200
         end if
         posval_ = jchar-long_
         posend_ = jchar-1
         if (long_.gt.0) then
            if (type_.eq.'numb') then
               call ctonum
               if(posdec_.gt.0) posdec_=posval_+posdec_-1
               jchar = posend_
            else 
              if (quote_.eq.' ') then
                 if (long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
                 if (long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
              end if
            end if
         end if
         goto 200

         
100      test_=.true.
         if(otestf.eq.'no ' .or. type_.eq.' ')  goto 120
         if(name(1:lname).eq.nametb(1:lnametb))   goto 200
120      call tbxxgitm(name(1:lname))
200      list_ =loopnl
         if(type_.eq.'null') test_=.false.
         if(type_.ne.'null'.and.type_.ne.'char'.and.                    &
     &     type_.ne.'text'.and.type_.ne.'numb') type_='char'
!DBG     print *,' leaving test_ ', type_, depth_, strg_(1:long_)
         return
         end

!
!
!
!
!
! >>>>>> Set or Reference a bookmark
!
         function bkmrk_(mark)
!
         logical   bkmrk_
         include   'ciftbx.sys'
!
         integer   mark,ii,nitem
         character*4 flag
         bkmrk_=.true.
         if(mark.eq.0) then
           do ii=1,MAXBOOK
             if(ibkmrk(1,ii).lt.0)      goto 100
           enddo
           bkmrk_=.false.
           call tbxxwarn(' More than MAXBOOK bookmarks requested')
           return
100        mark=ii
           ibkmrk(1,ii)=iname
           ibkmrk(2,ii)=irecd
           ibkmrk(3,ii)=jchar
           if(iname.gt.0) then
             ibkmrk(2,ii) = trecd(iname)
             ibkmrk(3,ii) = tchar(iname)
           endif
           ibkmrk(4,ii)=0
           if(iname.gt.0) then
             if(nloop(iname).ne.0.and.                                  &
     &         loopnl.eq.nloop(iname).and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               ibkmrk(2,ii)=looprd(1)
               ibkmrk(3,ii)=max(0,loopch(1)-1)
               ibkmrk(4,ii)=loopct
             endif
           endif
           ibkmrk(5,ii) = depth_
           ibkmrk(6,ii) = index_
         else
           if(ibkmrk(1,mark).lt.0) then
             bkmrk_=.false.
             return
           endif
           iname=ibkmrk(1,mark)
           irecd=ibkmrk(2,mark)
           loopct=ibkmrk(4,mark)
           loop_=.false.
           text_=.false.
           textfl = 'no '
           loopnl=-1
           testfl='no '
           if(iname.gt.0) then
            if(nloop(iname).ne.0.and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               looprd(nitem+1)=ibkmrk(2,mark)
               loopch(nitem+1)=ibkmrk(3,mark)
               do ii = 1,nitem
                 lloop(ii+iname-iloop(iname))=loopct-1
               enddo
               loopct=loopct-1
               if(lloop(iname).gt.0) then
                 loop_=.true.
                 loopnl=nloop(iname)
               endif
             endif
           endif
           jchar=MAXBUF
           if(irecd.gt.0) then
             irecd=irecd-1
             call getlin(flag)
             jchar=ibkmrk(3,mark)
           endif
           depth_=0
           index_=0
           if (ibkmrk(5,mark).gt.0) then
200           call getstr
              if (depth_ .lt. 1) then
                call tbxxwarn(                                          &
     &          ' Bookmark for list, array, tuple or table corrupted')
                go to 210
              end if
              if(ibkmrk(5,mark).ne.depth_                               &
     &          .or. ibkmrk(6,mark).ne.index_ ) go to 200
           endif
210        ibkmrk(1,mark)=-1
           mark=0
         endif
         return
         end
!
!
!
!
!
!
! >>>>>> Find the location of the requested item in the CIF
!        The argument "name" may be a data item name, blank
!        for the next such item.  The argument "type" may be
!        blank for unrestricted acceptance of any non-comment
!        string (use cmnt_ to see comments), including loop headers,
!        "name" to accept only the name itself and "valu"
!        to accept only the value, or "head" to position to the
!        head of the CIF.  Except when the "head" is requested,
!        the position is left after the data item provided.
!
         function find_(name,type,strg)
!
         logical   find_
         include   'ciftbx.sys'
         character  name*(*),type*(*),strg*(*),flag*4
         character  jjbuf*(MAXBUF)
         integer    jjchar,jjrecd,jjlast,jjlrec,jjjrec,jjdepth,jindex
!
!DBG     print *,' Entering find ', name, type
         find_  = .false.
         strg   = ' '
         long_  = 0
         jjchar = jchar
         jjrecd = lrecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjdepth = depth_
         jindex = index_
         jjbuf  = ' '
         if(lastch.gt.0) jjbuf(1:lastch)=buffer(1:lastch)
         if(type.eq.'head') then
           lrecd = min(nrecd,recend_)
           irecd = max(0,recbeg_-1)
           jchar=MAXBUF+1
           depth_=0
           call getlin(flag)
           if(flag.eq.'fini')       goto 300
           find_=.true.
           lrecd=max(0,recbeg_-1)
           return
         endif
         if(name.ne.' ') then
           testfl='no '
           call tbxxgitm(name)
           if(iname.eq.0) goto 300
           if(type.eq.'valu') then
             list_=loopnl
             strg=strg_(1:long_)
             find_=.true.
             return
           endif
           if(type.eq.'name'.or.loopnl.eq.0) then
             irecd=trecd(iname)-1
             call getlin(flag)
             jchar=tchar(iname)
             depth_=0
             posnam_=jchar+1
             call getstr
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           if(type.eq.' ') then
             irecd=loorec(loopnl)-1
             call getlin(flag)
             jchar=loopos(loopnl)
             depth_=0
             call getstr
             posval_=loopos(loopnl)
             if(tabx_) posval_=loopox(loopnl)
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           call tbxxerr(' Call to find_ with invalid arguments')
         endif
         if(name.eq.' ') then
           go to 200
190        if (text_.or.depth_.gt.0) then
              call getstr
              if (type_.eq.'fini')  goto 300
              if (type_.ne.'null')  goto 190  
           end if     
200        call getstr
           if(type_.eq.'fini')      goto 300
           if(type.ne.' '.and.                                          &
     &      (type_.eq.'data'.or.type_.eq.'save'.or.                     &
     &      type_.eq.'glob'))   goto 300
           if(type.eq.'name'.and.type_.ne.'name')  goto 190
           if(type.eq.'valu'.and.                                       &
     &       type_.ne.'numb'.and.type_.ne.'text'                        &
     &      .and.type_.ne.'char'.and.type_.ne.'null') goto 190
           find_=.true.
           strg=strg_(1:long_)
           if(type_.eq.'name') then
             posnam_=jchar-long_
           else
             posval_=jchar-long_
             if(quote_.ne.' '.and.quote_.ne.';')                        &
     &         posval_=posval_-1
             if(quote_.eq.'''''''' .or.quote_.eq.'"""')                 &
     &         posval_=posval_-2
           endif
           recn_=irecd
           return
         endif

!
!        Search failed, restore pointers
!
300      irecd  = jjrecd
         lastch = jjlast
         lrecd  = jjlrec
         jchar  = jjchar
         depth_ = jjdepth
         index_ = jindex

         buffer(1:1) = ' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd  = jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_  = irecd
!
         return
         end
!
!
!
!
!
!
! >>>>>> Get the next data name in the data block
!
         function name_(temp)
!
         logical    name_
         include   'ciftbx.sys'
         character  temp*(*)
!
         name_=.false.
         temp=' '
         iname=iname+1
         if(iname.gt.nname)  goto 100
         name_=.true.
         temp=dtag(iname)
         if(ddict(iname).ne.0) temp=dictag(ddict(iname))
100      return
         end
!
!
!
!
!
!
! >>>>>> Extract a number data item and its standard deviation
!        This version return single precision numbers
!
         function numb_(temp,numb,sdev)
!
         logical    numb_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         integer    lname
         real       numb,sdev

!DBG     print *,'***>>> Entering numb_ for ', temp

!
         call tbxxclc(name,lname,temp,len(temp))
         if(testfl.eq.'yes')    goto 100
         if(depth_.eq.0)        goto 120
         if(name(1:1).ne.' '.and.name(1:1).ne.char(0).and.              &
     &     name(1:lname).ne.nametb(1:lnametb))     goto 120
     
         numb_ = .false.
         call getstr
         if (type_.ne.'numb') go to 200 
         if (type_.eq.'numb') then
            call ctonum
            if(posdec_.gt.0) posdec_=posval_+posdec_-1
            numb_ = .true.
            if (depth_.gt.0) jchar=jchar-1
         end if
!DBG     print *,'***>>> In numb_ strg_ ', strg_(1:long_)
         go to 200

100      if(name(1:lname).eq.nametb(1:lnametb))   goto 150
!
120      call tbxxgitm(name(1:lname))
!
150      continue
!DBG     print *,'***>>> In numb_ strg_ ', strg_(1:long_)
         numb_=.false.
         if(type_.ne.'numb') goto 200
         numb_=.true.
         numb =sngl(numbtb)
         if(sdevtb.ge.0.0) sdev=sngl(sdevtb)
!
200      testfl='no '
         return
         end
!
!
!
!
!
!
! >>>>>> Extract a number data item and its standard deviation
!        This version returns double precision numbers
!
         function numd_(temp,numb,sdev)
!
         logical    numd_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         integer    lname
         double precision numb,sdev
!DBG     print *,'***>>> Entering numb_ for ', temp
!
         call tbxxclc(name,lname,temp,len(temp))
         if(testfl.eq.'yes')    goto 100
         if(depth_.eq.0)        goto 120
         if(name(1:1).ne.' '.and.name(1:1).ne.char(0).and.              &
     &     name(1:lname).ne.nametb(1:lnametb))     goto 120
     
         numd_ = .false.
         call getstr
         if (type_.ne.'numb') go to 200 
         if (type_.eq.'numb') then
            call ctonum
            if(posdec_.gt.0) posdec_=posval_+posdec_-1
            numd_ = .true.
            if (depth_.gt.0) jchar=jchar-1
         end if
!DBG     print *,'***>>> In numd_ strg_ ', strg_(1:long_)
         go to 200

100      if(name(1:lname).eq.nametb(1:lnametb))   goto 150
!
120      call tbxxgitm(name(1:lname))
!
150      numd_=.false.
!DBG     print *,'***>>> In numd_ strg_ ', strg_(1:long_)
         if(type_.ne.'numb') goto 200
         numd_=.true.
         numb =numbtb
         if(sdevtb.ge.0.0) sdev=sdevtb
!
200      testfl='no '
         return
         end
!
!
!
!
!
!
! >>>>>> Extract a character data item.
!
         function char_(temp,strg)
!
         logical    char_, charnp_
         include   'ciftbx.sys'
         character  temp*(*), strg*(*)
         integer    lstrg,nstrg

         nstrg = len(strg)
         char_ = charnp_(temp,strg,lstrg)

         if (lstrg.lt.len(strg)) strg(lstrg+1:nstrg) = ' '
         return
         end

!
!
!
!
!
!
! >>>>>> Extract a character data item, no padding.
!
         function charnp_(temp,strg,lstrg)
!
         logical    charnp_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  strg*(*),flag*4
         integer    lstrg
         character*1 slash
         character*4 otype
         integer    ltemp, lname, klow
         integer    lastnb
!
         slash = rsolidus(1:1)
         ltemp = lastnb(temp)
         otype = type_
         call tbxxclc(name,lname,temp,ltemp)
         if(testfl.eq.'yes')    goto 100
         if(.not.text_.and.depth_.eq.0)         goto 120
         if(name(1:1).ne.' '.and.name(1:1).ne.char(0).and.              &
     &     name(1:lname).ne.nametb(1:lnametb))     goto 120
         charnp_=.false.
         lstrg = 1
         strg=' '
         call getstr
         if (type_.eq.'fini') goto 200
         if (otype.eq.'text' .and. (.not. text_) .and.long_.eq.0) then
           quote_=' '
           textfl = 'no'
           charnp_=.false.
           type_ = 'null'
           goto 200
         end if
         posval_ = jchar-long_
         posend_ = jchar-1
         if (long_.gt.0) then
            strg=strg_(1:long_)
            lstrg=long_
            if (type_.eq.'numb') then
               call ctonum
               if(posdec_.gt.0) posdec_=posval_+posdec_-1
            else 
              if (quote_.eq.' ') then
                 if (long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
                 if (long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
              end if
            end if
         end if
         charnp_=.true.
         goto 200
!
100      if(name(1:lname).eq.nametb(1:lnametb))     goto 150
!
120      quote_=' '
         call tbxxgitm(name(1:lname))
         text_=.false.
         if(type_.eq.'null') then
           charnp_=.false.
           text_=.false.
           textfl = 'no '
           strg_=' '
           long_=0
           goto 200
         endif
!
!        strg_(1:long_) loaded with item
!
150      charnp_=.true.
         strg(1:1)=' '
         lstrg = 1
         if(long_.gt.0) then
           strg=strg_(1:long_)
           lstrg = long_
         endif
         if(type_.eq.'char' )   goto 200
         charnp_=.false.
         if(type_.ne.'text')   goto 200
         charnp_=.true.
         call getlin(flag)
         jchar=MAXBUF+1
         if(flag.eq.'fini')    goto 200
         if(buffer(1:1).eq.';')then
           jchar=2
           textfl = 'no '
           quote_=';'
           goto 200
         endif
         irecd=irecd-1
         text_=.true.
         if (long_.gt.0) then
         if (unfold_ .and. strg(long_:long_).eq.slash) then
170        klow = long_
           long_ = long_-1
           call getlin(flag)
           if(flag.eq.'fini')    goto 210
           if(buffer(1:1).eq.';') then
             jchar=2
             textfl = 'no '
           goto 210
           endif
           quote_=' '
           jchar=lastch+1
           long_=min(len(strg_),klow+max(1,lastch)-1)
           strg_(klow:long_)=buffer(1:max(1,lastch))
           strg(long_:long_)=' '
           if(lastch.gt.0) then
             long_=min(len(strg),klow+lastch-1)
             if(long_.ge.klow) strg(klow:long_)=buffer(1:lastch)
           endif
           if( strg(long_:long_).eq.slash ) go to 170
         endif
         endif
!
200      testfl='no '
         if(long_.eq.0) strg(1:1)=' '
         lstrg = max(1,long_)
!DBG     print *,' Leaving charnp_ text_, type_, quote_: ',
!DBG *     charnp_,text_,type_,quote_
!DBG     print *, ':>>>:'//strg_(1:lstrg)
         if (type_.eq.'char' .and. quote_.eq.' ') then 
           if (strg(1:lstrg).eq.'?'.or. strg(1:lstrg).eq.'.')           &
     &       type_='null'
         end if
         return
!
210      text_ = .false.
         go to 200
!
         end
!
!
!
!
!
! >>>>>> Extract a comment or terminal delimiter field
!        backing up to a prior delimiter, depth_ will not
!        be changed even when crossing a terminal delimiter
!
         function cotdb_(strg,lstrg,istd,posstart,recstart)
!
         logical   cotdb_
         logical   istd
         integer   posstart,recstart
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,                                 &
     &     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   lstrg
         integer   ipp,itpos
         integer   klow
         character*1 slash
!
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         istd = .false.
         slash = rsolidus(1:1)
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if (irecd.ne.recstart) then
           irecd = recstart-1
           call getlin(flag)
           if(flag.eq.'fini') then
             strg='fini'
             jchar=MAXBUF+1
             lstrg=4
             cotdb_=.false.
             posstart = jchar
             recstart = irecd
             go to 300
           endif
         end if
         jchar = posstart
         strg=' '
         lstrg=0
         cotdb_=.false.
100      jchar=jchar+1
         if (jchar.gt.jjchar.and.irecd.ge.jjrecd) go to 300
         if(jchar.le.lastch)     goto 140
!
!....... Read a new line
!
         call getlin(flag)
         if(flag.eq.'fini') then
           cotdb_=.false.
           posstart = jchar
           recstart = irecd
           strg='fini'
           lstrg=4
           go to 300
         endif
         jchar=0
         strg=char(0)
         lstrg=1
         posnam_=0
         quote_=' '
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
!
!....... Process this character in the line
!
         c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         if(depth_.gt.0 .and.                                           &
     &    ((c.eq.']'.and.rdbkt_)                                        &
     &      .or.(c.eq.')'.and.rdprn_)                                   &
     &      .or.(c.eq.'}'.and.rdbrc_))) go to 250
         goto 300
!
!        For a tab, when not expanding to blanks, accept
!        that single character as a comment
!
190      lstrg=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
!
!....... Accept the remainder of the line as a comment
!
200      lstrg=lastch-jchar
         quote_=buffer(jchar:jchar)
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
         posnam_=itpos
         if(lstrg.gt.0) then
           strg = buffer(jchar+1:lastch)
         endif
         if(lstrg.le.0) then
           strg=' '
           lstrg=1
         endif
         if (strg.eq.slash .and. unfold_) go to 390
         jchar=MAXBUF+1
220      cotdb_=.true.
         posstart = jchar
         recstart = irecd
         go to 300
         
!
!....... Accept the next character as a terminal delimiter
!        in a bracketed construct
!
250      lstrg=1
         quote_=' '
         posval_=jchar
         strg = buffer(jchar:jchar)
         istd =.true.
         cotdb_=.true.
         posstart = jchar
         recstart = irecd

!
!....... restore pointers and exit
!
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
!
!....... Got a comment with a folding flag
!
390      klow = 1
         lrecd=jjlrec
         cotdb_=.true.
         strg(1:1)=' '
400      jjchar = MAXBUF+1
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         lstrg=0
         go to 420
410      jchar=jchar+1
         if(jchar.le.lastch) go to 450
420      call getlin(flag)
         jchar = 1
         jjchar = 1
         if(flag.eq.'fini') then
           cotdb_=.false.
           strg='fini'
           jchar=MAXBUF+1
           lstrg=lastnb(strg)
           posstart = jchar
           recstart = irecd
           go to 300
         endif
         jchar=1
450      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 400
!
!....... Process this character in the line
!
         c=buffer(jchar:jchar)
         if(c.eq.' '.or.c.eq.tab)      goto 410
         if(c.eq.'#')       goto 470
         posstart = jchar
         recstart = irecd
         goto 300
!
!....... Accept the remainder of the line as part of the comment
!
470      lstrg=lastch-jchar
         itpos=jchar
         if(lastch.gt.jchar)                                            &
     &      strg(klow:min(len(strg),klow+lastch-2)) =                   &
     &      buffer(jchar+1:lastch)
         klow=lastnb(strg)
         if (strg(klow:klow).eq.slash) then
           strg(klow:klow)=' '
           go to 400
         endif
         jchar=MAXBUF+1
         lstrg = klow
         lrecd=jjlrec
         posstart = jchar
         recstart = irecd
         goto 300

         end
!
!
!
! >>>>>> Extract a comment or terminal delimiter field.
!
         function cotd_(strg,istd)
!
         logical   cotd_
         logical   istd
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,                                 &
     &     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   ipp,itpos
         integer   klow
         character*1 slash
!
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         istd = .false.
         slash = rsolidus(1:1)
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         strg=' '
         long_=0
         cotd_=.false.
         if (depth_.eq.0.and.jchar.gt.0) go to 105
100      jchar=jchar+1
105      if(jchar.le.lastch)     goto 140
!
!....... Read a new line
!
         call getlin(flag)
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=4
           cotd_=.false.
           return
         endif
         jchar=0
         strg=char(0)
         long_=1
         posnam_=0
         quote_=' '
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
!
!....... Process this character in the line
!
         c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         if(depth_.gt.0 .and.                                           &
     &    ((c.eq.']'.and.rdbkt_)                                        &
     &      .or.(c.eq.')'.and.rdprn_)                                   &
     &      .or.(c.eq.'}'.and.rdbrc_))) go to 250
         goto 300
!
!        For a tab, when not expanding to blanks, accept
!        that single character as a comment
!
190      long_=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
!
!....... Accept the remainder of the line as a comment
!
200      long_=lastch-jchar
         quote_=buffer(jchar:jchar)
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
         posnam_=itpos
         if(long_.gt.0) then
           strg = buffer(jchar+1:lastch)
         endif
         if(long_.le.0) then
           strg=' '
           long_=1
         endif
         if (strg.eq.slash .and. unfold_) go to 390
         jchar=MAXBUF+1
220      lrecd=jjlrec
         cotd_=.true.
         return
         
!
!....... Accept the next character as a terminal delimiter
!        in a bracketed construct
!
250      long_=1
         quote_=' '
         depth_ = depth_-1
         posval_=jchar
         strg = buffer(jchar:jchar)
         jchar=jchar+1
         lrecd=jjlrec
         istd =.true.
         cotd_=.true.
         return

!
!....... Found a non-comment field, restore pointers
!
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
!
!....... Got a comment with a folding flag
!
390      klow = 1
         lrecd=jjlrec
         cotd_=.true.
         strg(1:1)=' '
400      jjchar = MAXBUF+1
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         long_=0
         go to 420
410      jchar=jchar+1
         if(jchar.le.lastch) go to 450
420      call getlin(flag)
         jchar = 1
         jjchar = 1
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=lastnb(strg)
           return
         endif
         jchar=1
450      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 400
!
!....... Process this character in the line
!
         c=buffer(jchar:jchar)
         if(c.eq.' '.or.c.eq.tab)      goto 410
         if(c.eq.'#')       goto 470
         goto 500
!
!....... Accept the remainder of the line as part of the comment
!
470      long_=lastch-jchar
         itpos=jchar
         if(lastch.gt.jchar)                                            &
     &      strg(klow:min(len(strg),klow+lastch-2)) =                   &
     &      buffer(jchar+1:lastch)
         klow=lastnb(strg)
         if (strg(klow:klow).eq.slash) then
           strg(klow:klow)=' '
           go to 400
         endif
         jchar=MAXBUF+1
         long_ = klow
         lrecd=jjlrec
         return
!
!....... Found a non-comment field, restore pointers, but return the
!        comment found so far
!
500      jchar = jjchar
         return

         end
!         
         subroutine tbxxbtab
         include  'ciftbx.sys'
         if (jchar.gt.0 .and. jchar.le.lastch) then
           if (buffer(jchar:jchar).eq.tab                               &
     &      .and..not.tabx_)                                            &
     &      jchar=jchar-1
         end if
         return
         end
!         
         subroutine tbxxetab
         include  'ciftbx.sys'
         if (jchar.gt.1 .and. jchar.le.lastch) then
           if (buffer(jchar-1:jchar-1).eq.tab                           &
     &      .and..not.tabx_) then
            jchar = jchar-1 
            buffer(jchar:jchar) = ' '
            end if
         end if
         return
         end
!
!
!
!
!
!
! >>>>>> Extract a comment field.
!
         function cmnt_(strg)
!
         logical   cmnt_
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,                                 &
     &     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   ipp,itpos
         integer   klow
         character*1 slash
!
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         slash = rsolidus(1:1)
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         strg=' '
         long_=0
         cmnt_=.false.
         if (depth_.eq.0 .and. jchar.gt.0) go to 105
100      jchar=jchar+1
105      if(jchar.le.lastch)     goto 140
!
!....... Read a new line
!
         call getlin(flag)
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=4
           cmnt_=.false.
           return
         endif
         jchar=0
         strg=char(0)
         long_=1
         posnam_=0
         quote_=' '
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
!
!....... Process this character in the line
!
         c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         goto 300
!
!        For a tab, when not expanding to blanks, accept
!        that single character as a comment
!
190      long_=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
!
!....... Accept the remainder of the line as a comment
!
200      long_=lastch-jchar
         quote_=buffer(jchar:jchar)
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
         posnam_=itpos
         if(long_.gt.0) then
           strg = buffer(jchar+1:lastch)
         endif
         if(long_.le.0) then
           strg=' '
           long_=1
         endif
         if (strg.eq.slash .and. unfold_) go to 390
         jchar=MAXBUF+1
220      lrecd=jjlrec
         cmnt_=.true.
         return
!
!....... Found a non-comment field, restore pointers
!
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
!
!....... Got a comment with a folding flag
!
390      klow = 1
         lrecd=jjlrec
         cmnt_=.true.
         strg(1:1)=' '
400      jjchar = MAXBUF+1
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         long_=0
         go to 420
410      jchar=jchar+1
         if(jchar.le.lastch) go to 450
420      call getlin(flag)
         jchar = 1
         jjchar = 1
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=lastnb(strg)
           return
         endif
         jchar=1
450      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 400
!
!....... Process this character in the line
!
         c=buffer(jchar:jchar)
         if(c.eq.' '.or.c.eq.tab)      goto 410
         if(c.eq.'#')       goto 470
         goto 500
!
!....... Accept the remainder of the line as part of the comment
!
470      long_=lastch-jchar
         itpos=jchar
         if(lastch.gt.jchar)                                            &
     &      strg(klow:min(len(strg),klow+lastch-2)) =                   &
     &      buffer(jchar+1:lastch)
         klow=lastnb(strg)
         if (strg(klow:klow).eq.slash) then
           strg(klow:klow)=' '
           go to 400
         endif
         jchar=MAXBUF+1
         long_ = klow
         lrecd=jjlrec
         return
!
!....... Found a non-comment field, restore pointers, but return the
!        comment found so far
!
500      jchar = jjchar
         return

         end
!
!
!
!
!
!
! >>>>>> Return the delimiter prior to the most recently
!        examined value
!
         function delim_(depth,delim,posdlm,recdlm)
!
         logical   delim_
         integer   depth
         integer   posdlm
         integer   recdlm
         character*(*) delim

         include  'ciftbx.sys'
         delim = ' '
         delim_ = .false.
         posdlm = 0
         if (depth .ge.0 .and. depth .le.depth_) then
           delim = delimstack(depth+1)
           delim_ = .true.
           posdlm = posdlmstk(depth+1)
           recdlm = recdlmstk(depth+1)
         end if

         return

         end
!
!
!
!
!
! >>>>> Convert name string to lower case
!
         function tbxxlocs(name)
!
         include     'ciftbx.sys'
         character    tbxxlocs*(MAXBUF)
         character    temp*(MAXBUF),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
!
         temp=name
         kln = lastnb(name)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)
100      continue
         tbxxlocs=temp
         return
         end
!
!
!
!
!
! >>>>> Convert name string to lower case as subroutine
!
         subroutine tbxxnlc(loname, name)
!
         include     'ciftbx.sys'
         character    temp*(MAXBUF),loname*(*),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lolen,olen
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
!
         lolen = len(loname)
         olen = len(name)
         kln = min(MAXBUF,lolen,olen)
         kln = lastnb(name(1:kln))
         temp(1:kln)=name(1:kln)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)
100      continue
         loname=temp(1:kln)
         return
         end
!
!
!
!
!
! >>>>> Convert counted name string to lower case as subroutine
!       with counts
!
         subroutine tbxxclc(loname, lloname, name, lname)
!
         include     'ciftbx.sys'
         character    temp*(MAXBUF),loname*(*),name*(*)
         integer      lloname, lname
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lolen,olen
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
!
         lolen = len(loname)
         olen = min(len(name),lname)
         kln = min(MAXBUF,lolen,olen)
         kln = lastnb(name(1:kln))
         temp(1:kln)=name(1:kln)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)
100      continue
         loname(1:kln)=temp(1:kln)
         lloname = kln
         return
         end

!
!
!
!
!
! >>>>> Convert name string to upper case
!
         function tbxxupcs(name)
!
         include     'ciftbx.sys'
         character    tbxxupcs*(MAXBUF)
         character    temp*(MAXBUF),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
!
         temp=name
         kln = lastnb(name)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(low,c)
         if(j.ne.0) temp(i:i)=cap(j:j)
100      continue
         tbxxupcs=temp
         return
         end
!
!
!
!
!
! >>>>> Convert name string to upper case as subroutine
!
         subroutine tbxxnupc(upname, name)
!
         include     'ciftbx.sys'
         character    temp*(MAXBUF),upname*(*),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      olen,uplen
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
!
         uplen = len(upname)
         olen = len(name)
         kln = min(MAXBUF,uplen,olen)
         kln = lastnb(name(1:kln))
         temp(1:kln)=name(1:kln)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(low,c)
         if(j.ne.0) temp(i:i)=cap(j:j)
100      continue
         upname=temp(1:kln)
         return
         end

!
!
!
!
!
! >>>>>> Get the data item associated with the tag.
!
         subroutine tbxxgitm(name)
!
         include   'ciftbx.sys'
         SAVE
         character name*(*)
         character flag*4
         character*1 slash
         integer   iitem,nitem,npakt
         integer   kchar,loopi,i,jdict,itpos,ipp
         integer   lastnb
!
         slash = rsolidus(1:1)
!
!....... Find requested dataname in hash list
!
         lnametb=lastnb(name)
         nametb(1:lnametb)=name(1:lnametb)
!DBG     print *,' Entering tbxxgitm: ', name(1:lnametb),' ',
!DBG *     tcheck, vcheck
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         valid_ = .false.
         quote_=' '
         jdict = 0
         strg_= ' '
         long_=1
         if(name(1:1).eq.'_')       goto 100
         type_='null'
         ttype_='    '
         depth_=0
         index_=0
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         tagname_=' '
         goto 1000
100      call hash_find(nametb(1:lnametb),                              &
     &     dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,                   &
     &     iname)
         if(iname.gt.0)             goto 180
         if(dictfl.ne.'yes') then
         call hash_find(nametb(1:lnametb),                              &
     &     dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,jdict)
         if(jdict.ne.0) then
!DBG        print *,' found entry ', jdict, dicxtyp(jdict)
           dictype_=dicxtyp(jdict)
           if(dcindex(jdict).ne.0) diccat_=dcname(dcindex(jdict))
           dicname_=nametb(1:lnametb)
           if(aroot(jdict).ne.0) then
             dicname_=dictag(aroot(jdict))
             call hash_find(dicnam(aroot(jdict)),                       &
     &         dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,               &
     &         iname)
             if(iname.gt.0)      goto 180
           endif
           type_='null'
           ttype_='    '
           depth_=0
           index_=0
           tagname_=' '
           strg_=' '
           long_=1
           go to 1000
         endif
         endif

         continue
         type_='null'
         ttype_='    '
         depth_=0
         index_=0
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         long_=1
         goto 1000
!
!
180      tagname_=dtag(iname)
         if(ddict(iname).ne.0) tagname_=dictag(ddict(iname))
         posnam_=tchar(iname)
         if(tabx_)posnam_=xchar(iname)
         if(nloop(iname).le.0)      goto 500
!
!....... Process loop packet if first item request
!
         if(nloop(iname).ne.loopnl) goto 200
         if(lloop(iname).lt.loopct) goto 300
         if(loop_)                  goto 230
200      loop_=.true.
         depth_ = 0
         loopct=0
         loopnl=nloop(iname)
         nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=drecd(iname)-1
         call getlin(flag)
         jchar=max(0,dchar(iname)-1)
!DBG     if(jchar.lt.0) write(6,'(7H dchar ,i5)') jchar
         do 220 i=1,nitem
220      lloop(i+iname-iloop(iname))=0
         goto 240
!
!....... Read a packet of loop items
!
230      nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=looprd(nitem+1)-1
         call getlin(flag)
         jchar=loopch(nitem+1)
!DBG     if(jchar.lt.0) write(6,'(7H loopch,i5)') jchar
240      iitem=0
250      iitem=iitem+1
         quote_=' '
         text_=.false.
         if(iitem.le.nitem)     goto 255
         loopch(iitem)=jchar
         looprd(iitem)=irecd
         goto 270
255      call getstr
         loopch(iitem)=jchar-long_
         if(quote_.ne.' ') then
           if (quote_.eq.';') then
             loopch(iitem)=1
           else
             if (quote_.eq.''''''''.or.quote_.eq.'"""') then
               loopch(iitem)=jchar-long_-3
             else
               loopch(iitem)=jchar-long_-1
             end if
           end if
         end if
         loopln(iitem)=long_
         looprd(iitem)=irecd
         
         if (text_ .or.depth_ .gt. 0) then
         if (depth_.gt.0) then
           loopch(iitem)= posbrkstk(1)
           loopln(iitem)= 1
           looprd(iitem)= srecd
         end if
260        call getstr
           if (type_.eq.'fini') call tbxxerr(' Unexpected end of data')
           if (text_.or.depth_ .gt. 0) goto 260
         end if
         goto 250
270      loopct=loopct+1
         if(loopct.lt.npakt)    goto 300
         loop_=.false.
!
!....... Point to the loop data item
!
300      lloop(iname)=lloop(iname)+1
         loopi=iloop(iname)
         irecd=looprd(loopi)-1
         call getlin(flag)
         long_=loopln(loopi)
         kchar=loopch(loopi)
         if ((buffer(kchar:kchar).eq.'(' .and. rdprn_)                  &
     &     .or. (buffer(kchar:kchar).eq.'[' .and. rdbkt_)               &
     &     .or. (buffer(kchar:kchar).eq.'{' .and. rdbrc_)) then
           if (kchar.gt.1) then
             if (buffer(kchar-1:kchar-1).eq.''''                        &
     &          .or. buffer(kchar-1:kchar-1).eq.'"') goto 550
           end if
           jchar = kchar-1
           call getstr
!DBG         print *,' strg_ ', strg_(1:max(1,long_))
!DBG         print *,' depth_ ', depth_
           itpos=jchar-long_
           posval_=itpos
           posend_=itpos+long_-1
           jchar=kchar+long_
           if(jchar.le.MAXBUF) then
             if(buffer(jchar:jchar).ne.' ' .and.                        &
     &       buffer(jchar:jchar).ne.tab)                                &
     &       jchar=jchar+1
           endif
           if(type_.eq.'numb') then
             call ctonum
             if(posdec_.gt.0) posdec_=posval_+posdec_-1
           endif
           go to 1000
         end if
         goto 550
!
!....... Point to the non-loop data item
!
500      irecd=drecd(iname)-1
         call getlin(flag)
         kchar=dchar(iname)+1
         long_=iloop(iname)
         loop_=.false.
         loopct=0
         loopnl=0
!
!....... Place data item into variable string and make number
!
550      type_=dtype(iname)
         quote_=' '
         text_=.false.
         dictype_=dxtyp(iname)
         diccat_='(none)'
         if(cindex(iname).gt.0) diccat_=dcname(cindex(iname))
         if(cindex(iname).lt.0) diccat_=cname(-cindex(iname))
         if(diccat_.eq.' ') diccat_='(none)'
         dicname_=dtag(iname)
         if(ddict(iname).ne.0) then
           if (aroot(ddict(iname)).ne.0) then
             dicname_=dictag(aroot(ddict(iname)))
           endif
         endif
         strg_=' '
         if(long_.gt.0) then
!DBG       if (kchar.le.0)call tbxxwarn(' kchar le 0')
           strg_(1:long_)=buffer(kchar:kchar+long_-1)
           if ((buffer(kchar:kchar).eq.'('.and.rdprn_)                  &
     &       .or. (buffer(kchar:kchar).eq.'['.and.rdbkt_)               &
     &       .or. (buffer(kchar:kchar).eq.'{'.and.rdbrc_)) then
             if (kchar.gt.1) then
               if (buffer(kchar-1:kchar-1).eq.''''.or.                  &
     &           buffer(kchar-1:kchar-1).eq.'"') then
                 go to 555
               end if        
             end if
!DBG         print *,' getitm: kchar, irecd ',kchar,irecd
             jchar = kchar-1
!DBG         print *,' strg_ ', strg_(1:max(1,long_))
             call getstr
             if(type_.eq.'numb') then
               call ctonum
               if(posdec_.gt.0) posdec_=posval_+posdec_-1
             endif
!DBG         print *,' strg_ ', strg_(1:max(1,long_))
!DBG         print *,' depth_ ', depth_
             go to 1000
           end if
         endif
555      itpos=kchar
         posval_=itpos
         posend_=itpos+long_-1
         jchar=kchar+long_
         if(jchar.le.MAXBUF) then
           if(buffer(jchar:jchar).ne.' ' .and.                          &
     &       buffer(jchar:jchar).ne.tab) jchar=jchar+1
         endif
         quote_=' '
         if(kchar.gt.1) then
           if(buffer(kchar-1:kchar-1).ne.' ' .and.                      &
     &        buffer(kchar-1:kchar-1).ne.tab) then
             quote_=buffer(kchar-1:kchar-1)
             if (kchar.gt.3.and.rdtq_) then
               if (buffer(kchar-3:kchar-1).eq.                          &
     &           quote_//quote_//quote_) then
                 quote_ = buffer(kchar-3:kchar-1)
               end if
             end if
           endif
         endif
         if(type_.eq.'char' .and. kchar.eq.1 .and.                      &
     &     buffer(1:1).eq.';') then
           type_='text'
           fold_=.false.
           quote_=';'
         endif
         if(type_.eq.'text') then
           if(buffer(1:1).eq.';') then
             quote_=';'
             if (clipt_.or.long_.lt.2) then
             strg_(1:1)=' '
             if (strg_(1:long_).eq.(' '//slash) ) then
               fold_=.true.
               if(unfold_) then
                 strg_(1:long_)=slash
                   long_=1
                 endif
               endif
             else
               do ipp = 2,long_
                 strg_(ipp-1:ipp-1)=strg_(ipp:ipp)
               end do
               long_=long_-1
               if (strg_(1:long_).eq.slash) then
                 fold_=.true.
                 if (unfold_) then
                 long_=1
                 endif
               endif
             endif
           else
             type_='char'
             if (quote_.eq.';') quote_=' '
           endif
         endif
         if(type_.eq.'numb') then
           call ctonum
           if(posdec_.gt.0) posdec_=posval_+posdec_-1
         endif
         if(type_.eq.'char' .and. strg_.eq.' '.and.nblank_)             &
     &     type_='null'
         if (quote_.ne.' ') goto 1000
         if (long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
         if (long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
         if (tcheck.eq.'yes') then
           call hash_find(nametb(1:lnametb),                            &
     &       dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,jdict)
           if (jdict.gt.0)                                              &
     &     call tbxxckv(jdict)
         endif
!
1000     return
         end

!
!
!
!
!
!
!
! >>>>>> Convert string to integer, marking non-digit
!
!
         function tbxxsti(xstr,nondig)
         integer tbxxsti
         character *(*) xstr
         integer nondig, i
         integer sign, digits, kdv

         tbxxsti = 0
         digits = 0
         nondig = 0
         sign = 1
         do i = 1,len(xstr)
           kdv = ichar(xstr(i:i))-ichar('0')
           if (digits.eq.0) then
             if (xstr(i:i).eq.'-') then
               sign = -1
               digits = 1
             else
               if (xstr(i:i).eq.'+') then
                 sign = 1
                 digits = 1
               else
                 if (kdv.ge.0 .and. kdv.le.9) then
                   digits = 1
                   tbxxsti = kdv
                 else
                   if (xstr(i:i).ne.' ') then
                     nondig = i
                     return
                   endif
                 endif
               endif
             endif
           else
             if (kdv.ge.0 .and.kdv.le.9) then
               tbxxsti = tbxxsti*10+kdv
             else
               tbxxsti = sign*tbxxsti
               nondig = i
               return
             endif
           endif
         enddo
         tbxxsti = sign*tbxxsti
         return

         end



!
!
!
!
!
!
!
! >>>>>> Convert string to double, marking non-digit
!
!
         function tbxxstd(xstr,nondig)
         double precision tbxxstd
         integer tbxxsti
         character *(*) xstr
         integer nondig, i
         integer sign, digits, kdv
         integer idp, eval
         tbxxstd = 0.0
         digits = 0
         nondig = 0
         sign = 1
         idp = 0
         do i = 1,len(xstr)
           kdv = ichar(xstr(i:i))-ichar('0')
           if (i.lt.len(xstr)                                           &
     &       .and. (xstr(i:i).eq.'e'                                    &
     &       .or. xstr(i:i).eq.'E'                                      &
     &       .or. xstr(i:i).eq.'d'                                      &
     &       .or. xstr(i:i).eq.'D'                                      &
     &       .or. xstr(i:i).eq.'q'                                      &
     &       .or. xstr(i:i).eq.'Q')) then
             eval = tbxxsti(xstr(i+1:len(xstr)),nondig)
             tbxxstd = sign*tbxxstd*10.**eval
             if (nondig.ne.0) nondig=nondig+i+1
             return
           endif
           if (i.lt.len(xstr) .and. digits .ne.0                        &
     &       .and. (xstr(i:i).eq.'+'                                    &
     &       .or. xstr(i:i).eq.'-')) then
             eval = tbxxsti(xstr(i:len(xstr)),nondig)
             tbxxstd = sign*tbxxstd*10.**eval
             if (nondig.ne.0) nondig=nondig+i
             return
           endif
           if (xstr(i:i).eq.'.'.and.idp.eq.0) then
             idp = i
             digits = 1
           endif
           if (digits.eq.0) then
             if (xstr(i:i).eq.'-') then
               sign = -1
               digits = 1
             else
               if (xstr(i:i).eq.'+') then
                 sign = 1
                 digits = 1
               else
                 if (kdv.ge.0 .and. kdv.le.9) then
                   digits = 1
                   tbxxstd = kdv
                 else
                   if (xstr(i:i).ne.' ') then
                     nondig = i
                     return
                   endif
                 endif
               endif
             endif
           else
             if (kdv.ge.0 .and.kdv.le.9) then
               if (idp.eq.0) then
                 tbxxstd = tbxxstd*10.+kdv
               else
                 tbxxstd = tbxxstd+kdv*(10.**(idp-i))
               endif
             else
               if (i.ne.idp) then
                 tbxxstd = sign*tbxxstd
                 nondig = i
                 return
               endif
             endif
           endif
         enddo
         tbxxstd = sign*tbxxstd
         return

         end


!
!
!
!
!
!
!
! >>>>>> Validate the string in strg_(1:long_) of type type_
!        against the dictionary item at jdict
!
!
         subroutine tbxxckv(jdict)
         integer jdict
!
         include   'ciftbx.sys'

         character*(MAXBUF) temp, target, lcvalue
         integer tbxxfstb
         integer tbxxsti
         double precision tbxxstd
         integer tlen
         logical igood, isword, nolo, nohi
         integer fblank,ftab, symop, xlate
         integer yyyy, mm, dd, hr, mi, se, sf, tz
         integer nondig, prevdig, ldt, ldn
         logical enumflg
         integer lastnb
         integer kptr, ltarget, icptr
         integer ilolo, ilohi, ihilo, ihihi
         integer llcvalue

         valid_ = .false.
         igood = .false.
         enumflg = .false.
         kptr = 0
         if (long_ .lt. 1) return

         fblank = index(strg_(1:long_),' ')
         ftab = index(strg_(1:long_),tab)
         ldt = max(1,lastnb(dictype_))
         ldn = max(1,lastnb(dicnam(jdict)))
         isword = .true.
         if (fblank.ne.0 .or. ftab.ne.0) isword =.false.

         if (type_.eq.'null') igood = .true.

         if ((type_.eq.'char' .or. type_.eq. 'numb').and. isword) then
           if (dictype_.eq.'uchar3') then
             if (long_.eq.3.or.                                         &
     &            (long_.eq.4.and.strg_(1:1).eq.'+')                    &
     &            ) igood = .true.
             go to 90
           endif

           if (dictype_.eq.'uchar1') then
             if (long_.eq.1.or.                                         &
     &            (long_.eq.2.and.strg_(1:1).eq.'+')                    &
     &          ) igood = .true.
             go to 90
           endif

           if (dictype_(1:4).eq.'symo') then
             symop = tbxxsti(strg_(1:long_),nondig)
             xlate = 0
             if (nondig.ne.0.and.nondig.lt.long_) then
               if (strg_(nondig:nondig).eq.'_') then
                 xlate = tbxxsti(strg_(nondig+1:long_),nondig)
               endif
             endif

             if (nondig.eq.0 .and.                                      &
     &         symop .ge. 1 .and.                                       &
     &         symop .le. 192 .and.                                     &
     &         xlate .ge. 0 .and.                                       &
     &         xlate .le. 1000) igood =.true.
             go to 90
           endif

           if (dictype_(1:5).eq.'yyyy-') then

             mm=-1
             dd=-1
             hr=0
             mi =0
             se=0
             sf=0
             tz = 0

             yyyy = tbxxsti(strg_(1:long_),nondig)
             if (nondig.ne.0.and.nondig.lt.long_) then
             if (strg_(nondig:nondig).eq.'-') then
               prevdig = nondig
               mm = tbxxsti(strg_(nondig+1:long_),nondig)
               if (nondig.ne.0) nondig=prevdig+nondig
               if (nondig.ne.0.and.nondig.lt.long_) then
               if (strg_(nondig:nondig).eq.'-') then
                 prevdig = nondig
                 dd = tbxxsti(strg_(nondig+1:long_),nondig)
                 if (nondig.ne.0) nondig=prevdig+nondig
                 if (nondig.ne.0.and.nondig.lt.long_) then
                 if (strg_(nondig:nondig).eq.'T'                        &
     &             .or. strg_(nondig:nondig).eq.'t'                     &
     &             .or. strg_(nondig:nondig).eq.':') then
                   prevdig = nondig
                   hr = tbxxsti(strg_(nondig+1:long_),nondig)
                   if (nondig.ne.0) nondig=prevdig+nondig
                   if (nondig.ne.0.and.nondig.lt.long_) then
                   if (strg_(nondig:nondig).eq.':') then
                     prevdig = nondig
                     mi = tbxxsti(strg_(nondig+1:long_),nondig)
                     if (nondig.ne.0) nondig=prevdig+nondig
                     if (nondig.ne.0.and.nondig.lt.long_) then
                     if (strg_(nondig:nondig).eq.':') then
                       prevdig = nondig
                       se = tbxxsti(strg_(nondig+1:long_),nondig)
                       if (nondig.ne.0) nondig=prevdig+nondig
                       if (nondig.ne.0.and.nondig.lt.long_) then
                       if (strg_(nondig:nondig).eq.'.') then
                         prevdig = nondig
                         sf = tbxxsti(strg_(nondig+1:long_),nondig)
                         if (nondig.ne.0) nondig=prevdig+nondig
                       endif
                       endif
                     endif
                     endif
                   endif
                   endif
                 endif
                 endif
               endif
               endif
             endif
             endif
             if (nondig.ne.0) then
             if (strg_(nondig:nondig).eq.'-'                            &
     &         .or. strg_(nondig:nondig).eq.'+') then
               tz = tbxxsti(strg_(nondig+1:long_),nondig)
             endif
             endif
             if (nondig.eq.0                                            &
     &         .and. yyyy .ge. 0  .and. yyyy .lt. 10000                 &
     &         .and. mm .gt. 0 .and. mm .lt. 13                         &
     &         .and. dd .gt. 0 .and. dd .lt. 32                         &
     &         .and. hr .ge. 0 .and. hr .lt. 25                         &
     &         .and. mi .ge. 0 .and. mi .lt. 61                         &
     &         .and. se .ge. 0 .and. se .lt. 61                         &
     &         .and. sf .ge. 0                                          &
     &         .and. tz .ge. 0 .and. tz .lt. 25 ) igood =.true.
             go to 90
           endif

           if (dictype_(1:4).eq.'char'                                  &
     &       .or. dictype_(1:4).eq.'ucha'                               &
     &       .or. dictype_(1:4).eq.'code'                               &
     &       .or. dictype_(1:4).eq.'ucod'                               &
     &       .or. dictype_(1:4).eq.'line'                               &
     &       .or. dictype_(1:4).eq.'ulin'                               &
     &       .or. dictype_(1:3).eq.'any'                                &
     &       .or. dictype_(1:4).eq.'atco'                               &
     &       .or. dictype_(1:4).eq.'phon'                               &
     &       .or. dictype_(1:4).eq.'emai'                               &
     &       .or. dictype_(1:4).eq.'fax'                                &
     &       .or. dictype_(1:4).eq.'text')  then
             igood = .true.
             go to 90
           endif

           if (dictype_(1:4).eq.'numb'                                  &
     &       .or. dictype_(1:3).eq.'int'                                &
     &       .or. dictype_(1:4).eq.'floa') then
             tbxxintr = tbxxsti(strg_(1:long_),nondig)
             if (nondig.eq.0) then
               igood = .true.
               go to 90
             endif
             if (strg_(nondig:nondig).eq.'('                            &
     &         .and. nondig .lt. long_) then
               tbxxintr = tbxxsti(strg_(nondig+1:long_),nondig)
               if (nondig.gt.0) then
                 if (strg_(nondig:nondig).eq.')') then
                   igood = .true.
                   go to 90
                 endif
               endif
             endif
             if (dictype_(1:4).eq.'numb'                                &
     &         .or. dictype_(1:4).eq.'floa') then
               if (type_.eq.'numb') igood = .true.
             endif
             go to 90
           endif
           go to 90
         endif

         if (type_.eq.'char') then
           if (dictype_(1:4).eq.'text'                                  &
     &       .or. dictype_(1:3).eq.'any'                                &
     &       .or. dictype_(1:4).eq.'line'                               &
     &       .or. dictype_(1:4).eq.'ulin'                               &
     &       .or. dictype_(1:4).eq.'phon'                               &
     &       .or. dictype_(1:4).eq.'atco'                               &
     &       .or. dictype_(1:4).eq.'phon'                               &
     &       .or. dictype_(1:4).eq.'char'                               &
     &       .or. dictype_(1:4).eq.'ucha' ) igood = .true.
           go to 90
         endif

         if (type_.eq.'text') then
           if (dictype_(1:4).eq.'text'                                  &
     &       .or. dictype_(1:3).eq.'any'                                &
     &       .or. dictype_(1:4).eq.'char'                               &
     &       .or. dictype_(1:4).eq.'ucha' ) igood = .true.
           go to 90
         endif


 90      continue

         if (.not.igood) then
           call tbxxwarn(' Dictionary type '//dictype_(1:ldt)//         &
     &     ' for '//dicnam(jdict)(1:ldn)//                              &
     &     ' not matched by '//strg_(1:long_))
           return
         endif
         kptr = deindex(jdict)
         if (kptr.eq.0 .or. type_.eq.'null') then
           valid_ = .true.
           return
         endif



         call tbxxclc(lcvalue,llcvalue,strg_,long_)

100      if (kptr.ne.0) then

           tlen = tbxxfstb(temp,ivtsbp(kptr),.false.)
           if (tlen.gt.0) then
             call tbxxclc(target,ltarget,temp,tlen)
             if (ivtvet(kptr) .eq. 0) then
               enumflg = .true.
               if (target(1:ltarget).eq.lcvalue(1:llcvalue)) then
                 valid_ = .true.
                 return
               endif
               if(type_.eq.'numb'                                       &
     &           .and. (dictype_(1:4).eq.'numb'                         &
     &              .or. dictype_(1:3).eq.'int'                         &
     &              .or. dictype_(1:4).eq.'floa')) then
                 if (tbxxstd(target(1:ltarget),nondig)                  &
     &              .eq.numbtb) then
                    valid_= .true.
                    return
                 endif
               endif
             else
               enumflg = .false.
               icptr = index(target(1:ltarget),':')
               ilolo = 1
               ilohi = icptr-1
               ihilo = icptr+1
               ihihi = ltarget
               nolo = .true.
               if (ilohi.ge.ilolo) then
                 nolo = .false.
                 if (target(ilolo:ilohi).eq.'.')                        &
     &             nolo = .true.
               endif
               nohi = .true.
               if (ihihi.ge.ihilo) then
                 nohi = .false.
                 if (target(ihilo:ihihi).eq.'.')                        &
     &             nohi = .true.
               endif
              if (dictype_(1:4).eq.'numb'                               &
     &         .or. dictype_(1:3).eq.'int'                              &
     &         .or. dictype_(1:4).eq.'floa') then
              if (nolo.and.(.not.nohi)) then
                 if ((ivtvet(kptr).gt.0                                 &
     &             .and. numbtb .lt.                                    &
     &               tbxxstd(target(ihilo:ihihi),nondig)) .or.          &
     &               (ivtvet(kptr).lt.0                                 &
     &             .and. numbtb .le.                                    &
     &               tbxxstd(target(ihilo:ihihi),nondig))) then
                   valid_= .true.
                   return
                 endif
               endif
               if (nohi.and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0                                 &
     &             .and. numbtb .gt.                                    &
     &               tbxxstd(target(ilolo:ilohi),nondig)) .or.          &
     &               (ivtvet(kptr).lt.0                                 &
     &             .and. numbtb .ge.                                    &
     &               tbxxstd(target(ilolo:ilohi),nondig))) then
                   valid_= .true.
                   return
                 endif
               endif
               if ((.not.nohi).and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0                                 &
     &             .and. numbtb .lt.                                    &
     &               tbxxstd(target(ihilo:ihihi),nondig)                &
     &             .and. numbtb .gt.                                    &
     &               tbxxstd(target(ilolo:ilohi),nondig)) .or.          &
     &               (ivtvet(kptr).lt.0                                 &
     &             .and. numbtb .le.                                    &
     &               tbxxstd(target(ihilo:ihihi),nondig)                &
     &             .and. numbtb .ge.                                    &
     &               tbxxstd(target(ilolo:ilohi),nondig))) then
                   valid_= .true.
                   return
                 endif
               endif
               else
               if (nolo.and.(.not.nohi)) then
                 if ((ivtvet(kptr).gt.0                                 &
     &             .and. lcvalue(1:llcvalue) .lt.                       &
     &               target(ihilo:ihihi)) .or.                          &
     &               (ivtvet(kptr).lt.0                                 &
     &             .and. lcvalue(1:llcvalue) .le.                       &
     &               target(ihilo:ihihi))) then
                   valid_= .true.
                   return
                 endif
               endif
               if (nohi.and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0                                 &
     &             .and. lcvalue(1:llcvalue) .gt.                       &
     &               target(ilolo:ilohi)) .or.                          &
     &               (ivtvet(kptr).lt.0                                 &
     &             .and. lcvalue(1:llcvalue) .ge.                       &
     &               target(ilolo:ilohi))) then
                   valid_= .true.
                   return
                 endif
               endif
               if ((.not.nohi).and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0                                 &
     &             .and. lcvalue(1:llcvalue) .lt.                       &
     &               target(ihilo:ihihi)                                &
     &             .and. lcvalue(1:llcvalue) .gt.                       &
     &               target(ilolo:ilohi)) .or.                          &
     &               (ivtvet(kptr).lt.0                                 &
     &             .and. lcvalue(1:llcvalue) .le.                       &
     &               target(ihilo:ihihi)                                &
     &             .and. lcvalue(1:llcvalue) .ge.                       &
     &               target(ilolo:ilohi))) then
                   valid_= .true.
                   return
                 endif
               endif
              endif
             endif
           endif
           kptr = ivtnxt(kptr)
           go to 100
         endif
         continue
         if (enumflg) then
           call tbxxwarn(' Dictionary type '//dictype_(1:ldt)//         &
     &      ' for '//dicnam(jdict)(1:ldn)//', '//                       &
     &      strg_(1:long_)//                                            &
     &      ' not in dictionary list of values')
         else
         call tbxxwarn(' Dictionary type '//dictype_(1:ldt)//           &
     &      ' for '//dicnam(jdict)(1:ldn)//                             &
     &      ' range not matched by '//strg_(1:long_))
         endif

         return
         end
!
!
!
!
!
!
!
!
! >>>>>> Test for separator
!
!
       logical function tbxxtsts(c)
       include   'ciftbx.sys'
       character*1 c
       tbxxtsts = .true.
       if (rdrcqt_) return
       if (c.eq.' ') return
       if (c.eq.tab) return
       tbxxtsts = .false.
       if (depth_ .eq. 0) return
       if (c.eq.',' .or.                                                &
     &   c.eq.':' .or.                                                  &
     &   ((c.eq.')'.or.c.eq.'(') .and. rdprn_) .or.                     &
     &   ((c.eq.'}'.or.c.eq.'{') .and. rdbrc_) .or.                     &
     &   ((c.eq.']'.or.c.eq.'[') .and. rdbkt_) ) then
         tbxxtsts=.true.
         return
       end if
       return
       end
!
!
!
!
!
!
!
!
! >>>>>> Test for terminal treble quote
!
!      tests buffer(jchar:lastch) for a terminal treble quote
!      and returns the location in jtloc
!
!
       logical function tbxxtttq(jtloc)
       include   'ciftbx.sys'
       character*1 slash
       logical escaped
       logical tbxxtsts
       integer i,jtloc

       slash = rsolidus(1:1)
       tbxxtttq = .false.
       jtloc = jchar

       if (rdrcqt_) then

!        Process according to CIF2 rules on quotes
         escaped = .false.
         if (jchar .le. lastch-2) then
           do i = jchar,lastch-2
             if (.not.escaped) then
               if (buffer(i:i).eq.slash) then
                 escaped = .true.
               else
                 if (buffer(i:i+2).eq.quote_) then
                   jtloc = i
                   tbxxtttq = .true.
                   return
                 end if
               end if
             else
               escaped=.false.
             end if
           end do
         end if
         return
       else
!
!        Process according to CIF1 rules
         if (jchar .le. lastch-2) then
           do i = jchar,lastch-2
             if (buffer(i:i+2).eq.quote_) then
               if (i.lt.lastch-2) then
                 if (tbxxtsts(buffer(i+3:i+3))) then
                   jtloc = i
                   tbxxtttq = .true.
                   return
                 end if
               else
                 jtloc = i
                 tbxxtttq = .true.
                 return
               end if
             end if
           end do
         end if
         return
       end if
         end


!
!
!
!
!
!
!
! >>>>>> Read the next string from the file
!
!
       subroutine getstr
!
!      On entry, jchar is set to one less than the next character
!      to be read, on the line given by irecd, which is assumed
!      to have been loaded into buffer, with lastch set to the
!      position of the last character
!
!      if depth_ is greater than 0, then statestack(depth_),
!      brackstack(depth_) and indexstack(depth_) give the
!      state of the scan within a list, array, tuple of table
!
!      If the state is 0, we are starting a search for a token
!      if the state is 2, we had a token last pass and need
!      to find a comma, colon or a terminating ) ] or } before looking
!      for the next token.
!
!      On entry, the state of text_ is used to determine if we are continuing
!      a text field or a bracketed construct.  The case of a text field within
!      a bracketed construct is handled by checking quote_.  If it is any of
!      ';', "'''", or '"""', then the next line needs to be read if text_ is
!      true.  If it any any other string, then the next element of the
!      bracketed construct needs to be read
!
!      If the depth is zero, text_ will be cleared one getstr call
!      prior to the last, empty read.  This change is not made for treble
!      quoted strings
!
!      In a bracketed construct text_ will be left set after the read
!      of the last read with text and the next read will return a null
!      type_
!
       include   'ciftbx.sys'
       integer   i,j,jj(11),im,ip
       integer   state,innerdepth
       logical   quoted
       logical   escaped
       character c*1,num*21,flag*4
       character slash*1
       logical   tbxxtsts, tbxxtttq
       integer   jtloc

       data num/'0123456789+-.()EDQedq'/
       slash = rsolidus(1:1)
       im = 0
!DBG   print *,' entering getstr type, text_, quote_, depth_ jchar: '
!DBG   print *,  type_, text_, quote_, depth_, jchar


!      If text_ is true, we may be continuing a text field, a
!      treble-quoted string or a bracketed construct
!
!      We deal with the first 2 cases here
!
       if (text_) then
         if (quote_.eq.';'.or.quote_.eq.'"""'.or.quote_.eq."'''") then
!DBG     print *,' processing next line'
           call getlin(flag)
           if (flag.eq.'fini') then
             type_='fini'
             text_=.false.
             depth_=0
             ttype_='    '
             depth_=0
             quote_=' '
             quoted=.false.
             goto 500
           end if
!
!          Handle the case of a text field
!          This is terminated by \n; which
!          is unconditionally recognized if rdrcqt_ is true
!          or which requires a trailing separator if depth_
!          is greater than 0

!DBG   print *, 'read line:',buffer(1:lastch)
           if (quote_.eq.';') then
!DBG     print *, 'semicolon quote_ detected'
             if (buffer(1:1).eq.';') then
!DBG         print *, 'terminal semicolon detected'
               if (lastch.gt.1) then
                 if (.not.tbxxtsts(buffer(2:2))) goto 10
               end if
               jchar = 2
               long_ = 0
               strg_(1:1) = ' '
               text_=.false.
               quote_ = ' '
               type_='null'
               goto 500
             end if
!            Here is the line we have read is part of the text field
10           continue
!DBG         print *,'processing as text field'
             jchar = lastch+1
             strg_(1:lastch) = buffer(1:lastch)
!DBG         print *, ' buffer before backup ', buffer(1:lastch)
             long_ = lastch
             text_ = .true.
             
             if (depth_.eq.0) then
               call getlin(flag)
               if(flag.eq.'fini') then
                 text_ = .false.
               else
                 if (buffer(1:1).eq.';') then
                   text_ = .false.
                   jchar = 2
                 end if
               end if
               if (text_) then
                 irecd = irecd-1
                 jchar=MAXBUF+1
               else
                 go to 500
               endif
!DBG         print *, ' buffer after backup ', buffer(1:lastch)
             end if
             goto 500
           else
!
!          Handle the case of a treble-quoted string
!          This is terminated by an unescaped quote
!
             if (tbxxtttq(jtloc)) then
               long_ = jtloc-jchar
               if (long_.eq.0) then
                 long_=0
                 strg_(1:1) = ' '
               else
                 strg_(1:long_)=buffer(jchar:jtloc)
               end if
               jchar = jtloc+3
               text_=.false.
               goto 500
             else
               jchar = lastch+1
               strg_(1:lastch) = buffer(1:lastch)
               long_ = lastch
               text_ = .true.
               goto 500
             end if
           end if
         end if
       end if
!
!      Now we are sure we are done with the multiline quoted cases
       quoted=.false.
       quote_=' '
       ttype_=' '
       if (depth_ .gt. 0) then
         ttype_ = typestack(1)
         type_ = typestack(depth_)
         state = statestack(depth_)
         index_= indexstack(depth_)
         go to 3000
       end if

!      We are not in a bracketed construct

       if(irecd.gt.0.and.                                               &
     &   jchar.le.1.and.lastch.gt.0) then
         jchar=1
         goto 140
       end if
100    jchar=jchar+1
       if(jchar.le.lastch)     goto 150
!
!....... Read a new line
!
110    call getlin(flag)
       type_='fini'
       dictype_=type_
       diccat_='(none)'
       dicname_=' '
!DBG   write(6,'(/5i5,a)')
!DBG *   irecd,jrecd,lrecd,nrecd,lastch, buffer(1:lastch)
       if(flag.eq.'fini')  goto 500
!
!....... Test if the new line is the start of a text sequence
!
140    if(buffer(1:1).ne.';') goto 150
       type_='text'
       quote_=';'
       jchar=lastch+1
       long_=lastch
       if (clipt_) then
         strg_(1:long_)=buffer(1:long_)
         strg_(1:1)=' '
       else
         if (long_.eq.1) then
           strg_(1:long_) = ' '
           long_ = 0
         else
           long_ = long_-1
           strg_(1:long_)=buffer(2:long_+1)
         endif
       endif
       text_ = .true.
       goto 500
!
!....... Process this character in the line
!
150    c=buffer(jchar:jchar)
       ip = jchar
       if(c.eq.' ')       goto 100
       if(c.eq.tab)       goto 100
       if(c.eq.'#')       goto 110
       if(c.eq.'''')      goto 300
       if(c.eq.'"')       goto 300
       if(c.eq.'('.and.rdprn_)       goto 350
       if(c.eq.'['.and.rdbkt_)       goto 360
       if(c.eq.'{'.and.rdbrc_)       goto 370
       if(c.ne.'_')       goto 200

       type_='name'
       goto 210
!
!....... Span blank delimited token; test if a number or a character
!
200    type_='numb'
       im=0
       quoted=.false.
       quote_=' '
       do 205 i=1,11
205    jj(i)=0
210    ip = jchar
       do 250 i=jchar,lastch
         ip = i
         if(buffer(i:i).eq.' ')       goto 400
         if(buffer(i:i).eq.tab)       goto 400
         if(type_.ne.'numb')          goto 250
         j=index(num,buffer(i:i))
         if(j.eq.0)                 type_='char'
         if(j.le.10) then
           im=im+1
           goto 250
         endif
         if(j.gt.13.and.im.eq.0) type_='char'
         jj(j-10)=jj(j-10)+1
250    continue
       i=lastch+1
       ip = i
       if(type_.ne.'numb') goto 400
       do 270 j=1,5
         if((jj(j).gt.1.and.j.gt.2) .or.                                &
     &     jj(j).gt.2)             type_='char'
270    continue
       goto 400       
!
!....... Span quote delimited token; assume character
!
300    type_='char'
       quoted=.true.
       jchar=jchar+1
       if (rdtq_ .and. jchar+1 .le. lastch                              &
     &   .and. buffer(jchar:jchar+1).eq.c//c) then
         quote_ = c//c//c
         jchar = jchar+2
         if (tbxxtttq(jtloc)) then
           text_=.false.
         else
           jtloc = lastch
           text_=.true.
           type_='text'
         endif
         long_ = jtloc-jchar
         if (long_.eq.0) then
           strg_=' '
         else
           strg_(1:long_)=buffer(jchar:jtloc)
         endif
         goto 500
       end if
       escaped = .false.
       ip = jchar
       do 320 i=jchar,lastch
         ip = i
         if (.not.escaped.and.rdrcqt_) then
           if (c.eq.slash) then
             escaped = .true.
             go to 320
           end if
         end if
         if (escaped) then
           escaped = .false.
           go to 320
         end if
         if(buffer(i:i).ne.c)             goto 320
         if(rdrcqt_)                      goto 400
         if(i+1.ge.lastch)                goto 400
         if (tbxxtsts(buffer(i+1:i+1)))   goto 400
320    continue
!DBG   write(6,'(a,4i5,a)')
!DBG *     '**** ',irecd,lastch,i,jchar,buffer(jchar:i)
       call tbxxwarn(' Quoted string not closed')
       i = lastch+1
       goto 400

!
!...... Here to start a bracketed construct
!
350    type_='tupl'
       go to 390
360    type_='list'
       go to 390
370    type_='tabl'
       if (.not. rdbkt_) type_='list'
390    depth_=1
       srecd=irecd
       ttype_=type_
       typestack(depth_) = type_
       brackstack(depth_) = c
       posbrkstk(depth_) = ip
       if (c.eq.':') brackstack(depth_) = ' '
       delimstack(depth_+1) = c
       posdlmstk(depth_+1) = ip
       recdlmstk(depth_+1) = irecd
       indexstack(depth_) = 1
       statestack(depth_) = 0
       state = 0
       go to 3100

!
!....... Here to process within a bracketed construct
!
3000   continue
!DBG   print *,' Processing in backeted construct '
       if(irecd.gt.0.and.                                               &
     &   jchar.le.1.and.lastch.gt.0) then
         jchar=1
         goto 3140
       endif

3100   jchar = jchar+1
       if(jchar.le.lastch)     goto 3150
!
!....... Read a new line
!
3110   call getlin(flag)
       if(flag.eq.'fini')  goto 3500

!
!....... Test if the new line is the start of a text sequence
!
3140   continue
!DBG   print *,buffer(jchar:lastch)
       if (buffer(1:1).ne.';') goto 3150
         type_='text'
         quote_=';'
         jchar=lastch+1
         long_=lastch
         if (long_ .gt. 1) then
         if (clipt_) then
           strg_(2:long_)=buffer(2:long_)
           strg_(1:1) = ' '
         else
           strg_(1:long_-1)=buffer(2:long_)
         long_ = long_-1
         endif
       else
         strg_(1:1)=' '
         long_ = 0
       endif
       state = 2
       statestack(depth_) = 2
       goto 500
!
!..... Process this character in the line
!      within a bracket construct
!
3150   continue
!DBG   print *,buffer(jchar:lastch)
       c=buffer(jchar:jchar)
       ip = jchar
       if(c.eq.' ')       goto 3100
       if(c.eq.tab)       goto 3100
       if(c.eq.'#')       goto 3110
       if(c.eq.'''')      goto 3300
       if(c.eq.'"')       goto 3300
       if(c.eq.'('.and.rdprn_)       goto 3350
       if(c.eq.'['.and.rdbkt_)       goto 3360
       if(c.eq.'{'.and.rdbrc_)       goto 3370
       if(c.eq.':'.and.rdcolon_)     goto 3380
       if(c.eq.',')       goto 3160
       if((c.eq.')' .and. brackstack(depth_).eq.'(') .or.               &
     &   (c.eq.'}' .and. brackstack(depth_).eq.'{') .or.                &
     &   (c.eq.']' .and. brackstack(depth_).eq.'[')) goto 3160
       if(c.eq.')' .and.                                                &
     &   ((brackstack(depth_).eq.'{').or.                               &
     &   (brackstack(depth_).eq.'[')))                                  &
     &     call tbxxwarn(' Unbalanced ) treated as comma')
       if(c.eq.']' .and.                                                &
     &   ((brackstack(depth_).eq.'{').or.                               &
     &   (brackstack(depth_).eq.'(')))                                  &
     &     call tbxxwarn(' Unbalanced ] treated as comma')
       if(c.eq.'}' .and.                                                &
     &   ((brackstack(depth_).eq.'(').or.                               &
     &   (brackstack(depth_).eq.'[')))                                  &
     &     call tbxxwarn(' Unbalanced } treated as comma')
       go to 3200
!
!..... Process comma or close bracket found within bracketed construct
!
3160   continue
       if (state .eq. 2) then
         state = 0
         if (c.eq.',') then
           index_ = index_+1
           indexstack(depth_) = index_
           delimstack(depth_+1) = c
           posdlmstk(depth_+1) = ip
           recdlmstk(depth_+1) = irecd
           goto 3100
         endif
       endif
       depth_ = depth_-1
!DBG   print *,' decreasing depth ',depth_, recn_, 
!DBG *  buffer(jchar:lastch)

       type_='null'
       long_ = 1
       strg_(1:long_) = ' '
       if (depth_ .eq.0 ) goto 500
       state = 2
       statestack(depth_) = 2
       delimstack(depth_+1) = c
       posdlmstk(depth_+1) = ip
       recdlmstk(depth_+1) = irecd
       go to 500
!
!..... Process colon found within bracketed construct
!      treat as a comma if already started
!
3380   continue
       if (state .eq. 2) then
         state = 0
         index_ = index_+1
         indexstack(depth_) = index_
         delimstack(depth_) = c
         posdlmstk(depth_+1) = ip
         recdlmstk(depth_+1) = irecd
         goto 3100
       endif
       type_='null'
       long_ = 1
       strg_(1:long_) = ' '
       state = 2
       statestack(depth_) = 2
       delimstack(depth_) = c
       go to 500


!
!..... Span blank delimited token; test if a number or a character
!
3200   type_='numb'
       im=0
       innerdepth = depth_
       do 3205 i=1,11
3205     jj(i)=0
       ip = jchar
       do 3250 i=jchar,lastch
       ip = i
       if((buffer(i:i).eq.'('.and.rdprn_)                               &
     &   .or.(buffer(i:i).eq.'{'.and.rdbrc_)                            &
     &   .or.(buffer(i:i).eq.'['.and.rdbkt_)) then
         if (depth_ .ge. MAXDEPTH) then
           call tbxxerr(' Stack overflow, increase MAXDEPTH')
         end if
         depth_=depth_+1
!DBG   print *,' increasing depth ',depth_, recn_, 
!DBG *   buffer(jchar:lastch)
         typestack(depth_) = type_
         indexstack(depth_) = 1
         brackstack(depth_) = buffer(i:i)
       endif
       if((buffer(i:i).eq.')'.and.brackstack(depth_).eq.'(')            &
     &   .or.(buffer(i:i).eq.'}'.and.brackstack(depth_).eq.'{')         &
     &   .or.(buffer(i:i).eq.']'.and.brackstack(depth_).eq.'[')         &
     &   .or.(buffer(i:i).eq.' '.and.brackstack(depth_).eq.' ')) then
         if (depth_.gt.innerdepth) then
           depth_=innerdepth
           call tbxxwarn(                                               &
     & ' Failed to balance brackets in blank-delimited token')
!DBG   print *,' decreasing depth ',depth_, recn_, 
!DBG *  buffer(jchar:lastch)
         else
           go to 3395
         endif
       endif
       if(buffer(i:i).eq.' ')       goto 3400
       if(buffer(i:i).eq.tab)       goto 3400
       if(buffer(i:i).eq.',')       goto 3390
       if(buffer(i:i).eq.':'.and.rdcolon_) go to 3390
       if(buffer(i:i).eq.brackstack(depth_)) goto 3400
       if(type_.ne.'numb')          goto 3250
       j=index(num,buffer(i:i))
       if(j.eq.0)                 type_='char'
       if(j.le.10) then
         im=im+1
         goto 3250
       endif
       if(j.gt.13.and.im.eq.0) type_='char'
       jj(j-10)=jj(j-10)+1
3250   continue
       i=lastch+1
       ip = i
       if(type_.ne.'numb') goto 3400
       do 3270 j=1,5
         if((jj(j).gt.1.and.j.gt.2) .or.                                &
     &     jj(j).gt.2)             type_='char'
3270   continue
       go to 3400

!
!..... Span '\'' or '\"' quote delimited token; assume character
!
3300   type_='char'
       quoted=.true.
       jchar=jchar+1

       if (rdtq_ .and. jchar+1 .le. lastch                              &
     &   .and. buffer(jchar:jchar+1).eq.c//c) then
         quote_ = c//c//c
         jchar = jchar+2
         if (tbxxtttq(jtloc)) then
           text_=.false.
         else
           jtloc = lastch
           text_=.true.
           type_='text'
         endif
         long_ = jtloc-jchar
         if (long_.eq.0) then
           strg_=' '
         else
           strg_(1:long_)=buffer(jchar:jtloc)
         endif
         goto 500
       end if

       escaped = .false.
      
!DBG   print *,'Processing quoted string '
!DBG   print *,buffer(jchar:lastch)

       ip = jchar
       do 3320 i=jchar,lastch
         ip = i
         if (.not.escaped.and.rdrcqt_) then
           if (c.eq.slash) then
             escaped = .true.
             go to 3320
           end if
         end if
         if (escaped) then
           escaped = .false.
           go to 3320
         end if
!DBG     print *,'i,c,buffer(i:i)',i,c,buffer(i:i) 

         if(rdrcqt_)                      goto 3400
         if(buffer(i:i).ne.c)             goto 3320
         if(i+1.gt.lastch)                goto 3400
         if (tbxxtsts(buffer(i+1:i+1)))   goto 3400
3320    continue
       go to 3400
!.....  Span (-delimited tuple
3350    type_ = 'tupl'
       go to 3375
!.....  Span [-delimited list or array
3360    type_ = 'list'
       go to 3375
!.....  Span { delimited table
3370    type_ = 'tabl'
        if (.not. rdbkt_) type_='list'
3375    continue
       if (depth_ .ge. MAXDEPTH) then
         call tbxxerr(' Stack overflow, increase MAXDEPTH')
       end if
       depth_=depth_+1
!DBG   print *,' increasing depth ',depth_, recn_, 
!DBG *  buffer(jchar:lastch)
       typestack(depth_) = type_
       indexstack(depth_) = 1
       brackstack(depth_) = buffer(ip:ip)
       state = 0
       statestack(depth_) = state
       go to 3100


3390    if(depth_ .ne. innerdepth) then
         call tbxxwarn(' failed to close bracketed string')
         depth_ = innerdepth
        endif
        go to 3400


3395   ip = ip-1

3400   state = 2
       statestack(depth_) = state
       go to 400

3500   type_='fini'
       dictype_=type_
       diccat_='(none)'
       dicname_=' '
       if (depth_ .gt.0) then
         call tbxxwarn(                                                 &
     &     ' File ended in unterminated bracketed construct')
         depth_ = 0
       endif
       go to 500



!
!..... Store the string for the getter
!
400    long_=0
       strg_=' '
         if(ip.gt.jchar) then
         long_=ip-jchar
         strg_(1:long_)=buffer(jchar:ip-1)
       endif
       jchar=ip
       quote_=' '
       if(quoted) then
         quote_=buffer(jchar:jchar)
         if (depth_.eq.0) jchar =jchar+1
       endif
       if(type_.ne.'char'.or.quoted.or.depth_.gt.0) goto 500
       if(strg_(1:5).eq.'data_') then
          type_='data'
          depth_=0
       end if
       if(strg_(1:5).eq.'loop_') then
         type_='loop'
         depth_=0
       end if
!DBG   if (strg_(1:max(1,long_)).eq.'?') print *,long_,strg_(1:1)
       if(long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
       if(long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
       if(strg_(1:5).eq.'save_') then
         type_='save'
         depth_=0
       end if
       if(long_.eq.7.and. strg_(1:7).eq.'global_') then
         type_='glob'
         depth_=0
       end if
!
500    continue
!DBG   print *,' leaving getstr with strg: ', strg_(1:long_)
!DBG   print *,' leaving getstr type, text_, quote_, depth_, jchar: '
!DBG   print *, type_,', ',text_,', ',quote_,', ',depth_,', ',jchar
       return
       end
!
!
!
!
!
!
! >>>>>> Convert a character string into a number and its esd
!
!                                          Q
!                                          D+
!                                          E-
!                                +         +
!           number string        -xxxx.xxxx-xxx(x)
!           component count CCNT 11111222223333444
!           (with at least 1 digit in the mantissa)
!
         subroutine ctonum
!
         integer   lastnb
         include  'ciftbx.sys'
         character test*26,c*1
         integer*4 m,nchar
         integer*4 ccnt,expn,msin,esin,ndec,ids,nmd
         integer*4 nms,ned,nef,nes
         double precision numb,sdev,ntemp,mant
         data test /'0123456789+.-()EDQedq :,]}'/
!
         numbtb=0.D0
         sdevtb=-1.D0
         numb=1.D0
         sdev=0.D0
         ccnt=0
         mant=0.D0
         expn=0.
         msin=+1
         esin=+1
         ndec=0
         ids=0
         nmd=0
         nms=0
         ned=0
         nef=0
         nes=0
         type_='char'
         posdec_=0
         esddig_=0
         if(long_.eq.1.and.                                             &
     &     index('0123456789',strg_(1:1)).eq.0) goto 500
         lzero_=.false.
         decp_=.false.
!
!....... Loop over the string and identify components
!
!        The scan works in phases
!          ccnt = 0   processing looking for first digit
!          ccnt = 1   processing before decimal point
!          ccnt = 2   processing after decimal point
!          ccnt = 3   processing exponent
!          ccnt = 4   processing standard deviation
!
         do 400 nchar=1,long_
!
         c=strg_(nchar:nchar)
         m=index(test,c)
         if(m.eq.0)     goto 500
         if(m.gt.10)    goto 300
!
!....... Process the digits
!
         if(ccnt.eq.0)  ccnt=1
         if(ccnt.eq.2)  ndec=ndec+1
         if(ccnt.gt.2)  goto 220
         ntemp=m-1
         if (ndec.eq.0) then
           mant=mant*10.D0+ntemp
         else
           mant=mant+ntemp/10.D0**(ndec)
         endif
         nmd=nmd+1
         if(ccnt.eq.1.and.mant.ne.0.D0) ids=ids+1
         goto 400
220      if(ccnt.gt.3)  goto 240
         expn=expn*10+m-1
         goto 400
240      esddig_=esddig_+1
         ntemp=m-1
         sdev=sdev*10.D0+ntemp
         sdevtb=1.D0
         goto 400
!
!....... Process the characters    . + - ( ) E D Q
!
300      if(c.ne.'.')  goto 320
         decp_=.true.
         if(nchar.gt.1.and.mant.eq.0.d0) then
           if(strg_(nchar-1:nchar-1).eq.'0') lzero_=.true.
         endif
         if(ccnt.gt.1) goto 500
         posdec_=nchar
         ccnt=2
         goto 400
!
320      if(nmd.eq.0.and.m.gt.13) goto 500
         if(c.ne.'(')  goto 340
         if(posdec_.eq.0) posdec_=nchar
         ccnt=4
         goto 400
!
340      if(posdec_.eq.0.and.ccnt.gt.0) posdec_=nchar
         if(c.eq.')' .or. c.eq.' ')  goto 400
         if(ccnt.eq.3 .and. ned.gt.0) goto 500
         if(m.gt.13) then
           if (nef.gt.0) goto 500
           nef = nef+1
           ccnt = 3
           esin = 1
         else
           if(ccnt.gt.0) then
             if (nes.gt.0) goto 500
             nes = nes+1
             ccnt = 3
             esin = 12-m
           else
             if (nms.gt.0) goto 500
             nms = nms+1
             ccnt=1
             msin=12-m
           endif
         endif
!
400      continue
!
         if(posdec_.eq.0) posdec_=lastnb(strg_(1:long_))+1
!
!....... String parsed; construct the numbers
!
         expn=expn*esin
         if(expn+ids.gt.-minexp) then
           call tbxxwarn(' Exponent overflow in numeric input')
           expn=-minexp-ids
         endif
         if(expn.lt.minexp) then
           call tbxxwarn(' Exponent underflow in numeric input')
           expn=minexp
         endif
         if(expn-ndec.lt.0) numb=1./10.D0**abs(expn-ndec)
         if(expn-ndec.gt.0) numb=10.D0**(expn-ndec)
         if(sdevtb.gt.0.0) sdevtb=numb*sdev
         numb=1.D0
         if(expn.lt.0) numb=1./10.D0**abs(expn)
         if(expn.gt.0) numb=10.D0**(expn)
         ntemp=msin
         numbtb=numb*mant*ntemp
         type_='numb'
!
500      return
         end
!
!
!
!
!
!
! >>>>>> Read a new line from the direct access file
!
         subroutine getlin(flag)
!
         include  'ciftbx.sys'
         character flag*4
         integer   kpp,lpp,mpp,npp,ir
         integer   tbxxrld
         integer   lip,mp,kip,ip,mip,mis
         integer   icpos,itpos,ixpos,ixtpos
!
         irecd=irecd+1
         jchar=1
         kpp = 0
         if(irecd.eq.jrecd.and.                                         &
     &     irecd.gt.recbeg_.and.                                        &
     &     irecd.le.recend_)  goto 200
         if(irecd.le.min(lrecd,recend_))  goto 100
         irecd=min(lrecd,recend_)+1
         buffer(1:1)=' '
         lastch=0
         jchar=MAXBUF+1
         jrecd=-1
         flag='fini'
         goto 200
100      continue
         lpp=-1
         mpp=-1
         npp=kpp
         call tbxxflin(irecd,lip,kpp,mp,kip,ip,mip,mis)
         if (lip.eq.0) then
           buffer(1:1) = ' '
           lastch = 1
           go to 130
         endif
         do ir = 1,NUMPAGE
           if(iabs(mppoint(ir)).eq.kpp) then
             lpp = ir
             goto 120
           endif
           if(mppoint(ir).eq.0) then
             lpp=ir
           else
             if(iabs(iabs(mppoint(ir))-kpp)                             &
     &         .gt.iabs(npp-kpp)) then
               mpp=ir
               npp=iabs(mppoint(ir))
             endif
           endif
         enddo
!
!        failed to find page as resident
!        remove a target page
!
         if(lpp.eq.-1)lpp=mpp
         if(lpp.eq.-1)lpp=1
         if (mppoint(lpp).lt.0) then
           write(dirdev,'(a)',rec=-mppoint(lpp)) pagebuf(lpp)
         endif
         mppoint(lpp)=kpp
         read(dirdev,'(a)',rec=kpp) pagebuf(lpp)
120      lastch = tbxxrld(buffer,pagebuf(lpp)(mp:NUMCPP), .false.)
130      recn_=irecd
         jrecd=irecd
         flag=' '
         if (lastch.gt.0 .and. tabx_) then
           icpos=1
           itpos=1

140        ixpos=index(buffer(icpos:lastch),tab)
           ixtpos=ixpos+itpos-1
           if(ixpos.gt.0.and.ixtpos.le.MAXBUF) then
             ixtpos=((ixtpos+7)/8)*8
             if(ixpos.gt.1) then
               bufntb(itpos:ixtpos)=                                    &
     &           buffer(icpos:ixpos+icpos-2)
             else
             bufntb(itpos:ixtpos)=' '
             endif
             itpos=ixtpos+1
             icpos=ixpos+icpos
             goto 140
           else
             bufntb(itpos:min(MAXBUF,itpos+lastch-icpos))=              &
     &         buffer(icpos:lastch)
           endif
           buffer(1:min(MAXBUF,itpos+lastch-icpos))=                    &
     &        bufntb(1:min(MAXBUF,itpos+lastch-icpos))
           lastch = min(MAXBUF,itpos+lastch-icpos)
         endif
200      return
         end
!
!
!
!
!
!
! >>>>>> Write error message and exit.
!
         subroutine tbxxerr(mess)
         character*(*) mess
         call tbxxcmsg('error',mess)
         stop
         end
!
!
!
!
!
!
! >>>>>> Write warning message and continue.
!
         subroutine tbxxwarn(mess)
         character*(*) mess
         call tbxxcmsg('warning',mess)
         return
         end
!
!
!
!
!
!
! >>>>>> Write a message to the error device
!
         subroutine tbxxcmsg(flag,mess)
!
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) flag
         character*(*) mess
         character*(MAXBUF)  tline
         character*5   btype
         integer       ll,ls,ltry,ii,i
!
         btype = 'data_'
         if(save_) btype = 'save_'
         if(.not.glob_) then
         tline= ' ciftbx '//flag//': '                                  &
     &   //file_(1:longf_)//' '//btype                                  &
     &   //bloc_(1:max(1,lastnb(bloc_)))//' line:'
         else
         tline= ' ciftbx '//flag//': '                                  &
     &   //file_(1:longf_)//' global_'//' line:'
         endif
         ll = max(1,lastnb(tline))
         write(errdev,'(a,i7)')tline(1:ll),irecd
         ll=len(mess)
         ls=1
100      if(ll-ls.le.79) then
           write(errdev,'(1X,a)') mess(ls:ll)
           return
         else
           ltry = min(ll,ls+79)
           do ii = ls+1,ltry
           i = ltry-ii+ls+1
           if(mess(i:i).eq.' ') then
             write(errdev,'(1X,a)') mess(ls:i-1)
             ls=i+1
             if(ls.le.ll) go to 100
             return
           endif
           enddo
           write(errdev,'(1X,a)') mess(ls:ltry)
           ls=ltry+1
           if(ls.le.ll) go to 100
           return
         endif
         end
!
!
!
!
! >>>>>> Create a named file.
!
         function pfile_(fname)
!
         logical   pfile_
         include   'ciftbx.sys'
         logical   test
         integer   lfname
         integer   i
         character fname*(*)
!
!....... Test if a file by this name is already open.
!
         if(pfilef.eq.'yes') call close_
         pfilef='no '
         file_(1:longf_) = ' '
         lfname = len(fname)
         file_(1:lfname)=fname
         do 120 i=1,lfname
         if(file_(i:i).eq.' ') goto 140
120      continue
         i = lfname+1
140      if (i.gt.1) then
           inquire(file=file_(1:i-1),exist=test)
           pfile_=.false.
           longf_ = i-1
           if(test)            goto 200
         else
           file_ = ' '
           pfile_ = .true.
           longf_ = 1
         endif
!
!....... Open up a new CIF
!
         if (file_(1:1) .ne. ' ')  then
         open(unit=outdev,file=file_(1:longf_),status='NEW',            &
     &                    access='SEQUENTIAL',                          &
     &                    form='FORMATTED')
         precn_=0
         endif
         pfile_=.true.
         pfilef='yes'
         nbloc=0
         pchar=1+lprefx
         pcharl=0
         obuf=prefx
         obuf(pchar:MAXBUF)=' '
200      ploopn = 0
         ploopc = 0
         ploopf = 'no '
         ptextf = 'no '
         pdepth_ = 0
         pdelimstack(1) = ' '
         pposdlmstk(1) = 0
         plcat = ' '
         pdblok = ' '
         plhead(1) = ' '
         if (xmlout_) then
           call tbxxpstr('<?xml version="1.0"?>')
         endif
         return
         end


!
!
!
!
!
! <<<<<< Substitute item in data block XML translation
!
         function tbxxxsub(oblok,xstring)
         include    'ciftbx.sys'
         character   oblok*(*)
         character   xstring*(*)
         character   tbxxxsub*(MAXBUF)
         integer     ii, jj, kk
         integer     lastnb
         jj = 1
         tbxxxsub = ' '
         do ii = 1,lastnb(xstring)
           if(xstring(ii:ii).ne.'%') then
             tbxxxsub(jj:jj) = xstring(ii:ii)
             jj = jj+1
           else
             do kk = 1,lastnb(oblok)
               tbxxxsub(jj:jj) = oblok(kk:kk)
               jj = jj+1
             enddo
           endif
         enddo
         return
         end

!
!
!
!
!
! >>>>>> Store a data block command in the CIF
!        Call with blank name to close current block only
!
         function pdata_(name)
!
         logical   pdata_
         SAVE
         include  'ciftbx.sys'
         character name*(*),temp*(MAXBUF)
         character dbloc(100)*(NUMCHAR)
         character tbxxxsub*(MAXBUF)
         integer   i
         integer   lastnb
!
         pdata_=.true.
         if(ptextf.eq.'yes') call tbxxeot
         if(pdepth_ .gt.0)   call tbxxebkt
         if(ploopn.ne.0)     call tbxxelp

         if(psaveo) then
           pchar=-1
           if(pposval_.ne.0) then
             pchar=lprefx+1
             call tbxxpstr(' ')
             pchar=lprefx+pposval_
             pposval_=0
           endif
           if (xmlout_) then
             call tbxxpxct('save_',' ')
           else
             call tbxxpstr('save_')
           endif
           psaveo=.false.
         endif
         if (pdblok(1:1).ne.' ') then
           if (xmlout_) then
             if (xmdata.eq.0) then
               call tbxxpxct(pdblok,' ')
             else
               call tbxxpxct(tbxxxsub(pdblok,xmlate(xmdata)),' ')
             endif
           endif
           pdblok=' '
         endif
         if(globo_) then
           pchar=-1
           temp='global_'
           pdblok='global_'
           psaveo=.false.
           goto 135
         endif
!
!....... Check for duplicate data name
!
         temp=name
         if(temp.eq.' ')        goto 200
         if(saveo_)             goto 130
         pdata_=.false.
         do 110 i=1,nbloc
         if(temp.eq.dbloc(i))   goto 130
110      continue
         pdata_ = .true.
         goto 125
!
!....... Save block name and put data_ statement
!
125      nbloc=nbloc+1
         if(nbloc.le.100) dbloc(nbloc)=temp(1:min(NUMCHAR,MAXBUF))
         pdblok = temp(1:min(NUMCHAR,MAXBUF))
130      pchar=-1
         temp='data_'//name
         if(saveo_) temp='save_'//name
         if(globo_) temp='global_'
         psaveo=saveo_
135      if(pposnam_.gt.0) then
           pchar=lprefx+1
           call tbxxpstr(' ')
           pchar=lprefx+pposnam_
           pposnam_=0
         endif
         if (xmlout_) then
           if (globo_) then
             call tbxxpxot('global_',' ')
           else
             if (xmdata.eq.0) then
               call tbxxpxot(pdblok,' ')
             else
               call tbxxpxot(tbxxxsub(pdblok,xmlate(xmdata)),' ')
             endif
             if (saveo_) then
               call tbxxpxot('save_',' ')
             endif
           endif
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         pchar=lprefx
         plcat = ' '
         ploopn = 0
!
200      return
         end

!
!
!
!
!
!
! >>>>>> Process a name to extract the category and item
!
         subroutine tbxxgcat(name,type,flag,tflag,mycat,myxcat,         &
     &     item,xitem,nroot)
!
         character  name*(*),mycat*(*),item*(*),nroot*(*),type*4
         character  myxcat*(*),xitem*(*)
         include   'ciftbx.sys'
         character  xxxtemp*(NUMCHAR)
         logical    flag,tflag
         integer    lastnb,kpl,npl
         character  str1*(NUMCHAR), str2*(NUMCHAR)
         integer    kdc, lmycat


         item  = name
         xitem = ' '
         nroot = name
         mycat = ' '
         myxcat = ' '
         flag = .true.
         tflag = .true.
         if(vcheck.eq.'yes') then
           kdc = 0
           call tbxxdck(name,type,flag,tflag)
           if (xdchk.ne.0) then
             kdc = dcindex(xdchk)
             if (xmindex(xdchk).ne.0) xitem = xmlate(xmindex(xdchk))
           endif
           if (aliaso_.and.xdchk.ne.0) then
             if (aroot(xdchk).ne.0) then
               nroot = dictag(aroot(xdchk))
               kdc = dcindex(aroot(xdchk))
             endif
           endif
           if (kdc.ne.0) then
             mycat = dcname(kdc)
             myxcat = ' '
             if (xmcind(kdc).ne.0) myxcat = xmlate(xmcind(kdc))
           endif
         else
           call tbxxcat(name,mycat,lmycat)
         endif
         kpl = lastnb(mycat)
         npl = lastnb(name)
         call tbxxnlc(str1, mycat)
         call tbxxnlc(str2, name)
         if (mycat(1:1).ne.' ' .and. name(1:1).eq.'_') then
           if(str1(1:kpl).eq.str2(2:kpl+1) .and. npl .gt. kpl+2 .and.   &
     &       (name(kpl+2:kpl+2).eq.'.' .or.                             &
     &       name(kpl+2:kpl+2).eq.'_') ) then
             item = name(kpl+3:npl)
           else
             item = name(2:npl)
           endif
         else
           if (mycat(1:1).eq.' ' .and. plcat(1:1).ne.' '                &
     &       .and. name(1:1).eq.'_') then
             call tbxxnlc(str1, plcat)
             kpl = lastnb(plcat)
             if(str1(1:kpl).eq.str2(2:kpl+1) .and. npl .gt. kpl+2 .and. &
     &         (name(kpl+2:kpl+2).eq.'.' .or.                           &
     &         name(kpl+2:kpl+2).eq.'_') ) then
               mycat = plcat
               item = name(kpl+3:npl)
             else
               item = name(2:npl)
             endif
           else
             item = name
             if (item(1:1).eq.'_') item = name(2:npl)
           endif
         endif
         if (xmlong_) then
           item = name
           if (item(1:1).eq.'_') item = name(2:npl)
         endif
         call tbxxnupc(xxxtemp,mycat)
         mycat = xxxtemp
         return
         end

!
!
!
!
!
!
! >>>>>> Put a number into the CIF, perhaps with an esd appended
!
         function pnumb_(name,numb,sdev)
!
         logical    pnumb_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         real       numb,sdev
         double precision dnumb,dsdev,dprec
         integer    kmn
         integer    lastnb
!
         pnumb_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call tbxxeot
         if (pdepth_ .gt.0.and.name(1:1).ne.char(0))                    &
     &     call tbxxebkt
!
         if(name(1:1).eq.' '.or.name(1:1).eq.char(0))                   &
     &     goto 110
         call tbxxgcat(name,'numb',flag,tflag,mycat,myxcat,             &
     &     item,xitem,temp)
         pnumb_=flag
         if(ploopn.ne.0)        call tbxxelp
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)            &
     &       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,plxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         go to 120
!
110      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
!
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         dprec=decprc
         dnumb=numb
         dsdev=sdev
         call tbxxpnum(dnumb,dsdev,dprec)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not correct type')
             call tbxxpstr(char(0))
           endif
         endif
         if (xmlout_) then
           if (ploopn.gt.1 .and.ploopc.gt.0) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
         endif
!
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         pesddig_=0
         return
         end
!
!
!
!
!
!
! >>>>>> Put a double precision number into the CIF, perhaps
!        with an esd appended
!
         function pnumd_(name,numb,sdev)
!
         logical    pnumd_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         double precision numb,sdev
         integer    kmn
         integer    lastnb
!
         pnumd_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call tbxxeot
         if (pdepth_ .gt.0.and.name(1:1).ne.char(0))                    &
     &      call tbxxebkt
!
         if(name(1:1).eq.' '.or.name(1:1).eq.char(0))                   &
     &      goto 110
         call tbxxgcat(name,'numb',flag,tflag,mycat,myxcat,             &
     &     item,xitem,temp)
         pnumd_=flag
         if(ploopn.ne.0)        call tbxxelp
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)            &
     &       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,myxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         go to 120
!
110      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
!
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         call tbxxpnum(numb,sdev,dpprc)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not correct type')
             call tbxxpstr(char(0))
           endif
         endif
         if (xmlout_) then
           if (ploopn.gt.1 .and.ploopc.gt.0) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
         endif
!
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         pesddig_=0
         return
         end
!
!
!
!
!
!
! >>>>>> Put a character string into the CIF.
!
         function pchar_(name,string)
!
         logical    pchar_
         include   'ciftbx.sys'
         logical    flag,tflag
         logical    pdelim_
         character  name*(*),temp*(NUMCHAR),string*(*)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         character  line*(MAXBUF),strg*(MAXBUF)
         character*3 tsq,tdq,pqt
         integer    i, j, kfold, pql
         integer    lstring
         integer    lastnb
         integer    kmn, ic
         character*1 slash
!
         slash = rsolidus(1:1)
         pchar_=.true.
         flag  =.true.
         tflag =.true.
         tsq = ''''''''
         tdq = '"""'
         temp  =name
         lstring = lastnb(string)
         if (lstring .gt. MAXBUF) then
           call tbxxwarn(                                               &
     &       'Output CIF line longer than MAXBUF, truncated')
           lstring = MAXBUF
         endif
         pqt = pquote_
         pql = lastnb(pqt)
         if(ptextf.eq.'yes') call tbxxeot
         if(pdepth_ .gt.0.and.name(1:1).ne.char(0))                     &
     &     call tbxxebkt
!
         if(name(1:1).eq.' '.or.name(1:1).eq.char(0))                   &
     &     goto 110
         call tbxxgcat(name,'char',flag,tflag,mycat,myxcat,             &
     &     item,xitem,temp)
         pchar_=flag
         if(ploopn.ne.0)        call tbxxelp
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.gt.0) pchar=posnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)            &
     &       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,plxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         go to 120
!
110      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
!
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         i=1
         if (string(1:1).eq.char(0)) go to 210
         if (xmlout_) then
           do ic = 1,lstring
            if ( string(ic:ic).eq.'&'                                   &
     &        .or. string(ic:ic).eq.'<'                                 &
     &        .or. string(ic:ic).eq.'>' ) then
              if(i.lt.MAXBUF) then
                line(i:i) = '&'
              endif
              if (i.lt.MAXBUF) then
                if( string(ic:ic).eq.'&' ) then
                  line(i:MAXBUF)='amp;'
                  i = i+4
                endif
                if( string(ic:ic).eq.'<' ) then
                  line(i:MAXBUF)='lt;'
                  i = i+3
                endif
                if( string(ic:ic).eq.'>' ) then
                  line(i:MAXBUF)='gt;'
                  i = i+3
                endif
              endif
              if (i.gt.MAXBUF+1) then
                i = MAXBUF+1
              endif
            else
              if(i.lt.MAXBUF) then
                line(i:i) = string(ic:ic)
                i = i+1
              endif
            endif
           enddo
           if (i.gt.1) i = i-1
           if (i.lt.MAXBUF) line(i+1:MAXBUF) = ' '
         else
           line=string
           i = lstring
         endif
         if(pposval_.ne.0.and.pposend_.ge.pposval_)                     &
     &      i=max(i,pposend_-pposval_+1)
         if(pfold_ .ne. 0 .and. lstring .gt. min(pfold_,line_) )        &
     &      go to 290
         if (i .gt. MAXBUF) then
           call tbxxwarn(                                               &
     &       'Output CIF line longer than MAXBUF, truncated')
           i = MAXBUF
         endif
         if(pquote_.ne.' ')   go to 150
         do 140 j=i,1,-1
         if(line(j:j).eq.' ') go to 150
140      continue
         if((line(1:1).eq.'_'                                           &
     &     .or. line(i:i).eq.'_'                                        &
     &     .or. line(1:1).eq.''''                                       &
     &     .or. line(1:1).eq.'"'                                        &
     &     .or. line(1:1).eq.';'                                        &
     &     .or. line(1:1).eq.'('                                        &
     &     .or. line(1:1).eq.'['                                        &
     &     .or. line(1:1).eq.'(')                                       &
     &     .and.line(1:i).ne.'''.'''                                    &
     &     .and.line(1:i).ne.'''?'''                                    &
     &     .and.line(1:i).ne.'"."'                                      &
     &     .and.line(1:i).ne.'"?"') go to 150
         strg=line(1:i)
         goto 200
150      if(pqt.eq.';'                                                  &
     &     .or. pqt.eq. tsq                                             &
     &     .or. pqt.eq. tdq                                             &
     &     .or. pqt.eq. '('                                             &
     &     .or. pqt.eq. '{'                                             &
     &     .or. pqt.eq. '[')       go to 190
         if(line(1:i).eq.' '.and.nblanko_) then
           strg = '.'
           i = 1
           if(pposval_.ne.0) then
             pchar=pposval_+lprefx
           endif
           call tbxxpstr(strg(1:i))
           go to 210
         endif
         if(pqt.eq.'"' .or. pqt.eq.'"""')       go to 170
         do 160 j=1,i-1
         if(line(j:j).eq.''''.and.                                      &
     &     (line(j+1:j+1).eq.' '.or.line(j+1:j+1).eq.tab))              &
     &     goto 170
160      continue
         strg=''''//line(1:i)//''''
         i=i+2
         pqt = ''''
         pql = 1
         goto 200
170      do 180 j=1,i-1
         if(line(j:j).eq.'"'.and.                                       &
     &     (line(j+1:j+1).eq.' '.or.line(j+1:j+1).eq.tab))              &
     &     goto 190
180      continue
         strg='"'//line(1:i)//'"'
         i=i+2
         pqt = '"'
         pql = 1
         if(pfold_ .gt. 1 .and. i .gt. min(pfold_,line_) ) go to 290
         goto 200
190      pchar=-1
         if (xmlout_) then
           if (pqt.eq.';') then
           strg = '<![CDATA['//line(1:i)
           i = i+9
         else
             strg = pqt(1:pql)//'<![CDATA['//line(1:i)
             i = i+9+pql
           endif
         else
           if (pclipt_ .and. pqt.eq.';' .and. line(1:1).eq.' ') then
             strg =pqt(1:pql)//line(2:i)
             i = i-1+pql
           else
             strg =pqt(1:pql)//line(1:i)
             i = i+pql
           endif
         endif
         ptextf='yes'
         call tbxxpstr(strg(1:i))
         pchar=-1
         ptextf='no '
         if (pqt.eq.'('                                                 &
     &     .or.pqt.eq.'{'                                               &
     &     .or.pqt.eq.'[') then
           if (pqt.eq.'(') call tbxxpstr(')')
           if (pqt.eq.'{') call tbxxpstr('}')
           if (pqt.eq.'[') call tbxxpstr(']')
         else
           call tbxxpstr(pqt(1:pql))
         endif
         if (xmlout_) then
           call tbxxpstr(']]>')
         endif
         pchar=lprefx
         call tbxxpstr(' ')
         strg =                                                         &
     &    ' Converted pchar_ output to text for: '//string(1:lstring)
         call tbxxwarn(strg)
         goto 210
!
200      if(pposval_.ne.0) then
           pchar=pposval_+lprefx
           if(pqt.ne.' ') pchar=pchar-pql
         endif
         if (pdepth_.gt.0) then
           if(pstatestack(pdepth_).eq.2)                                &
     &     tbxxrslt = pdelim_(',',.false.,0)
         end if
         call tbxxpstr(strg(1:i))
         if (pdepth_.gt.0) pstatestack(pdepth_) = 2
210      if(.not.flag) then
           pchar = pcharl+4
           if (.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if((.not.tflag).and.line(1:i).ne.'.'.and.                      &
     &     line(1:i).ne.'?'.and.pqt.eq.' ') then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not correct type')
             call tbxxpstr(char(0))
           endif
         endif
         if (xmlout_) then
           if (ploopn.gt.1 .and. ploopc.gt.0) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
         endif


         pposval_=0
         pposdec_=0
         pposnam_=0
         pposend_=0
         pquote_=' '
         return

!
!        fold a string to min(pfold_,line_)
!
290      pchar=-1
         pqt = ';'
         pql = 1
         if (xmlout_) then
           call tbxxpstr('<![CDATA['//slash)
         else
           call tbxxpstr(';'//slash)
         endif
         kfold = min(pfold_,line_)
         call tbxxpfs(string(1:lstring),' ',kfold)
         pchar=-1
         ptextf='no '
         if (xmlout_) then
           call tbxxpstr(']]>')
         else
           call tbxxpstr(';')
         endif
         pchar=lprefx
         call tbxxpstr(' ')
         if (pdepth_.gt.0) pstatestack(pdepth_) = 2
         goto 210

         end
!
!
!
!
!
! >>>>>> Put a comment in the output CIF
!
         function pcmnt_(string)
!
         logical    pcmnt_
         include   'ciftbx.sys'
         character  string*(*), temp*(MAXBUF)
         character*3 pqt
         character*1 slash
         integer lstring, kfold, pql, ic, ik
         integer lastnb
!
         slash = rsolidus(1:1)
         lstring = lastnb(string)
         pqt = pquote_
         pql = lastnb(pquote_)
         kfold = min(pfold_,line_)
         if(ptextf.eq.'yes') call tbxxeot
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if(string.eq.' '.or.                                           &
     &     (string.eq.char(0)) .or.                                     &
     &     (string.eq.tab.and.(.not.ptabx_))) then
           if(string.eq.' ') pchar=-1
           if (pqt.eq.'#') then
             temp(1:1+lstring) = pqt(1:pql)//string(1:lstring)
             call tbxxpstr(temp(1:1+lstring))
           else
             call tbxxpstr(string)
           endif
           if(string.eq.' ') call tbxxpstr(char(0))
         else
           if ((kfold .ne. 0) .and.                                     &
     &        ((xmlout_ .and. (max(pchar,1)+8+lstring.gt.kfold))        &
     &        .or.((.not.xmlout_) .and.                                 &
     &        ((max(pchar,1)+lstring).gt.kfold)))) then
           if (xmlout_) then
             call tbxxpstr('<!-- '//slash)
             call tbxxpfs(string(1:lstring),' ',kfold)
           else
             call tbxxpstr('#'//slash)
             call tbxxpfs(string(1:lstring),"#",kfold)
           endif
           call tbxxpstr(char(0))
           else
           if (xmlout_) then
             temp(1:5) = '<!-- '
             ik = 6
             do ic = 1,lastnb(string)
               if (string(ic:ic).eq.'-'.and.                            &
     &           temp(ik-1:ik-1).eq.'-'.and.                            &
     &           ik.lt.MAXBUF-4) then
                 temp(ik:ik) = slash
                 ik = ik+1
               endif
               if (ik.lt.MAXBUF-4) then
                 temp(ik:ik) = string(ic:ic)
                 ik = ik+1
               endif
             enddo
             temp(ik:ik+3) = ' -->'
             ik = ik+4
             if (ik.lt.MAXBUF) temp(ik:MAXBUF) = ' '
           else
             temp='#'//string
           endif
           call tbxxpstr(temp(1:lastnb(temp)))
           call tbxxpstr(char(0))
           endif
         endif
         pcmnt_=.true.
         pposnam_=0
         if(string.ne.tab)pchar=lprefx+1
         return
         end
!
!
!
!
!
!
!
! >>>>>> Put a text sequence into the CIF.
!
         function ptext_(name,string)
!
         logical    ptext_
         integer    lastnb
         include   'ciftbx.sys'
         logical    flag,tflag
         integer    ll
         character  name*(*),temp*(NUMCHAR),string*(*),store*(NUMCHAR)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         character  temp2*(MAXBUF)
         character  slash*1
         character*3 pqt
         integer    pql
         integer    kmn
         integer    kfold
         data       store/'                                '/
!
!DBG     print *,' ptext_, pclipt_ ', pclipt_
         ptext_=.true.
         flag  =.true.
         tflag =.true.
         slash = rsolidus(1:1)
         pqt = pquote_
         if (pqt .eq. ' ') pqt = ';'
         pql = lastnb(pqt)
         ll=lastnb(string)
         temp=name
         if(ptextf.eq.'no ')    goto 100
         if(temp.eq.store)      goto 150
         call tbxxeot
         if (pdepth_ .gt.0)  call tbxxebkt
!
100      if(name(1:1).ne.' ')   goto 110
         if(ptextf.eq.'yes')    goto 150
         goto 120
!
110      if(ploopn.ne.0)        call tbxxelp
         call tbxxgcat(name,'char',flag,tflag,mycat,myxcat,             &
     &     item,xitem,temp)
         ptext_=flag
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)            &
     &       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,plxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp)
         endif
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
           else
             call tbxxpstr('#< not correct type')
           endif
         endif
         go to 130
!
120      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
!
130      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         ptextf='yes'
         store=temp
         if (pfold_.eq.0) then
           pchar=-1
           if (xmlout_) then
             if (pclipt_ .and. string(1:1).eq.' '.and.ll.gt.1) then
               if (pqt.eq.';') then
             temp2 = '<![CDATA['//string(2:ll)
               else
                 temp2 = '<![CDATA['//pqt(1:pql)//string(2:ll)
               endif
             else
               if (pqt.eq.';') then
                 temp2 = '<![CDATA['//string(1:ll)
           else
                 temp2 = '<![CDATA['//pqt(1:pql)//string(1:ll)
               endif
             endif
           else
             if (pclipt_ .and. string(1:1).eq.' '.and.ll.gt.1) then
               temp2 = pqt(1:pql)//string(2:ll)
             else
               temp2 = pqt(1:pql)//string(1:ll)
             endif
           endif
           call tbxxpstr(temp2(1:lastnb(temp2)))
           pchar=-1
           return
         endif
         pchar=-1
         if (xmlout_) then
           if (pqt.eq.';') then
             call tbxxpstr('<![CDATA['//slash)
           else
             call tbxxpstr('pqt(1:pql)//<![CDATA['//slash)
           endif
         else
           call tbxxpstr(pqt(1:pql)//slash)
         endif
150      pchar=-1
         kfold = 0
         if (pfold_ .gt. 0 ) kfold = pfold_
         if (line_ .gt. 0 .and. pfold_ .gt. line_) then
            kfold = line_
         endif
         if (kfold .gt. 0) then
           call tbxxpfs(string,' ',kfold)
         else
           call tbxxpstr(string(1:max(1,ll)))
         endif
         continue
         pchar=-1
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         return
         end
!
!
!
!
!
!
! >>>>>> Put a folded string to the output CIF
!
         subroutine tbxxpfs(string,prefix,kfold)
!
         include   'ciftbx.sys'
         character *(*) string,prefix
         character *(MAXBUF) temp
         character *1 slash
         logical stabl
         integer kfold
         integer sploopn
         integer i, klow, khi, kpref, klen
         integer lastnb

         slash = rsolidus(1:1)
         sploopn = ploopn
         ploopn = -1
         stabl = tabl_
         tabl_ = .false.
         if (kfold .lt. 4) then
           call                                                         &
     &     tbxxwarn(                                                    &
     &       'Invalid attempt to fold output line, limit reset to 4')
           pfold_ = 4
           kfold = 4
         endif
         klen = lastnb(string)
         kpref = len(prefix)
         if (prefix.eq.' ') kpref=0
         klow = 1
 100     khi = klen
         if (khi.gt.klow+kfold-1-kpref) then
           khi = klow+kfold-1-kpref-1
           do i = khi,klow+1,-1
             if(string(i:i).eq.' ') then
               khi = i
               go to 120
             endif
           enddo
 120       if (kpref.gt.0) then
             temp(1:kpref+khi-klow+2) = prefix//string(klow:khi)//slash
           else
             temp(1:kpref+khi-klow+2) = string(klow:khi)//slash
           endif
           pchar = -1
           call tbxxpstr(temp(1:kpref+khi-klow+2))
           call tbxxpstr(char(0))
           klow = khi+1
           go to 100
         else
           if (string(khi:khi).eq.slash) then
             if (khi.lt.klow+kfold-1-kpref) then
                if (kpref.gt.0) then
                  temp(1:kpref+khi-klow+2) =                            &
     &              prefix//string(klow:khi)//slash
                  pchar = -1
                  call tbxxpstr(temp(1:kpref+khi-klow+2))
                  call tbxxpstr(char(0))
                  pchar = -1
                  call tbxxpstr(prefix)
                else
                  temp(1:khi-klow+2) = string(klow:khi)//slash
                  pchar = -1
                  call tbxxpstr(temp(1:khi-klow+2))
                  pchar = -1
                  call tbxxpstr(' ')
                endif
                call tbxxpstr(char(0))
             else
                if (kpref.gt.0) then
                  temp(1:kpref+khi-klow+1) = prefix//string(klow:khi)
                  pchar = -1
                  call tbxxpstr(temp(1:kpref+khi-klow+1))
                  call tbxxpstr(char(0))
                  temp(1:kpref+2) = prefix//slash//slash
                  pchar = -1
                  call tbxxpstr(temp(1:kpref+2))
                  call tbxxpstr(char(0))
                  call tbxxpstr(prefix)
                else
                  pchar = -1
                  call tbxxpstr(string(klow:khi))
                  call tbxxpstr(char(0))
                  pchar = -1
                  call tbxxpstr(slash//slash)
                  call tbxxpstr(char(0))
                  call tbxxpstr(' ')
                endif
                call tbxxpstr(char(0))
                pchar = -1
             endif
           else
             pchar = -1
             if (kpref.gt.0) then
               temp(1:kpref+khi-klow+1)=prefix//string(klow:khi)
               call tbxxpstr(temp(1:kpref+khi-klow+1))
             else
               call tbxxpstr(string(klow:khi))
             endif
             call tbxxpstr(char(0))
           endif
         endif
         pchar = -1
         ploopn = sploopn
         tabl_ = stabl
         return
         end
!
!
!
!
!
!
! >>>>>> Put a delimiter symbol into the CIF.
!
         function pdelim_(delim,force,posdlm)
!
         logical    pdelim_
         character*(*) delim
         logical    force
         integer    posdlm
         
         include   'ciftbx.sys'
         pdelim_ = .true.
         if (ptextf.eq.'yes') call tbxxeot
         if (delim.eq.'('                                               &
     &     .or. delim.eq.'['                                            &
     &     .or. delim.eq.'{') then
           if (pdepth_.gt.0) then
             if (pstatestack(pdepth_).eq.2) then
               pposdlmstk(pdepth_+1) = pchar
               call tbxxpstr(',')
               pdelimstack(pdepth_+1) = ','
               pstatestack(pdepth_) = 1
             end if
           else
             if (ploopf.eq.'yes') ploopc=0
             ploopf='no '
             ploopc = ploopc+1
             if (ploopc.gt.ploopn) ploopc=ploopc-ploopn
           end if
           pdepth_ = pdepth_+1
           pbrackstack(pdepth_) = delim
           pdelimstack(pdepth_+1) = delim
           pstatestack(pdepth_) = 1
           pposbrkstk(pdepth_) = posdlm
           pposdlmstk(pdepth_+1) = posdlm 
           go to 100
         end if
         if (delim.eq.'}') then
           if (pdepth_.eq.0) then
             if (.not.force) pdelim_ = .false.
           else
             if (pbrackstack(pdepth_) .ne.'{'                           &
     &         .and. .not.force) pdelim_ = .false.
           end if
         end if
         if (delim.eq.')') then
           if (pdepth_.eq.0) then
             if (.not.force) pdelim_ = .false.
           else
             if (pbrackstack(pdepth_) .ne.'('                           &
     &         .and. .not.force) pdelim_ = .false.
           end if
         end if
         if (delim.eq.']') then
           if (pdepth_.eq.0) then
             if (.not.force) pdelim_ = .false.
           else
             if (pbrackstack(pdepth_) .ne.'['                           &
     &         .and. .not.force) pdelim_ = .false.
           end if
         end if
         if (delim.eq.':' .and. pdelimstack(pdepth_+1).eq.':') then
            if (.not.force) pdelim_= .false.
         end if
         if (.not.pdelim_) return
         if (delim.eq.'}' .or. delim.eq.']' .or. delim.eq.')') then
           pdepth_ = pdepth_-1
!DBG       print *,' decreasing pdepth ',pdepth_, precn_
           if (pdepth_ .gt.0) pstatestack(pdepth_) = 2
           go to 100
         end if
         if (delim.eq.',' .or. delim.eq.':') then
           pstatestack(pdepth_) = 1
           pdelimstack(pdepth_+1) = delim
           pposdlmstk(pdepth_+1) = posdlm
         end if
100      continue 
         if (posdlm.ne.0) pchar = lprefx+posdlm
         call tbxxpstr(delim)
         return
         end
         
!
!
!
!
!
!
! >>>>>> Put a loop_ data name into the CIF.
!
         function ploop_(name)
!
         logical    ploop_
         integer    kpc,kpl,npl
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*), temp*(NUMCHAR), mycat*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         character  item*(NUMCHAR), str*(NUMCHAR)
         character  shead*(NUMCHAR),xshead*(XMLCHAR)
         integer    lastnb
!
         ploop_=.true.
         flag  =.true.
         if(ptextf.eq.'yes')    call tbxxeot
         if (pdepth_ .gt.0)     call tbxxebkt
         if(ploopn.ne.0 .and. ploopf.ne.'yes' .and. name(1:1).eq.' ') then
           call tbxxelp
         endif
         temp = ' '
         mycat = ' '
         item = ' '
         shead = plhead(1)
         xshead = plxhead(1)(1:min(XMLCHAR,NUMCHAR))
         str = ' '
         if(name(1:1).eq.' ')   goto 100
!
         call tbxxgcat(name,'    ',flag,tflag,mycat,myxcat,             &
     &     item,xitem,str)
         ploop_ = flag
         if (ploopn.ne.0 .and. ploopf.ne.'yes') then
           if (plcat.eq.mycat) then
             plcat = ' '
             call tbxxelp
             plcat = mycat
             plxcat = myxcat
           else
             call tbxxelp
           endif
         endif
         if (xmlout_) then
           if (plcat(1:1).ne.' '.and.ploopn.eq.0) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             shead = ' '
             xshead = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = ' '
               plxcat = ' '
             endif
           endif
         endif
         if(tabl_.and.pposnam_.eq.0) then
           temp='    '//str(1:NUMCHAR-4)
         else
           temp=str
         endif
         plhead(max(ploopn,0)+2) = item
         plxhead(max(ploopn,0)+2) = xitem
100      if(ploopn.ne.0)        goto 120
         ploopf='yes'
         pchar=-1
         if(pposval_.ne.0) then
           pchar=lprefx+1
!           call tbxxpstr(' ')
           pchar=pposval_+lprefx
         else
           if(pposnam_.ne.0) then
             pchar=lprefx+1
             call tbxxpstr(' ')
             pchar=pposnam_+lprefx+1
           endif
         endif
         if (xmlout_) then
           if (shead(1:1).ne.' ') then
             call tbxxpxct (shead,xshead)
           endif
         else
           call tbxxpstr('loop_')
         endif
         pchar=-1
         if(name(1:1).eq.' ') then
           ploopn=-1
           plhead(1) = ' '
           plxhead(1) = ' '
           return
         endif
120      if(ploopn.le.0) then
           if (xmlout_.and.plcat.ne.mycat) then
             call tbxxpxct(plcat,plxcat)
             plcat = mycat
             plxcat = myxcat
             call tbxxpxot(mycat,myxcat)
           endif
         else
           if(xmlout_ .and. plcat.ne.mycat) then
             kpl = lastnb(plcat)
             if(mycat(1:1).eq.' ') then
                mycat = '(none)'
                myxcat = '_NONE_ '
             endif
             npl = lastnb(mycat)
             kpc = pchar
             call tbxxpstr('<!--  mixed categories '//plcat(1:kpl)      &
     &                      //' and '//mycat(1:npl)//' -->')
             pchar = kpc
           endif
         endif
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if (.not. xmlout_) then
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         if(flag)               goto 130
         if(.not.tabl_) pchar=lprefx+57
         if (xmlout_) then
           call tbxxpstr('<!-- '//temp(1:lastnb(temp))//                &
     &     ' not in dictionary -->')
           call tbxxpstr(char(0))
         else
           call tbxxpstr('#< not in dictionary')
           call tbxxpstr(char(0))
         endif
130      pchar=lprefx+1
         ploopn=max(ploopn,0)+1
         ploopc = 0
!
         return
         end
!
!
!
!
!
! >>>>>> Create or clear a prefix string
!        Any change in the length of the prefix string flushes
!        pending text, if any,  loops and partial output lines
!
         function prefx_(strg,lstrg)
!
         logical    prefx_
         include   'ciftbx.sys'
         character  strg*(*)
         integer    lstrg,mxline
!
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         if(lstrg.ne.lprefx.and.pcharl.gt.0) then
           pchar=-1
           call tbxxpstr(' ')
         endif
         if (lstrg.le.0) then
           prefx=' '
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx
           lprefx=0
         else
           if(lstrg.gt.mxline) then
             call tbxxwarn(' Prefix string truncated')
           endif
           prefx=strg
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx+lstrg
           obuf(1:min(mxline,lstrg))=prefx
           lprefx=lstrg
           if(mxline-lprefx.lt.NUMCHAR) then
             call tbxxwarn(' Output prefix may force line overflow')
           endif
         endif
         prefx_=.true.
         return
         end
!
!
!
!
!
!
! >>>>>> Close the CIF
!
         subroutine close_
!
         include   'ciftbx.sys'
         character  tbxxxsub*(MAXBUF)
!
         if(ptextf.eq.'yes') call tbxxeot
         if (pdepth_ .gt.0)  call tbxxebkt
         if(ploopn.ne.0)     call tbxxelp
         if (xmlout_) then
           if (plhead(1)(1:1).ne.' ')                                   &
     &       call tbxxpxct(plhead(1),plxhead(1))
           if (plcat(1:1).ne.' ') call tbxxpxct(plcat,plxcat)
           if (pdblok(1:1).ne.' ')  then
             if (xmdata.eq.0) then
               call tbxxpxct(pdblok,' ')
             else
               call tbxxpxct(tbxxxsub(pdblok,xmlate(xmdata)),' ')
             endif
           endif
         endif
         pdblok = ' '
         plcat = ' '
         plxcat = ' '
         plhead(1) = ' '
         plxhead(1) = ' '
         if(pcharl.ge.lprefx+1) then
           pchar=-1
           call tbxxpstr(' ')
         endif
         if (file_(1:1) .ne. ' ') then
           file_(1:longf_) = ' '
           longf_ = 1
           close(outdev)
           precn_=0
         endif
         return
         end
!
!
! >>>>>  Clean out characters not valid for an XML name
!
!        An XML name may begin with a letter, '_' or ':'
!        and may contain letters, digits, '_', ':', '.' or '-'
!
!        Note that the full Unicode character set would also permit
!        combining characters and extender characters, but these
!        have no representation in a 128 character ASCII set
!
!
         function tbxxxcln(xstring,lstr)
         logical tbxxxcln
         character*(*) xstring
         integer lstr, ii, ix
         character*10 chkstr1
         character*28 chkstr2
         character*28 chkstr3
         character*1 c
         data chkstr1/'0123456789'/
         data chkstr2/'_:ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data chkstr3/'abcdefghijklmnopqrstuvwxyz.-'/
         tbxxxcln = .true.
         do ii = 1,lstr
           c = xstring(ii:ii)
           if (c.eq.' ') return
           ix = index(chkstr2,c)
           if (ix.eq.0) then
             ix = index(chkstr3,c)
             if (ix.eq.0.and.ii.gt.1) then
               ix = index(chkstr1,c)
             endif
             if(ix.eq.0.or.(ix.gt.26.and.ii.eq.1)) then
               xstring(ii:ii) = '_'
               tbxxxcln = .false.
             endif
           endif
         enddo
         return
         end




!
!
! >>>>>> Put out the given string as an xml open tag
!
!        Note that the string may have embedded blanks and
!        parameters.  The second argument is an optional
!        translation to be used if non-blank.
!
         subroutine tbxxpxot(string,xstring)
!
         integer lastnb
         include 'ciftbx.sys'
         character sbuf*(MAXBUF)
         character*(*) string, xstring
         integer ik
         logical tbxxxcln

         if (string(1:1).eq.' ') return
         sbuf(1:1) = '<'
         if (xstring(1:1).eq.' ') then
           ik = lastnb(string)
           sbuf(2:ik+1)=string(1:ik)
         else
           ik = lastnb(xstring)
           sbuf(2:ik+1)=xstring(1:ik)
         endif
         sbuf(ik+2:ik+2) = '>'
         pchar = -1
         if (.not.tbxxxcln(sbuf(2:ik+1),ik)) then
           call tbxxwarn(' XML required remapping for '//sbuf(2:ik+1))
         endif
         call tbxxpstr(sbuf(1:ik+2))
         return
         end


!
!
! >>>>>> Put out the given string as an xml close tag
!
!        Note that the string may have embedded blanks and
!        parameters.  Only the first token will be used for close.
!        The second argument is an optional translation to be
!        used if non-blank
!
         subroutine tbxxpxct(string, xstring)
!
         include 'ciftbx.sys'
         character sbuf*(MAXBUF)
         character*(*) string, xstring
         integer ik
         logical tbxxxcln

         if (string(1:1).eq.' ') return
         sbuf(1:2) = '</'
         if (xstring(1:1).eq.' ') then
           do ik = 1,len(string)
             if (string(ik:ik).eq.' ') go to 100
           enddo
           ik = len(string)+1
 100       ik = ik-1
           sbuf(3:ik+2)=string(1:ik)
         else
           do ik = 1,len(xstring)
             if (xstring(ik:ik).eq.' ') go to 200
           enddo
           ik = len(xstring)+1
 200       ik = ik-1
           sbuf(3:ik+2)=xstring(1:ik)
         endif
         sbuf(ik+3:ik+3) = '>'
         if (.not.tbxxxcln(sbuf(3:ik+2),ik)) then
           call tbxxwarn(' XML required remapping for '//sbuf(3:ik+2))
         endif
         pchar = -1
         call tbxxpstr(sbuf(1:ik+3))
         return
         end

!
!
!
!
!
! >>>>>> Put the string into the output CIF buffer
!
         subroutine tbxxpstr(string)
!
         integer    lastnb
         include   'ciftbx.sys'
         SAVE
         character  string*(*),temp*(MAXBUF),bfill*(MAXBUF)
         character  temp2*(MAXBUF)
         integer    i,ii,mxline,ioffst,ifree,icpos,itpos
         integer    ixpos,ixtpos,it,im,kbin,kpass
         integer    lstring
         logical    pflush,waslop
         data       waslop /.false./

!
!DBG     print *,' entry to tbxxpstr, pchar, pcharl, string'
!DBG     print *, pchar, ', ',pcharl,', ', string
         bfill = ' '
         lstring = min(MAXBUF,lastnb(string))
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         if(pfold_.gt.0) then
           if (pfold_ .lt. lprefx+lstring) then
             call tbxxwarn('Invalid value of pfold_, reset')
             pfold_ = min(line_,lprefx+lstring)
           endif
           mxline=min(mxline,pfold_)
         endif
         temp(1:lstring)=string
         temp2=temp
         pflush=.false.
         if(pchar.lt.0) pflush=.true.
!
         do 100 i=lstring,1,-1
         if(temp(i:i).eq.' ')              goto 100
         if(ptabx_.and.temp(i:i).eq.tab) goto 100
         goto 110
100      continue
         i=0
         it=i
!
!....... Organise the output of loop_ items
!
110      if(i.eq.0)             goto 130
         if(i.eq.1.and.string.eq.tab) goto 130
         if(i.eq.1.and.string.eq.char(0)) then
           pcharl=MAXBUF
           goto 200
         endif
         if((.not.xmlout_).and.temp(1:1).eq.'#')   goto 130
         if(xmlout_.and.temp(1:1).eq.'<') go to 130
         if(ploopf.eq.'yes')    goto 130
         if(ptextf.eq.'yes')    goto 130
         if(ploopn.le.0)        goto 130
         if(pdepth_.gt.0)       goto 130
         ploopc=ploopc+1
         if((align_.or.tabl_).and.pdepth_.eq.0 )then
           if(ploopc.gt.ploopn) then
             if(pcharl.gt.lprefx) pflush=.true.
             ploopc=1
             if(pchar.gt.0) pchar=lprefx+1
           endif
           if(pchar.lt.0)    goto 130
           if(tabl_) then
           kbin=(mxline-lprefx)/8
           if(ploopn.lt.kbin) then
             if(kbin/(ploopn+1).gt.1) then
             pchar=9+lprefx+                                            &
     &         (ploopc-1)*8*(kbin/(ploopn+1))
             else
             pchar=1+lprefx+                                            &
     &         (ploopc-1)*8*(kbin/ploopn)
             endif
           else
             if(ploopc.le.kbin) then
               pchar=1+lprefx+(ploopc-1)*8
             else
               kpass=(ploopc-kbin-1)/(kbin-1)+1
               pchar=2*kpass+1+lprefx+                                  &
     &           mod(ploopc-kbin-1,kbin-1)*8
             endif
           endif
           else
             if(ptabx_) then
             icpos=1
             itpos=1
120          ixpos = 0
             if (icpos.le.i) ixpos=index(temp(icpos:i),tab)
             ixtpos=(pchar+itpos-1+ixpos)
             ixtpos=((ixtpos+7)/8)*8
             if(ixpos.gt.0) then
               if(ixpos.gt.1) then
                 temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
               else
                 temp2(itpos:ixtpos-pchar+1)=' '
               endif
               icpos=ixpos+1
               itpos=ixtpos+2-pchar
               if(icpos.le.i) goto 120
               it=itpos-1
             else
               if(icpos.le.i) then
                 temp2(itpos:itpos+i-icpos)=temp(icpos:i)
                 it=itpos+i-icpos
               endif
             endif
             endif
             if((pchar+i).gt.mxline+1.or.                               &
     &          (ptabx_.and.pchar+it.gt.mxline+1)) then
               if(pcharl.gt.lprefx)pflush=.true.
               pchar=lprefx+1
             endif
           endif
         else
           if(ploopc.le.ploopn)   goto 130
           ploopc=mod(ploopc-1,ploopn)+1
         endif
!
!....... Is the buffer full and needs flushing?
!
130      if(i.eq.1.and.string.eq.tab) then
           if(pcharl.gt.lprefx) then
             if(obuf(pcharl:pcharl).eq.' ') pcharl=pcharl-1
           endif
         endif
         if(pdepth_.gt.0                                                &
     &     .and.string(1:1).eq.'#'                                      &
     &     .and. pcharl.gt.lprefx) then
           if (obuf(pcharl:pcharl).ne.' ') then
             pcharl = pcharl+1
             obuf(pcharl:pcharl) = ' '
           end if
         end if
         if(pchar.le.pcharl.and.pcharl.gt.lprefx) pflush=.true.
         pchar=max(lprefx+1,pchar)
         if (string.ne.'('                                              &
     &     .and.string.ne.'['                                           &
     &     .and.string.ne.'{'                                           &
     &     .and.string.ne.')'                                           &
     &     .and.string.ne.']'                                           &
     &     .and.string.ne.'}'                                           &
     &     .and.string.ne.','                                           &
     &     .and.string.ne.':'                                           &
     &     .and.(ploopf.eq.'yes'.or.ploopn.le.0).and.tabl_)             &
     &     pchar=((pchar-lprefx+6)/8)*8+1+lprefx
         if(ptabx_) then
           icpos=1
           itpos=1
135        ixpos=0
           if(icpos.le.i) ixpos=index(temp(icpos:i),tab)
           ixtpos=(pchar+itpos-1+ixpos)
           ixtpos=((ixtpos+7)/8)*8
           if(ixpos.gt.0) then
             if(ixpos.gt.1) then
               temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
             else
               temp2(itpos:ixtpos-pchar+1)=' '
             endif
             icpos=ixpos+1
             itpos=ixtpos+2-pchar
             if(icpos.le.i) goto 135
             it=itpos-1
           else
             if(icpos.le.i) then
               temp2(itpos:itpos+i-icpos)=temp(icpos:i)
               it=itpos+i-icpos
             endif
           endif
         endif
         if((pchar+i).gt.mxline+1.or.                                   &
     &     (ptabx_.and.pchar+it.gt.mxline+1)) then
            pflush=.true.
            pchar=mxline+1-i
            if (xmlout_) pchar = 1
            if (pdepth_.gt.0) then
              pchar = 1
              if (pposbrkstk(pdepth_)+i.lt.mxline)                      &
     &          pchar = pposbrkstk(pdepth_)+1
            end if
            pchar=max(lprefx+1,pchar)
         endif
         if(.not.pflush)  goto 150
         if(pcharl.gt.lprefx) then
           if(waslop.or.(.not.tabl_)) goto 145
           ioffst=0
           pcharl=max(lastnb(obuf(1:pcharl)),lprefx+1)
           ifree=mxline-pcharl
           if(ifree.gt.0) then
           im=numtab+2
           if(numtab.gt.0.and.numtab.le.MAXTAB) then
             if(obuf(itabp(numtab):itabp(numtab)).eq.'#')               &
     &         im=im-1
           endif
           if(ifree.ge.16.and.im.lt.4.and.                              &
     &       (obuf(1+lprefx:1+lprefx).ne.'#'                            &
     &        .and.((.not.xmlout_).or.(                                 &
     &             obuf(1+lprefx:1+lprefx).ne.'<'                       &
     &        .and.obuf(1+lprefx:1+lprefx).ne.']'))                     &
     &        .and.obuf(1+lprefx:1+lprefx).ne.';'                       &
     &        .and.obuf(1+lprefx:1+lprefx).ne.'_'                       &
     &        .and.obuf(1+lprefx:1+lprefx).ne.' '                       &
     &        .and.obuf(1+lprefx:5+lprefx).ne.'data_'                   &
     &        .and.obuf(1+lprefx:5+lprefx).ne.'save_'                   &
     &        .and.obuf(1+lprefx:5).ne.'loop_')) then
             temp(1+lprefx:pcharl)=obuf(1+lprefx:pcharl)
             obuf(1+lprefx:pcharl+8)=                                   &
     &         bfill(1:8)//temp(1+lprefx:pcharl)
             ioffst = 8
             ifree=ifree-8
             pcharl=pcharl+8
           endif
           do ii=1,min(MAXTAB,numtab)
             icpos=itabp(ii)+ioffst
             if(icpos.gt.pcharl)   goto 145
             if(im.lt.4) then
             itpos=(max(icpos-lprefx,                                   &
     &         ii*(mxline-lprefx)/im)+6)/8
             itpos=itpos*8+1+lprefx
             else
             itpos=(max(icpos-lprefx,                                   &
     &         ii*(mxline-lprefx)/im)+4)/6
             itpos=itpos*6+1+lprefx
             endif
             if((obuf(icpos:icpos).eq.''''.or.                          &
     &          obuf(icpos:icpos).eq.'"').and.                          &
     &          itpos.gt.icpos) itpos=itpos-1
             if(itpos-icpos.gt.ifree) itpos=icpos+ifree
             if(itpos.gt.icpos) then
               temp(1:pcharl-icpos+1)=                                  &
     &           obuf(icpos:pcharl)
               if(i.lt.numtab) then
                 ixpos=itabp(ii+1)+ioffst
                 if(ixpos.gt.icpos+itpos-icpos+1) then
                   if(obuf(ixpos-(itpos-icpos+1):ixpos-1).eq.           &
     &               bfill(1:itpos-icpos+1)) then
                     temp(ixpos-itpos+1:pcharl-itpos+1)=                &
     &               obuf(ixpos:pcharl)
                     pcharl=pcharl-(itpos-icpos)
                   endif
                 endif
               endif
               obuf(icpos:pcharl+itpos-icpos)=                          &
     &           bfill(1:itpos-icpos)//temp(1:pcharl-icpos+1)
               ifree=ifree-(itpos-icpos)
               ioffst=ioffst+itpos-icpos
               pcharl=pcharl+itpos-icpos
             endif
             if(ifree.le.0)      goto 145
           enddo
           endif
145        pcharl=max(1,lastnb(obuf))
           write(outdev,'(a)') obuf(1:pcharl)
         else
           if(precn_.gt.0) then
           if(lprefx.gt.0) then
           write(outdev,'(a)') obuf(1:lprefx)
           else
           write(outdev,'(a)')
           endif
           else
           precn_=precn_-1
           endif
         endif
         waslop=.false.
         precn_=precn_+1
         do ii = 1,MAXTAB
           itabp(ii)=0
         enddo
         numtab=0
         if(lprefx.gt.0) then
           obuf=prefx(1:lprefx)
         else
           obuf=' '
         endif
!
!....... Load the next item into the buffer
!
150      pcharl=pchar+i
         if(ptabx_) pcharl=pchar+it
         waslop= ploopf.eq.'no '.and.ploopn.gt.0.and.align_
         if(i.eq.0) then
           if(pcharl.eq.lprefx+1.and.                                   &
     &       obuf(lprefx+1:lprefx+1).eq.' ') pcharl=pcharl-1
             pchar=pcharl+1
             if (pdepth_.gt.0) then
               if(pcharl.gt.lprefx+1                                    &
     &           .and.obuf(pcharl:pcharl).eq.' ') then
                 pchar = pcharl
                 pcharl = pcharl-1
               end if
             endif
           goto 200
         endif
         if(ptabx_) then
           obuf(pchar:pcharl)=temp2(1:it)
         else
           if(string.eq.tab) pcharl=pcharl-1
           obuf(pchar:pcharl)=string(1:i)
         endif
         if(pchar.gt.1+lprefx) then
           numtab=numtab+1
           if(numtab.le.MAXTAB) itabp(numtab)=pchar
         endif
         pchar=pcharl+1
         if (pdepth_.gt.0.and.pcharl.gt.lprefx+1                        &
     &     .and.obuf(pcharl:pcharl).eq.' ') then
           pchar = pcharl
           pcharl = pcharl-1
         endif
         if(pchar.gt.mxline+2) then
           if (pfold_.eq.0) then
             call tbxxwarn(' Output CIF line longer than line_')
           else
             call tbxxwarn(                                             &
     &         ' Output CIF line longer than line_ or pfold_')
           endif
         endif
!
200      continue
!
!DBG     print *,' exit from tbxxpstr, pchar, pcharl, obuf'
!DBG     print *, pchar, ', ',pcharl,', ', obuf(1:pcharl)

         return
         end
!
!
!
!
!
! >>>>>> Convert the number and esd to string nnnn(m), limited
!        by relative precision prec
!
         subroutine tbxxpnum(numb,sdev,prec)
!
         include   'ciftbx.sys'
         character  string*30,temp*30,c*1,sfmt*8
         double precision numb,sdev,prec,xxnumb,xsdev,slog
         integer    i,iexp,ifp,ii,jj,j,jlnz,jn,kexp,m,ixsdev,islog
         integer    kdecp,ibexp,lexp
!
         kdecp=0
         lexp = 0
         jn = 0
         if (sdev.gt.abs(numb)*prec) then
           if (iabs(esdlim_).ne.esdcac) then
!
!            determine the number of digits set by esdlim_
!
             if (iabs(esdlim_).lt.9 .or.iabs(esdlim_).gt.99999) then
               call tbxxwarn(' Invalid value of esdlim_ reset to 19')
               esdlim_ = 19
             endif
!
!            determine the number of esd digits
!
             esddigx = int(1.+alog10(float(iabs(esdlim_))))
             esdcac = iabs(esdlim_)
           endif
!
!          if esdlim_ < 0, validate pesddig_
!
           if (esdlim_.lt. 0 )then
             if (pesddig_.lt.0 .or. pesddig_.gt.5) then
               call tbxxwarn(' Invalid value of pesddig_ reset to 0')
               pesddig_ = 0
             endif
           endif
!
!          determine kexp, the power of 10 necessary
!          to present sdev as an integer in the range
!          (esdlim_/10,esdlim_] or [1,-esdlim_] if esdlim_ < 0
!
           slog = dlog10(sdev)
           islog = int(slog+1000.)
           islog = islog-1000
           kexp = -islog+esddigx
!
!          Adjust exponent kexp, so that sdev*10**kexp
!          is in the interval (esdlim_/10,esdlim_] or [1,-esdlim_]
!
 20        if (kexp.lt.minexp) then
             call tbxxwarn(' Underflow of esd')
             ixsdev = 0
             go to 30
           endif
           if (kexp.gt.-minexp) then
             call tbxxwarn(' Overflow of esd')
             ixsdev = 99999
             go to 30
           endif
           xsdev = sdev*10.D0**kexp
           ixsdev = int(xsdev+.5)
           if (ixsdev.gt.iabs(esdlim_)) then
             kexp = kexp -1
             go to 20
           endif
           if (ixsdev.lt.(iabs(esdlim_)+5)/10) then
             kexp = kexp+1
             go to 20
           endif
!
!          lexp holds the number of trailing zeros which may be
!          sacrificed in the esd if the number itself has
!          trailing zeros in the fraction which is permitted if
!          esdlim_ is negative
!
!          If esdlim_ is negative and pesddig_ is .gt.0,
!          pesddig_ will be used to force the number of digits
!          in which case lexp has the number of digits that
!          must be sacrificed (lexp > 0) or zeros to add (lexp < 0)
!
           lexp=0
           if(esdlim_.lt.0) then
             if(pesddig_.gt.0) then
25             continue
               if(ixsdev*10**(-lexp).ge.10**(pesddig_))then
                 if(lexp.gt.0)                                          &
     &             ixsdev=ixsdev-5*10**(lexp-1)
                 ixsdev=ixsdev+5*10**lexp
                 lexp=lexp+1
                 goto 25
               endif
               if(ixsdev.lt.10**(pesddig_-1+lexp)                       &
     &           .and.lexp.gt.0) then
                 if(ixsdev*10**(-lexp).le.iabs(esdlim_))then
                   lexp =lexp-1
                   if(lexp.ge.0) then
                     ixsdev=ixsdev-5*10**lexp
                   endif
                   if(lexp.gt.0) then
                     ixsdev=ixsdev+5*10**(lexp-1)
                   endif
                   goto 25
                 endif
               endif
               kexp=kexp-lexp
               ixsdev = ixsdev/(10**lexp)
               lexp=0
             else
             do ii = 1,4
               if(mod(ixsdev,10**ii).ne.0) go to 30
               lexp = ii
             enddo
             endif
           endif
!
!          We need to present the number to the same scaling
!          at first, but will adjust to avoid Ennn notation
!          if possible
!
 30        xxnumb = dabs(numb)*10.d0**kexp+.5
           if(xxnumb*prec .gt.1.D0) then
             call tbxxwarn(' ESD less than precision of machine')
             ixsdev=0
           endif
           if(numb.lt.0.d0) xxnumb = -xxnumb
           write(string,ndpfmt)xxnumb
           if(xxnumb.lt.1.d0 .and. xxnumb.ge.0.d0)                      &
     &        string='                         0.0E0'
           if(xxnumb.gt.-1.d0 .and. xxnumb.lt.0.d0)                     &
     &        string='                        -0.0E0'
!
!          Extract the power of 10
!
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = string(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 40
               endif
             endif
           enddo
           call tbxxerr(' Internal error in tbxxpnum')
!
!          Scan the rest of the string shifting the
!          decimal point to get an integer
!
40         ifp = 0
           j=1
           do ii = 1,i-1
           c = string(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               temp(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.temp(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 50
               endif
             else
               if(c.eq.'.') then
                 ifp=1
                 if(iexp.le.0) goto 50
               endif
             endif
           endif
           enddo
!
!          The string from 1 to j-1 has an integer
!          If iexp < 0, we present a 0.  If iexp > 0
!          we pad with zeros
!
50         if(j.eq.2 .and. temp(1:1).eq.'-') then
              temp(1:2)='-0'
              j=3
              iexp=0
           endif
           if(j.eq.1 .or. iexp.lt.0) then
             temp(1:1)='0'
             j=2
             iexp = 0
             if(xxnumb.lt.0.d0) then
               temp(1:2)='-0'
               j=3
             endif
           endif
           if (iexp.gt.0) then
             do ii = 1,iexp
             temp(j:j)='0'
             j=j+1
             enddo
             iexp=0
           endif
           string=temp(1:j-1)
!
!          We have the number for which the presentation
!          would be nnnnnE-kexp.  If kexp is gt 0, we can
!          decrease it and introduce a decimal point
!
           jj=0
           if(index('0123456789',temp(1:1)).eq.0) jj=1
           if(kexp.gt.0.and.kexp.lt.j-jj+8) then
             if(kexp.lt.j-1) then
               if(plzero_ .and.                                         &
     &           j-1-kexp.eq.1.and.temp(1:1).eq.'-') then
                 string=temp(1:j-1-kexp)//'0.'//                        &
     &             temp(j-kexp:j-1)
                 j=j+2
               else
                 string=temp(1:j-1-kexp)//'.'//                         &
     &           temp(j-kexp:j-1)
                 j=j+1
               endif
               kexp = 0
             else
               if(jj.ne.0)string(1:1)=temp(1:1)
               if(plzero_) then
                 string(1+jj:2+jj)='0.'
                 do ii=1,kexp-(j-1-jj)
                   string(2+jj+ii:2+jj+ii)='0'
                 enddo
                 string(3+jj+(kexp-(j-1-jj)):30)=                       &
     &             temp(1+jj:j-1)
                 j=j+2+kexp-(j-1-jj)
               else
                 string(1+jj:1+jj)='.'
                 do ii=1,kexp-(j-1-jj)
                   string(1+jj+ii:1+jj+ii)='0'
                 enddo
                 string(2+jj+(kexp-(j-1-jj)):30)=                       &
     &             temp(1+jj:j-1)
                 j=j+1+kexp-(j-1-jj)
               endif
               kexp=0
             endif
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.gt.0.and.kdecp.lt.j-1.and.lexp.gt.0) then
             jj=0
             do ii = 1,min(lexp,j-1-kdecp)
               c = string(j-ii:j-ii)
               if(c.ne.'0') goto 60
               jj=jj+1
             enddo
60           j=j-jj
             ixsdev=ixsdev/10**jj
             if(.not.pdecp_.and.string(j-1:j-1).eq.'.') then
               j=j-1
               kdecp=0
             endif
           endif
           if(kdecp.eq.0) then
             kdecp=j
             if(pdecp_) then
               if(plzero_.and.                                          &
     &           (j.eq.1 .or. (j.eq.2.and.string(1:1).eq.'-'))) then
                 string(j:j)='0'
                 j=j+1
               endif
               string(j:j)='.'
               j=j+1
             endif
           endif
           if(kexp.ne.0) then
             write(temp(1:5),'(i5)') -kexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
!
!          if there is a standard deviation
!          append it in parentheses
!
           if(ixsdev.ne.0) then
             write(temp(1:5),'(i5)') ixsdev
             string(j:j)='('
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
             string(j:j)=')'
             j=j+1
           endif
         else
!
!          There is no standard deviation, just write numb
!          But limit to the digits implied by prec
!
           slog = dlog10(min(.1D0,max(prec,dpprc)))
           islog = int(slog+1000.5)
           islog = islog-1000
           kexp = -islog
           write(sfmt,'(5h(D30.,i2,1h))') kexp
           write(temp,sfmt)numb
!
!          Now have the number in the form
!          [sign][0].nnnnnnnnDeee
!          which, while sufficient, is not neat
!          we reformat for the case 0<=eee<=kexp
!
!
!          Extract the power of 10
!
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = temp(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 140
               endif
             endif
           enddo
           call tbxxerr(' Internal error in tbxxpnum')
!
!          Scan the rest of the string shifting the
!          decimal point to get a number with exponent 0,
!          if possible
!
140        ifp = 0
           j=1
           do ii = 1,i-1
           jn=ii
           c = temp(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               string(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.string(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 150
               endif
             else
               if(c.eq.'.') then
                 ifp = -1
                 if(iexp.le.0) goto 150
               endif
             endif
           endif
           enddo
150        if(plzero_ .and.                                             &
     &       (j.eq.1 .or.(j.eq.2.and.string(1:1).eq.'-'))) then
             string(j:j)='0'
             j=j+1
           endif
           string(j:j)='.'
           ifp = j
           j = j+1
           jlnz = j-1
           do ii = jn+1,i-1
             c = temp(ii:ii)
             if (c.ne.' ')then
               m=index('0123456789',c)
               if(m.ne.0) then
                 string(j:j)=c
                 j=j+1
                 if(m.ne.1)jlnz=j
                 if(m.eq.1.and.ifp.ge.1.and.                            &
     &             pposdec_.ne.0.and.pposend_.ne.0) then
                   if(j-1-ifp-min(iexp,0).le.pposend_-pposdec_)         &
     &               jlnz=j
                 endif
               else
                 goto 160
               endif
             endif
           enddo
160        j=jlnz
           if(j.eq.1) then
            string(1:1)='0'
            j=2
           endif
           if(iexp.lt.0.and.iexp.gt.-7.and.ifp.lt.j-1.and.              &
     &       ifp.ne.0.and.j-ifp-iexp.le.kexp) then
             temp(1:ifp)=string(1:ifp)
             do ii = 1,-iexp
               temp(ifp+ii:ifp+ii) = '0'
             enddo
             temp(ifp-iexp+1:j-iexp-1) = string(ifp+1:j-1)
             j = j-iexp
             iexp=0
             string(1:j-1) = temp(1:j-1)
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.eq.0) then
             kdecp=j
             if(pdecp_) then
               string(kdecp:kdecp)='.'
               j=j+1
             endif
           endif
           if(iexp.ne.0) then
             write(temp(1:5),'(i5)')iexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
         endif
!
         if(j.lt.1) then
           string(1:1)='0'
           j=2
         endif
         if(kdecp.lt.1)kdecp=j
         if(pposdec_.ne.0) then
           pchar=lprefx+pposdec_-kdecp+1
         else
           if(pposval_.ne.0)pchar=lprefx+pposval_
         endif
         call tbxxpstr(string(1:j-1))
         return
         end
!
!
!
!
!
! >>>>>> Check dictionary for data name validation
!
         subroutine tbxxdck(name,type,flag,tflag)
!
         include   'ciftbx.sys'
         logical    flag,tflag
         integer    nln
         character  name*(*),temp*(NUMCHAR),                            &
     &              type*4
!
         flag=.true.
         tflag=.true.
         nln = min(len(name),len(temp))
         call tbxxnlc(temp(1:nln),name)
         call hash_find(temp(1:nln),                                    &
     &     dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,xdchk)
         if(xdchk.eq.0) goto 150
         if(tcheck.eq.'no ')          goto 200
         if(type.eq.dictyp(xdchk))    goto 200
         if(type.eq.'    ')           goto 200
         if(dictyp(xdchk).eq.'text' .and. type.eq.'char') goto 200
         if(dictyp(xdchk).eq.'char' .and. type.eq.'numb') goto 200
         tflag=.false.
         goto 200
150      flag=.false.
200      continue
         return
         end
!
!
!
!
!
! >>>>>> End of text string
!
         subroutine tbxxeot
!
         include   'ciftbx.sys'
!
         character*3 pqt
         integer pql
         integer lastnb

         pqt = pquote_
         pql = lastnb(pqt)
         if (pqt.eq.' ') pqt = ';'
         if (pqt.eq.'(') pqt = ')'
         if (pqt.eq.'{') pqt = '}'
         if (pqt.eq.'[') pqt = ']'

         if(ptextf.ne.'yes') then
           call tbxxwarn(' Out-of-sequence call to end text block')
           return
         endif
         ptextf='no '
         pchar=-1

         if (xmlout_) then
           if (pqt.eq.';') then
           call tbxxpstr(']]>')
           else
             call tbxxpstr(']]>'//pqt(1:pql))
           endif
           if (ploopn.gt.1) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
           if (ploopn.le.0) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
           endif
         else
           call tbxxpstr(pqt(1:pql))
         endif
         if (pqt.eq.';') then
         call tbxxpstr(char(0))
         else
           call tbxxpstr(' ')
         endif
         return
         end
!
!
!
!
!
! >>>>>> End of bracketed structure detected;
!        close all open levels
!
         subroutine tbxxebkt
!
         include   'ciftbx.sys'
         integer   i
         character*1 cd
         if (pdepth_ .eq. 0) return
         do i = 1,pdepth_
           cd = '}'
           if (pbrackstack(1+pdepth_-i).eq.'(' ) cd = ')'
           if (pbrackstack(1+pdepth_-i).eq.'[' ) cd = ']'
           pchar = max(pcharl,lprefx+pposbrkstk(1+pdepth_-i))
           call tbxxpstr(cd)
         end do
         pdepth_ = 0
         return
         end
!
!
!
!
!
! >>>>>> End of loop detected; check integrity and tidy up pointers
!
         subroutine tbxxelp
!
         include   'ciftbx.sys'
         integer   i
!
         if(ploopn.eq.0)          goto 200
         if(ploopn.eq.-1) then
           if (xmlout_) then
             plcat = ' '
             plxcat = ' '
             plhead(1) = 'DUMMY'
             plxhead(1) = ' '
           else
             call tbxxpstr('_DUMMY')
           endif
           ploopn=1
           ploopc=0
           call tbxxwarn(                                               &
     &       ' Missing: missing loop_ name set as _DUMMY')
         endif
         if (xmlout_ .and. ploopn.eq.1 .and.                            &
     &     ploopf.ne.'yes') then
           call tbxxpxct(plhead(2),plxhead(2))
         endif
         if(ploopn.eq.ploopc)     goto 200
         do i=ploopc+1,ploopn
         if (xmlout_) then
           call tbxxpxot(plhead(i+1),plxhead(1+1))
           call tbxxpstr('DUMMY')
           call tbxxpxct(plhead(i+1),plxhead(i+1))
         else
           call tbxxpstr('DUMMY')
         endif
         enddo
         call tbxxwarn(                                                 &
     &         ' Missing: missing loop_ items set as DUMMY')
         plhead(1) = ' '
         plxhead(1) = ' '
!
200      ploopc=0
         ploopn=0
         if (xmlout_) then
           call tbxxpxct(plhead(1),plxhead(1))
           plhead(1) = ' '
           call tbxxpxct(plcat,plxcat)
           plcat = ' '
         endif
         return
         end
!
!
!
!
!
!
! >>>>>> Set common default values
!
         block data
!
         include   'ciftbx.sys'
         data cifdev     /1/
         data outdev     /2/
         data dirdev     /3/
         data errdev     /6/
         data recbeg_    /1/
         data recend_    /0/
         data loopct     /0/
         data nhash      /0/
         data ndict      /0/
         data nname      /0/
         data nbloc      /0/
         data ploopn     /0/
         data ploopc     /0/
         data xmnxlat    /0/
         data xmdata     /0/
         data rsolidus   /'\\'/
         data ploopf     /'no '/
         data ptextf     /'no '/
         data pfilef     /'no '/
         data testfl     /'no '/
         data textfl     /'no '/
         data vcheck     /'no '/
         data tcheck     /'no '/
         data catchk     /'yes'/
         data parchk     /'yes'/
         data align_     /.true./
         data append_    /.false./
         data tabl_      /.true./
         data tabx_      /.true./
         data ptabx_     /.true./
         data text_      /.false./
         data loop_      /.false./
         data ndcname    /0/
         data ncname     /0/
         data rdprn_     /.false./
         data rdbrc_     /.false./
         data rdbkt_     /.false./
         data rdtq_      /.false./
         data rdrcqt_    /.false./
         data rdcolon_   /.false./
         data save_      /.false./
         data saveo_     /.false./
         data psaveo     /.false./
         data glob_      /.false./
         data globo_     /.false./
         data alias_     /.true./
         data aliaso_    /.false./
         data nblank_    /.false./
         data nblanko_   /.false./
         data decp_      /.false./
         data pdecp_     /.false./
         data lzero_     /.false./
         data plzero_    /.false./
         data xmlout_    /.false./
         data catkey     /NUMDICT*.false./
         data xmlong_    /.true./
         data dchash     /NUMHASH*0/
         data dichash    /NUMHASH*0/
         data dhash      /NUMHASH*0/
         data dcchain    /NUMDICT*0/
         data aroot      /NUMDICT*0/
         data keychain   /NUMDICT*0/
         data ccatkey    /NUMDICT*0/
         data cindex     /NUMBLOCK*0/
         data deindex    /NUMDICT*0/
         data dcindex    /NUMDICT*0/
         data line_      /80/
         data lastch     /0/
         data dictype_   /' '/
         data dicname_   /' '/
         data dicver_    /' '/
         data diccat_    /' '/
         data tagname_   /' '/
         data plcat      /' '/
         data plhead     /NUMLP1*' '/
         data prefx      /' '/
         data file_      /' '/
         data longf_     /1/
         data tbxver_    /'CIFtbx version 4.1.0 29 Nov 2009'/
         data lprefx     /0/
         data esdlim_    /19/
         data esddig_    /0/
         data pesddig_   /0/
         data esdcac     /19/
         data esddigx    /2/
         data esdfmt     /'(e12.2)'/
         data edpfmt     /'(d12.2)'/
         data ndpfmt     /'(d30.14)'/
         data decprc     /1.e-6/
         data dpprc      /1.d-14/
         data decmin     /1.e-37/
         data dpmin      /1.d-307/
         data minexp     /-307/
         data itabp      /MAXTAB*0/
         data jrect      /-1/
         data numtab     /0/
         data recn_      /0/
         data precn_     /0/
         data posnam_    /0/
         data posval_    /0/
         data posdec_    /0/
         data posend_    /0/
         data pposnam_   /0/
         data pposval_   /0/
         data pposdec_   /0/
         data pposend_   /0/
         data quote_     /' '/
         data pquote_    /' '/
         data unfold_    /.false./
         data fold_      /.false./
         data clipt_     /.true./
         data pclipt_    /.true./
         data pfold_     /0/
         data ibkmrk     /MAXBOOK*-1,MAXBOOK*-1,                        &
     &                    MAXBOOK*-1,MAXBOOK*-1,                        &
     &                    MAXBOOK*-1,MAXBOOK*-1/
         data lnametb    /1/
         data nametb     /' '/

         end
!
!
!       change the following include to include 'clearfp_sun.f'
!       for use on a SUN
!
!        include 'clearfp.f90'

