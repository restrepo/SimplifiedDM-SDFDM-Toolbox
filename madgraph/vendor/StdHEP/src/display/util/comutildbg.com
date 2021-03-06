$ !
$ ! VMS procedure to compile and link modules for Histoscope tool.
$ !
$ ! SCCS ID: @(#)comutildbg.com	1.3     2/24/94
$ ON ERROR THEN EXIT
$ COMPILE := CC/DEBUG/NOOPT/OBJ=[.DBGOBJ]
$ @[-.HISTO]COMSYM
$ SET VERIFY
$ COMPILE DIALOGF.C
$ COMPILE FILEUTILS.C
$ COMPILE GETFILES.c
$ COMPILE MISC.C  
$ COMPILE STRINGUTILS.C
$ COMPILE HELP.C
$ COMPILE PSUTILS.C
$ COMPILE PREFFILE.C
$ COMPILE FONTSEL.C
$ COMPILE PRINTUTILS.C
$ !
$ LIBR/CREATE/OBJ [.DBGOBJ]LIBUTIL [.DBGOBJ]DIALOGF, FILEUTILS, GETFILES, -
	MISC, STRINGUTILS, HELP, PSUTILS, PRINTUTILS, PREFFILE, FONTSEL
$ !
$ COMPILE VMSUTILS.C
$ !
$ LIBR/CREATE/OBJ [.DBGOBJ]VMSUTILS [.DBGOBJ]VMSUTILS
$ !
$ COMPILE TESTDIALOGF.C
$ LINK/DEBUG/EXE=[.DBGOBJ] [.DBGOBJ]TESTDIALOGF, LIBUTIL.OLB/LIB, -
	[-.HISTO]HISTO_OPTIONS_FILE/OPT
$ COMPILE TESTGETFILES.C
$ LINK/DEBUG/EXE=[.DBGOBJ] [.DBGOBJ]TESTGETFILES, LIBUTIL.OLB/LIB, -
	[-.HISTO]HISTO_OPTIONS_FILE/OPT
$ COMPILE TESTGREEK.C
$ LINK/DEBUG/EXE=[.DBGOBJ] [.DBGOBJ]TESTGREEK, LIBUTIL.OLB/LIB, -
	[-.HISTO]HISTO_OPTIONS_FILE/OPT
$ COMPILE TESTPRINT.C
$ LINK/DEBUG/EXE=[.DBGOBJ] [.DBGOBJ]TESTPRINT, LIBUTIL.OLB/LIB, -
	[-.HISTO]HISTO_OPTIONS_FILE/OPT
