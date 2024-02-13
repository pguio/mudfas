
ACLOCAL_AMFLAGS  = -I m4

SUBDIRS = src lib tests matlab doc share

EXTRA_DIST = ChangeLog.0 m4/ChangeLog

## Clean up amongst other template repositories
clean-local:
	-rm -rf cxx_repository ti_files ii_files
	-rm -rf so_locations

update-cl:
	cvs2cl.pl -P --no-times --no-wrap --no-common-dir -U $(HOME)/CVSROOT/users --accum \
	  -I m4 -l "-d>2006-08-28"
