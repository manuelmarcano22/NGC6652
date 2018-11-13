#!/bin/sh

usage () {
    cat <<__EOF__
usage: $(basename $0) [-hlp] [-u user] [-X args] [-d args]
  -h        print this help text
  -l        print list of files to download
  -p        prompt for password
  -u user   download as a different user
  -X args   extra arguments to pass to xargs
  -d args   extra arguments to pass to the download program

__EOF__
}

username=mmarcano22
xargsopts=
prompt=
list=
while getopts hlpu:xX:d: option
do
    case $option in
    h)  usage; exit ;;
    l)  list=yes ;;
    p)  prompt=yes ;;
    u)  prompt=yes; username="$OPTARG" ;;
    X)  xargsopts="$OPTARG" ;;
    d)  download_opts="$OPTARG";;
    ?)  usage; exit 2 ;;
    esac
done

if test -z "$xargsopts"
then
   #no xargs option speficied, we ensure that only one url
   #after the other will be used
   xargsopts='-L 1'
fi

if [ "$prompt" != "yes" ]; then
   # take password (and user) from netrc if no -p option
   if test -f "$HOME/.netrc" -a -r "$HOME/.netrc"
   then
      grep -ir "dataportal.eso.org" "$HOME/.netrc" > /dev/null
      if [ $? -ne 0 ]; then
         #no entry for dataportal.eso.org, user is prompted for password
         echo "A .netrc is available but there is no entry for dataportal.eso.org, add an entry as follows if you want to use it:"
         echo "machine dataportal.eso.org login mmarcano22 password _yourpassword_"
         prompt="yes"
      fi
   else
      prompt="yes"
   fi
fi

if test -n "$prompt" -a -z "$list"
then
    trap 'stty echo 2>/dev/null; echo "Cancelled."; exit 1' INT HUP TERM
    stty -echo 2>/dev/null
    printf 'Password: '
    read password
    echo ''
    stty echo 2>/dev/null
fi

# use a tempfile to which only user has access 
tempfile=`mktemp /tmp/dl.XXXXXXXX 2>/dev/null`
test "$tempfile" -a -f $tempfile || {
    tempfile=/tmp/dl.$$
    ( umask 077 && : >$tempfile )
}
trap 'rm -f $tempfile' EXIT INT HUP TERM

echo "auth_no_challenge=on" > $tempfile
# older OSs do not seem to include the required CA certificates for ESO
echo "check-certificate=off"  >> $tempfile
if [ -n "$prompt" ]; then
   echo "--http-user=$username" >> $tempfile
   echo "--http-password=$password" >> $tempfile

fi
WGETRC=$tempfile; export WGETRC

unset password

if test -n "$list"
then cat
else xargs $xargsopts wget $download_opts 
fi <<'__EOF__'
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:32:34.698/VIMOS.2003-06-04T04:32:34.698.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:35:15.634.AT/VIMOS.2003-06-04T04:35:15.634.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:32:34.699/VIMOS.2003-06-04T04:32:34.699.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:32:34.699.AT/VIMOS.2003-06-04T04:32:34.699.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:33:04.204.AT/VIMOS.2003-06-04T04:33:04.204.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:33:04.203/VIMOS.2003-06-04T04:33:04.203.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:35:15.635.AT/VIMOS.2003-06-04T04:35:15.635.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:33:04.204/VIMOS.2003-06-04T04:33:04.204.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:34:46.144/VIMOS.2003-06-04T04:34:46.144.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:34:46.144.AT/VIMOS.2003-06-04T04:34:46.144.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:32:34.698.AT/VIMOS.2003-06-04T04:32:34.698.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:34:46.145/VIMOS.2003-06-04T04:34:46.145.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:34:46.145.AT/VIMOS.2003-06-04T04:34:46.145.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:35:15.634/VIMOS.2003-06-04T04:35:15.634.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/STAGING/VIMOS.2003-06-04T04:33:04.203.AT/VIMOS.2003-06-04T04:33:04.203.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/375473/SAF/VIMOS.2003-06-04T04:35:15.635/VIMOS.2003-06-04T04:35:15.635.fits.Z" -P data_with_raw_calibs

__EOF__