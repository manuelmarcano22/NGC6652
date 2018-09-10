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
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:13:50.004.NL/VIMOS.2003-06-22T03:13:50.004.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:35:34.035/VIMOS.2003-06-22T03:35:34.035.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:57.404/VIMOS.2003-06-22T02:40:57.404.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:13:50.004/VIMOS.2003-06-22T03:13:50.004.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:19:13.485/VIMOS.2003-06-22T02:19:13.485.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:35:34.036.NL/VIMOS.2003-06-22T03:35:34.036.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:19:13.485.NL/VIMOS.2003-06-22T02:19:13.485.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:36:05.268/VIMOS.2003-06-22T03:36:05.268.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:35:34.035.NL/VIMOS.2003-06-22T03:35:34.035.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:36:05.267.NL/VIMOS.2003-06-22T03:36:05.267.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:26.085/VIMOS.2003-06-22T02:40:26.085.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:19:13.484.NL/VIMOS.2003-06-22T02:19:13.484.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:14:21.278.NL/VIMOS.2003-06-22T03:14:21.278.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:14:21.277/VIMOS.2003-06-22T03:14:21.277.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:18:42.137/VIMOS.2003-06-22T02:18:42.137.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:57.404.NL/VIMOS.2003-06-22T02:40:57.404.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:18:42.136/VIMOS.2003-06-22T02:18:42.136.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:35:34.036/VIMOS.2003-06-22T03:35:34.036.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:18:42.136.NL/VIMOS.2003-06-22T02:18:42.136.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:13:50.005/VIMOS.2003-06-22T03:13:50.005.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:57.405/VIMOS.2003-06-22T02:40:57.405.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:18:42.137.NL/VIMOS.2003-06-22T02:18:42.137.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:19:13.484/VIMOS.2003-06-22T02:19:13.484.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:14:21.277.NL/VIMOS.2003-06-22T03:14:21.277.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:36:05.267/VIMOS.2003-06-22T03:36:05.267.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:26.086.NL/VIMOS.2003-06-22T02:40:26.086.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:26.086/VIMOS.2003-06-22T02:40:26.086.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:36:05.268.NL/VIMOS.2003-06-22T03:36:05.268.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:13:50.005.NL/VIMOS.2003-06-22T03:13:50.005.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T03:14:21.278/VIMOS.2003-06-22T03:14:21.278.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:57.405.NL/VIMOS.2003-06-22T02:40:57.405.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369791/SAF/VIMOS.2003-06-22T02:40:26.085.NL/VIMOS.2003-06-22T02:40:26.085.NL.txt"

__EOF__
