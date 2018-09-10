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
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:29:15.295/VIMOS.2003-06-05T15:29:15.295.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:07.255.NL/VIMOS.2003-06-04T06:01:07.255.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:31:00.308/VIMOS.2003-06-05T15:31:00.308.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:16:21.529/M.VIMOS.2008-06-17T18:16:21.529.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:26:36.197/VIMOS.2003-06-06T11:26:36.197.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:16:06.327/M.VIMOS.2008-06-17T18:16:06.327.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:11:52.230/M.VIMOS.2008-06-17T18:11:52.230.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T16:02:06.241/VIMOS.2003-06-05T16:02:06.241.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:26:36.196/VIMOS.2003-06-06T11:26:36.196.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:59:29.115/VIMOS.2003-06-05T15:59:29.115.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:25:54.584/VIMOS.2003-06-04T06:25:54.584.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:20:44.067/M.VIMOS.2008-06-17T18:20:44.067.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T23:28:45.412/VIMOS.2003-06-04T23:28:45.412.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:12:53.176/M.VIMOS.2008-06-17T18:12:53.176.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T16:01:13.720/VIMOS.2003-06-05T16:01:13.720.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:32:45.291/VIMOS.2003-06-05T15:32:45.291.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:30:07.745/VIMOS.2003-06-05T15:30:07.745.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T16:02:06.242/VIMOS.2003-06-05T16:02:06.242.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:59:29.116/VIMOS.2003-06-05T15:59:29.116.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:07.255/VIMOS.2003-06-04T06:01:07.255.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:16:09.125/M.VIMOS.2008-06-17T18:16:09.125.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/STAGING/VIMOS.2003-06-04T06:01:07.256.AT/VIMOS.2003-06-04T06:01:07.256.xml"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:22.789.NL/VIMOS.2003-06-04T05:39:22.789.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/STAGING/VIMOS.2003-06-04T05:39:22.788.AT/VIMOS.2003-06-04T05:39:22.788.xml"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:25:54.585/VIMOS.2003-06-04T06:25:54.585.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:31:52.806/VIMOS.2003-06-05T15:31:52.806.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:29:13.671/VIMOS.2003-06-06T11:29:13.671.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:20:47.901/M.VIMOS.2008-06-17T18:20:47.901.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:30:07.746/VIMOS.2003-06-05T15:30:07.746.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:32:45.292/VIMOS.2003-06-05T15:32:45.292.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T01:06:56.200/VIMOS.2003-06-05T01:06:56.200.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:16:18.569/M.VIMOS.2008-06-17T18:16:18.569.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:31:00.309/VIMOS.2003-06-05T15:31:00.309.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:52.286.NL/VIMOS.2003-06-04T05:39:52.286.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:07.256/VIMOS.2003-06-04T06:01:07.256.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:11:28.741/M.VIMOS.2008-06-17T18:11:28.741.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:31:52.807/VIMOS.2003-06-05T15:31:52.807.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:36.627.NL/VIMOS.2003-06-04T06:01:36.627.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:13:17.157/M.VIMOS.2008-06-17T18:13:17.157.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:22.788/VIMOS.2003-06-04T05:39:22.788.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T01:06:56.201/VIMOS.2003-06-05T01:06:56.201.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:30:07.163/VIMOS.2003-06-06T11:30:07.163.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:52.286/VIMOS.2003-06-04T05:39:52.286.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:24:37.708/VIMOS.2003-06-04T06:24:37.708.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:07.256.NL/VIMOS.2003-06-04T06:01:07.256.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:11:23.749/M.VIMOS.2008-06-17T18:11:23.749.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T16:00:21.440/VIMOS.2003-06-05T16:00:21.440.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T23:30:31.990/VIMOS.2003-06-04T23:30:31.990.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T16:01:13.719/VIMOS.2003-06-05T16:01:13.719.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:30:07.162/VIMOS.2003-06-06T11:30:07.162.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:27:28.545/VIMOS.2003-06-06T11:27:28.545.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:36.626/VIMOS.2003-06-04T06:01:36.626.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:28:21.089/VIMOS.2003-06-06T11:28:21.089.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:24:37.707/VIMOS.2003-06-04T06:24:37.707.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:29:13.672/VIMOS.2003-06-06T11:29:13.672.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:58:36.343/VIMOS.2003-06-05T15:58:36.343.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:36.627/VIMOS.2003-06-04T06:01:36.627.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:27:28.546/VIMOS.2003-06-06T11:27:28.546.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/STAGING/VIMOS.2003-06-04T05:39:22.789.AT/VIMOS.2003-06-04T05:39:22.789.xml"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:22.788.NL/VIMOS.2003-06-04T05:39:22.788.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T01:08:21.944/VIMOS.2003-06-05T01:08:21.944.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:52.285.NL/VIMOS.2003-06-04T05:39:52.285.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:29:15.296/VIMOS.2003-06-05T15:29:15.296.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:22.789/VIMOS.2003-06-04T05:39:22.789.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T15:58:36.344/VIMOS.2003-06-05T15:58:36.344.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T06:01:36.626.NL/VIMOS.2003-06-04T06:01:36.626.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/M.VIMOS.2008-06-17T18:11:46.360/M.VIMOS.2008-06-17T18:11:46.360.fits"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-04T05:39:52.285/VIMOS.2003-06-04T05:39:52.285.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T16:00:21.441/VIMOS.2003-06-05T16:00:21.441.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-06T11:28:21.090/VIMOS.2003-06-06T11:28:21.090.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/369584/SAF/VIMOS.2003-06-05T01:08:21.943/VIMOS.2003-06-05T01:08:21.943.fits.Z"

__EOF__
