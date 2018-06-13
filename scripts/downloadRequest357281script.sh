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
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:41:51.285/VIMOS.2003-06-24T05:41:51.285.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:34:05.212/VIMOS.2003-06-04T05:34:05.212.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:28:01.922/VIMOS.2003-06-22T07:28:01.922.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.847/VIMOS.2003-06-24T04:24:10.847.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:37:12.448.NL/VIMOS.2003-06-02T03:37:12.448.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:10.005/VIMOS.2003-06-02T03:08:10.005.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:32:48.145/VIMOS.2003-06-04T05:32:48.145.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:33:17.640/VIMOS.2003-06-04T05:33:17.640.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:14.569/VIMOS.2003-06-04T05:16:14.569.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.180.NL/VIMOS.2003-06-24T06:08:54.180.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:09:33.600/VIMOS.2003-06-06T07:09:33.600.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:43.862.NL/VIMOS.2003-06-04T05:16:43.862.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:54:19.488/VIMOS.2003-06-03T08:54:19.488.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:34:46.145/VIMOS.2003-06-04T04:34:46.145.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:19.816.NL/VIMOS.2003-06-06T06:52:19.816.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:27:27.728/VIMOS.2003-06-25T03:27:27.728.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:26:24.085/VIMOS.2003-06-04T06:26:24.085.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:53:50.086.NL/VIMOS.2003-06-03T08:53:50.086.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.028.NL/VIMOS.2003-06-25T02:14:35.028.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.221.NL/VIMOS.2003-06-25T02:41:33.221.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:27.082.NL/VIMOS.2003-06-04T04:00:27.082.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:24:11.398.NL/VIMOS.2003-06-06T06:24:11.398.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.242/VIMOS.2003-06-25T03:03:17.242.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:14:43.459.NL/VIMOS.2003-06-04T04:14:43.459.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.496.NL/VIMOS.2003-06-24T05:47:09.496.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:34.593/VIMOS.2003-06-03T09:08:34.593.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:44.311.NL/VIMOS.2003-06-04T03:46:44.311.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:30.780/VIMOS.2003-06-02T03:22:30.780.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:12.788.NL/VIMOS.2003-06-04T03:46:12.788.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:25:54.584/VIMOS.2003-06-04T06:25:54.584.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:54:19.487/VIMOS.2003-06-03T08:54:19.487.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:25.774.NL/VIMOS.2003-06-06T06:38:25.774.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:34:46.144/VIMOS.2003-06-04T04:34:46.144.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:41:25.313.NL/VIMOS.2003-06-22T06:41:25.313.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:27:27.727/VIMOS.2003-06-25T03:27:27.727.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:13:50.005.NL/VIMOS.2003-06-22T03:13:50.005.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:58.044.NL/VIMOS.2003-06-02T03:22:58.044.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.986.NL/VIMOS.2003-06-24T04:24:10.986.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.025.NL/VIMOS.2003-06-25T02:14:35.025.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:14:09.879/VIMOS.2003-06-06T07:14:09.879.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T04:43:29.247.NL/VIMOS.2003-06-06T04:43:29.247.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:40:29.120/VIMOS.2003-06-04T04:40:29.120.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:27:30.997/VIMOS.2003-06-22T07:27:30.997.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.241/VIMOS.2003-06-25T03:03:17.241.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:01.681.NL/VIMOS.2003-06-04T06:33:01.681.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:40.126.NL/VIMOS.2003-06-06T06:52:40.126.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:07.255/VIMOS.2003-06-04T06:01:07.255.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:22.789.NL/VIMOS.2003-06-04T05:39:22.789.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:25:54.585/VIMOS.2003-06-04T06:25:54.585.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:26.085/VIMOS.2003-06-22T02:40:26.085.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:44.311/VIMOS.2003-06-04T03:46:44.311.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:19:13.484.NL/VIMOS.2003-06-22T02:19:13.484.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:25:10.257/VIMOS.2003-06-03T09:25:10.257.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:40:04.763.NL/VIMOS.2003-06-03T08:40:04.763.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:34:05.213/VIMOS.2003-06-04T05:34:05.213.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:26:24.086/VIMOS.2003-06-04T06:26:24.086.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.103.NL/VIMOS.2003-06-25T03:03:17.103.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.848/VIMOS.2003-06-24T04:24:10.848.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:32:48.144/VIMOS.2003-06-04T05:32:48.144.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:33:17.641/VIMOS.2003-06-04T05:33:17.641.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:14.568/VIMOS.2003-06-04T05:16:14.568.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:31.174.NL/VIMOS.2003-06-04T06:33:31.174.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:07.256/VIMOS.2003-06-04T06:01:07.256.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:34.594.NL/VIMOS.2003-06-03T09:08:34.594.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:03:15.214.NL/VIMOS.2003-06-22T07:03:15.214.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:36.627.NL/VIMOS.2003-06-04T06:01:36.627.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:05.558/VIMOS.2003-06-06T06:38:05.558.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:44.310/VIMOS.2003-06-04T03:46:44.310.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:37.264.NL/VIMOS.2003-06-02T03:08:37.264.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.242.NL/VIMOS.2003-06-25T03:03:17.242.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:13:50.004.NL/VIMOS.2003-06-22T03:13:50.004.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:41:25.314.NL/VIMOS.2003-06-22T06:41:25.314.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:19.817/VIMOS.2003-06-06T06:52:19.817.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:43.863.NL/VIMOS.2003-06-04T05:16:43.863.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.102/VIMOS.2003-06-25T03:03:17.102.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:25:10.256/VIMOS.2003-06-03T09:25:10.256.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:59:47.203/VIMOS.2003-06-22T03:59:47.203.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:12.788/VIMOS.2003-06-04T03:46:12.788.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:07.256.NL/VIMOS.2003-06-04T06:01:07.256.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:24:37.708/VIMOS.2003-06-04T06:24:37.708.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:03:15.213/VIMOS.2003-06-22T07:03:15.213.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:05.559/VIMOS.2003-06-06T06:38:05.559.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:24:11.398/VIMOS.2003-06-06T06:24:11.398.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:37:12.448/VIMOS.2003-06-02T03:37:12.448.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:27:30.998/VIMOS.2003-06-22T07:27:30.998.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:58.575/VIMOS.2003-06-04T04:00:58.575.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:18:42.137.NL/VIMOS.2003-06-22T02:18:42.137.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:58.045.NL/VIMOS.2003-06-02T03:22:58.045.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:40:29.119.NL/VIMOS.2003-06-04T04:40:29.119.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:36.626/VIMOS.2003-06-04T06:01:36.626.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:53:50.085.NL/VIMOS.2003-06-03T08:53:50.085.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:24:11.397.NL/VIMOS.2003-06-06T06:24:11.397.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:59:47.204/VIMOS.2003-06-22T03:59:47.204.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:24:37.707/VIMOS.2003-06-04T06:24:37.707.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.220.NL/VIMOS.2003-06-25T02:41:33.220.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:14:43.458.NL/VIMOS.2003-06-04T04:14:43.458.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:12.787.NL/VIMOS.2003-06-04T03:46:12.787.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:37:12.449/VIMOS.2003-06-02T03:37:12.449.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:40:04.764/VIMOS.2003-06-03T08:40:04.764.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:37:12.449.NL/VIMOS.2003-06-02T03:37:12.449.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:30.779/VIMOS.2003-06-02T03:22:30.779.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.707/VIMOS.2003-06-24T04:55:50.707.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.055/VIMOS.2003-06-24T06:08:54.055.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:19:13.485/VIMOS.2003-06-22T02:19:13.485.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:35:34.036.NL/VIMOS.2003-06-22T03:35:34.036.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:39:35.356/VIMOS.2003-06-03T08:39:35.356.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:36.627/VIMOS.2003-06-04T06:01:36.627.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:53:50.085/VIMOS.2003-06-03T08:53:50.085.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:02:44.256/VIMOS.2003-06-22T07:02:44.256.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:33:04.204/VIMOS.2003-06-04T04:33:04.204.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:27.083.NL/VIMOS.2003-06-04T04:00:27.083.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:05.185.NL/VIMOS.2003-06-03T09:08:05.185.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:35:15.634/VIMOS.2003-06-04T04:35:15.634.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:40:54.318/VIMOS.2003-06-22T06:40:54.318.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.495.NL/VIMOS.2003-06-24T05:47:09.495.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:52.285.NL/VIMOS.2003-06-04T05:39:52.285.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:40:04.763/VIMOS.2003-06-03T08:40:04.763.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:19.816/VIMOS.2003-06-06T06:52:19.816.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.054/VIMOS.2003-06-24T06:08:54.054.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:40:54.318.NL/VIMOS.2003-06-22T06:40:54.318.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.706/VIMOS.2003-06-24T04:55:50.706.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:28:01.923/VIMOS.2003-06-22T07:28:01.923.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:19:13.484/VIMOS.2003-06-22T02:19:13.484.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:39:35.355/VIMOS.2003-06-03T08:39:35.355.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:53:50.086/VIMOS.2003-06-03T08:53:50.086.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:33:04.203/VIMOS.2003-06-04T04:33:04.203.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:14.568.NL/VIMOS.2003-06-04T05:16:14.568.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:12.787/VIMOS.2003-06-04T03:46:12.787.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.055.NL/VIMOS.2003-06-24T06:08:54.055.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:02:44.255/VIMOS.2003-06-22T07:02:44.255.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:03:15.214/VIMOS.2003-06-22T07:03:15.214.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:40:29.120.NL/VIMOS.2003-06-04T04:40:29.120.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:24:11.397/VIMOS.2003-06-06T06:24:11.397.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:40:54.319/VIMOS.2003-06-22T06:40:54.319.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:35:15.635/VIMOS.2003-06-04T04:35:15.635.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.027/VIMOS.2003-06-25T02:14:35.027.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:07.255.NL/VIMOS.2003-06-04T06:01:07.255.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:11:40.954/VIMOS.2003-06-06T07:11:40.954.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:13:50.004/VIMOS.2003-06-22T03:13:50.004.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:36:45.182.NL/VIMOS.2003-06-02T03:36:45.182.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.567/VIMOS.2003-06-24T04:55:50.567.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.633.NL/VIMOS.2003-06-24T05:47:09.633.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:30.780.NL/VIMOS.2003-06-02T03:22:30.780.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.241.NL/VIMOS.2003-06-25T03:03:17.241.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:39:35.355.NL/VIMOS.2003-06-03T08:39:35.355.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:05.185/VIMOS.2003-06-03T09:08:05.185.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:01.681/VIMOS.2003-06-04T06:33:01.681.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:27:27.582/VIMOS.2003-06-25T03:27:27.582.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:27.082/VIMOS.2003-06-04T04:00:27.082.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:15:12.954.NL/VIMOS.2003-06-04T04:15:12.954.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:58.045/VIMOS.2003-06-02T03:22:58.045.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:30.779.NL/VIMOS.2003-06-02T03:22:30.779.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:14:21.277/VIMOS.2003-06-22T03:14:21.277.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:40:29.119/VIMOS.2003-06-04T04:40:29.119.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:31.175/VIMOS.2003-06-04T06:33:31.175.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:15:12.954/VIMOS.2003-06-04T04:15:12.954.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:18:42.136.NL/VIMOS.2003-06-22T02:18:42.136.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:23:51.163.NL/VIMOS.2003-06-06T06:23:51.163.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:11:40.953/VIMOS.2003-06-06T07:11:40.953.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:13:50.005/VIMOS.2003-06-22T03:13:50.005.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:58.574.NL/VIMOS.2003-06-04T04:00:58.574.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.566/VIMOS.2003-06-24T04:55:50.566.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:05.186/VIMOS.2003-06-03T09:08:05.186.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:56:02.396/VIMOS.2003-06-02T03:56:02.396.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.567.NL/VIMOS.2003-06-24T04:55:50.567.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.987.NL/VIMOS.2003-06-24T04:24:10.987.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:27:27.581/VIMOS.2003-06-25T03:27:27.581.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:36:05.268.NL/VIMOS.2003-06-22T03:36:05.268.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:14:21.278/VIMOS.2003-06-22T03:14:21.278.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.351/VIMOS.2003-06-25T02:41:33.351.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:25.774/VIMOS.2003-06-06T06:38:25.774.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:31.174/VIMOS.2003-06-04T06:33:31.174.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:05.186.NL/VIMOS.2003-06-03T09:08:05.186.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:15:12.953/VIMOS.2003-06-04T04:15:12.953.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.054.NL/VIMOS.2003-06-24T06:08:54.054.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:58.574/VIMOS.2003-06-04T04:00:58.574.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.706.NL/VIMOS.2003-06-24T04:55:50.706.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.986/VIMOS.2003-06-24T04:24:10.986.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T03:46:44.310.NL/VIMOS.2003-06-04T03:46:44.310.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:56:02.397/VIMOS.2003-06-02T03:56:02.397.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T04:00:18.414/VIMOS.2003-06-22T04:00:18.414.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.388/VIMOS.2003-06-24T05:17:34.388.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:35:34.035.NL/VIMOS.2003-06-22T03:35:34.035.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:05:51.793/VIMOS.2003-06-22T03:05:51.793.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:41:25.313/VIMOS.2003-06-22T06:41:25.313.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.252.NL/VIMOS.2003-06-24T05:17:34.252.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:14:21.278.NL/VIMOS.2003-06-22T03:14:21.278.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.848.NL/VIMOS.2003-06-24T04:24:10.848.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.028/VIMOS.2003-06-25T02:14:35.028.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:14.569.NL/VIMOS.2003-06-04T05:16:14.569.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:35:34.036/VIMOS.2003-06-22T03:35:34.036.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:25.773/VIMOS.2003-06-06T06:38:25.773.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T04:00:18.413/VIMOS.2003-06-22T04:00:18.413.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.103/VIMOS.2003-06-25T03:03:17.103.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.027.NL/VIMOS.2003-06-25T02:14:35.027.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:52.286.NL/VIMOS.2003-06-04T05:39:52.286.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.387.NL/VIMOS.2003-06-24T05:17:34.387.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:26.086.NL/VIMOS.2003-06-22T02:40:26.086.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:40:54.319.NL/VIMOS.2003-06-22T06:40:54.319.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T06:41:25.314/VIMOS.2003-06-22T06:41:25.314.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.387/VIMOS.2003-06-24T05:17:34.387.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.352.NL/VIMOS.2003-06-25T02:41:33.352.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:39:59.620.NL/VIMOS.2003-06-04T04:39:59.620.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:55:35.127/VIMOS.2003-06-02T03:55:35.127.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:01.682/VIMOS.2003-06-04T06:33:01.682.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:27.083/VIMOS.2003-06-04T04:00:27.083.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:57.405.NL/VIMOS.2003-06-22T02:40:57.405.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:22:58.044/VIMOS.2003-06-02T03:22:58.044.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:54:19.488.NL/VIMOS.2003-06-03T08:54:19.488.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:10.004.NL/VIMOS.2003-06-02T03:08:10.004.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:49:18.387/VIMOS.2003-06-04T06:49:18.387.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:12:01.120/VIMOS.2003-06-06T07:12:01.120.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.566.NL/VIMOS.2003-06-24T04:55:50.566.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:34.593.NL/VIMOS.2003-06-03T09:08:34.593.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:22.788/VIMOS.2003-06-04T05:39:22.788.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:25:07.209/VIMOS.2003-06-04T06:25:07.209.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:43.863/VIMOS.2003-06-04T05:16:43.863.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:35:34.035/VIMOS.2003-06-22T03:35:34.035.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:57.404/VIMOS.2003-06-22T02:40:57.404.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:06:23.075/VIMOS.2003-06-22T03:06:23.075.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T04:43:29.247/VIMOS.2003-06-06T04:43:29.247.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:08:16.803/VIMOS.2003-06-06T07:08:16.803.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:40.127/VIMOS.2003-06-06T06:52:40.127.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.634/VIMOS.2003-06-24T05:47:09.634.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:52.286/VIMOS.2003-06-04T05:39:52.286.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:32:58.132/VIMOS.2003-06-24T06:32:58.132.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:49:18.386/VIMOS.2003-06-04T06:49:18.386.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:36:05.268/VIMOS.2003-06-22T03:36:05.268.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:25:39.665/VIMOS.2003-06-03T09:25:39.665.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:09:33.599/VIMOS.2003-06-06T07:09:33.599.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:39:59.620/VIMOS.2003-06-04T04:39:59.620.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:55:35.128/VIMOS.2003-06-02T03:55:35.128.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:23:51.163/VIMOS.2003-06-06T06:23:51.163.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.495/VIMOS.2003-06-24T05:47:09.495.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:36:05.267.NL/VIMOS.2003-06-22T03:36:05.267.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.181.NL/VIMOS.2003-06-24T06:08:54.181.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:14:09.880/VIMOS.2003-06-06T07:14:09.880.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:36:45.181.NL/VIMOS.2003-06-02T03:36:45.181.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.221/VIMOS.2003-06-25T02:41:33.221.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.634.NL/VIMOS.2003-06-24T05:47:09.634.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:36:45.181/VIMOS.2003-06-02T03:36:45.181.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:18:42.137/VIMOS.2003-06-22T02:18:42.137.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.388.NL/VIMOS.2003-06-24T05:17:34.388.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:12:01.121/VIMOS.2003-06-06T07:12:01.121.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:18:42.136/VIMOS.2003-06-22T02:18:42.136.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:32:34.699/VIMOS.2003-06-04T04:32:34.699.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:25:07.208/VIMOS.2003-06-04T06:25:07.208.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:57.405/VIMOS.2003-06-22T02:40:57.405.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:23:51.162.NL/VIMOS.2003-06-06T06:23:51.162.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:08:16.804/VIMOS.2003-06-06T07:08:16.804.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:15:12.953.NL/VIMOS.2003-06-04T04:15:12.953.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:40.126/VIMOS.2003-06-06T06:52:40.126.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.252/VIMOS.2003-06-24T05:17:34.252.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.987/VIMOS.2003-06-24T04:24:10.987.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.181/VIMOS.2003-06-24T06:08:54.181.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:32:58.278/VIMOS.2003-06-24T06:32:58.278.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:32:58.133/VIMOS.2003-06-24T06:32:58.133.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.633/VIMOS.2003-06-24T05:47:09.633.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:49:47.885/VIMOS.2003-06-04T06:49:47.885.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:39:59.621/VIMOS.2003-06-04T04:39:59.621.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:36:05.267/VIMOS.2003-06-22T03:36:05.267.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:25:39.666/VIMOS.2003-06-03T09:25:39.666.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:23:51.162/VIMOS.2003-06-06T06:23:51.162.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:05:51.792/VIMOS.2003-06-22T03:05:51.792.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T04:43:29.246.NL/VIMOS.2003-06-06T04:43:29.246.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:26.086/VIMOS.2003-06-22T02:40:26.086.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:02:44.256.NL/VIMOS.2003-06-22T07:02:44.256.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:25.773.NL/VIMOS.2003-06-06T06:38:25.773.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:39:35.356.NL/VIMOS.2003-06-03T08:39:35.356.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:01.682.NL/VIMOS.2003-06-04T06:33:01.682.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:32:34.698/VIMOS.2003-06-04T04:32:34.698.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:05.558.NL/VIMOS.2003-06-06T06:38:05.558.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:02:44.255.NL/VIMOS.2003-06-22T07:02:44.255.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.025/VIMOS.2003-06-25T02:14:35.025.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.351.NL/VIMOS.2003-06-25T02:41:33.351.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:14:43.458/VIMOS.2003-06-04T04:14:43.458.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:37.264/VIMOS.2003-06-02T03:08:37.264.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:41:51.425/VIMOS.2003-06-24T05:41:51.425.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:38:05.559.NL/VIMOS.2003-06-06T06:38:05.559.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.251/VIMOS.2003-06-24T05:17:34.251.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:19:13.485.NL/VIMOS.2003-06-22T02:19:13.485.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:08:36.980/VIMOS.2003-06-06T07:08:36.980.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:32:58.279/VIMOS.2003-06-24T06:32:58.279.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:13:49.717/VIMOS.2003-06-06T07:13:49.717.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T06:08:54.180/VIMOS.2003-06-24T06:08:54.180.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T09:08:34.594/VIMOS.2003-06-03T09:08:34.594.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:49:47.884/VIMOS.2003-06-04T06:49:47.884.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:24:10.847.NL/VIMOS.2003-06-24T04:24:10.847.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:00:58.575.NL/VIMOS.2003-06-04T04:00:58.575.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:22.788.NL/VIMOS.2003-06-04T05:39:22.788.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:09:53.781/VIMOS.2003-06-06T07:09:53.781.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T04:55:50.707.NL/VIMOS.2003-06-24T04:55:50.707.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.352/VIMOS.2003-06-25T02:41:33.352.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.024.NL/VIMOS.2003-06-25T02:14:35.024.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T07:03:15.213.NL/VIMOS.2003-06-22T07:03:15.213.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:34:34.705/VIMOS.2003-06-04T05:34:34.705.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:39:59.621.NL/VIMOS.2003-06-04T04:39:59.621.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:57.404.NL/VIMOS.2003-06-22T02:40:57.404.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:22.789/VIMOS.2003-06-04T05:39:22.789.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T04:14:43.459/VIMOS.2003-06-04T04:14:43.459.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:33:31.175.NL/VIMOS.2003-06-04T06:33:31.175.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:41:51.286/VIMOS.2003-06-24T05:41:51.286.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:14:35.024/VIMOS.2003-06-25T02:14:35.024.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:16:43.862/VIMOS.2003-06-04T05:16:43.862.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T06:01:36.626.NL/VIMOS.2003-06-04T06:01:36.626.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:13:49.718/VIMOS.2003-06-06T07:13:49.718.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T04:43:29.246/VIMOS.2003-06-06T04:43:29.246.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:06:23.076/VIMOS.2003-06-22T03:06:23.076.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:10.005.NL/VIMOS.2003-06-02T03:08:10.005.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:41:51.424/VIMOS.2003-06-24T05:41:51.424.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:10.004/VIMOS.2003-06-02T03:08:10.004.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:54:19.487.NL/VIMOS.2003-06-03T08:54:19.487.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:37.265.NL/VIMOS.2003-06-02T03:08:37.265.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:17:34.251.NL/VIMOS.2003-06-24T05:17:34.251.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T03:14:21.277.NL/VIMOS.2003-06-22T03:14:21.277.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:39:52.285/VIMOS.2003-06-04T05:39:52.285.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-24T05:47:09.496/VIMOS.2003-06-24T05:47:09.496.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-03T08:40:04.764.NL/VIMOS.2003-06-03T08:40:04.764.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:08:36.981/VIMOS.2003-06-06T07:08:36.981.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:08:37.265/VIMOS.2003-06-02T03:08:37.265.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:19.817.NL/VIMOS.2003-06-06T06:52:19.817.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T03:03:17.102.NL/VIMOS.2003-06-25T03:03:17.102.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T06:52:40.127.NL/VIMOS.2003-06-06T06:52:40.127.NL.txt"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-02T03:36:45.182/VIMOS.2003-06-02T03:36:45.182.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-25T02:41:33.220/VIMOS.2003-06-25T02:41:33.220.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-06T07:09:53.782/VIMOS.2003-06-06T07:09:53.782.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-04T05:34:34.706/VIMOS.2003-06-04T05:34:34.706.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/357281/SAF/VIMOS.2003-06-22T02:40:26.085.NL/VIMOS.2003-06-22T02:40:26.085.NL.txt"

__EOF__
