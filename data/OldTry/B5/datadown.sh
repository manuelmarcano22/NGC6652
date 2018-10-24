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
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:03.627/VIMOS.2003-06-24T11:28:03.627.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:45:54.300/VIMOS.2003-06-24T10:45:54.300.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.633.NL/VIMOS.2003-06-24T05:47:09.633.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:50:19.638/VIMOS.2003-06-23T14:50:19.638.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:52:04.864/VIMOS.2003-06-23T14:52:04.864.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.180.NL/VIMOS.2003-06-24T06:08:54.180.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:47:11.252/VIMOS.2003-06-24T10:47:11.252.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:30:43.186/VIMOS.2003-06-24T11:30:43.186.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:16:06.327/M.VIMOS.2008-06-17T18:16:06.327.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T05:47:09.496.AT/VIMOS.2003-06-24T05:47:09.496.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:45:54.169/VIMOS.2003-06-24T10:45:54.169.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:16:03.380/M.VIMOS.2008-06-17T18:16:03.380.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.496.NL/VIMOS.2003-06-24T05:47:09.496.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:50:19.637/VIMOS.2003-06-23T14:50:19.637.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:03.763/VIMOS.2003-06-24T11:28:03.763.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:47:11.251/VIMOS.2003-06-24T10:47:11.251.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:52:04.865/VIMOS.2003-06-23T14:52:04.865.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:16:12.090/M.VIMOS.2008-06-17T18:16:12.090.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:30:43.187/VIMOS.2003-06-24T11:30:43.187.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:48:34.459/VIMOS.2003-06-23T14:48:34.459.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:20:44.067/M.VIMOS.2008-06-17T18:20:44.067.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T06:08:54.055.AT/VIMOS.2003-06-24T06:08:54.055.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:12:53.176/M.VIMOS.2008-06-17T18:12:53.176.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:49:27.030/VIMOS.2003-06-23T14:49:27.030.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:43:31.778/VIMOS.2003-06-24T10:43:31.778.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:34:16.581/VIMOS.2003-06-24T10:34:16.581.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.054.NL/VIMOS.2003-06-24T06:08:54.054.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:56.625/VIMOS.2003-06-24T11:28:56.625.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:03.762/VIMOS.2003-06-24T11:28:03.762.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:16:09.125/M.VIMOS.2008-06-17T18:16:09.125.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:13:13.956/M.VIMOS.2008-06-17T18:13:13.956.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:11:18.672/M.VIMOS.2008-06-17T18:11:18.672.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:49:27.031/VIMOS.2003-06-23T14:49:27.031.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:43:31.636/VIMOS.2003-06-24T10:43:31.636.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:44:42.923/VIMOS.2003-06-24T10:44:42.923.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:20:47.901/M.VIMOS.2008-06-17T18:20:47.901.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:11:31.782/M.VIMOS.2008-06-17T18:11:31.782.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:38:01.522/VIMOS.2003-06-24T10:38:01.522.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:45:54.301/VIMOS.2003-06-24T10:45:54.301.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:27:11.267/VIMOS.2003-06-24T11:27:11.267.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T05:47:09.633.AT/VIMOS.2003-06-24T05:47:09.633.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:50:19.473/VIMOS.2003-06-23T14:50:19.473.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:56.626/VIMOS.2003-06-24T11:28:56.626.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T06:08:54.180.AT/VIMOS.2003-06-24T06:08:54.180.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:44:42.924/VIMOS.2003-06-24T10:44:42.924.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:48:34.615/VIMOS.2003-06-23T14:48:34.615.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:11:28.741/M.VIMOS.2008-06-17T18:11:28.741.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:51:12.021/VIMOS.2003-06-23T14:51:12.021.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T05:47:09.495.AT/VIMOS.2003-06-24T05:47:09.495.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:27:11.268/VIMOS.2003-06-24T11:27:11.268.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:50:19.472/VIMOS.2003-06-23T14:50:19.472.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.634/VIMOS.2003-06-24T05:47:09.634.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:51:12.166/VIMOS.2003-06-23T14:51:12.166.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:32:58.132/VIMOS.2003-06-24T06:32:58.132.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:11:23.749/M.VIMOS.2008-06-17T18:11:23.749.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.495/VIMOS.2003-06-24T05:47:09.495.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:47:11.389/VIMOS.2003-06-24T10:47:11.389.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:48:34.616/VIMOS.2003-06-23T14:48:34.616.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:49:26.872/VIMOS.2003-06-23T14:49:26.872.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:52:04.693/VIMOS.2003-06-23T14:52:04.693.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T06:08:54.054.AT/VIMOS.2003-06-24T06:08:54.054.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.181.NL/VIMOS.2003-06-24T06:08:54.181.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:20:50.827/M.VIMOS.2008-06-17T18:20:50.827.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.634.NL/VIMOS.2003-06-24T05:47:09.634.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:51:12.022/VIMOS.2003-06-23T14:51:12.022.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:27:11.124/VIMOS.2003-06-24T11:27:11.124.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:32:58.278/VIMOS.2003-06-24T06:32:58.278.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.181/VIMOS.2003-06-24T06:08:54.181.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:56.487/VIMOS.2003-06-24T11:28:56.487.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.633/VIMOS.2003-06-24T05:47:09.633.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:32:58.133/VIMOS.2003-06-24T06:32:58.133.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:52:04.694/VIMOS.2003-06-23T14:52:04.694.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:43:31.637/VIMOS.2003-06-24T10:43:31.637.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:49:26.873/VIMOS.2003-06-23T14:49:26.873.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:44:42.776/VIMOS.2003-06-24T10:44:42.776.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:29:50.113/VIMOS.2003-06-24T11:29:50.113.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:48:34.460/VIMOS.2003-06-23T14:48:34.460.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:45:54.168/VIMOS.2003-06-24T10:45:54.168.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.055/VIMOS.2003-06-24T06:08:54.055.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:27:11.125/VIMOS.2003-06-24T11:27:11.125.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:32:58.279/VIMOS.2003-06-24T06:32:58.279.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:56.486/VIMOS.2003-06-24T11:28:56.486.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.180/VIMOS.2003-06-24T06:08:54.180.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:29:50.114/VIMOS.2003-06-24T11:29:50.114.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:26:46.402/VIMOS.2003-06-24T10:26:46.402.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:29:49.974/VIMOS.2003-06-24T11:29:49.974.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:30:43.334/VIMOS.2003-06-24T11:30:43.334.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/M.VIMOS.2008-06-17T18:20:40.981/M.VIMOS.2008-06-17T18:20:40.981.fits" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.495.NL/VIMOS.2003-06-24T05:47:09.495.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:30:31.616/VIMOS.2003-06-24T10:30:31.616.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:43:31.777/VIMOS.2003-06-24T10:43:31.777.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:44:42.777/VIMOS.2003-06-24T10:44:42.777.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T06:08:54.181.AT/VIMOS.2003-06-24T06:08:54.181.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/STAGING/VIMOS.2003-06-24T05:47:09.634.AT/VIMOS.2003-06-24T05:47:09.634.xml" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T10:47:11.390/VIMOS.2003-06-24T10:47:11.390.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:28:03.628/VIMOS.2003-06-24T11:28:03.628.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.054/VIMOS.2003-06-24T06:08:54.054.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-23T14:51:12.165/VIMOS.2003-06-23T14:51:12.165.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:29:49.973/VIMOS.2003-06-24T11:29:49.973.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T05:47:09.496/VIMOS.2003-06-24T05:47:09.496.fits.Z" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T06:08:54.055.NL/VIMOS.2003-06-24T06:08:54.055.NL.txt" -P data_with_raw_calibs
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/370380/SAF/VIMOS.2003-06-24T11:30:43.335/VIMOS.2003-06-24T11:30:43.335.fits.Z" -P data_with_raw_calibs

__EOF__
