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
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:17.818/VIMOS.2003-06-03T11:21:17.818.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:31:06.794/VIMOS.2003-06-02T12:31:06.794.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:47.235/VIMOS.2003-06-03T11:21:47.235.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:37:13.007/VIMOS.2003-06-02T12:37:13.007.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:27:06.395/VIMOS.2003-06-02T19:27:06.395.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:36:21.589/VIMOS.2003-06-02T12:36:21.589.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:22:26.545/VIMOS.2003-06-03T11:22:26.545.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:31:43.663/VIMOS.2003-06-02T12:31:43.663.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:19:58.071/VIMOS.2003-06-03T11:19:58.071.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:25:50.846/VIMOS.2003-06-03T11:25:50.846.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:23:43.086/VIMOS.2003-06-02T19:23:43.086.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:27:16.647/VIMOS.2003-06-03T11:27:16.647.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:17.817/VIMOS.2003-06-03T11:21:17.817.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:27:06.394/VIMOS.2003-06-02T19:27:06.394.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:47.234/VIMOS.2003-06-03T11:21:47.234.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:23:04.866/VIMOS.2003-06-02T19:23:04.866.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:30:29.892/VIMOS.2003-06-02T12:30:29.892.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:37:13.006/VIMOS.2003-06-02T12:37:13.006.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:22:26.544/VIMOS.2003-06-03T11:22:26.544.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:31:43.662/VIMOS.2003-06-02T12:31:43.662.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:25:50.845/VIMOS.2003-06-03T11:25:50.845.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:19:20.059/VIMOS.2003-06-03T11:19:20.059.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:23:34.177/VIMOS.2003-06-02T19:23:34.177.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:24:50.854/VIMOS.2003-06-02T19:24:50.854.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:24:21.516/VIMOS.2003-06-02T19:24:21.516.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:19:49.472/VIMOS.2003-06-03T11:19:49.472.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:23:09.654/VIMOS.2003-06-03T11:23:09.654.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:23:04.865/VIMOS.2003-06-02T19:23:04.865.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:06.381/VIMOS.2003-06-03T11:21:06.381.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:25:00.010/VIMOS.2003-06-02T19:25:00.010.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:30:29.891/VIMOS.2003-06-02T12:30:29.891.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:25:29.353/VIMOS.2003-06-02T19:25:29.353.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:28:51.439/VIMOS.2003-06-02T19:28:51.439.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:29:44.361/VIMOS.2003-06-02T19:29:44.361.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:30:08.120/VIMOS.2003-06-02T19:30:08.120.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:26:20.261/VIMOS.2003-06-03T11:26:20.261.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:28:51.438/VIMOS.2003-06-02T19:28:51.438.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:24:21.515/VIMOS.2003-06-02T19:24:21.515.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:26:37.073/VIMOS.2003-06-02T19:26:37.073.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:32:57.518/VIMOS.2003-06-02T12:32:57.518.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:31:06.793/VIMOS.2003-06-02T12:31:06.793.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:24:50.853/VIMOS.2003-06-02T19:24:50.853.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:25:29.354/VIMOS.2003-06-02T19:25:29.354.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:06.380/VIMOS.2003-06-03T11:21:06.380.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:36:21.588/VIMOS.2003-06-02T12:36:21.588.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:20:36.973/VIMOS.2003-06-03T11:20:36.973.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:23:43.085/VIMOS.2003-06-02T19:23:43.085.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:32:57.517/VIMOS.2003-06-02T12:32:57.517.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:32:20.583/VIMOS.2003-06-02T12:32:20.583.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:25:38.518/VIMOS.2003-06-02T19:25:38.518.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:30:08.121/VIMOS.2003-06-02T19:30:08.121.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:26:20.262/VIMOS.2003-06-03T11:26:20.262.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:24:12.380/VIMOS.2003-06-02T19:24:12.380.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:30:37.451/VIMOS.2003-06-02T19:30:37.451.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:28:22.123/VIMOS.2003-06-02T19:28:22.123.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:26:37.074/VIMOS.2003-06-02T19:26:37.074.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:25:26.548/VIMOS.2003-06-03T11:25:26.548.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:27:29.546/VIMOS.2003-06-02T19:27:29.546.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:33:47.805/VIMOS.2003-06-02T12:33:47.805.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:26:47.234/VIMOS.2003-06-03T11:26:47.234.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:24:57.138/VIMOS.2003-06-03T11:24:57.138.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:29:15.041/VIMOS.2003-06-02T19:29:15.041.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:24:12.379/VIMOS.2003-06-02T19:24:12.379.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:20:36.974/VIMOS.2003-06-03T11:20:36.974.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:35:30.201/VIMOS.2003-06-02T12:35:30.201.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:57.131/VIMOS.2003-06-03T11:21:57.131.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:25:38.519/VIMOS.2003-06-02T19:25:38.519.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:32:20.584/VIMOS.2003-06-02T12:32:20.584.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:34:39.026/VIMOS.2003-06-02T12:34:39.026.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:24:02.318/VIMOS.2003-06-03T11:24:02.318.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:30:37.452/VIMOS.2003-06-02T19:30:37.452.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:19:49.473/VIMOS.2003-06-03T11:19:49.473.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:23:34.178/VIMOS.2003-06-02T19:23:34.178.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:33:47.806/VIMOS.2003-06-02T12:33:47.806.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:25:26.549/VIMOS.2003-06-03T11:25:26.549.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:28:22.122/VIMOS.2003-06-02T19:28:22.122.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:29:15.040/VIMOS.2003-06-02T19:29:15.040.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:35:30.202/VIMOS.2003-06-02T12:35:30.202.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:26:47.235/VIMOS.2003-06-03T11:26:47.235.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:29:44.362/VIMOS.2003-06-02T19:29:44.362.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:24:57.137/VIMOS.2003-06-03T11:24:57.137.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:25:00.011/VIMOS.2003-06-02T19:25:00.011.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:23:09.653/VIMOS.2003-06-03T11:23:09.653.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T12:34:39.027/VIMOS.2003-06-02T12:34:39.027.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:24:02.319/VIMOS.2003-06-03T11:24:02.319.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:23:39.091/VIMOS.2003-06-03T11:23:39.091.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:24:31.674/VIMOS.2003-06-03T11:24:31.674.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:27:16.646/VIMOS.2003-06-03T11:27:16.646.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:20:27.486/VIMOS.2003-06-03T11:20:27.486.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:26:07.846/VIMOS.2003-06-02T19:26:07.846.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:27:58.893/VIMOS.2003-06-02T19:27:58.893.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:27:29.545/VIMOS.2003-06-02T19:27:29.545.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:19:20.060/VIMOS.2003-06-03T11:19:20.060.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:24:31.675/VIMOS.2003-06-03T11:24:31.675.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:23:39.090/VIMOS.2003-06-03T11:23:39.090.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:19:58.072/VIMOS.2003-06-03T11:19:58.072.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:20:27.487/VIMOS.2003-06-03T11:20:27.487.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:26:07.845/VIMOS.2003-06-02T19:26:07.845.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-03T11:21:57.130/VIMOS.2003-06-03T11:21:57.130.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361651/SAF/VIMOS.2003-06-02T19:27:58.894/VIMOS.2003-06-02T19:27:58.894.fits.Z"

__EOF__
