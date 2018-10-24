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
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:31:25.407/VIMOS.2003-06-03T01:31:25.407.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:24:15.394/VIMOS.2003-06-02T23:24:15.394.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:33:53.194/VIMOS.2003-06-03T01:33:53.194.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:38:18.702/VIMOS.2003-06-02T23:38:18.702.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:32:08.623/VIMOS.2003-06-03T01:32:08.623.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:54:26.441/VIMOS.2003-06-02T23:54:26.441.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:31:07.301/VIMOS.2003-06-02T23:31:07.301.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:44:10.439/VIMOS.2003-06-03T09:44:10.439.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:22:41.860/VIMOS.2003-06-02T23:22:41.860.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:53:32.527/VIMOS.2003-06-03T09:53:32.527.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:37:33.630/VIMOS.2003-06-03T09:37:33.630.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:14:35.763/VIMOS.2003-06-03T01:14:35.763.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:25:45.574/VIMOS.2003-06-02T23:25:45.574.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:59:14.081/VIMOS.2003-06-02T23:59:14.081.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:24:42.760/VIMOS.2003-06-02T23:24:42.760.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:19:44.929/VIMOS.2003-06-03T10:19:44.929.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T00:44:29.325/VIMOS.2003-06-03T00:44:29.325.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:57:39.135/VIMOS.2003-06-02T23:57:39.135.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:23:09.220/VIMOS.2003-06-02T23:23:09.220.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:46:42.366/VIMOS.2003-06-02T23:46:42.366.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:33:50.366/VIMOS.2003-06-03T09:33:50.366.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:21:37.088/VIMOS.2003-06-02T23:21:37.088.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:16:49.603/VIMOS.2003-06-03T10:16:49.603.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:29:42.774/VIMOS.2003-06-03T01:29:42.774.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:31:25.406/VIMOS.2003-06-03T01:31:25.406.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:33:23.817/VIMOS.2003-06-03T01:33:23.817.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:38:18.701/VIMOS.2003-06-02T23:38:18.701.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:51:45.192/VIMOS.2003-06-03T09:51:45.192.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T00:57:23.075/VIMOS.2003-06-03T00:57:23.075.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:33:53.195/VIMOS.2003-06-03T01:33:53.195.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:24:15.393/VIMOS.2003-06-02T23:24:15.393.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:37:04.214/VIMOS.2003-06-03T09:37:04.214.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:14:06.387/VIMOS.2003-06-03T01:14:06.387.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:32:08.624/VIMOS.2003-06-03T01:32:08.624.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:25:45.573/VIMOS.2003-06-02T23:25:45.573.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:22:41.861/VIMOS.2003-06-02T23:22:41.861.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:53:32.526/VIMOS.2003-06-03T09:53:32.526.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:32:38.001/VIMOS.2003-06-03T01:32:38.001.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:24:42.761/VIMOS.2003-06-02T23:24:42.761.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:16:20.187/VIMOS.2003-06-03T10:16:20.187.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T00:44:29.324/VIMOS.2003-06-03T00:44:29.324.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:43:41.029/VIMOS.2003-06-03T09:43:41.029.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:46:42.367/VIMOS.2003-06-02T23:46:42.367.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:56:33.504/VIMOS.2003-06-02T23:56:33.504.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:30:39.932/VIMOS.2003-06-02T23:30:39.932.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:47:49.782/VIMOS.2003-06-02T23:47:49.782.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:21:09.719/VIMOS.2003-06-02T23:21:09.719.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:29:42.775/VIMOS.2003-06-03T01:29:42.775.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:21:13.245/VIMOS.2003-06-03T10:21:13.245.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:48:17.159/VIMOS.2003-06-02T23:48:17.159.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:32:38.002/VIMOS.2003-06-03T01:32:38.002.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:33:23.816/VIMOS.2003-06-03T01:33:23.816.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:51:45.191/VIMOS.2003-06-03T09:51:45.191.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T00:57:23.076/VIMOS.2003-06-03T00:57:23.076.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:17:47.525/VIMOS.2003-06-03T10:17:47.525.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:37:04.213/VIMOS.2003-06-03T09:37:04.213.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:14:06.386/VIMOS.2003-06-03T01:14:06.386.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:18:16.939/VIMOS.2003-06-03T10:18:16.939.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:30:56.031/VIMOS.2003-06-03T01:30:56.031.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:23:09.219/VIMOS.2003-06-02T23:23:09.219.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:48:17.158/VIMOS.2003-06-02T23:48:17.158.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:41:18.472/VIMOS.2003-06-03T09:41:18.472.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:05:49.426/VIMOS.2003-06-03T01:05:49.426.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:59:43.450/VIMOS.2003-06-02T23:59:43.450.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:21:09.720/VIMOS.2003-06-02T23:21:09.720.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:58:08.503/VIMOS.2003-06-02T23:58:08.503.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:56:33.503/VIMOS.2003-06-02T23:56:33.503.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:45:49.853/VIMOS.2003-06-03T09:45:49.853.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:30:39.933/VIMOS.2003-06-02T23:30:39.933.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:47:49.783/VIMOS.2003-06-02T23:47:49.783.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:21:13.244/VIMOS.2003-06-03T10:21:13.244.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:43:41.030/VIMOS.2003-06-03T09:43:41.030.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:16:20.188/VIMOS.2003-06-03T10:16:20.188.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:44:10.440/VIMOS.2003-06-03T09:44:10.440.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:46:14.992/VIMOS.2003-06-02T23:46:14.992.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:54:55.819/VIMOS.2003-06-02T23:54:55.819.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:37:33.631/VIMOS.2003-06-03T09:37:33.631.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:17:47.526/VIMOS.2003-06-03T10:17:47.526.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:54:26.442/VIMOS.2003-06-02T23:54:26.442.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T00:34:56.369/VIMOS.2003-06-03T00:34:56.369.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:31:07.302/VIMOS.2003-06-02T23:31:07.302.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:14:35.764/VIMOS.2003-06-03T01:14:35.764.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:18:16.938/VIMOS.2003-06-03T10:18:16.938.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:59:43.449/VIMOS.2003-06-02T23:59:43.449.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:30:56.030/VIMOS.2003-06-03T01:30:56.030.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:58:08.504/VIMOS.2003-06-02T23:58:08.504.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:59:14.082/VIMOS.2003-06-02T23:59:14.082.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T10:16:49.604/VIMOS.2003-06-03T10:16:49.604.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T01:05:49.427/VIMOS.2003-06-03T01:05:49.427.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:41:18.473/VIMOS.2003-06-03T09:41:18.473.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:45:49.854/VIMOS.2003-06-03T09:45:49.854.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-03T09:33:50.367/VIMOS.2003-06-03T09:33:50.367.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:21:37.087/VIMOS.2003-06-02T23:21:37.087.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:54:55.818/VIMOS.2003-06-02T23:54:55.818.fits.Z"
"https://dataportal.eso.org/dataPortal/api/requests/mmarcano22/361656/SAF/VIMOS.2003-06-02T23:46:14.993/VIMOS.2003-06-02T23:46:14.993.fits.Z"

__EOF__
