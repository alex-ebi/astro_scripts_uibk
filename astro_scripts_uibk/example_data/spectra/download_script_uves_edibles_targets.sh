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

hostname=dataportal.eso.org
username=anonymous
anonymous=
xargsopts=
prompt=
list=
while getopts hlpu:xX:d: option
do
    case $option in
	h) usage; exit ;;
	l) list=yes ;;
	p) prompt=yes ;;
	u) prompt=yes; username="$OPTARG" ;;
	X) xargsopts="$OPTARG" ;;
	d) download_opts="$OPTARG";;
	?) usage; exit 2 ;;
    esac
done

if [ "$username" = "anonymous" ]; then
    anonymous=yes
fi

if [ -z "$xargsopts" ]; then
    #no xargs option specified, we ensure that only one url
    #after the other will be used
    xargsopts='-L 1'
fi

netrc=$HOME/.netrc
if [ -z "$anonymous" -a -z "$prompt" ]; then
    # take password (and user) from netrc if no -p option
    if [ -f "$netrc" -a -r "$netrc" ]; then
	grep -ir "$hostname" "$netrc" > /dev/null
	if [ $? -ne 0 ]; then
            #no entry for $hostname, user is prompted for password
            echo "A .netrc is available but there is no entry for $hostname, add an entry as follows if you want to use it:"
            echo "machine $hostname login anonymous password _yourpassword_"
            prompt="yes"
	fi
    else
	prompt="yes"
    fi
fi

if [ -n "$prompt" -a -z "$list" ]; then
    trap 'stty echo 2>/dev/null; echo "Cancelled."; exit 1' INT HUP TERM
    stty -echo 2>/dev/null
    printf 'Password: '
    read password
    echo ''
    stty echo 2>/dev/null
    escaped_password=${password//\%/\%25}
    auth_check=$(wget -O - --post-data "username=$username&password=$escaped_password" --server-response --no-check-certificate "https://www.eso.org/sso/oidc/accessToken?grant_type=password&client_id=clientid" 2>&1 | awk '/^  HTTP/{print $2}')
    if [ ! $auth_check -eq 200 ]
    then
        echo 'Invalid password!'
        exit 1
    fi
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
echo "check_certificate=off" >> $tempfile
echo "content_disposition=on" >> $tempfile
echo "continue=on" >> $tempfile
if [ -z "$anonymous" -a -n "$prompt" ]; then
    echo "http_user=$username" >> $tempfile
    echo "http_password=$password" >> $tempfile
fi
WGETRC=$tempfile; export WGETRC

unset password

if [ -n "$list" ]; then
    cat
else
    xargs $xargsopts wget $download_opts 
fi <<'__EOF__'
https://archive.eso.org/downloadportalapi/readme/d0d2d8f8-a56c-4751-9c30-0b5c5a5b20ca
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:13.629
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.427
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.267
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:49.804
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.458
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.693
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.572
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:49.127
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.967
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T18:20:31.641
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:29.969
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:47:53.640
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T18:20:31.885
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:42.814
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:01.540
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:43.071
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:52:30.016
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T13:46:31.412
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:48.892
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T18:20:32.220
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.580
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T18:10:38.901
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.831
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:31.314
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.835
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T08:41:12.198
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:57.893
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:49.969
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:42.822
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.899
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.779
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:58.215
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:48.659
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:45:30.676
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:20.475
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T12:20:58.489
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.350
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:49.309
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.128
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.598
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:01.327
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.765
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:42.839
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T16:59:59.248
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:47:53.623
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T09:10:02.290
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:59.016
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.432
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:33:10.476
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T17:00:01.313
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.179
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:27:50.720
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T18:58:41.265
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T15:54:10.759
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:52:30.717
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:01.576
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:58.128
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:02.256
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:00.086
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.509
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:49.749
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.856
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:50.043
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:45:30.657
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:58.436
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:48.704
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.569
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:01.793
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:48.639
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:29.638
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:42.601
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.004
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T08:41:12.941
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.328
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.824
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:01.744
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:58.054
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:49.042
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:48.725
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:48.805
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:15.392
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.336
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.113
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.879
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:33.117
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:57.214
https://dataportal.eso.org/dataportal_new/file/ADP.2021-09-18T18:19:12.174
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:14.363
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:33:10.022
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:15:18.419
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T09:10:02.761
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T17:47:32.115
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:11.665
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T12:20:59.095
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:02.804
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.694
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:19.812
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.570
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:58.372
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.168
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.200
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:14.497
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.947
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:58.497
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.247
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:19.727
https://dataportal.eso.org/dataportal_new/file/ADP.2020-10-30T17:20:51.520
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.939
https://dataportal.eso.org/dataportal_new/file/ADP.2020-10-30T17:20:51.404
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:33:09.756
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T09:10:03.228
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.898
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.284
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.611
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:58.845
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T17:47:32.129
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:14.268
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T14:26:52.995
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:49.038
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:58.363
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T08:41:12.180
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.362
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.559
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:15:18.315
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T19:20:22.026
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:26.716
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:48.914
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:59.286
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:32.890
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:13.592
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:01.933
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:27.257
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T19:13:23.440
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:42.657
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:00.329
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:50.330
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:58.503
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T19:20:22.650
https://dataportal.eso.org/dataportal_new/file/ADP.2020-11-02T08:37:23.461
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T12:20:59.227
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:19.052
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.181
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.558
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T14:26:53.024
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:27.363
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:57.496
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:00.220
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:50.453
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:35.996
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T19:20:22.522
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:48.737
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.389
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T08:41:12.801
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.949
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:01.038
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:31.249
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:32.979
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T08:41:12.247
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.687
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:20.607
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:49.069
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:59.414
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:00.587
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:27.283
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T14:26:53.003
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:15.410
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:58.761
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.717
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:15:18.192
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.439
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.713
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:02.073
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:57.594
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:11.729
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:02.096
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.311
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:02.862
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T18:10:39.558
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:59.762
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T14:26:49.999
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T17:00:01.161
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:10.609
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:47:46.037
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:27:50.685
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:52:30.685
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:19.300
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:25:21.100
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.906
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:31.623
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:18.855
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:50.404
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T09:10:02.357
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.772
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T12:20:59.145
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:31.180
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:47:46.088
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:35.840
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:20.305
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:26.882
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.734
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:33:10.777
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.910
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:36.873
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:58.225
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:31.515
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.780
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:43.013
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:58.748
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:43.373
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:57.853
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:15:18.135
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.918
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.302
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:00.814
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:01.290
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T16:02:58.410
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.062
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:20.437
https://dataportal.eso.org/dataportal_new/file/ADP.2020-10-30T17:20:51.237
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:18.189
https://dataportal.eso.org/dataportal_new/file/ADP.2021-09-18T18:19:12.326
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.181
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:33.070
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T17:44:19.416
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:45:30.595
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-12T12:20:59.137
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:35.954
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:33:10.784
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T16:26:09.034
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:58.797
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T09:18:22.922
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:32:13.649
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:26.670
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T14:26:53.057
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:05:51.768
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-08T14:42:27.472
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:03.208
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:57.748
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:43.602
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:02.317
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T09:10:02.251
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:57.633
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:24:57.077
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:48.989
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T06:48:43.293
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:52:30.115
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:34.631
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:03:48.861
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:45:30.690
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:52:29.946
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:48.791
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:15:17.921
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:06:58.882
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:02.174
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:52:09.708
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:33.298
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T14:14:35.443
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T09:10:02.385
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-16T09:07:02.785
https://dataportal.eso.org/dataportal_new/file/ADP.2020-11-02T08:44:55.725
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:47:46.193
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T08:35:49.790
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T17:47:33.328
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T17:26:00.229
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T15:05:30.136
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T07:47:46.997
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T15:54:10.958
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T10:47:54.994
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-09T11:08:35.979
https://dataportal.eso.org/dataportal_new/file/ADP.2020-06-10T18:00:48.679
__EOF__
