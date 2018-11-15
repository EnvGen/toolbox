# -n option disables auto-logon
FILE_SUFFIX="fastq.gz"
for upload_dir in "$@"
do
ftp -n webin.ebi.ac.uk <<End-Of-Session
user $ENA_USER $ENA_PASS
binary
prompt
lcd "$upload_dir"
!ls *"$FILE_SUFFIX"
mput *"$FILE_SUFFIX"
bye
End-Of-Session
done
