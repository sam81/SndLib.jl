
echo -n "FTP password: "
read passwd
ncftpput -R -v -u "samcarcagno" -p "$passwd" ftp.samcarcagno.altervista.org SndLib/site/ ../site/* 
