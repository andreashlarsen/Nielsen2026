## prepare
#sudo apt-get install keyutils
#sudo apt-get install cifs-utils

## mount personal drive
if [ -d "mnt" ]
then
  echo trying to make mount folder "mnt", but it already exists - all content of mnt would be removed, so no mounting taking place
  echo edit this script: mnt_personal_ku_drive.sh, or rename your existing folder "mnt" 
else
  mkdir -p mnt
  echo did not find directory mnt
  sudo mount -t cifs -o username=plt572,uid=andreas,sec=ntlmssp,iocharset=utf8,domain=unicph.domain,dir_mode=0700 //unicph.domain/users/plt572 /home/andreas/mnt/
fi

## unmount (recommended before shutting down system)
#sudo unmount /home/andreas/mnt/
