
cd UMBRELLA

tpr=tpr_files_umbrella.dat
pullx=pullx_files_umbrella.dat
pullf=pullf_files_umbrella.dat
rm $tpr $pullx $pullf

for i in  27  40  42  63  69  80  94 100 117 121 138 145 148 154 162 178 181 186 197 204 206 212 221 224 226 229 240 243 246 248 251 261 266 269 271 272 276 288 289 291 294 301 302 303 308 312 317 320 321 322 323 326 327 339 340 349 354 356 357 366 370 371 374 377 
do
    echo "umbrella${i}.tpr" >> $tpr
    echo "umbrella${i}_pullf.xvg" >> $pullf
done


gmx wham -it $tpr -if $pullf -hist -temp 310 -b 1000 -quiet 
#gmx wham -it $tpr -if $pullf -hist -temp 310 -b 50000 -nBootstrap 200 -quiet

xmgrace -nxy histo.xvg
xmgrace profile.xvg
#xmgrace bsResult.xvg

cd ..

