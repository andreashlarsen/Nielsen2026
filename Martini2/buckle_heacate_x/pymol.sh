i=303 # 27  40  42  63  69  80  94 100 117 121 138 145 148 154 162 178 181 186 197 204 206 212 221 224 226 229 240 243 246 248 251 261 266 269 271 272 276 288 289 291 294 301 302 303 308 312 317 320 321 322 323 326 327 339 340 349 354 356 357 366 370 371 374 377

cd UMBRELLA

cat << EOF > script.pml
load conf$i.gro
load_traj umbrella$i.xtc
remove resname W
remove resname WF
remove resname ION
show spheres
extract LIP, resname POPC resname POPR resname POPS resname CHOL
color orange, conf$i 
show cell
set orthoscopic, on
EOF

pymol script.pml
rm script.pml

cd ..
