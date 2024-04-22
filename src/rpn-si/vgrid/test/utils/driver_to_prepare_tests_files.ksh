#!/bin/bash

set -e

cat > directives <<EOF
 SORTIE(STD,1000,A)
 GRILLE(PS,30,30,-140.,310.,15000.,20.,NORD)
 HEURE(-1)
 CHAMP(TT)
 CHAMP(P0)
 END
EOF

DATE=20100608

CHEMIN[1]=lam/east.eta; HEURE[1]=12; SUFFIX[1]=_2.5km
CHEMIN[2]=reghyb;       HEURE[2]=12; SUFFIX[2]=
CHEMIN[3]=regeta;       HEURE[3]=12; SUFFIX[3]=
CHEMIN[4]=glbhyb;       HEURE[4]=12; SUFFIX[4]=
CHEMIN[5]=glbeta;       HEURE[5]=12; SUFFIX[5]=

for I in $(seq 1 5)
do

   rm -f o
   pgsm -iment ${CMCGRIDF}/prog/${CHEMIN[I]}/${DATE}${HEURE[I]}_000${SUFFIX[I]}\
        -ozsrt o -i directives

   rm -f px
   r.hy2pres -iment o -ozsrt px -var TT 
   rm -f data/${CHEMIN[I]##*/}
   editfst -s o px -d data/${CHEMIN[I]##*/} -i 0

done

 
