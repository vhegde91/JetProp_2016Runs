#!/bin/sh    

for i in B C D E F G H
do
    echo $i
#    grep 'Run2016'$i'.*03Feb2017.*HTMHT' ~/inputFileList_main.txt > Run2016$i\_03Feb2017_HTMHT.txt
    mv Run2016$i\_03Feb2017_HTMHT.root Run2016$i\_03Feb2017_HTMHT_Incl.root
done