#!/bin/sh 
#change file to nexus format
cd ./output/birds/trees/
awk '{print " = [&R] " $0;}' ./BirdzillaEricson_100Trees.tree > ./BirdzillaEricson_100Trees.tree2
awk '{print NR,$0}' ./BirdzillaEricson_100Trees.tree2 > ./BirdzillaEricson_100Trees.tree3
awk '{print "\ttree tree_",$0}' ./BirdzillaEricson_100Trees.tree3 > ./BirdzillaEricson_100Trees.tree4
echo 'begin trees;' | cat  - ./BirdzillaEricson_100Trees.tree4 > temp && mv temp ./BirdzillaEricson_100Trees.tree5
echo '#NEXUS' | cat  - ./BirdzillaEricson_100Trees.tree5 > temp && mv temp ./BirdzillaEricson_100Trees.tree6
echo "end;" >> ./BirdzillaEricson_100Trees.tree6
mv ./BirdzillaEricson_100Trees.tree6 ./BirdzillaEricson_100Trees.nex
sed -i -e 's/tree_ /tree_/' ./BirdzillaEricson_100Trees.nex
#delete intermediate files
rm ./BirdzillaEricson_100Trees.nex-e
rm ./BirdzillaEricson_100Trees.tree[0-9]
#compressing 
