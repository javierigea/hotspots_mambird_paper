#!/bin/sh 
cd ./raw_data/bird_trees_posterior/
#get 1000 random lines (trees)
gshuf -n 1000 ./BirdzillaEricsonAllTrees.tre > ./BirdzillaEricson_1000Trees.tree
#change file to nexus format
awk '{print " = [&R] " $0;}' ./BirdzillaEricson_1000Trees.tree > ./BirdzillaEricson_1000Trees.tree2
awk '{print NR,$0}' ./BirdzillaEricson_1000Trees.tree2 > ./BirdzillaEricson_1000Trees.tree3
awk '{print "\ttree tree_",$0}' ./BirdzillaEricson_1000Trees.tree3 > ./BirdzillaEricson_1000Trees.tree4
echo 'begin trees;' | cat  - ./BirdzillaEricson_1000Trees.tree4 > temp && mv temp ./BirdzillaEricson_1000Trees.tree5
echo '#NEXUS' | cat  - ./BirdzillaEricson_1000Trees.tree5 > temp && mv temp ./BirdzillaEricson_1000Trees.tree6
echo "end;" >> ./BirdzillaEricson_1000Trees.tree6
mv ./BirdzillaEricson_1000Trees.tree6 ./BirdzillaEricson_1000Trees.nex
sed -i -e 's/tree_ /tree_/' ./BirdzillaEricson_1000Trees.nex
#delete intermediate files
rm ./BirdzillaEricson_1000Trees.nex-e
rm ./BirdzillaEricson_1000Trees.tree[0-9]
#compressing 
gzip ./BirdzillaEricson_1000Trees.tree
gzip ./BirdzillaEricsonAllTrees.tre
