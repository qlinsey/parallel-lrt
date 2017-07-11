#!/usr/bin/bash
#./lrt grid500 -m 4 -x 8 -y 32  -p &>o.txt
for g in {100,200,400,600,800,1000} #grid
        do
                for x in {8,16,32,64,128,256} #thread
                do
                        for y in {32,64,256,1024}
    do
          for i in {1..6}
           do
                   echo "start $i th measuring  grid $g, blocks="$x ",$y "
                 #echo "start $i th measuring multi-core-lrt  grid"$g", threads="$t
              ./lrt grid$g -m 4 -x $x -y $y -p #&>o.txt
         done
 done
    done
 done
