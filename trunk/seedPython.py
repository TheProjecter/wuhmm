#!/usr/bin/env python

################################################################################
#
#       Copyright (C) Patrick Cahan 2007-2009
#
#       Contact: pcahan@wustl.edu
#
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#################################################################################


import sys

def overlap(str1, stp1, str2, stp2):
    if(( (str1<=stp2) & (stp1>=str2)) |
      ( (stp1>=str2) & (stp1<=stp2))):
        ans = 1
    else:
        ans=0

    return ans
   

def matchRank2(val,match,seedType):
    if seedType==2:
        if match==1:
            match=2
        else:
            match=match-1
    
    if val==match:
        return 1
    else:
        return 0


def seedSearch2(positions,ranks,minSeed,start,x):
    out = list()
    cont = 0
    stp = start
    num_probes = 1
    i = start+1;

    if i < len(ranks):
        if matchRank2(ranks[i],x,1):
            stp = i
            num_probes=num_probes+1
            i=i+1
            cont = 1

        elif i+1 < len(ranks):
            if matchRank2(ranks[i],x,2) & matchRank2(ranks[i+1],x,1):
                stp = i+1
                num_probes=num_probes+2
                i=i+2
                cont = 1

        if cont==1:
            while i < len(ranks):
                if matchRank2(ranks[i],x,1):
                    stp = i
                    num_probes=num_probes+1
                    i=i+1
                else:
                    break
            

    out.append(stp)
    out.append(positions[start])
    out.append(positions[stp])
    out.append(num_probes)
    return out


def write_seed(filename,start,stp,num_probes):
    filename.write(str(start)+"\t"+str(stp)+"\t"+str(num_probes)+"\n")


def expandSeeds(fname, positions, starts, stps, num_probes, fact):    
    # Expand seeds
    i = 0;
    while i < len(starts):
        marg = num_probes[i]*fact
        newStart = positions.index(starts[i]);
        
        if(newStart<=marg):
            starts[i] = positions[0]
        else:
            starts[i] = positions[newStart-marg]

        print str(newStart)+": "+str(starts[i])
        
        newStop = positions.index(stps[i])
        
        
        if(newStop >= (len(positions)-marg)):
            stps[i] = positions[len(positions)-1]
        else:
            stps[i] = positions[newStop+marg]

        print str(newStop)+": "+str(stps[i])

        i = i+1
    
    # merge seeds
    i = 0;
    while i < len(starts):
        if (i+1) >= len(starts):
            break;
        else:
            if(overlap(starts[i], stps[i], starts[i+1], stps[i+1])):
                starts.pop(i+1);
                stps[i]=stps[i+1];
                stps.pop(i+1);
                print str(starts[i])+" "+str(stps[i])
            else:
                i=i+1

    # write new seeds
    write_seed(fname, "str", "stp", "num_probes")
    i = 0;
    while i < len(starts):
        write_seed(fname, starts[i], stps[i], 3)
        i = i + 1

nclusts = int(sys.argv[1]);
minSeed = int(sys.argv[2]);

ranks = list();
positions = list();
clustsfile = open("./ngForSeeds.csv",'r')

for line in clustsfile:
    l=line.split("\t");
    ranks.append(int(l[2]));
    positions.append(float(l[0]));
    
fname_out = open("mySeeds.csv",'w')
seedStarts = list()
seedStops  = list()
seedNProbes = list()

i = 0
while i < len(positions):
    if ranks[i]==1:
        x = 1
        out = seedSearch2(positions,ranks,minSeed,i,x)
        i = out[0]+1
        if(out[3]>=minSeed):
            seedStarts.append(out[1])
            seedStops.append(out[2])
            seedNProbes.append(out[3])
            #print str(out[0])+" "+str(out[1])+" "+str(out[2])+" "+str(out[3])+"\n" 
        
    elif ranks[i]==nclusts:        
        x = nclusts
#        out = seedSearch(chrs[i],positions,ranks,minSeed,i,x, nclusts)
        out = seedSearch2(positions,ranks,minSeed,i,x)
        i = out[0]+1
        if(out[3]>=minSeed):
            seedStarts.append(out[1])
            seedStops.append(out[2])
            seedNProbes.append(out[3])
            #print str(out[0])+" "+str(out[1])+" "+str(out[2])+" "+str(out[3])+"\n" 
        
    else:
        i = i + 1

#expandSeeds(fname_out,positions, seedStarts, seedStops, seedNProbes,10)     

# write new seeds
write_seed(fname_out, "str", "stp", "num_probes")
i = 0;
while i < len(seedStarts):
    write_seed(fname_out, seedStarts[i], seedStops[i], seedNProbes[i])
    i = i + 1


