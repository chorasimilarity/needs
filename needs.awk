BEGIN { 

# this variant 04.08.2017
#
# reduces genes.perm, permutation style, generalization of  http://imar.ro/~mbuliga/chemlambda_strains.html
#
# owner: chorasimilarity (Marius Buliga, http://chorasimilarity.wordpress.com/,  marius.buliga@imar.ro , marius.buliga@gmail.com)
#

  colors_node_types(); 
  markers();

       }

{

     if (NF == 7) {  if ($2 != voidc) {succ[$1]=$2;} 
                     if ($3 != voidc ) {ccus[$1]=$3;} 
                     if ($4 != voidc) {gamma[$1]=$4;} 
                     if ($5 != voidc) {ammag[$1]=$5;} 
                     if ($6 != voidc) {edge[$1]++;}
                     if ($7 != voidc) {col[$1]=$7;} 
                     tag[$1]=""; } 
     else { numobj=$1; voidc=$2; edgec=$3;}

}

END { 

#### How many cycles? 

   cyclecount=100;

#### Things we need

   patterns(vecsep, vecint);
   colpatterns(vecsep, vecint);
   color_id(supersep, vecint);
   move_needs(vecint);
   reactions(vecint);

#### Random seed

   random_seed();

#### Definition of "tag". Mind that "col" array is made of names, with the exception of edges which do not have "col"

  for (a in tag) { if (a in ammag) { tag[a]=tag[a] "1";} else {tag[a]=tag[a] "0";}   
                   if (a in gamma) { tag[a]=tag[a] "1";} else {tag[a]=tag[a] "0";}
                   if (a in succ)  { tag[a]=tag[a] "1";} else {tag[a]=tag[a] "0";}
                   if (a in ccus)  { tag[a]=tag[a] "1";} else {tag[a]=tag[a] "0";}
                   if (a in col)   { tag[a]=tag[a] "1";} else {tag[a]=tag[a] "0";}
                 }

#### Stacks of ingredients

   for (a in need) { split(need[a], tempo, vecint); for (b in tempo) { c=tempo[b]; stack[c]=""; } }

#### Stack of moves

   movestack="";

#### Unseen edges initialization

   for (a in edge) {unseen[a]=0;}

####    - Main cycle: 

      cntactive=1;

      for (cyclec = 1; cyclec <= cyclecount; cyclec++) {

####             - movestack is split into a vector move_stack 

      lenmov=split(movestack, move_stack, vecsep); 
      la=length(missing); if (la!=0) {for (a = 1; a <= la; a++) {delete missing[a];}}
      cntmiss=0;
      cntreact=0;

####  add a break if movestack is empty and there's no missing ingredients 

      if (cntactive==0) {break}

####  if there are moves to perform

      if (lenmov!=0) { 


      for (nmov = 1; nmov <= lenmov; nmov++) {
                  split(move_stack[nmov],npair,vecint);
                  nod=npair[1];
                  it=npair[2];

                  lneed=split(need[it],stuff,vecint); cntneed=0; permpat="";

                  if (lneed!=0) { for (foru = 1; foru <= lneed; foru++) {kuku=stuff[foru]; 
                                                                         gugu=stack[kuku];
                                                                         if (gugu!=""){cntneed++;} else { cntmiss++; missing[cntmiss]=kuku; }
                                                                        }
                                }

####             - the reaction is attempted by completing what's needed from the top of stacks of ingredients. If the completion is successful then ...

                  if (lneed==cntneed) { if (lneed!=0) { for (foru = 1; foru <= lneed; foru++) {kuku=stuff[foru]; 
                                                                         gugu=stack[kuku];
                                                                         gaga=popleft(gugu, supersep, vecsep);
                                                                         split(gaga, fifi, supersep);
                                                                         nodo=fifi[1];
                                                                         stack[kuku]=fifi[1]; 
                                                                         split(seen_node[nodo],dpair,supersep);
                                                                         permpat=permpat dpair[2] vecsep;
                                                                                               }
                                                      } else {permpat=vecsep;}

                                       split(seen_node[nod],dpair,supersep);
                                       permpat=dpair[2] vecsep permpat;
                                       permpat=substr(permpat,1,length(permpat)-1);
####             - check if the complete pattern is a permutation, if so ... 

                                       du=isperm(permpat, vecsep, vecint);
                                       if (du==0) {

####             - the  block_up and block_down unseen and seen_node are updated

                                       lpat=split(permpat,permut,vecsep);
                                       for (ael = 1; ael <= lpat; ael++) { split(permut[ael],succsee,vecint); 
                                                                           succ_up[ael]=succsee[1]; delete block_up[succ_up[ael]]; 
                                                                           if (succ_up[ael] in edge) { unseen[succ_up[ael]]=0; 
                                                                                                       if (succ_up[ael] in seen_node) 
                                                                                                          {delete seen_node[succ_up[ael]]; }}
                                                                           succ_down[ael]=succsee[2]; delete block_down[succ_down[ael]];
                                                                           if (succ_down[ael] in edge) { unseen[succ_down[ael]]=0; 
                                                                                                       if (succ_down[ael] in seen_node) 
                                                                                                          {delete seen_node[succ_down[ael]]; }}
                                                                         }

####             - the reaction is performed 

                   split(react[it],dothemove,vecint); 
                   for (ael = 1; ael <= lpat; ael++) { fu=succ_up[ael]; fd=succ_down[dothemove[ael]]; succ[fu]=fd; ccus[fd]=fu; }
                   cntreact++;
                                                  }
                                      }
                  else { unseen[nod]=0; delete seen_node[nod]; delete move_stack[nmov]; }
             
#end  for (nmov = 1; nmov <= lenmov; nmov++) 
                                                             }  
          movestack="";
         
# end if (lenmov!=0)
                     } 

        cntactive=cntmiss+cntreact;

####    - add the missing ingredients with add_ingr

        if (cntmiss!=0) { for (ael = 1; ael <= cntmiss; ael++) { s=missing[ael]; numobj=add_ingr(numobj,s, vecsep, vecint);  } }

####    - see jumps randomly for a while and builds the seen_node, the stacks of ingredients and moves

####                - some method for exploring the edges which gives a nod node

   unseenstack=shuffle(unseen, vecint); if (unseenstack!="") {split(unseenstack, shuffled_edge,vecint);}

####                - populate stacks with function see
  
   for (nodi in shuffled_edge) { nod=shuffled_edge[nodi]; 
                                    if (nod in unseen) { for (patte in pattern) { a=see(nod, patte, supersep, vecsep, vecint); 
                                                                    if (a!="error") { break}
                                                                                }
                                                       }
                                  }

# end       for (cyclec = 1; cyclec <= cyclecount; cyclec++)
}

# end END
}
