#---------------------- Insert ID -------------------------#
awk 'BEGIN{FS=","; OFS=","}
function custom_sort(i1,v1,i2,v2){
         l1=length(v1); l2=length(v2);
         if (l1 == l2) { 
             return (v1 > v2)? 1:-1
         } else { 
             return l2 - l1
         }   
     } NR==1 { $0=$0 FS "ID" }
       NR>1 { a[1]=$11; a[2]=$12; asort(a,b,"custom_sort");
       $(NF+1) = sprintf("%s:%s_%s_%s",$1,$2,b[1],b[2])
     }1' T2D_european_noID.txt > T2D_european.txt | column -t


#------------------ Insert ID in multiple files ---------------------#
for i in basalganglia-EUR-30 cortex-AFR-40 cortex-EAS-30 cortex-EUR-80 hippocampus-EUR-30;
do awk 'BEGIN{FS=","; OFS=","}
function custom_sort(i1,v1,i2,v2){
         l1=length(v1); l2=length(v2);
         if (l1 == l2) { 
             return (v1 > v2)? 1:-1
         } else { 
             return l2 - l1
         }   
     } NR==1 { $0=$0 FS "ID" }
       NR>1 { a[1]=$8; a[2]=$9; asort(a,b,"custom_sort");
       $(NF+1) = sprintf("%s:%s_%s_%s",$10,$11,b[1],b[2])
     }1' eQTL_${i}_hg19.csv > eQTL_${i}.csv;
done



