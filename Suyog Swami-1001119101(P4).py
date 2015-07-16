#-----------------------------------------------------------------
#Name: Suyog Swami 
#hierarchical agglomerative clustering (HAC) Implementation
#Data Mining- Programming Assignment 4
#------------------------------------------------------------------
import csv
import math
from matplotlib.pyplot import show,plot
from scipy.cluster.hierarchy import dendrogram
from operator import  itemgetter
       
legen=[]
ln={}
init_cluster=[]
init=[]
def main():
    
    Z,SSE,label_array,K=hac("P4csegrad.csv","complete")   
    draw_dendogram(Z,legen,None,True)
    mention_x=5      # Specify k min for sses
    mention_y=39     # Specify k max for sse
    if SSE==[] and K==[]:
        print("Sorry SSE chart cannot be calculated for single linakge or complete linkage method as we dont cancluate centroid")
    else:
        sse(SSE,K,mention_x,mention_y)
    
def add_one(init_cluster):
    for a_in in init_cluster:
        a_in.append(1)
    
def hac(filename,method):
    with open(filename,'r') as csvfile:
        readers=csv.reader(csvfile)
        i=0 
        for r in readers:
            if i==0:
                i=i+1
                pass
            else:
                legen.append(r[0])
                ln[r[0]]=1
                init_cluster.append(r)
    le=legen.copy()
    legend=legen.copy()
    init=init_cluster.copy()
    add_one(init_cluster)
    if method=='complete':
        Z,SSE,K=complete_linkage(init_cluster,init,legend,le,ln)
        return(Z,SSE,legen,K)
    if method=='single':
        Z,SSE,K=single_linkage(init_cluster,init,legend,le,ln)
        return(Z,SSE,legen,K)        
    if method=='centroid':
        Z,SSE,K=centroid(init_cluster,init,legend,le)
        return(Z,SSE,legen,K)
    if method=='average':
        Z,SSE,K=centroid(init_cluster,init,legend,le)
        return(Z,SSE,legen,K)

def single_linkage(init_cluster,init,legend,le,ln):
    l=len(init_cluster)
    euclid_list=[]
    euclid_dict={}
    for ini in init_cluster:
        euclidean_dist=0
        ini_index=init_cluster.index(ini)
        for j in range(init_cluster.index(ini)+1,l):
            inj=init_cluster[j]
            inj_index=init_cluster.index(inj)
            euclidean_dist=math.sqrt((float(ini[1])-float(inj[1]))**2+\
            (float(ini[2])-float(inj[2]))**2+\
            (float(ini[3])-float(inj[3]))**2+\
            (float(ini[4])-float(inj[4]))**2+\
            (float(ini[5])-float(inj[5]))**2)
            euclid_list.append([legend[ini_index],legend[inj_index],euclidean_dist,ln.get(inj[0])+ln.get(ini[0])])
            euclid_dict[legend[ini_index],legend[inj_index]]=euclidean_dist
    Z=[]
    o=0
    x=0
    while o<l-1:
        k=2
        todelete=[]
        new_club_cluster=[]
        el=[]
        new_club_cluster=(sorted(euclid_list,key=itemgetter(2),reverse=False)[0])
        z=[legend.index(new_club_cluster[0]),legend.index(new_club_cluster[1]),new_club_cluster[2],ln.get(new_club_cluster[0])+ln.get(new_club_cluster[1])]
        if z[0]>l+1 :
            z[0]=new_club_cluster[0]
        if z[1]>l+1:
            z[1]=new_club_cluster[1]
        ln[l+o]=ln.get(new_club_cluster[0])+ln.get(new_club_cluster[1])        
        todelete.append(new_club_cluster) 
        for i in euclid_list:
            if i[1]==new_club_cluster[0] :
                legend.append(l+x)                
                el.append([i[0],l+x,(min(euclid_dict[i[0],i[1]],euclid_dict[i[0],new_club_cluster[1]])),ln.get(i[0])+ln.get(l+o)])
                euclid_dict[i[0],l+x]=min(euclid_dict[i[0],i[1]],euclid_dict[i[0],new_club_cluster[1]])
                todelete.append(i)
                todelete.append([i[0],new_club_cluster[1],euclid_dict[i[0],new_club_cluster[1]],ln.get(i[0])+ln.get(new_club_cluster[1])])

            elif i[0]==new_club_cluster[0] and legend.index(i[1])>legend.index(new_club_cluster[1]):
                legend.append(l+x)                
                el.append([i[1],l+x,(min(euclid_dict[i[0],i[1]],euclid_dict[new_club_cluster[1],i[1]])),ln.get(i[1])+ln.get(l+o)])
                euclid_dict[i[1],l+x]=min(euclid_dict[i[0],i[1]],euclid_dict[new_club_cluster[1],i[1]])                
                todelete.append(i)            
                todelete.append([new_club_cluster[1],i[1],euclid_dict[new_club_cluster[1],i[1]],ln.get(i[1])+ln.get(new_club_cluster[1])])            
                            
            elif i[0]==new_club_cluster[0] and legend.index(new_club_cluster[0])<legend.index(i[1])<legend.index(new_club_cluster[1]):                     
                legend.append(l+o)                
                el.append([i[1],l+x,(min(euclid_dict[i[0],i[1]],euclid_dict[i[1],new_club_cluster[1]])),ln.get(i[1])+ln.get(l+o)])
                euclid_dict[i[1],l+x]=(min(euclid_dict[i[0],i[1]],euclid_dict[i[1],new_club_cluster[1]]))
                todelete.append(i)
                todelete.append([i[1],new_club_cluster[1],euclid_dict[i[1],new_club_cluster[1]],ln.get(i[1])+ln.get(new_club_cluster[1])])
                                
        for e in el:
            euclid_list.append(e)
        euclid_list = [n for n in euclid_list if n not in todelete]
        
        Z.append(z)
        o=o+1
        x=x+1
        k=k+1
    return(Z,[],[])
    
def complete_linkage(init_cluster,init,legend,le,ln):
    l=len(init_cluster)
    euclid_list=[]
    euclid_dict={}
    for ini in init_cluster:
        euclidean_dist=0
        ini_index=init_cluster.index(ini)
        for j in range(init_cluster.index(ini)+1,l):
            inj=init_cluster[j]
            inj_index=init_cluster.index(inj)
            euclidean_dist=math.sqrt((float(ini[1])-float(inj[1]))**2+\
            (float(ini[2])-float(inj[2]))**2+\
            (float(ini[3])-float(inj[3]))**2+\
            (float(ini[4])-float(inj[4]))**2+\
            (float(ini[5])-float(inj[5]))**2)
            euclid_list.append([legend[ini_index],legend[inj_index],euclidean_dist,ln.get(inj[0])+ln.get(ini[0])])
            euclid_dict[legend[ini_index],legend[inj_index]]=euclidean_dist
    Z=[]
    o=0
    x=0
    while o<l-1:
        todelete=[]
        new_club_cluster=[]
        el=[]
        new_club_cluster=(sorted(euclid_list,key=itemgetter(2),reverse=True)[0])
        z=[legend.index(new_club_cluster[0]),legend.index(new_club_cluster[1]),new_club_cluster[2],ln.get(new_club_cluster[0])+ln.get(new_club_cluster[1])]
        if z[0]>l+1 :
            z[0]=new_club_cluster[0]
        if z[1]>l+1:
            z[1]=new_club_cluster[1]
        ln[l+o]=ln.get(new_club_cluster[0])+ln.get(new_club_cluster[1])        
        todelete.append(new_club_cluster) 
        for i in euclid_list:
            if i[1]==new_club_cluster[0] :
                legend.append(l+x)                
                el.append([i[0],l+x,(max(euclid_dict[i[0],i[1]],euclid_dict[i[0],new_club_cluster[1]])),ln.get(i[0])+ln.get(l+o)])
                euclid_dict[i[0],l+x]=max(euclid_dict[i[0],i[1]],euclid_dict[i[0],new_club_cluster[1]])
                todelete.append(i)
                todelete.append([i[0],new_club_cluster[1],euclid_dict[i[0],new_club_cluster[1]],ln.get(i[0])+ln.get(new_club_cluster[1])])

            elif i[0]==new_club_cluster[0] and legend.index(i[1])>legend.index(new_club_cluster[1]):
                legend.append(l+x)                
                el.append([i[1],l+x,(max(euclid_dict[i[0],i[1]],euclid_dict[new_club_cluster[1],i[1]])),ln.get(i[1])+ln.get(l+o)])
                euclid_dict[i[1],l+x]=max(euclid_dict[i[0],i[1]],euclid_dict[new_club_cluster[1],i[1]])                
                todelete.append(i)            
                todelete.append([new_club_cluster[1],i[1],euclid_dict[new_club_cluster[1],i[1]],ln.get(i[1])+ln.get(new_club_cluster[1])])            
                            
            elif i[0]==new_club_cluster[0] and legend.index(new_club_cluster[0])<legend.index(i[1])<legend.index(new_club_cluster[1]):                     
                legend.append(l+o)                
                el.append([i[1],l+x,(max(euclid_dict[i[0],i[1]],euclid_dict[i[1],new_club_cluster[1]])),ln.get(i[1])+ln.get(l+o)])
                euclid_dict[i[1],l+x]=(max(euclid_dict[i[0],i[1]],euclid_dict[i[1],new_club_cluster[1]]))
                todelete.append(i)
                todelete.append([i[1],new_club_cluster[1],euclid_dict[i[1],new_club_cluster[1]],ln.get(i[1])+ln.get(new_club_cluster[1])])
                                
        for e in el:
            euclid_list.append(e)
        euclid_list = [n for n in euclid_list if n not in todelete]
        
        Z.append(z)
        
        o=o+1
        x=x+1
    return(Z,[],[])
    
def centroid(init_cluster,init,legend,leg):
    i=1
    k=2
    K=[]
    SSE=[]
    li=len(init)
    #print(li)
    Z=[]
    while i<li:
        l=len(init_cluster)
        euclid_list=[]
        for ini in init_cluster:
            euclidean_dist=0
            ini_index=init_cluster.index(ini)
            for j in range(init_cluster.index(ini)+1,l):
                inj=init_cluster[j]
                inj_index=init_cluster.index(inj)
                euclidean_dist=math.sqrt((float(ini[1])-float(inj[1]))**2+\
                                    (float(ini[2])-float(inj[2]))**2+\
                                    (float(ini[3])-float(inj[3]))**2+\
                                    (float(ini[4])-float(inj[4]))**2+\
                                    (float(ini[5])-float(inj[5]))**2)
                euclid_list.append([legend[ini_index],legend[inj_index],euclidean_dist])
        new_club_cluster=[]
        z=[]
        new_club_cluster=(sorted(euclid_list,key=itemgetter(2),reverse=True)[0])
        cnt=init_cluster[legend.index(new_club_cluster[0])][6]+init_cluster[legend.index(new_club_cluster[1])][6]
        z=[leg.index(new_club_cluster[0]),leg.index(new_club_cluster[1]),new_club_cluster[2],cnt]        
        l1=legend.index(new_club_cluster[0])
        l2=legend.index(new_club_cluster[1])
        indices=l1,l2        
        cluster=[li+i,(float(init_cluster[l1][1])+float(init_cluster[l2][1]))/2,\
                        (float(init_cluster[l1][2])+float(init_cluster[l2][2]))/2,\
                        (float(init_cluster[l1][3])+float(init_cluster[l2][3]))/2,\
                        (float(init_cluster[l1][4])+float(init_cluster[l2][4]))/2,\
                        (float(init_cluster[l1][5])+float(init_cluster[l2][5]))/2,
                        cnt]
        init_cluster.append(cluster)
        sse=0        
        for s in indices:
            sse=sse+(math.sqrt((float(cluster[1])-float(init_cluster[s][1]))**2+\
                (float(cluster[2])-float(init_cluster[s][2]))**2+\
                (float(cluster[3])-float(init_cluster[s][3]))**2+\
                (float(cluster[4])-float(init_cluster[s][4]))**2+\
                (float(cluster[5])-float(init_cluster[s][5]))**2))**2
                    
        SSE.append(sse)
        K.append(k)
        legend.append(li+i)
        leg.append(li+i)
        Z.append(z)               
        indices = l1,l2
        init_cluster = [n for j, n in enumerate(init_cluster) if j not in indices]
        legend = [m for j, m in enumerate(legend) if j not in indices]
        i=i+1
        k=k+1
    return(Z,SSE,K)
def draw_dendogram(Z,labels,K,displaylabels):
    if K==None:
        truncate_mode=None
        p=0
    else:
        truncate_mode='lastp'
        p=K
        
    if displaylabels==False:
        no_labels=True
    else:
        no_labels=False
        
    orientation='left'
# Display a dendrogram without leaf labels
    dendrogram(Z, labels=labels, no_labels=no_labels, orientation=orientation,truncate_mode=truncate_mode,p=p)
    show()

def sse(SSE,K,mention_x,mention_y):
    K_new=[]
    SSE_new=[]
    for i in range(K.index(mention_x),K.index(mention_y)):
        K_new.append(K[i])
        SSE_new.append(SSE[i])
        
    plot(K_new, SSE_new)
    show()


if __name__ == '__main__':
    main()
