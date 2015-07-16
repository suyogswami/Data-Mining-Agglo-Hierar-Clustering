from matplotlib.pyplot import show, plot
from scipy.cluster.hierarchy import dendrogram
Z=[[3, 5, 0.159463005114039, 2], [4, 6, 0.2752299947316794, 3], [1, 2, 0.4970557815778824, 2], [7, 8, 0.6034587392688916, 5], [0, 9, 1.052588613846834, 6]]

#Z = [[1, 4, 1.41421356, 2],  [0, 6, 2, 3],  [2, 7, 2.82842712, 4],  [3, 5, 3.16227766, 2],  [8, 9, 4, 6]]
label_array = ['a','b','c','d','e','f'] 

# Display a dendrogram  
dendrogram(Z, labels = label_array)
show()

#Display a dendrogram, from left to right
dendrogram(Z, labels = label_array, orientation = 'left')
show()

# Display a dendrogram without leaf labels
dendrogram(Z, no_labels = True, orientation = 'left')
show()

#Display top 3 levels of the dendrogram
dendrogram(Z, labels = label_array, orientation = 'left', truncate_mode = 'lastp', p =5 )
show()

#Draw a Line plot using Matplotlib
#You can use Line plot to show SSE vs. # of clusters
kvalues = [2,3,4,5,6,7,8,9,10]
SSEs = [100, 65, 45, 30, 20, 15, 12, 10, 9]
plot(kvalues, SSEs)
show()
