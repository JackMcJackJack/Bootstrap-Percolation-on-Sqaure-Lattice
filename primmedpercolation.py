import numpy as np
import matplotlib.pyplot as plt

class BootstrapPerc2D():
    #BLOCKED, OPEN SPACE/site that can be percolated  = 0, 1
    
    def __init__(self,x,y,p,m):
        self.x_size = x
        self.y_size = y
        
        self.p = p #probability a site can be percolated
        self.m = m
        
        self.safe_sites = np.zeros((self.x_size,self.y_size))
    
        self.n=0 #used to show the number of required sweeps before stable
        
        self.site_labels = np.ones((self.x_size,self.y_size))
        #this needs to be removed
        self.periodic_site_labels =np.ones((self.x_size,self.y_size))
        self.geometry_number = 4


    
    def get_initial(self):
        self.field = np.random.choice([0,1],
                                      size = (self.x_size,self.y_size),
                                      p =[1-self.p,self.p])
    
    def cull_once(self,field, m):
        neighbour_values = np.zeros((self.x_size,self.y_size))
        new_field = np.copy(field)
        
        for i in range(self.x_size):
            for j in range(self.y_size):
                if self.safe_sites[i][j] == 1:
                #if we know this site is safely unoccupied, that will never change & we dont do anything
                    continue
                else:
                    pass
                #x = field[i][j]
                occupied_neighbours = 0
                
                #brute force check every neighbour
                #this doesnt actually represent the directions, but doesnt matter a whole lot
                directions = ['right','left','up','down'] 
                for direction, position in zip(directions,[(i+1,j),(i-1,j),(i,j-1),(i,j+1)]):
                    try:
                        directions[directions.index(direction)] = field[position]
                    except:
                        pass
             
                for direction in directions:
                    if direction == 1:
                        occupied_neighbours += 1
                    #else: If it's a 0 I don't care, and if its not a valid direction
                    #it will remain a string and I don't care.
                    #if its 1 i do care
                neighbour_values[i][j] = occupied_neighbours
        for i in range(self.x_size):
            for j in range(self.y_size):
                if neighbour_values[i][j] < m:#if it doesnt meet the threshold
                    new_field[i][j] = 0
                    try:
                        #once we know its empty, we dont need to check it again
                        self.safe_sites[i][j] = 1
                        #everything around it is a danger now though
                        self.safe_sites[i-1][j] = 0
                        self.safe_sites[i+1][j] = 0
                        self.safe_sites[i][j-1] = 0
                        self.safe_sites[i][j+1] = 0
                    except:
                        pass
        return new_field

    def bootstrap_perc(self,starting_field):
        field = starting_field
        new_field = self.cull_once(field,self.m)
                
        while np.array_equal(field,new_field) != True:
            self.n+=1
            print('sweep =',self.n)
            #print('new_field=\n',new_field)
            field = new_field
            self.field = field
            self.bootstrap_perc(field)
    #modified version of the hoshen-kopelman algorithm
    def clusters(self,F,periodic = False):
        largest_ID = 1
        x_size , y_size = np.shape(F)
        
        if not periodic:
            directions = self.directions_physical
            labels = self.site_labels
        else:
            directions = self.directions_periodic
            labels = self.periodic_site_labels
            
        for i in range(x_size):
            for j in range(y_size):
                if F[i][j]:
                    up , left = directions(F,i,j)
                    #if both are empty, let this be a new cluster (it could be unified later)
                    if up == 0 and left == 0:
                        labels[i][j] = largest_ID
                        largest_ID += 1
                    if up == 0 and left == 1:
                        #if theyre connected, they have the same ID
                        labels[i][j] = labels[i][j-1]
                    if up == 1 and left == 0:
                        labels[i][j] = labels[i-1][j]
                    #if both are filled, connect these clusters and give them the smallest ID of the two
                    if up == 1 and left == 1:
                        self.unify(i,j,labels)
                else:
                    #if there isnt anything there, tell the labels that (helps with periodic)
                    labels[i][j] = 0
                        
                        
    def unify(self,i,j,labels):
        smaller_label = min(labels[i-1][j],labels[i][j-1])
        bigger_label = max(labels[i-1][j],labels[i][j-1])

        #must remember to also connect these clusters by the i,j site
        labels[i][j] = smaller_label            
                    
        swallowed_cluster_indices = np.argwhere(labels == bigger_label)
        for i in swallowed_cluster_indices:
            labels[tuple(i)] = smaller_label
                        

    def directions_physical(self,F,i,j):
        if i > 0:
            if j > 0:
                up , left = F[i-1][j], F[i][j-1]
            elif j == 0:        
                up , left = F[i-1][j], 0
        elif i == 0:
            if j > 0:
                up , left = 0, F[i][j-1]
            elif j == 0:
                up , left = 0, 0
        return up, left
        
    def isolate_clusters(self,labels):
        labels = labels.astype(int)
        IDs = range(1,np.amax(labels)+1)
        all_clusters_isolated = []
        for cluster_ID in IDs:
            isolated_cluster = np.zeros((self.x_size,self.y_size))
            isolated_cluster[np.where(labels == cluster_ID)] = cluster_ID
            all_clusters_isolated += [isolated_cluster.astype(int)]
        return zip(all_clusters_isolated, IDs)

    def percolates(self,labels,data = False):
        percolation = []
        isolated_clusters = self.isolate_clusters(labels)
        percolating_cluster, percolating_cluster_size = [],0
        for cluster,i in isolated_clusters:
            touch_left = np.hsplit(cluster,self.y_size)[0]
            touch_right  = np.hsplit(cluster,self.y_size)[-1]
            touch_up  = np.vsplit(cluster,self.x_size)[0]
            touch_down  = np.vsplit(cluster,self.x_size)[-1]
            
            #print(i,'u',touch_down)
            if sum(touch_left) != 0 and sum(touch_right) != 0:
                percolation += ['Cluster %d percolates horizontally'%i]
                percolating_cluster += [i] #the cluster that percolates
                percolating_cluster_size += sum(sum(cluster))* i**(-1)#how many points make up this cluster
            if sum(sum(touch_up)) != 0 and sum(sum(touch_down)) != 0:
                percolation += ['Cluster %d percolates vertically'%i]
                percolating_cluster += [i] #the cluster that percolates
                percolating_cluster_size += sum(sum(cluster))* i**(-1)#how many points make up this cluster

        if percolation == []:
            percolation = 'No Percolation'
            
        if data == False:
            print(percolation)
        else:
            return(percolating_cluster,percolating_cluster_size)
            
            

def main():
    perc = BootstrapPerc2D(100,100,0.7,2)
    
    perc.get_initial()
    print('initial field = \n',perc.field)
    
    perc.bootstrap_perc(perc.field)
    print('Required every site to have at least %d percolatable adjacent sites in order to survive (except on the borders).' %(5-perc.m))
    print('Finished with %d iterations' %perc.n)
    print('final field = \n',perc.field)
    perc.clusters(perc.field)
    print('labelled field = \n',perc.site_labels.astype(int))
    print(perc.percolates(perc.site_labels.astype(int)))
    
#main()

def data(p_min,p_max,step,grid_size = 100):
    p_axis = [np.linspace(p_min,p_max,step),np.linspace(p_min,p_max,step),np.linspace(p_min,p_max,step),np.linspace(p_min,p_max,step)]
    cluster_probability = [np.linspace(0,0,step),np.linspace(0,0,step),np.linspace(0,0,step),np.linspace(0,0,step)]
    for m in [1,2,3,4]:
       for p , i in zip(p_axis[m-1],range(len(p_axis[m-1]))):
            print('i = ',i, ' m = ', m)
            if p > 0.75:
                grid_size = 50
            perc = BootstrapPerc2D(grid_size,grid_size,p,m)
            perc.get_initial()
            perc.bootstrap_perc(perc.field)
            perc.clusters(perc.field)
            print('labelled field = \n',perc.site_labels.astype(int))
            
            percolating_cluster_ID, percolating_cluster_size = perc.percolates(perc.site_labels,data = True)
            print(percolating_cluster_size)
            prob = percolating_cluster_size / (perc.x_size*perc.y_size)
            cluster_probability[m-1][i] = prob
    fig , axes = plt.subplots(1,1)
    for m in [1,2,3,4]:
        axes.plot(p_axis[m-1],cluster_probability[m-1])
        axes.set_title('Plot of probability of site being in spanning cluster against occupation probability with bootstrap thresholds of m = 1,2')
    fig.set_figheight(10)
    fig.set_figwidth(15)
    fig.legend(['m = 1','m=2','m=3','m=4'])
    plt.show()
    return cluster_probability

a = data(0.5,1,10)