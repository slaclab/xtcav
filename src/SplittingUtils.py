import numpy as np
from Utils import *

def SplitImage(image, n, islandsplitmethod,par1,par2):
    """
    Split an XTCAV image depending of different bunches, this function is still to be programmed properly
    Arguments:
      image: 3d numpy array with the image where the first index always has one dimension (it will become the bunch index), the second index correspond to y, and the third index corresponds to x
      n: number of bunches expected to find
    Output:
      outimage: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """

    if n==1:    #For one bunch, just the same image
        outimage=np.zeros((n,image.shape[0],image.shape[1]))
        outimage[0,:,:]=image    
    elif n==2:  #For two bunches,  the image on the top and on the bottom of the center of mass
        
        
        #outimage=HorizontalLineSplitting(image[0,:,:]) 
        #outimage=OptimalSplittingLine(image[0,:,:])       

        if islandsplitmethod == 'contourLabel':   
            outimage = IslandSplittingContour(image,par1,par2)
        elif islandsplitmethod == 'autothreshold':              
            outimage = IslandSplittingAutoThreshold(image,2)
        else:
            outimage = IslandSplitting(image,2)

    else:       #In any other case just copies of the image, for debugging purposes
        outimage=IslandSplitting(image,n)
    
     
    
    return outimage
    
    
def HorizontalLineSplitting(image):
    """
    Divides the image with a horizontal line at the center of mass
    Arguments:
      image: 2d numpy array with the image where the first index correspond to y, and the second index corresponds to x
    Output:
      outimage: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """
    
    outimage=np.zeros((2,image.shape[0],image.shape[1]))
    x=range(image.shape[1])
    y=range(image.shape[0])
    x0,y0=GetCenterOfMass(image[:,:],x,y)
    outimage[0,0:round(y0),:]=image[0:round(y0),:]
    outimage[1,round(y0):,:]=image[round(y0):,:]  
    
    return outimage
    
def RotatingLineSplitting(image):
    """
    Divides the image with a straight line crossing the center of mass at the optimal angle
    Arguments:
      image: 2d numpy array with the image where the first index correspond to y, and the second index corresponds to x
    Output:
      outimage: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """
    
    outimage=np.zeros((2,image.shape[0],image.shape[1]))
    x=range(image.shape[1])
    y=range(image.shape[0])
    x0,y0=GetCenterOfMass(image[:,:],x,y)
    
    angles=np.linspace(-np.pi/2, np.pi/2, num=50, endpoint=False)
    
    areas=np.zeros(len(angles))
    cutinds = np.zeros((len(angles),len(x)))
        
    for i in range(len(angles)): 
        dx=np.cos(angles[i])
        dy=np.sin(angles[i])
                        
        if (angles[i]>=0):
            lengthleft=np.amin([y0/np.sin(angles[i]),x0/np.cos(angles[i])])
            lengthright=np.amin([(len(y)-y0-1)/np.sin(angles[i]),(len(x)-x0-1)/np.cos(angles[i])])
        else:
            lengthleft=np.amin([(len(y)-y0-1)/np.sin(-angles[i]),x0/np.cos(angles[i])])
            lengthright=np.amin([y0/np.sin(-angles[i]),(len(x)-x0-1)/np.cos(angles[i])])
            
            
        for j in range(-int(lengthleft),int(lengthright)):        
            areas[i]=areas[i]+image[int(np.round(j*dy+y0)),int(np.round(j*dx+x0))]

        
    

    
    ind=np.where(areas == np.amin(areas))

    #In case of more than one minimum value, better to take the one in the middle
    optimalangle=angles[int((ind[0][-1]+ind[0][0])/2)]            
            
    mask=np.zeros((len(y),len(x)))
    splitline=(x-x0)*np.tan(optimalangle)+y0 
    
    for i in range(len(x)): 
        mask[:,i]=y<splitline[i]
    
    outimage[0,:,:]=image*mask
    outimage[1,:,:]=image*(1-mask)
    
    return outimage
    
    
def OptimalSplittingLine(image):
    """
    Find the optimal line to separate the two bunches and separates the image in two. Assuming they have different energy for each x value obtains the optimal y value between two peaks.
    Arguments:
      image: 2d numpy array with the image where the first index correspond to y, and the second index corresponds to x
    Output:
      outimage: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """
    #Use the center of mass as a starting point
    Nx=image.shape[1]
    Ny=image.shape[0]
    
    x=range(Nx)
    y=range(Ny)
    x0,y0=GetCenterOfMass(image,x,y)
    
    splitline=np.zeros(Nx, dtype=np.int32)
   
    
    x0=int(round(x0))
    y0=int(round(y0))
    
    #Filter the image along y axis
    fimage=im.filters.uniform_filter1d(image, 20, axis=0)
    
    #Each vertical slice of the right half of the image (left to right) from the x center of masses and then 
    #each vertical slice of the left half of the image (right to left) from the x center of masses  
    startPoint=y0
    for i in (range(x0,Nx)+range(x0-1,-1,-1)):
        #Find the point towards higher indices when the value starts to increase
        top=startPoint;    
        while (top<(Ny-1) and fimage[top,i]>=fimage[top+1,i]):
            top=top+1;    
        #Find the point towards lower indices when the value starts to increase
        bottom=startPoint;    
        while (bottom>0 and fimage[bottom,i]>=fimage[bottom-1,i]):
            bottom=bottom-1;
    
        #Use the middle value
        splitline[i]=int(round((top+bottom)/2))
        #The initial point for the next iteration will be the calculated point from he previous slice
        if (i!=(Nx-1)):
            startPoint=splitline[i];
        else: # except when we jump to the right half, in which case we we the value assigned for the initial point
            startPoint=splitline[x0];
        
        
    #For each x axis we split the slice
    outimage=np.zeros((2,image.shape[0],image.shape[1]))
    for i in range(image.shape[1]):
        outimage[0,0:splitline[i],i]=image[0:splitline[i],i]
        outimage[1,splitline[i]:,i]=image[splitline[i]:,i] 
        
    return outimage

    
def IslandSplitting(image,N):
    """
    Find islands in the picture and order them by area, returning the image in N images, ordered by area. Teh total area is one.
    Arguments:
      image: 2d numpy array with the image where the first index correspond to y, and thesecond index corresponds to x
      N: number of islands to return
    Output:
      outimages: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """
    #Use the center of mass as a starting point
    Nx=image.shape[1]
    Ny=image.shape[0]
    x=range(Nx)
    y=range(Ny)
    
    #Get a bool image with just zeros and ones
    imgbool=image>0
    
    #Calculate the groups
    groups, n_groups =im.measurements.label(imgbool);
    
    #Structure for the areas and the images
    areas=np.zeros(n_groups,dtype=np.float64)
    
    
    #Obtain the areas
    for i in range(0,n_groups):    
        areas[i]=np.sum(image*(groups==(i+1)))
            
    #Get the indices in descending area order
    orderareaind=np.argsort(areas)  
    orderareaind=np.flipud(orderareaind)

    #Check that the area of the second bunch is not smaller than a fraction of the first bunch (otherwise the split probably did not succeed)
    n_area_valid=1
    for i in range(1,n_groups): 
        if areas[orderareaind[i]]<1.0/20*areas[orderareaind[0]]:
            break
        else:
            n_area_valid+=1

    #Number of valid images for the output

    n_valid=np.amin([N,n_groups,n_area_valid])    
    
    #Obtain the separated images
    images=[]  
    for i in range(0,n_valid):    
        images.append(image*(groups==(orderareaind[i]+1)))
    
    #Calculate the angle of each large area island with respect to the center of mass
    x0,y0=GetCenterOfMass(image,x,y)        
    angles=np.zeros(n_valid,dtype=np.float64)
    xi=np.zeros(n_valid,dtype=np.float64)
    yi=np.zeros(n_valid,dtype=np.float64)
    for i in range(n_valid):
        xi[i],yi[i]=GetCenterOfMass(images[i],x,y)
        angles[i]=math.degrees(np.arctan2(yi[i]-y0,xi[i]-x0))
                
    #And we order the output based on angular distribution
    
    #If the distance of one of the islands to -180/180 angle is smaller than a certain fraction, we add an angle to make sure that nothing is close to the zero angle

    
    #Ordering in angles (counterclockwise from 3 oclock)
    #dist=180-abs(angles)
    #if np.amin(dist)<30.0/n_valid:
    #    angles=angles+180.0/n_valid
    #orderangleind=np.argsort(angles)  

    #Ordering in height (higher energy first)
    orderangleind=np.argsort(-yi)  

    #Structure for the output
    outimages=np.zeros((n_valid,Ny,Nx))        

    #Assign the proper images to the output
    for i in range(n_valid):
        outimages[i,:,:]=images[orderangleind[i]]
        
    #Renormalize to total area of 1
    outimages=outimages/np.sum(outimages)
    
    return outimages

def IslandSplittingAutoThreshold(image,N):
    """
    Find islands in the picture and order them by area, returning the image in N images, ordered by area. Teh total area is one.
    Arguments:
      image: 2d numpy array with the image where the first index correspond to y, and thesecond index corresponds to x
      N: number of islands to return
    Output:
      outimages: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """
    
    Nx=image.shape[1]
    Ny=image.shape[0]
    x=range(Nx)
    y=range(Ny)
    
    #Minimum and maximum edge for the threshold to be set
    thres0=image[image!=0].min()
    thres1=image[image!=0].max()
    
    n_valid=0
    iternum=0
    maxiter=17
    prevGood=False
    while (iternum<maxiter):
        #On each iteration we set the threshold to the middle value between the edges 
        if iternum==maxiter-1:
            thres=thres1
        #Except on the last iteration that we set it to the highest one
        else:
            thres=(thres0+thres1)/2
        #print thres
        imageaux=np.array(image)
        imageaux[imageaux<thres]=0
            
        #Get a bool image with just zeros and ones
        imgbool=imageaux>0
        
        #Calculate the groups
        groups, n_groups =im.measurements.label(imgbool);
   
        #Obtain the areas
        areas=np.zeros(n_groups,dtype=np.float64)
        for i in range(0,n_groups):    
            areas[i]=np.sum(image*(groups==(i+1)))
            
        #Get the indices in descending area order
        orderareaind=np.argsort(areas)  
        orderareaind=np.flipud(orderareaind)

        #Check that the area of the second bunch is not smaller than a fraction of the first bunch (otherwise the split probably did not succeed)
        n_area_valid=1
        for i in range(1,n_groups): 
            if areas[orderareaind[i]]<1.0/20*areas[orderareaind[0]]:
                break
            else:
                n_area_valid+=1

        #Number of valid images for the output
        n_valid=np.amin([N,n_groups,n_area_valid])    
        #print [thres0,thres1],thres,n_valid,n_groups,areas[orderareaind]        
        #If there are too few islands we decrease the upper limit, because we thresholded too much
        if n_valid<N:
            if prevGood==False: #If we have never seen a value that works, we keep going down
                thres1=thres
            else: #If we have seen a value that works we go up, because we went down too much from that value
                thres0=thres
                
        #If there are the right number of islands, we decrease the upper limit to see if we could get the same with a smaller threshold
        elif n_valid==N:
            prevGood=True
            thres1=thres
        #In any other case, we have to threshold more
        else:
            thres0=thres

        iternum+=1
    
    #Obtain the separated images
    images=[]  
    for i in range(0,n_valid):    
        images.append(image*(groups==(orderareaind[i]+1)))
    
    #Calculate the angle of each large area island with respect to the center of mass
    x0,y0=GetCenterOfMass(image,x,y)        
    angles=np.zeros(n_valid,dtype=np.float64)
    xi=np.zeros(n_valid,dtype=np.float64)
    yi=np.zeros(n_valid,dtype=np.float64)
    for i in range(n_valid):
        xi[i],yi[i]=GetCenterOfMass(images[i],x,y)
        angles[i]=math.degrees(np.arctan2(yi[i]-y0,xi[i]-x0))
                
    #And we order the output based on angular distribution
    
    #If the distance of one of the islands to -180/180 angle is smaller than a certain fraction, we add an angle to make sure that nothing is close to the zero angle

    
    #Ordering in angles (counterclockwise from 3 oclock)
    #dist=180-abs(angles)
    #if np.amin(dist)<30.0/n_valid:
    #    angles=angles+180.0/n_valid
    #orderangleind=np.argsort(angles)  

    #Ordering in height (higher energy first)
    orderangleind=np.argsort(-yi)  

    #Structure for the output
    outimages=np.zeros((n_valid,Ny,Nx))        

    #Assign the proper images to the output
    for i in range(n_valid):
        outimages[i,:,:]=images[orderangleind[i]]
        
    #Renormalize to total area of 1
    outimages=outimages/np.sum(outimages)
    
    return outimages
    
def IslandSplittingContour(image,ratio1,ratio2):
    """
    Find islands using the contour method in the picture and order them by area, returning two islands, ordered by area.
      1. Iteratively threshold the image (increase the the threshold level) up until the point that we have 2 big bunches
      2. Then use the labeling function in opencv to find how many contiguous groups exist in the data. 
      3. Then using the contours function in opencv, we find the pixels that correspond to the contours of each of the 2 big bunches.
      4. For each contiguous object D that is not one of the 2 big bunches, we take one representative pixel from D and calculate which contour it is closest to. We then give D the same label as the contour it is closest to.
      5. After this, we have still not labelled the non-zero pixels that were excluded due to thresholding(picture attached)
      6. We then repeat steps 2 and 4 on these points.
    Arguments:
      image: 2d numpy array with the image where the first index correspond to y, and thesecond index corresponds to x
    Output:
      outimages: 3d numpy array with the split image image where the first index is the bunch index, the second index correspond to y, and the third index corresponds to x
    """
    
    data = image
    k = 0  # initial threshold
    indicator = 0 # indicating that we do not yet have a satisfactory threshhold value
    while indicator == 0:  
        h = data>k # threshhold image
        labelled_array , num_features = im.measurements.label(h) #label blobs in images
        count =np.zeros(num_features)

        for g in range(1,num_features):
            temp = sum(sum(labelled_array==g)) #count number of pixels associated with each blob 
            count[g] = temp

        idx = (-count).argsort()[:3] #calculates 3 largest blobs.blob1,blob2, blob3
        var1 = sum(sum(labelled_array==idx[0]))
        var2 = sum(sum(labelled_array==idx[1]))
        var3 = sum(sum(labelled_array==idx[2]))

        if var1/var2 < ratio1 and var2/var3 > ratio2:#if blob1 and blob2 are approximately same size and blob2 is significantly bigger than blob3, then proceed. blob1>blob2>blob3 
            indicator =1
        k = k + .000002
   
    j1 = (labelled_array == idx[0])
    j2 = (labelled_array == idx[1])
    j1 = j1.astype(np.uint8)#convert from false/true to integer values
    j2 = j2.astype(np.uint8)
    contours1, hierarchy1 = cv2.findContours(j1,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)#find contours of 2 biggest blobs
    contours2, hierarchy2 = cv2.findContours(j2,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)

    j1 = (j1>1)
    j2 = (j2>1)

    t1 = np.nonzero(j1)
    t2 = np.nonzero(j2)
    labelled_array = np.matrix(labelled_array)
    copy = labelled_array

    t1 = np.matrix(t1)
    t2 = np.matrix(t2)

    t1 = t1.T
    t2 = t2.T

    # this is code to associate the smaller blobs with the two biggest blobs, based on distance to the contour
    # we take one representative pixel at random from each of the small blobs for this

    innerTrial1 =np.zeros(shape = (t1.shape[0],1,))#below is matrix manipulations to calculate distance

    #  the key idea is that dist(x,y)^2 = x'x - 2x'y + y'y where ' means transpose 
    #  (x is the location of the random pixel in the unassigned blob and y
    #  is the location of any pixel in the big-blob contour)
    #  x'x (the coordinate of the random pixel in the small blob) is constant and is ignored in this computation
    #  so we only need to compute 2x'y and y'y
    for k in range(0,t1.shape[0]):
        innerTrial1[k] = t1[k]*t1[k].T # y'y for blob1
         
    innerTrial2 = np.zeros(shape = (t2.shape[0],1,))
    for k in range(0,t2.shape[0]):
        innerTrial2[k] = t2[k]*t2[k].T # y'y for blob2

    for p in range(1,num_features):
        if p!=idx[0]:
            if p!=idx[1]:
                t = np.nonzero((copy == p))
                t = [t[0][0,0], t[1][0,0]]
                t= np.matrix(t)
                t = t.T
                
                temp1 = innerTrial1 - 2 * t1 * t # y'y - 2x'y for blob1
                temp2 = innerTrial2 - 2 * t2 * t # y'y - 2x'y for blob2
                temp1 =np.amin(temp1) # compute the minimum distance between large blob1 and small blob
                temp2 =np.amin(temp2) # compute the minimum distance between large blob2 and small blob
                if temp1 < temp2:#here we find which contour each unlabeled blob is closer too and associate the blob with that specific bunch
                    labelled_array[labelled_array==p]=idx[0]
                else:
                    labelled_array[labelled_array==p]=idx[1]

    if k>0:
        g = np.matrix(data)# Now we repeat the same process for pixels that were initially zeroed out because of the thresholding  
        temp = (labelled_array<1)
        temp = temp.astype(int)
        dent = np.multiply(g,temp)
        labelled_array2 , num_features2 = im.measurements.label(dent)
        copy = labelled_array2


        for p in range( 1, num_features2):
            t = np.nonzero((copy ==p))
            t = [t[0][0], t[1][0]]
            t = np.matrix(t)
            t = t.T
           
            temp1 = innerTrial1 - 2 * t1 * t
            temp2 = innerTrial2 - 2 * t2 * t
            temp1 = np.amin(temp1)
            temp2 = np.amin(temp2)
            if temp1 < temp2:
                labelled_array[labelled_array2==p]=idx[0]
            else:
                labelled_array[labelled_array2==p]=idx[1]

   
    # calculate the final output arrays
    Nx=image.shape[1]#here we split the orginal image into 2 bunches 
    Ny=image.shape[0]
    outimages=np.zeros((2,Ny,Nx))
    temp = (labelled_array == idx[0])
    temp = temp.astype(int)
    outimages[0,:,:]=np.multiply(g,temp)
    
    temp = (labelled_array == idx[1])
    temp = temp.astype(int)
    outimages[1,:,:]=np.multiply(g,temp)
    
    
    return outimages