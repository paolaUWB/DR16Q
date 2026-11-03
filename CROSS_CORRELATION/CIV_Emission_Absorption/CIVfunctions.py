import numpy as np
import joblib

def project(data, fit, perp_dist=False):
    #Use: Project scattered data (data) onto a line of best fit (fit)
    #Returns: 2-D array of (x,y) locations on the line of your orthogonal projection (and optionally the orthogonal distances from said line)
    
    #data: 2-D (N by 2) array of your data in some x-y space
    #fit: 2-D (N by 2) array of coordinates along your line of best fit
    #perp_dist: Bool; if True, also return the data points' perpendicular distances from the best-fit curve    

    locs = [] #this will contain xy locations along the best-fit curve corresponding to each data point's orthogonal projection onto said curve - python list faster than numpy array
    perp_distances = [] #to save perp distances from curve in    

    #Loop through each data point and compare with locations along fit
    for scat in data: 
        r2 = (scat[0]-fit[:,0])**2 + (scat[1]-fit[:,1])**2      #dist^2 of scat from each point along fit
        delta = [fit[np.argmin(r2), 0], fit[np.argmin(r2), 1]]  #save point along fit where dist^2 was minimum
        locs.append(delta)
        
        #Which side of the curve is this data point on?
        if (scat[1] >= fit[np.argmin(r2), 1]):
            side = 1.
        else:
            side = -1.
        perp_distances.append(np.sqrt(r2.min()) * side)  
        
    if perp_dist:
        return np.array(locs), np.array(perp_distances)
    return np.array(locs)


def CIV_distance(data_original, fit_original, step=1, logEW=True, path="./", perp_dist=False):
    #fit: N-by-2 array containing coordinates of points along best fit line 
    #data: N-by-2 [[x,y]] array of data 
    #NOTE: This really just caters to this situation (assumes monotonically decreasing fit_original)
    
    #Since we'll be manipulating the input arrays, work with copies
    data = data_original.copy()
    fit = fit_original.copy()

    #1) Scale all the data equally
    if logEW:
        scaler = joblib.load(path+"scalers/scaler_logEW.save")
        data = scaler.transform(data)
        fit = scaler.transform(fit)
    else:
        fit[:,1] = 10.**fit[:,1] #The saved curve's EW is on a logarithmic scale
        scaler = joblib.load(path+"scalers/scaler_linEW.save")
        data = scaler.transform(data)
        fit = scaler.transform(fit)

    #2) data is now each point's orthogonal projection onto fit
    ind_projcut = -((fit.shape[0]-1)%step)
    if not perp_dist:
        if ind_projcut!=0: data = project(data, fit[:ind_projcut])
        else: data = project(data, fit)
    else:
        if ind_projcut!=0: data, perpendicular_dist = project(data, fit[:ind_projcut], perp_dist=True)
        else: data, perpendicular_dist = project(data, fit, perp_dist=True)

    darr = [] #list to fill with distances along best-fit line for each point

    #3) Loop through each data point- start at tip of line and sum dist traveled until passing data point
    for scat in data:
        d = 0 #start at beginning of the line
        for i in range(0, fit.shape[0]-step, step):
            xp, x = fit[i,0], fit[i+step,0]
            yp, y = fit[i,1], fit[i+step,1] 
            dp = d
            d += np.sqrt((x-xp)**2 + (y-yp)**2)
            if yp >= scat[1] >= y: #if we pass the projected y-coord, save the distance traveled
                darr.append((d+dp)/2) #save the average between dprevious and d
                break
    
    if perp_dist:
        return np.array(darr), perpendicular_dist
    return np.array(darr)
