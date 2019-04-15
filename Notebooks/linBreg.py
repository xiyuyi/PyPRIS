import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt


# Implementation of shrink(x, mu)

def shrink(x, mu):
    x[np.where(x > mu)] -= mu
    x[np.where(x < -mu)] += mu
    return x

# Implementation of Stopping Criteria
def ifcont(res, crit):
    if len(res) >= 2:
        drop_perc = abs((res[len(res)-1] - res[len(res)-2])/res[len(res)-2])
        #print (res[len(res)-1], (res[len(res)-2], drop_perc)
        if drop_perc < crit:
            return False
        else: return True
    else: return True

# Implementation of Bregman Iteration

def linbreg(im, sensmx, threshold=800, stepsize=1e-5, d=0.000012, minstep = 1000, maxstep = 1e6, crit = 1e-8, dim = 2, plot = True):
    f = im
    
    A = sensmx
    
    A_t = np.transpose(A)
    
    n = A.shape[1]
    
    u = np.zeros(n)
    v = np.zeros(n)

    delta = abs(2 / np.linalg.norm(np.dot(A, A_t)) - d)
    mu = threshold  # shrink threshold
    step_size = mu * stepsize
    count = 0
    res = np.array([])

    while (count < minstep or ifcont(res, crit)):
        count += 1
        if count % 1e4 == 0: print(count)
        res = np.append(res, np.linalg.norm(np.dot(A, u) - f))
        v += np.dot(A_t, f - np.dot(A, u)) * step_size
        u = 1 * shrink(v, mu)
        u[np.where(u < 0)] = 0 #positivity constraint
        if count >= maxstep: break
            
    if (plot):
        
        if (dim == 2): 
            u = np.reshape(u, (int(np.sqrt(len(u))),int(np.sqrt(len(u)))))
            plt.subplot(121)
            plt.imshow(u)

            plt.subplot(122)
            plt.scatter(x=np.linspace(0, len(res), len(res)), y=np.log10(res), color='green')
            
            
    else:
        plt.scatter(x=np.linspace(0, len(res), len(res)), y=np.log10(res), color='green')
                                    
                                    

    return u
