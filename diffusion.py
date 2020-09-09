
plot_graphs=False

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct
import scipy
import scipy.integrate

from geometry import *



m=n=50
grid_shape=(m,n)




def mdct2(v,shape=None):
    '''
    ####CHANGE THIS SO THAT DIMENSION OF TRANSFORM LIMITED####

    v is an array of values defined on an (m_1,m_2,...m_k) grid, the Cosine Fourier transform of the
    function defined on [0,1]^len(v.shape), v[i_1,i_2,...,i_k] is value of
    function at (2*(i_1,i_2,...,i_k)+1)/v.shape (pointwise division).

    see page 106 of Strauss
    '''
    #temp=dct2(v)/np.prod(v.shape)
    #return temp

    return nmdct2(v,shape=shape)



def nmdct2(v,shape=None):
    '''
    new way of doing the transform (quicker)

    v is an array of values defined on an (m_1,m_2,...m_k) grid, the Cosine Fourier transform of the
    function defined on [0,1]^len(v.shape), v[i_1,i_2,...,i_k] is value of
    function at (2*(i_1,i_2,...,i_k)+1)/v.shape (pointwise division).

    see page 106 of Strauss
    '''
    #temp=dct2(v)/np.prod(v.shape)
    #return temp


    if shape is None:
        return (scipy.fftpack.dctn(v)/np.prod(v.shape))
        return dct2(v)/np.prod(v.shape)#dct returns totally new array
    else:
        return (scipy.fftpack.dctn(v)[:shape[0],:shape[1]]/(np.prod(shape)**1.5)*(np.prod(v.shape)**0.5))
        return (scipy.fftpack.dctn(v,shape=shape)/(np.prod(shape)**1.5)*(np.prod(v.shape)**0.5))
        return dct2(v)[:shape[0],:shape[1]]/np.prod(v.shape)#dct returns totally new array


def omdct2(v,shape=None):
    '''
    old way of doing the transform

    v is an array of values defined on an (m_1,m_2,...m_k) grid, the Cosine Fourier transform of the
    function defined on [0,1]^len(v.shape), v[i_1,i_2,...,i_k] is value of
    function at (2*(i_1,i_2,...,i_k)+1)/v.shape (pointwise division).

    see page 106 of Strauss
    '''

    if shape is None:
        return dct2(v)/np.prod(v.shape)#dct returns totally new array
    else:
        return dct2(v)[:shape[0],:shape[1]]/np.prod(v.shape)#dct returns totally new array

def fcos(A,x,box=[[0,1],[0,1]]):
    '''
    A - Cosine series of function f
    x - point at which we want to

    Estimate f at x using Cosine series given by A

    see page 106 of Strauss
    '''
    x=np.copy(x)

    l=np.array([box[i][1]-box[i][0] for i in range(len(x))])
    ret=np.copy(A)
    for i in range(x.shape[0]):

        temp=np.cos(np.arange(A.shape[i])*(np.pi*(x[i]-box[i][0])/(l[i])))
        temp[0]=temp[0]/2
        ret=np.tensordot(ret,temp,axes=(0,0))

    assert np.prod(ret.shape)==1
    return np.sum(ret)

def fcosd(A,x,axis,box=[[0,1],[0,1]]):
    '''
    A - Cosine series of function f
    x - point at which we want to
    axis - what axis of x is the derivative with respect to

    Estimate f_axis at x using Cosine series given by A

    see page 106 of Strauss
    '''

    l=np.array([box[i][1]-box[i][0] for i in range(len(x))])
    ret=np.copy(A)
    for i in range(x.shape[0]):
        temp=0
        if i==axis:
            temp=np.arange(A.shape[i])*(np.pi/l[i])*np.sin(np.arange(A.shape[i])*(np.pi*(x[i]-box[i][0])/l[i]))
            temp[0]=0#redundant
        else:
            temp=np.cos(np.arange(A.shape[i])*(np.pi*(x[i]-box[i][0])/l[i]))
            temp[0]=temp[0]/2
        ret=np.tensordot(ret,temp,axes=(0,0))

    assert np.prod(ret.shape)==1
    return np.sum(ret)

def heat_fcos_t(A,x,t,box=[[0,1],[0,1]]):
    '''
    Currently only for 2D
    '''

    l=np.array([box[i][1]-box[i][0] for i in range(len(x))])
    B=np.zeros(A.shape)

    for i in range(l.shape[0]):
        temp=(np.arange(A.shape[i])*np.pi/l[i])*(np.arange(A.shape[i])*np.pi/l[i])


        B=np.swapaxes(np.swapaxes(B,i,-1)+temp,i,-1)
    return fcos(np.exp(-B*t)*A,x,box=box)

def heat_fcosd_t(A,x,axis,t,box=[[0,1],[0,1]]):
    '''
    Currently only for 2D
    '''

    l=np.array([box[i][1]-box[i][0] for i in range(len(x))])
    #B=np.outer(np.arange(A.shape[0])*np.pi/l[0],np.arange(A.shape[1])*np.pi/l[1])
    B=np.zeros(A.shape)

    for i in range(l.shape[0]):
        temp=(np.arange(A.shape[i])*np.pi/l[i])*(np.arange(A.shape[i])*np.pi/l[i])

        #mul=lambda y:y+temp
        ##assert temp[0]==0
        #B=np.apply_along_axis(mul,i,B)


        B=np.swapaxes(np.swapaxes(B,i,-1)+temp,i,-1)
    return -fcosd(np.exp(-B*t)*A,x,axis,box=box)




class DiffusionSystem:
    '''
    Implemented in 2-D with Neumann Boundary conditions (0 flow in and out)
    Want to impllement various degrees of precision later (not hard).
    All diffusion system pretend coordinate system of the density grid is unit square.
    '''

    def __init__(self, initial_density_grid,shape=None,atol=1e-6,rtol=1e-6,A=None,box=[[0,1],[0,1]]):

        self.atol=atol
        self.rtol=rtol
        self.idg=initial_density_grid
        self.box=box
        if A is None:
            self.A=mdct2(self.idg,shape=shape) #cosine series of idg
        else:
            self.A=A

    def get_density(self,x,t):
        return heat_fcos_t(self.A,x,t,box=self.box)

    def get_density_gradient(self,x,t):
        return np.array([heat_fcosd_t(self.A,x,axis,t,box=self.box) for axis in range(x.shape[0])])

    def move_point_rk(self,point,T,t0=0,rtol=None,atol=None):
        '''
        point is a Point object or array

        Returns:
        Coordinate where point will end
        '''
        if rtol is None:
            rtol=self.rtol

        if atol is None:
            atol=self.atol

        c=None
        if not ( type(point)==list or type(point)==np.ndarray):
            c=np.copy(point.c)
        else:
            c=np.copy(point)

        integration_obj=scipy.integrate.RK45(lambda t,y:
                -self.get_density_gradient(y,t)/
                self.get_density(y,t)
                            ,t0,c,T,rtol=rtol,atol=atol)
        steps=0
        while True:
            crtol=rtol
            catol=atol
            if steps>=1000:

                print("stopped early")
                print(T)
                print(integration_obj.t)
                print(f"Trying again with rtol={rtol*2} and atol={atol*2}")
                return self.move_point_rk(point,T,t0=0,rtol=rtol*2,atol=atol*2)
            steps+=1
            try:
                integration_obj.step()
            except RuntimeError:
                break

        if plot_graphs:
            print(f"{steps} steps taken.")


        return integration_obj.y


    def move_point(self,point,T,n=None,t0=0):
        '''
        point is a Point object or array

        Returns:
        Coordinate where point will end
        '''


        return  self.move_point_rk(point,T,t0=t0)


    def move_points(self,points,T,n=None,t0=0):
        move_point_l=lambda point:self.move_point(point,T,n=n,t0=t0)
        return np.array([move_point_l(point) for point in points])


    def move_point_obj(self,point,T,n=None,t0=0):
        '''
        point is a Point object

        '''
        point.c=self.move_point_rk(point,T,t0=t0)

    def move_points_obj(self,points,T,n=None,t0=0):
        move_point_l=lambda point:self.move_point(point,T,n=n,t0=t0)

        print(f"There are {len(points)} points to evaluate.")
        for c,point in enumerate(points):
            if c%100==0:
                print(f"Point number {c} ({point}) being moved.")
            self.move_point_obj(point,T,n=n,t0=t0)

    def move_shape(self,shape,T,n=None,t0=0):
        self.move_points_obj(shape.get_points(),T,n=None,t0=0)

    def get_density_grid(self,mesh,t):
        '''
        mesh should be ij - of shape (#dimension,#dimension 1 elements, ...)
        '''
        trans=lambda x: heat_fcos_t(self.A,np.array([x[0],x[1]]),t)

        return np.apply_along_axis(trans,0,mesh)

    def get_density_grid(self,dims,t=0):
        return self.get_density_grid(np.array(np.meshgrid(*[crange(dims[i],0,1)
                                         for i in range(len(dims))
                                            ],indexing='ij')),0.0)

    def make_diffusion_system_mass_dict(mass_dict,cart_map,dims,shape,atol=1e-6,rtol=1e-6,box=[[0,1],[0,1]]):
        '''
        mass_dict should map from cart_maps region objects to a number, it may also map the string 'default' to a default value.
        '''

        mdg=MapDensityGrid.make_map_density_grid(mass_dict,cart_map,dims)

        if plot_graphs:
            plt.imshow(mdg.density_grid[:,::-1].T)
            plt.show(block=False)
        return DiffusionSystem(mdg.density_grid
                               ,shape=shape,atol=atol,rtol=rtol,box=box)

    def make_diffusion_system_map_density_grid(mdg,dims,shape,atol=1e-6,rtol=1e-6,box=[[0,1],[0,1]]):
        return DiffusionSystem(mdg.density_grid,shape=shape,atol=atol,rtol=rtol,box=box)


    def plot(self,shape=(100,100)):
        plt.figure(figsize=(8,8))
        grrd=self.get_density_grid(shape)
        plt.imshow(grrd.T[::-1,:])
