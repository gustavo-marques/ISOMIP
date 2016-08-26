import numpy as np
from pyvtk import *
import os

def VTKgen(lat,lon,mask,depth=None,h=None,temp=None,salt=None,rho=None,u=None,v=None,w=None,seaice=None,shelf_base=None,shelf_thick=None,writebottom=False,fname='test',dirname='VTK',date=None, t=0):
    """ A function that transform the output of couple models (ocean, sea-ice for now) into VTK files """

    NY,NX=lat.shape
    os.system('mkdir ' + dirname)
    NY,NX=lat.shape

    if depth is not None:
       newlat=np.resize(lat,(2,NY,NX))
       newlon=np.resize(lon,(2,NY,NX))
       newdepth,bottom=get_depth(h,depth,mask)
       pp = f3(newlon,newlat,newdepth)
       structure=StructuredGrid([2,NY,NX],pp)
       path_to_file = str('%s/%s-bathymetry.vtk' % (dirname,fname))
       if os.path.isfile(path_to_file):
           print ' \n' + '==> ' + 'Bathymetry has already been written, moving on ...\n' + ''
       else:
          # create bottom/shape and depths
          newdepth=f1(newdepth)
          bottom=f1(bottom)
          pointdata = PointData(Scalars(newdepth,name='Depth'), Scalars(bottom,name='Bottom9999'))
          # saving the data
          vtk = VtkData(structure,pointdata)
          vtk.tofile(dirname+'/'+fname+'-bathymetry','binary')

       if writebottom == True:
          print ' \n' + '==> ' + 'Writting tracers/vel. just at the bottom layer ...\n' + ''
          data=[]
          if temp is not None:
             tmp=np.zeros((2,NY,NX)) 
             if len(temp.shape)==2: # in case the user provides 2D array with bottom data
                 tmp[:,:,:]=temp[:,:]
             else:
                 tmp[:,:,:]=temp[-1,:,:]
                 
             temp=f1(tmp)
             data.append("Scalars(temp,name='Temp')")

          if salt is not None:
             tmp=np.zeros((2,NY,NX))
             if len(salt.shape)==2:
                 tmp[:,:,:]=salt[:,:] 
             else:
                tmp[:,:,:]=salt[-1,:,:]

             salt=f1(tmp)
             data.append("Scalars(salt,name='Salt')")

          if rho is not None:
             tmp=np.zeros((2,NY,NX))
             if len(rho.shape)==2:
                 tmp[:,:,:]=rho[:,:]
             else:
                tmp[:,:,:]=rho[-1,:,:]

             rho=f1(tmp)
             data.append("Scalars(rho,name='Neutral_density')")
 
          if u is not None and v is not None:
                w=np.zeros((2,NY,NX)) # no w vel for now
                tmpu=np.zeros((2,NY,NX))
                tmpv=np.zeros((2,NY,NX))
                if len(u.shape)==2:
                    tmpu[:,:,:]=u[:,:]
                    tmpv[:,:,:]=v[:,:]
                else:
                    tmpu[:,:,:]=u[-1,:,:]
                    tmpv[:,:,:]=v[-1,:,:]

                vel=f3(tmpu,tmpv,w)
                data.append("Vectors(vel,name='Velocity')")

          if temp is not None or salt is not None or rho is not None or u is not None:

              for d in range(len(data)):
                 if d==0:
                    tmp=data[d]
                 else:
                    tmp=tmp+','+data[d]

              s = str("PointData(%s)" % (tmp))
              pointdata = eval(s)
              # saving the data
              vtk = VtkData(structure,pointdata)
              if date is not None:
                 s = str("vtk.tofile('%s/%s-%s-bottom-%05d','binary')" % (dirname,fname,date,t))
                 eval(s)

              else:
                 s = str("vtk.tofile('%s/%s-bottom-%05d','binary')" % (dirname,fname,t))
                 eval(s)
   
            
 
    if shelf_base is not None and shelf_thick is not None:
       newlat=np.resize(lat,(2,NY,NX))
       newlon=np.resize(lon,(2,NY,NX))
       dum,z = get_iceshelf(shelf_base,shelf_thick)
       iceshelf = f1(dum)
       pp = f3(newlon,newlat,z)
       structure=StructuredGrid([2,NY,NX],pp)
       pointdata = PointData(Scalars(iceshelf,name='IceShelf9999'))
       vtk = VtkData(structure,pointdata)
       vtk.tofile(dirname+'/'+fname+'-ice-shelf','binary')

    if writebottom==False:
       data=[]
       if temp is not None:
         temp=f1(temp)
         data.append("Scalars(temp,name='Temp')")
         
       if salt is not None:
         salt=f1(salt)
         data.append("Scalars(salt,name='Salt')")

       if rho is not None:
         rho=f1(rho)
         data.append("Scalars(rho,name='Neutral_density')")

       if u is not None and v is not None:
         if w is not None:
            vel=f3(u,v,w)
         else:
            w=np.zeros(u.shape)
            vel=f3(u,v,w)

         data.append("Vectors(vel,name='Velocity')")
    
       if seaice is not None:
         NZ,NY,NX=h.shape
         sice1=np.zeros((NZ,NY,NX))
         sice1[0,:,:]=seaice[:,:]
         sice2=np.zeros((NZ,NY,NX))
         seaice[seaice>=0.15]=1.0 # all values >=15% are unit
         sice2[0,:,:]=seaice[:,:]
         seaice1=f1(sice1) 
         seaice2=f1(sice2) 
         data.append("Scalars(seaice1,name='Sea-ice')")
         data.append("Scalars(seaice2,name='Sea-ice-binary')")
       

       if temp is not None or salt is not None or rho is not None or u is not None or seaice is not None:
         NZ,NY,NX=h.shape
         # resize lon lat for real mesh
         newlat=np.resize(lat,(NZ,NY,NX))
         newlon=np.resize(lon,(NZ,NY,NX))
         pp = f3(newlon,newlat,h)
         structure=StructuredGrid([NZ,NY,NX],pp)

         for d in range(len(data)):
            if d==0:
               tmp=data[d]
            else:
               tmp=tmp+','+data[d]

         s = str("PointData(%s)" % (tmp))
         pointdata = eval(s)
         # saving the data
         vtk = VtkData(structure,pointdata)
         if date is not None:
             s = str("vtk.tofile('%s/%s-%s-%05d','binary')" % (dirname,fname,date,t))
             eval(s)

         else:
             s = str("vtk.tofile('%s/%s-%05d','binary')" % (dirname,fname,t))
             eval(s)


def get_iceshelf(bottom,thick):
    NY,NX=bottom.shape
    tmp=np.resize(thick,(2,NY,NX))
    ice_shelf=np.zeros((2,NY,NX))
    z_shelf=np.zeros((2,NY,NX))
    z_shelf[0,:,:]=bottom[:,:]
    z_shelf[1,:,:]=bottom[:,:] + thick[:,:]
    z_shelf = np.ma.masked_where(tmp<100,z_shelf)
    z_shelf.fill_value=np.nan
    ice_shelf[:]=9999
    ice_shelf = np.ma.masked_where(tmp<100,ice_shelf)
    ice_shelf.fill_value=np.nan
    return ice_shelf,z_shelf

def get_depth(h,D,mask):
    NZ,NY,NX=h.shape
    dep=np.zeros((2,NY,NX))
    bot=np.zeros((2,NY,NX))
    bot[0,:,:]=9999
    for i in range(NX):
      for j in range(NY):
         if mask[j,i]==False:
                 dep[:,j,i]=-D[j,i] # ocean/ice-shelf
         else:
                 dep[:,j,i]=0 # land

    return dep,bot

def f1(X):
    """function to sort one 3D array"""

    a,b,c=X.shape
    return [(X[k,i,j]) for j in range(c) for i in range(b) for k in range(a)]

def f3(X,Y,Z):
    """ function to sort three 3D matrices"""
    a,b,c=X.shape
    return [(X[k,i,j],Y[k,i,j],Z[k,i,j]) for j in range(c) for i in range(b) for k in range(a)]


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]
