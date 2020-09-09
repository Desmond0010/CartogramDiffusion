import numpy as np
import json
import matplotlib.pyplot as plt

do_verbose=False#alters verbosity of print statements.

class Point():
    counter=0
    meta_point_dict={}
    def __init__(self,x,y):
        self.c=np.array([x,y])
        self.id=Point.counter
        Point.meta_point_dict[self.id]=self
        Point.counter+=1



    def a_coords(points):
        '''
        Nx2 array of points, anticlockwise order
        '''
        return np.array([p.c for p in points])

    def distinct(p1,p2,ep=0):
        '''
        Test whether two points are at different coordinates

        ep- if nonzero, two close points will be considered the same
        '''
        return dist_sq(p1,p2)>ep

    def distinct_simple(p1,p2,ep=0):
        '''
        Test whether two points are at different coordinates
        '''
        return np.min(p1.c!=p2.c)

    def dist_sq(p1,p2):
        c1,c2=None,None
        if type(p1)==np.ndarray or type(p1)==list:
            c1=p1
            c2=p2
        else:
            c1=p1.c
            c2=p2.c

        return ((c1[0]-c2[0])*(c1[0]-c2[0])+(c1[1]-c2[1])*(c1[1]-c2[1]))

    def dist(p1,p2):
        return (Point.dist_sq(p1,p2))**(1/2)

    def get_c(self):
        return self.c

    def set_c(self,c):
        self.c=c

class PolygonPoint(Point):
    '''
    Extension of Point class to allow Point to have metadata giving its relation to polygons.
    '''

    def __init__(self,x,y):
        '''
        Point should be created and then it should be imputed to polygons using init_polygon_point()

        polygons_data = {'polygonobjecthash':{'prev':prev_point,'next':next_point,
            'prev_length':Point,
            'next_length':Point,
            'skip_length':length of line connecting two neighbouring points,
            'area_change':signed area change of polygon if this point is removed}

        prev and next should traverse polygon anticlockwise

        '''


        Point.__init__(self,x,y)
        #self.polygons_data=polygons_data
        self.polygons_data={}


    def insert_polygon_point(coordinate,prev_point,next_point):
        '''
        Insert a polygon point at the midpoint of two points and include in every polygon they share.
        '''

        polygons=[(polygon,prev_point.polygons_data[polygon]['next']==next_point)
                  for polygon in prev_point.polygons_data if
                   (prev_point.polygons_data[polygon]['prev']==next_point or
                   prev_point.polygons_data[polygon]['next']==next_point)]

        new_point=PolygonPoint(coordinate[0],coordinate[1])

        temp=[[multi_polygon for multi_polygon in polygon_entry[0].multi_polygons if multi_polygon is not None]
             for polygon_entry in polygons   ]
        mps=[]
        for multi_polygons in temp:
            mps+=multi_polygons

        mps=list(set(mps))
        for mp in mps:
            mp.N+=1


        for polygon_entry in polygons:
            polygon=polygon_entry[0]
            after=polygon_entry[1]
            polygon.points.insert(polygon.points.index(prev_point)+1*after,new_point)

            if after:#new_point after prev_point

                new_point.redo_polygons_data(polygon,prev_point=prev_point,next_point=next_point)
                prev_point.redo_polygons_data(polygon,next_point=new_point)
                next_point.redo_polygons_data(polygon,prev_point=new_point)

                #should also recalculate areas but I think that is depracated

            else:#new_point before prev_point
                new_point.redo_polygons_data(polygon,prev_point=next_point,next_point=prev_point)
                prev_point.redo_polygons_data(polygon,prev_point=new_point)
                next_point.redo_polygons_data(polygon,next_point=new_point)

            polygon.N+=1#point will never appear twice in same polygon


    def init_polygon_point(self, polygon,prev_point,next_point):
        '''
        Initialiases a PolygonPoint with respect to a polygon
        '''

        self.polygons_data[polygon]={
            'prev':prev_point,'next':next_point,
            'prev_length':Point.dist(self,prev_point),
            'next_length':Point.dist(self,next_point),
            'skip_length':Point.dist(prev_point,next_point),
            'area_change':-Polygon([
            prev_point,
                       self,
                       next_point],simple=True).calculate_area(),
        }


    def redo_polygons_data(self, polygon,prev_point=None,next_point=None):
        '''
        Initialiases a PolygonPoint with respect to a polygon
        '''
        if prev_point is None:
            prev_point=self.polygons_data[polygon]['prev']
        if next_point is None:
            next_point=self.polygons_data[polygon]['next']

        self.polygons_data[polygon]={
            'prev':prev_point,'next':next_point,
            'prev_length':Point.dist(self,prev_point),
            'next_length':Point.dist(self,next_point),
            'skip_length':Point.dist(prev_point,next_point),
            'area_change':-Polygon([
            prev_point,
                       self,
                       next_point],simple=True).calculate_area(),
        }



    def absorb(self, point):
        '''
        integrates polygon_data from other point and the merges points, objects which is self remains alive

        returns True if point successfully absorbed, else returns False
        '''


        if len(list(set(list(self.polygons_data.keys()))&set(list(point.polygons_data.keys()))))>0:
            '''
            If the two points are included in the same polygon do not merge
            '''
            return False

        multi_polygons_counted=[[],[]]

        for polygon in self.polygons_data:
            if type(polygon.multi_polygons)==list:
                for multi_polygon in polygon.multi_polygons:
                    if multi_polygon not in multi_polygons_counted[0]:
                        multi_polygons_counted[0].append(multi_polygon)

        for polygon in point.polygons_data:
            if type(polygon.multi_polygons)==list:
                for multi_polygon in polygon.multi_polygons:
                    if multi_polygon not in multi_polygons_counted[1]:
                        multi_polygons_counted[1].append(multi_polygon)

        for multi_polygon in (set(multi_polygons_counted[0]) & set(multi_polygons_counted[1])):
            multi_polygon.N-=1

        for polygon in point.polygons_data:
            #replace this with method __absorb_in
            assert point.getPrev(polygon).getNext(polygon)==point
            assert point.getNext(polygon).getPrev(polygon)==point
            point.getPrev(polygon).setNext(polygon,self)
            point.getNext(polygon).setPrev(polygon,self)

            assert polygon not in self.polygons_data
            self.polygons_data[polygon]=point.polygons_data[polygon]
        point.polygons_data={}

        return True


    def remove(self):
        '''
        does not delete object or remove from polygon point list as that would be too resource intensive, that must be done elsewhere
        '''
        assert len(list(self.polygons_data.keys()))>0
        multi_polygons_counted=[]#multi polygon will be here if we already reduced its point tally.

        for polygon in list(self.polygons_data.keys()):
            self.__remove_from(polygon)
            polygon.N-=1

            if type(polygon.multi_polygons)==list:#it may have no super structure
                for multi_polygon in polygon.multi_polygons:
                    if multi_polygon not in multi_polygons_counted:
                        multi_polygon.N-=1
                        multi_polygons_counted.append(multi_polygon)

        self.polygons_data={}
        #self.polygons_data={"message":"set to zero",'old_stuff':self.polygons_data}


    def __remove_from(self,polygon):
        '''
        Should only be called from self.remove()
        '''
        assert not ((self.polygons_data[polygon]['prev']==self) ^ (self.polygons_data[polygon]['next']==self))



        polygon.area   += self.polygons_data[polygon]['area_change']
        polygon.length = (polygon.length-self.polygons_data[polygon]['prev_length']
                          -self.polygons_data[polygon]['next_length']+self.polygons_data[polygon]['skip_length'])


        self.getPrev(polygon).init_polygon_point(
            polygon,self.getPrev(polygon).getPrev(polygon),self.getNext(polygon))


        self.getNext(polygon).init_polygon_point(
            polygon,self.getPrev(polygon),self.getNext(polygon).getNext(polygon))


    def removal_cost(self):
        '''
        The cost to remove this point, higher means more costly to remove.
        '''
        cost=0
        for polygon in self.polygons_data:
            mps=polygon.multi_polygons

            multi_polygon_modifier=1
            if len(mps)==2:
                multi_polygon_modifier=2#want to keep points on the borders more
            elif len(mps)>=3:#could just return a super high number
                return 1e100
                multi_polygon_modifier=10000

            for mp in mps:
                cost=max(cost,self.removal_cost_formula(mp,polygon)*multi_polygon_modifier)
        return cost

    def removal_cost_formula(self,mp,polygon):
        '''
        mp-multipolygon object which polygon is in.
        '''

        LENGTH_AREA=0.05
        return ((abs(self.polygons_data[polygon]['area_change']) +
                self.polygons_data[polygon]['prev_length']*self.polygons_data[polygon]['next_length']*LENGTH_AREA)
                *(
                self.polygons_data[polygon]['prev_length']+
                 self.polygons_data[polygon]['next_length']+
                  self.polygons_data[polygon]['skip_length']))

        area_w,length_w=1/2,1
        return (
            area_w*abs(self.polygons_data[polygon]['area_change'])/mp.area+
            length_w*4*(
                (
                    Point.dist(self,self.getPrev(polygon))*Point.dist(self,self.getNext(polygon))
                )
                /
                (
                    Point.dist(self.getNext(polygon),self.getPrev(polygon))
                    #*Point.dist(self.getNext(polygon),self.getPrev(polygon))
                )#*Point.dist(self.getNext(polygon),self.getPrev(polygon))
                /mp.length
            )
               )



    def getPrev(self,polygon):
        return self.polygons_data[polygon]['prev']

    def getNext(self,polygon):
        return self.polygons_data[polygon]['next']

    def setPrev(self,polygon,point):
        self.polygons_data[polygon]['prev']=point

    def setNext(self,polygon,point):
        self.polygons_data[polygon]['next']=point




class Shape():
    def __init__(self):
        self.extremities=None
        self.N=None
        self.points=None


    def contains_point(self,point):
        pass

    def contains_points(self,points):
        pass

    def random_point(self):

        point=Point(
            (self.extremities[0][1]-self.extremities[0][0])*np.random.unif(1)+self.extremities[0][0],
            (self.extremities[1][1]-self.extremities[1][0])*np.random.unif(1)+self.extremities[1][0]
        )
        while not self.contains_point(point):
            point=Point(
            (self.extremities[0][1]-self.extremities[0][0])*np.random.unif(1)+self.extremities[0][0],
            (self.extremities[1][1]-self.extremities[1][0])*np.random.unif(1)+self.extremities[1][0]
            )

        return point

    def get_points(self):
        return self.points

    def create_point_dict(self):
        point_dict={}
        for point in self.get_points():
            point_dict[point]=np.copy(point.c)

        return point_dict

    def set_points(self, point_dict):
        for point in self.get_points():
            point.set_c(point_dict[point])

class Polygon(Shape):
    def __init__(self,points,multi_polygons=None,simple=False):
        '''
        points should traverse perimeter in anticlockwise order
        multipolygons is a list of all multipolygons which contain this polygon (generally list of one element)

        multi_polygons could be replaced by more generic super shape object later
        simple should only be used for temporary polygons used for calculating area and such
        '''
        if type(points)==np.ndarray:
            self.points=[PolygonPoint(*points[i]) for i in range(points.shape[0])]
        else:
            if len(points)>0 and type(points[0])==list:
                self.points=[]
                for point in points:
                    self.points.append(PolygonPoint(point[0],point[1]))
            else:#list of Point objects
                self.points=points

        self.N=len(points)#always true as it is not allowed for one polygon to feature the same point object twise
        #this is dealt with by having two point objects share a coordinate.

        self.multi_polygons=None
        self.extremities=None

        if not simple:
            self.init_polygon_points()


        self.area=self.calculate_area()
        self.length=self.calculate_length()

        self.multi_polygons=multi_polygons

        self.extremities=[]#box format


        self.recalculate_extremities()

    def init_polygon_points(self):

        for i in range(len(self.points)):
            self.points[i].init_polygon_point(
                self,#polygon
                self.points[(i-1)%len(self.points)],
                self.points[(i+1)%len(self.points)])

    def recalculate_extremities(self):
        '''
        Changes extremities attribute.
        '''
        for point in self.points:
            if len(self.extremities)==0:
                for i in range(point.c.shape[0]):
                    self.extremities.append([point.c[i],point.c[i]])

            for i in range(point.c.shape[0]):
                if self.extremities[i][0]>point.c[i]:
                    self.extremities[i][0]=point.c[i]
                elif self.extremities[i][1]<point.c[i]:
                    self.extremities[i][1]=point.c[i]

    def update_points(self,main=True):
        '''
        Points removed from the polygon are not automatically removed from the points list. This redoes the point list.
        '''
        if len(self.points)==0:
            self.points=[]
            return
        new_p=[]
        for point in self.points:
            if self in point.polygons_data:
                new_p.append(point)
                break

        for i in range(self.N-1):
            new_p.append(new_p[-1].polygons_data[self]['next'])

        self.points=new_p
        if main:
            self.area=self.calculate_area()
            self.length=self.calculate_length()

    def plot_obj(self,color=None):
        '''
        Returns a matplotlib object of the polygon
        '''
        a=Point.a_coords(self.points)
        if a.shape[0]==0:
            return plt.plot()
        return plt.Polygon(a,color=color)

    def plot(self,color=None):
        Polygon.plot_polygons([self],color=color)


    def calculate_area_matrix(self):
        '''
        returns A, Symmetic (2N)x(2N) matrix s.t.
        Area = |z.T @ A @z|
        Where z is a list of points coordinates, alternating from x and y coordinate.
        '''
        Point.a_coords(self.points).flatten()

        A=np.zeros((2*self.N,2*self.N))
        for i in range(0,self.N):
            A[i*2,((i+1)*2+1)%(2*self.N)]=1/4
            A[((i+1)*2+1)%(2*self.N),i*2]=1/4#make 1/2, need matrix to be symmetric

            A[i*2,((i-1)*2+1)%(2*self.N)]=-1/4
            A[((i-1)*2+1)%(2*self.N),i*2]=-1/4#make 1/2, need matrix to be symmetric

        return A


    def calculate_area(self,scale=1):
        '''
        returns A, Symmetic (2N)x(2N) matrix s.t.
        Area = |z.T @ A @z|
        Where z is a list of points coordinates, alternating from x and y coordinate.
        '''
        area=0
        for i in range(self.N):
            area+=self.points[(i)%self.N].c[0]*self.points[(i+1)%self.N].c[1]*scale
            area-=self.points[(i)%self.N].c[1]*self.points[(i+1)%self.N].c[0]*scale

        area=abs(area)/2
        self.area=area
        return area

    def calculate_length(self,scale=1):
        '''
        Caclulates length of polygon.
        '''
        length=0
        for i in range(self.N):
            temp=0
            temp+=((self.points[(i)%self.N].c[0]-self.points[(i+1)%self.N].c[0])*scale*
                  (self.points[(i)%self.N].c[0]-self.points[(i+1)%self.N].c[0])*scale)
            temp+=((self.points[(i)%self.N].c[1]-self.points[(i+1)%self.N].c[1])*scale*
                  (self.points[(i)%self.N].c[1]-self.points[(i+1)%self.N].c[1])*scale)
            length+=temp**(1/2)

        self.length=length
        return length


    def plot_polygons(polygons,color=None):
        '''
        polygons is an object of matplotlib.Polygon not my Polygon class
        '''
        ax=plt.gca()
        for p in polygons:
            if len(p.points)==0:
                continue

            ax.add_patch(p.plot_obj(color=color))



    def distinctPoints(polygons):
        '''
        Returns of distinct points (as defined by object, not coordinate) in list of polygons.
        '''
        l=[]
        for p in polygons:
            l+=p.points

        l=list(set(l))#get rid of duplicates
        return l

    def a_coords(self,as_list=False):
        '''
        Nx2 array of points, anticlockwise order
        '''
        if as_list:
            return [list(p.c) for p in self.points]
        else:
            return np.array([p.c for p in self.points])

    def contains_point(self, point):
        '''
        returns simple true/false
        '''

        return ray_tracing_numpy(np.array([point.c[0]]),np.array([point.c[1]]),
                         self.a_coords)[0]

    def contains_points(self, points)->np.ndarray:
        '''
        points can be a list of Point objects or a numpy array

        returns true/false numpy array
        '''
        points_arr=None
        if type(points)==np.ndarray:
            points_arr=points
        else:
            points_arr=Point.a_coords(points)

        return ray_tracing_numpy(points_arr[:,0],points_arr[:,1],
                         self.a_coords())

    def get_points(self):
        return self.points

    def get_area(self):
        return self.area

    def get_length(self):
        return self.length

    def set_area(self,area):
        self.area=area

    def set_length(self,length):
        self.length=length



class MultiPolygon(Shape):
    #Multi-Polygon - number of polygons, 1, number of points, 2
    def __init__(self,polygons,simple=False):
        self.polygons=[]
        for polygon in polygons:
            if isinstance(polygon,Polygon):
                self.polygons.append(polygon)
                if polygon.multi_polygons is None:
                    polygon.multi_polygons=[self]
                else:
                    if self not in polygon.multi_polygons:
                        polygon.multi_polygons.append(self)
            else:#if list of single polygon- included due to data source
                assert type(polygon[0])==list and len(polygon)==1
                self.polygons.append(Polygon(polygon[0],multi_polygons=[self],simple=simple))

        self.extremities=[]
        self.recalculate_extremities(update_subshapes=False)


        self.N=None
        self.recalculate_N(update_subshapes=False)


    def get_points(self):
        points=[]
        for polygon in self.polygons:
            points+=polygon.get_points()
        return list(set(points))


    def calculate_area(self):
        area=0
        for polygon in self.polygons:
            area+=polygon.calculate_area()
        return area


    def calculate_length(self):
        length=0
        for polygon in self.polygons:
            length+=polygon.calculate_length()
        return length

    def get_area(self):
        area=0
        for polygon in self.polygons:
            area+=polygon.area
        return area

    def get_length(self):
        length=0
        for polygon in self.polygons:
            #length+=polygon.calculate_length(scale=scale)
            length+=polygon.length
        return length



    def recalculate_extremities(self,update_subshapes=True):
        if update_subshapes:
            for polygon in self.polygons:
                polygon.recalculate_extremities()

        self.extremities=[]
        for polygon in self.polygons:
            if len(self.extremities)==0:
                for i in range(len(polygon.extremities)):
                    self.extremities.append(
                        [polygon.extremities[i][0],polygon.extremities[i][1]])

            for i in range(len(polygon.extremities)):
                if self.extremities[i][0]>polygon.extremities[i][0]:
                    self.extremities[i][0]=polygon.extremities[i][0]

                if self.extremities[i][1]<polygon.extremities[i][1]:
                    self.extremities[i][1]=polygon.extremities[i][1]


    def update_points(self):
        '''
        Points removed from the polygon are not automatically removed from the points list. This redoes the point list.
        '''
        for polygon in self.polygons:
            polygon.update_points(main=False)

        self.calculate_area()
        self.calculate_length()


    def recalculate_N(self,update_subshapes=True):
        if update_subshapes:
            for polygon in self.polygons:
                polygon.update_points()

        point_list=[]
        for polygon in self.polygons:
            point_list+=polygon.points
        self.N=len(set(point_list))#this is necessary in case duplicate points are added


    def contains_point(self, point):
        '''
        returns simple true/false
        '''
        for polygon in self.polygons:
            if polygon.contains_point(point):
                return True
        return False

    def clean_polygons(self):
        ret=[]
        for polygon in self.polygons:
            if len(polygon.points)!=0:
                ret.append(polygon)
        i=0
        while i < len(self.polygons):
            polygon=self.polygons[i]
            if (polygon not in ret):
                self.polygons.remove(polygon)
            else:
                i+=1

    def contains_points(self, points)->np.ndarray:
        '''
        returns true/false numpy array
        '''
        points_arr=None
        if type(points)==np.ndarray:
            points_arr=points
        else:
            points_arr=Point.a_coords(points)

        self.clean_polygons()#this may not be necessary

        ret=np.zeros(points.shape[0])
        for i in range(len(self.polygons)):
            ret[:]=ret[:]+self.polygons[i].contains_points(points_arr)

        return ret!=0

    def reduce_points(self,n):
        '''
        Gets rid of points so only have n ledft
        '''

        '''
        points=[]
        for polygon in self.polygons:
            points+=[point for point in polygon.points]

        points=list(set(points))
        '''
        self.remove_points(self.N-n)

    def remove_points(self,n):
        '''
        Removes n lowest cost point in way which keeps cost lowest
        '''
        n_removed=0

        while n_removed<n:
            costs=[]
            removal_points=[]
            m=min(max((n-n_removed)//20,10),n-n_removed)
            points=[]
            for polygon in self.polygons:
                points+=[point for point in polygon.points  if polygon in point.polygons_data ]

            points=list(set(points))

            for point in points:
                cost=point.removal_cost()
                i=min(m,len(costs))
                while i>=1 and cost<costs[i-1]:
                    i-=1

                if i<m:
                    costs.insert (i,cost )
                    removal_points.insert(i,point)

                    if len(costs)==m+1:
                        costs.pop()
                        removal_points.pop()

            highest_cost=costs[-1]
            for point in removal_points:#change to <=
                if point.removal_cost()<=highest_cost:#as removing a point may change dynamic
                    point.remove()
                    n_removed+=1

                    if n_removed>=n:
                        return

    def make_geojson_geometry(self):
        '''
        return dictionary of this MultiPolygon in GeoJSON format.
        '''
        coordinates=[]
        for polygon in self.polygons:
            coordinates.append([polygon.a_coords(as_list=True)])
            #coordinates dimensions - polygon - list with one entry  - points list - point coordinates list

        return {'type':'MultiPolygon','coordinates':coordinates}

    def plot(self,color=None):
        Polygon.plot_polygons(self.polygons,color=color)

class Region(MultiPolygon):
    '''
    Extends MultiPolygon to include metadata
    '''
    def __init__(self,polygons,mass,color=None,id=None):
        if isinstance(polygons,MultiPolygon):
            MultiPolygon.__init__(self,polygons.polygons)
        else:
            MultiPolygon.__init__(self,polygons)

        self.mass, self.density=None,None
        self.resolve_mass(mass)

        if color is None:
            color=random_color()
        self.color=color
        self.id=id

    def resolve_mass(self,mass):
        '''
        Uses mass to determine density.
        '''
        self.mass=mass
        self.density=None
        area=self.calculate_area()

        if self.mass is not None and area>0:
            self.density=self.mass/area

    def get_mass(self):
        return self.mass


    def plot(self,color=None):
        if color is None:
            color=self.color
        Polygon.plot_polygons(self.polygons,color=color)

class RegionCollection(Shape):

    def __init__(self,regions):
        self.regions=regions

        self.extremities=[]
        for region in self.regions:
            if len(self.extremities)==0:
                for i in range(len(region.extremities)):
                    self.extremities.append(
                        [region.extremities[i][0],region.extremities[i][1]])

            for i in range(len(region.extremities)):
                if self.extremities[i][0]>region.extremities[i][0]:
                    self.extremities[i][0]=region.extremities[i][0]

                if self.extremities[i][1]<region.extremities[i][1]:
                    self.extremities[i][1]=region.extremities[i][1]

        self.mass=sum([region.mass for region in self.regions])


    def get_area(self):
        return sum([region.get_area() for region in self.regions])

    def get_length(self):
        return sum([region.get_length() for region in self.regions])

    def get_mass(self):
        return sum([region.get_mass() for region in self.regions])

    def get_points(self):
        points=[]
        for region in self.regions:
            points+=region.get_points()
        return list(set(points))

    def set_points(self, point_dict):
        for point in self.get_points():
            point.set_c(point_dict[point])

    def merge_points(self):
        point_dict={}#maps value to point object
        for region in self.regions:
            for polygon in region.polygons:
                for point in polygon.points:#must be a point object
                    assert point.getPrev(polygon).getNext(polygon)==point
                    assert point.getNext(polygon).getPrev(polygon)==point


        for region in self.regions:
            for polygon in region.polygons:
                for point in polygon.points : #if not point.absorbed#must be a point object
                    if polygon not in point.polygons_data:
                        continue

                    p_data=point.c.tobytes()
                    if p_data in point_dict:


                        if point_dict[p_data].absorb(point):#successfully merged
                            #assert point.getPrev(polygon).getNext(polygon)==point
                            #assert point.getNext(polygon).getPrev(polygon)==point
                            #point=point_dict[p_data]
                            point.absorbed=True

                        #assert len(list(point.polygons_data.keys()))>0

                        #if point cannot be absorbed it is just left out
                    else:
                        point_dict[p_data]=point
                        #assert len(list(point.polygons_data.keys()))>0

        '''
        for region in self.regions:
            for polygon in region.polygons:
                for point in polygon.points:#must be a point object
                    assert point.getPrev(polygon).getNext(polygon)==point
                    assert point.getNext(polygon).getPrev(polygon)==point

        '''

    def clean_polygons(self):
        for region in self.regions:
            region.clean_polygons()

    def calculate_area(self):
        area=0
        for region in self.regions:
            area+=region.calculate_area()
        return area

    def calculate_length(self):
        length=0
        for region in self.regions:
            length+=region.calculate_length()
        return length

    def calculate_average_density(self):
        return self.get_mass()/self.get_area()


    def contains_point(self, point):
        '''
        returns simple true/false
        '''
        for region in self.regions:
            if region.contains_point(point):
                return True
        return False

    def contains_points(self, points)->np.ndarray:
        '''
        returns true/false numpy array
        '''
        points_arr=None
        if type(points)==np.ndarray:
            points_arr=points
        else:
            points_arr=Point.a_coords(points)

        ret=np.zeros((points.shape[0],len(self.regions)))
        for i in range(len(self.regions)):
            ret[:,i]=self.regions[i].contains_points(points_arr)

        return np.max(ret,axis=1)


    def reduce_points(self,n,evenly=True):
        '''
        Gets rid of points so only have n left
        '''

        '''
        points=[]
        for region in self.regions:
            for polygon in region.polygons:
                points+=[point for point in polygon.points]

        points=list(set(points))
        '''

        N=sum([region.N for region in self.regions])
        #N=len(points)

        if evenly:
            region_n=n//len(self.regions)
            i=0
            counter=0
            while counter<len(self.regions):
                print(i)
                print(self.regions[i].N)
                counter+=1
                i=(i+1)%len(self.regions)
                if self.regions[i].N>region_n:
                    counter=0
                    m=min(
                        max(int((self.regions[i].N-region_n)//1.5),5)
                          ,self.regions[i].N-region_n)

                    self.regions[i].remove_points(m)


        else:
            self.remove_points(N-n)

    def remove_points(self,n):
        '''
        Removes n lowest cost point in way which keeps cost lowest
        '''
        if n<=0:
            return

        n_removed=0



        while n_removed<n:
            removal_points=[]
            costs=[]
            m=min(max((n-n_removed)//100,10),n-n_removed)
            points=[]
            for region in self.regions:
                for polygon in region.polygons:
                    points+=[point for point in polygon.points if polygon in point.polygons_data ]

            points=list(set(points))

            for point in points:
                cost=point.removal_cost()
                i=min(m,len(costs))
                while i>=1 and cost<costs[i-1]:
                    i-=1

                if i<m:
                    costs.insert (i,cost)
                    removal_points.insert(i,point)

                    if len(costs)==m+1:
                        costs.pop()
                        removal_points.pop()

            highest_cost=costs[-1]

            for point in removal_points:
                if point.removal_cost()<=highest_cost:#as removing a point may change dynamic

                    point.remove()
                    n_removed+=1



    def resolve_mass(self,mass_dict):
        '''
        mass dict maps from region object to number
        Does not change mass of wrapping CartMap
        '''
        for region in self.regions:
            if region in mass_dict:
                region.resolve_mass(mass_dict[region])
            else:
                if 'default' in mass_dict:
                    region.resolve_mass(mass_dict['default'])
                else:
                    raise KeyError("One of the regions is not given a mass and there is no default mass.")

    def plot(self):
        for region in self.regions:
            region.plot()


    def insert_points(self,cutoff_percent=0.05):
        '''
        To counter the lack of points in certain long straight borders
        '''
        def polygon_super_area(polygon):
            return sum([multi_polygon.get_area() for multi_polygon in polygon.multi_polygons])
        def point_length_score(point,polygon=None):
            if polygon==None:
                return max([point.polygons_data[polygon]['next_length']/polygon_super_area(polygon)
                            for polygon in point.polygons_data])
            else:
                return point.polygons_data[polygon]['next_length']/polygon_super_area(polygon)

        lengths=[point_length_score(point,polygon=None)
                 for point in self.get_points()]

        lengths=list(np.sort(np.array(lengths)))

        cutoff=lengths[-int(len(lengths)*cutoff_percent)]

        no_change=False
        n_inserted=0
        while not no_change:
            print(n_inserted)
            no_change=True
            for point in self.get_points():
                max_length=point_length_score(point,polygon=None)
                if max_length>cutoff:
                    no_change=False
                    for polygon in point.polygons_data:
                         if point_length_score(point,polygon=polygon)==max_length:
                            PolygonPoint.insert_polygon_point(
                                (
                                point.c+
                                point.polygons_data[polygon]['next'].c)/2,
                                point,
                                point.polygons_data[polygon]['next'])
                            n_inserted+=1
                            continue
        return n_inserted

class CartMap():

    def __init__(self,region_collection,frame,mass_choice='average'):

        self.region_collection=region_collection
        self.frame=frame
        self.regions=self.region_collection.regions
        self.extremities=self.region_collection.extremities


        self.background_mass, self.background_density, self.area=None,None,None

        self.resolve_background_mass()



    def resolve_background_mass(self,mass_choice='average'):

        self.area=self.calculate_area()
        self.background_density=self.region_collection.calculate_average_density()
        self.background_mass=None

        if self.background_density is not None:
            self.background_mass=self.background_density*(self.area-self.region_collection.get_area())





    def calculateOccupiedArea(self):
        area=0
        for region in self.regions:
            area+=region.calculate_area()
        return area

    def calculate_area(self):
        self.region_collection.calculate_area()
        return np.prod(np.array(self.frame.original_proxy_box)[:,1]-np.array(self.frame.original_proxy_box)[:,0])

    def contains_point(self, point):
        '''
        contains point in region collection
        returns simple true/false
        '''
        return self.region_collection.contains_point(point)

    def contains_points(self, points)->np.ndarray:
        '''
        contains point in region collection
        returns true/false numpy array
        '''
        return self.region_collection.contains_points(points)

    def get_points(self):
        return self.region_collection.get_points()

    def plot_map(self):
        for region in self.region_collection.regions:
            for polygon in region.polygons:
                if len(polygon.points)==0:
                    continue

                plt.plot(*np.array(
                    [point.c for point in polygon.points]+[polygon.points[0].c]
                ).T)

    def plot(self):
        self.region_collection.plot()

    def create_point_dict(self):
        return self.region_collection.create_point_dict()

    def set_points(self, point_dict):
        self.region_collection.set_points(point_dict)

    def resolve_mass(self,mass_dict):
        self.region_collection.resolve_mass(mass_dict)
        self.resolve_background_mass()

    def transform_points(self,points=None):
        '''
        points must be an interable
        '''
        if points is None:
            points=self.get_points()

        for point in points:
            point.c=self.frame.transform(point.c)

    def inverse_transform_points(self,points=None):
        '''
        points must be an iterable
        '''
        if points is None:
            points=self.get_points()

        for point in points:
            point.c=self.frame.inverse_transform(point.c)

    def insert_geometry_into_geojson(self,geojson_data,make_copy=True):
        '''
        This assumes that the regions in our region_collection and the regions in the geojson are in the same order
        '''
        if make_copy:
            geojson_data=geojson_data.copy()

        for i,region in enumerate(self.region_collection.regions):

            geojson_data['features'][i]['geometry']=region.make_geojson_geometry()

        return geojson_data


    def init_cart_map_geojson(geojson_data,cart_map_args={"FRAME_MULTIPLIER":1,"COUNTY_DETAIL":50},mass_dict={'default':1},color_dict={}):
        '''
        geojson_data
        Initialise a cart_map from a geojson dictionary.
        '''
        my_regions=[]
        for i in range(len(geojson_data['features'])):

            my_regions.append(Region(geojson_data['features'][i]['geometry']['coordinates'],
                                     mass_dict.get(geojson_data['features'][i]['id'],mass_dict['default']),
                                     color=color_dict.get(geojson_data['features'][i]['id'],None),
                                     id=geojson_data['features'][i]['id']))

        rc=RegionCollection(my_regions)

        cm=CartMap(rc,Frame(rc.extremities,cart_map_args.get("FRAME_MULTIPLIER",1)))
        rc.merge_points()#if the same point occurs twice in the same polygon it must not be be merged into one object
        update_points(cm)

        rc.reduce_points(len(rc.regions)*cart_map_args.get("COUNTY_DETAIL",50))
        update_points(cm)
        return cm

    def insert_points(self,cutoff_percent=0.05):
        self.calculate_area()#to update areas
        return self.region_collection.insert_points(cutoff_percent=cutoff_percent)

    def clean_polygons(self):
        self.region_collection.clean_polygons()

class Frame():

    def __init__(self,original_box,overshoot_multiplier, target_box=[[0,1],[0,1]]):

        self.original_box            = original_box
        self.target_box              = target_box
        self.overshoot_multiplier    = overshoot_multiplier

        center=[(self.original_box[0][1]+self.original_box[0][0])/2,
                (self.original_box[1][1]+self.original_box[1][0])/2]
        dims=[(self.original_box[0][1]-self.original_box[0][0]),
              (self.original_box[1][1]-self.original_box[1][0])]

        self.original_proxy_box=[
            [center[0]-dims[0]/2*self.overshoot_multiplier,center[0]+dims[0]/2*self.overshoot_multiplier],
            [center[1]-dims[1]/2*self.overshoot_multiplier,center[1]+dims[1]/2*self.overshoot_multiplier]
                                ]


    def transform(self,xy):

        return np.array([((xy[0]-self.original_proxy_box[0][0])/
                         (self.original_proxy_box[0][1]-self.original_proxy_box[0][0]))
                         *(self.target_box[0][1]-self.target_box[0][0])+self.target_box[0][0],
                         (xy[1]-self.original_proxy_box[1][0])/
                         (self.original_proxy_box[1][1]-self.original_proxy_box[1][0])
                         *(self.target_box[1][1]-self.target_box[1][0])+self.target_box[1][0]])


    def inverse_transform(self,xy):

        return np.array([(xy[0]-self.target_box[0][0])/(self.target_box[0][1]-self.target_box[0][0])
                         *(self.original_proxy_box[0][1]-self.original_proxy_box[0][0])
                         +self.original_proxy_box[0][0],
                         (xy[1]-self.target_box[1][0])/(self.target_box[1][1]-self.target_box[1][0])
                         *(self.original_proxy_box[1][1]-self.original_proxy_box[1][0])
                         +self.original_proxy_box[1][0]])


class MapDensityGrid():

    def __init__(self,cart_map,dims=[100,100],occupier_grid=None):
        '''
        dims gives dimension of the cart_map
        '''


        self.cart_map=cart_map
        self.dims=dims

        self.mesh=np.array(np.meshgrid(*[crange(dims[i],
                                                cart_map.frame.original_proxy_box[i][0],
                                                cart_map.frame.original_proxy_box[i][1])
                                         for i in range(len(cart_map.frame.original_proxy_box))],indexing='ij'))
        if occupier_grid is None:
            self.occupier_grid=self.determine_occupiers_grid()
        else:
            self.occupier_grid=occupier_grid
        self.density_grid=self.grid_occupier_density()
        #first [::-1]second [::-1] necessary as row column vs xy


    '''
    #currently not used but there may be a way to use this to make a more efficient algorithm.
    def determine_occupier(self,x,hint=None):
        #determines density of map at a coordinate
        if hint is not None:
            if hint.contains(x):
                return hint
        else:
            for region in self.cart_map.region_collection.regions:
                if region.contains(x):
                    return region
        return None
    '''

    def determine_occupiers_grid(self):

        grid=ij_mesh_to_point_grid(self.mesh)
        points=grid.reshape((grid.shape[0]*grid.shape[1],2))

        points_arr=None
        if type(points)==np.ndarray:
            points_arr=points#don`t need to copy
        else:
            points_arr=Point.a_coords(points)

        ret=np.zeros((points.shape[0]))
        for i in range(len(self.cart_map.regions)):
            ret[:]=(ret[:]*(self.cart_map.regions[i].contains_points(points_arr)==0)
                    +self.cart_map.regions[i].contains_points(points_arr)*(i+1))


        return ((ret-1).reshape((grid.shape[0],grid.shape[1]))).astype(int)
        #-1 if not occupied, ==occupiers if occupied


    def grid_occupier_density(self):
        counter=-1
        occ_den={-1:self.cart_map.background_density}
        for region in self.cart_map.regions:
            counter+=1
            occ_den[counter]=region.density
        f=np.vectorize(lambda x:occ_den[x])#to apply pointwise

        return f(self.occupier_grid)

    def make_map_density_grid(mass_dict,cart_map,dims,occupier_grid=None):
        '''
        This will change the mass structure of the cart_map
        '''
        cart_map.plot()
        cart_map.resolve_mass(mass_dict)

        if do_verbose:
            print(occupier_grid)
            print([region.density for region in cart_map.regions])

        return MapDensityGrid(cart_map,dims=dims,occupier_grid=occupier_grid)

def load_json(path):

    with open(path,encoding='utf-8') as file:
        print(f"Loading: {file}")
        data=(json.load(file))
    return data

def save_json(path,data):

    with open(path,'w', encoding='utf-8') as file:
        print(f"Saving: {file}")
        data=json.dump(data, file, ensure_ascii=False, indent=4)
    return data

def geojson_standardise(geojson_data):
    geojson_data_copy=geojson_data.copy()
    for i in range(len(geojson_data_copy['features'])):
        if geojson_data_copy['features'][i]['geometry']['type']=='Polygon':
            geojson_data_copy['features'][i]['geometry']['type']='MultiPolygon'
            geojson_data_copy['features'][i]['geometry']['coordinates']=[geojson_data_copy['features'][i]['geometry']['coordinates']]

    return geojson_data_copy


def update_points(cm):
    for i in range(len(cm.regions)):
        cm.regions[i].update_points()

    for region in cm.regions:
        for polygon in region.polygons:
            polygon.set_area(polygon.calculate_area())


def random_color():
    return (np.random.random(size=3))

def ij_mesh_to_point_grid(mesh):

    my_mesh=np.array(mesh)#in case mesh is a list
    return np.swapaxes(np.swapaxes(my_mesh,0,2),0,1)

def ray_tracing_numpy(x :np.ndarray,y:np.ndarray,poly:np.ndarray)->np.ndarray:
    '''
    x,y coordinates of point
    poly - list of polygon points

    This method counts methods on
    Source:
    https://gis.stackexchange.com/questions/62925/shapely-not-installing-correctly
    '''
    if len(poly)==0:
        return np.zeros(len(x),np.bool_)

    n = len(poly)
    inside = np.zeros(len(x),np.bool_)
    p2x = 0.0
    p2y = 0.0
    xints = 0.0
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        idx = np.nonzero((y > min(p1y,p2y)) & (y <= max(p1y,p2y)) & (x <= max(p1x,p2x)))[0]
        if p1y != p2y:
            xints = (y[idx]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
        if p1x == p2x:
            inside[idx] = ~inside[idx]
        else:
            idxx = idx[x[idx] <= xints]
            inside[idxx] = ~inside[idxx]

        p1x,p1y = p2x,p2y
    return inside

def crange(n,a,b):
    '''
    Takes line segment [a,b] and breaks it into n equal length subsegments. Returns an array of the midpoints of these subsegments.
    '''
    return ( (np.arange(n)/n+1/(2*n)))*(b-a)+a

