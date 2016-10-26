from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
from math import sin, cos, sqrt, pi, acos 
from numpy import array, add, dot, cross, arccos
from sys import exit

class pycell:
    def __init__(self):

	# pymol object
	self._obj = None

	# pymol cell name
	self._name = None

        # cell parameters 
        self._a = None
	self._b = None
	self._c = None
	self._alpha = None
	self._beta = None
	self._gamma = None

	# cell alignment
	self._cell_alignment = "GMX"

	# unit cell vectors
        self._A = None
	self._B = None
	self._C = None

        # symmetry list for trajectories
        self._symmetrylist = []

    def find_cell(self):
	if self._obj != None:
            cell_param = cmd.get_symmetry(self._obj) 
            self._a = cell_param[0]
	    self._b = cell_param[1]
	    self._c = cell_param[2]
	    self._alpha = cell_param[3]
	    self._beta = cell_param[4]
	    self._gamma = cell_param[5]

    def set_cell(self, a, b, c, alpha, beta, gamma):
        self._a = a
	self._b = b
	self._c = c
	self._alpha = alpha
	self._beta = beta
	self._gamma = gamma
	
	# init cell vectors with new parameters
	self.init_cell()

    def draw_trajectory(self, filepath, radius=0.1, color_type = 'white', name = "cell"): 

        self._name = name
	self.read_trajctorysymmetry(filepath)
	N = 1
	for sym in self._symmetrylist:
            self.set_cell(sym[0],sym[1],sym[2],sym[3],sym[4],sym[5])
            self.draw(radius, color_type, N)
            N += 1

    def read_trajctorysymmetry(self, filepath):
        self._symmetrylist = []	    
        with open(filepath) as infile:
            for line in infile:
                if len(line)>0:
		    line = line.split()
		    #print line
		    if line[0] == 'CRYST1':
		        a = float(line[1])
		        b = float(line[2])
		        c = float(line[3])
		        alpha = float(line[4])
		        beta = float(line[5])
		        gamma = float(line[6])
                        self._symmetrylist.append([a,b,c,alpha,beta,gamma])

    def get_unit_cell_vectors_MS6(self, collinear_vector_length, plane_vector_length, free_vector_length, alpha, beta, gamma):
    
        # collinear vector
        collinear_vector = array([0,0,collinear_vector_length])
    
        # find plane_vector, rotate around (1,0,0)
        axis = array([1,0,0])
        plane_vector = array([0,1,0])
        angle = self.vector_angle(collinear_vector,array([0,0,0]),plane_vector) - alpha*(pi/180)
        plane_vector = self.rotate(plane_vector, axis, angle)
        plane_vector = self.norm_vector(plane_vector)*plane_vector_length  
    
        # find free_vector
        free_vector = array([1,0,0])
        angle = self.vector_angle(plane_vector, array([0,0,0]), free_vector) - gamma*(pi/180)
        free_vector = self.rotate(free_vector, collinear_vector, angle )
    
        angle = self.vector_angle(collinear_vector, array([0,0,0]), free_vector) - beta*(pi/180) 
        axis = self.perpvector(collinear_vector, array([0,0,0]), free_vector)
        free_vector = self.rotate(free_vector, axis, angle )
        free_vector = self.norm_vector(free_vector)*free_vector_length

        return collinear_vector, plane_vector, free_vector

    def get_unit_cell_vectors_GMX(self, collinear_vector_length, plane_vector_length, free_vector_length, alpha, beta, gamma):
    
        # collinear vector
        collinear_vector = array([collinear_vector_length,0,0])
    
        # find plane_vector, rotate around (1,0,0)
        axis = array([0,0,1])
        plane_vector = array([0,1,0])
        angle = self.vector_angle(collinear_vector,array([0,0,0]),plane_vector) - gamma*(pi/180)
        plane_vector = self.rotate(plane_vector, -axis, angle)
        plane_vector = self.norm_vector(plane_vector)*plane_vector_length  
    
        # find free_vector
        free_vector = array([0,0,1])
        angle = self.vector_angle(plane_vector, array([0,0,0]), free_vector) - alpha*(pi/180)
        free_vector = self.rotate(free_vector, -collinear_vector, angle )
    
        angle = self.vector_angle(collinear_vector, array([0,0,0]), free_vector) - beta*(pi/180) 
        axis = self.perpvector(collinear_vector, array([0,0,0]), free_vector)
        free_vector = self.rotate(free_vector, axis, angle )
        free_vector = self.norm_vector(free_vector)*free_vector_length

        return collinear_vector, plane_vector, free_vector

    def set_cell_aligment(self, cell_aligment = "MS6" ):
        self._cell_alignment = cell_aligment 

    def init_cell(self):

        if (self._a == None):
	    print "cell parameters not set"
	    exit()

        if self._cell_alignment == "MS6":  
            # material_studio 6, C along Z, B in YZ plane
            self._A, self._B, self._C = self.get_unit_cell_vectors_MS6(self._c, self._b, self._a, self._alpha, self._beta, self._gamma)
        elif self._cell_alignment == "GMX":  
            # Gromacs, A along X, B in YX plane
            self._A, self._B, self._C = self.get_unit_cell_vectors_GMX(self._a, self._b,self._c, self._alpha, self._beta, self._gamma)
        else:
	    print "ERROR in cell alignment"
	    exit()

    def draw(self, radius=0.1, color_type = 'white', state = 1):
	
        size = radius*15.
        origin_offset = radius * -25.
      
        vert_000 = [0.,0.,0] 
        vert_100 = self._A.tolist()  
        vert_010 = self._B.tolist()
        vert_001 = self._C.tolist()
        vert_110 = add(self._A,self._B).tolist() 
        vert_011 = add(self._B,self._C).tolist()
        vert_101 = add(self._A,self._C).tolist()
        vert_111 = add( add(self._A,self._B), self._C).tolist()
     
        color_list = self.get_colorlist(color_type)

        cell = []
        cell.append(CYLINDER)
        cell = cell + vert_000 + vert_100 + [radius] + color_list[0][0] + color_list[0][1] 
        cell.append(CYLINDER)
        cell = cell + vert_000 + vert_010 + [radius] + color_list[1][0] + color_list[1][1]
        cell.append(CYLINDER)
        cell = cell + vert_000 + vert_001 + [radius] + color_list[2][0] + color_list[2][1]
        cell.append(CYLINDER)
        cell = cell + vert_100 + vert_110 + [radius] + color_list[3][0] + color_list[3][1]
        cell.append(CYLINDER)
        cell = cell + vert_100 + vert_101 + [radius] + color_list[4][0] + color_list[4][1]
        cell.append(CYLINDER)
        cell = cell + vert_010 + vert_110 + [radius] + color_list[5][0] + color_list[5][1]
        cell.append(CYLINDER)
        cell = cell + vert_010 + vert_011 + [radius] + color_list[6][0] + color_list[6][1]
        cell.append(CYLINDER)
        cell = cell + vert_001 + vert_101 + [radius] + color_list[7][0] + color_list[7][1]
        cell.append(CYLINDER)
        cell = cell + vert_001 + vert_011 + [radius] + color_list[8][0] + color_list[8][1]
        cell.append(CYLINDER)
        cell = cell + vert_110 + vert_111 + [radius] + color_list[9][0] + color_list[9][1]
        cell.append(CYLINDER)
        cell = cell + vert_101 + vert_111 + [radius] + color_list[10][0] + color_list[10][1]
        cell.append(CYLINDER)
        cell = cell + vert_011 + vert_111 + [radius] + color_list[11][0] + color_list[11][1]
      
        cmd.load_cgo(cell, self._name, state)
      
      #  text = [COLOR, 1.0, 0.0, 1.0,]
      
      #  cyl_text(text,plain,[origin_offset,origin_offset,-1],'Origin',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]],color=[1.0,0.0,1.0])
      #  cyl_text(text,plain,map(None,U.orthogonalize((1.05,0.0,0.0))),'A',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]],color=[1.0,0.0,0.0])
      #  cyl_text(text,plain,map(None,U.orthogonalize((0.0,1.05,0.0))),'B',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]],color=[0.0,1.0,0.0])
      #  cyl_text(text,plain,map(None,U.orthogonalize((0.0,0.0,1.05))),'C',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]],color=[0.0,0.0,1.0])
      
      #  cmd.load_cgo(text,name+'_labels')


    def get_colorlist(self, color_type):

	color_list = [],[],[],[],[],[],[],[],[],[],[],[]

	if color_type == 'white':
            for i in range(0, len(color_list)):
                color_list[i].append([1,1,1])
                color_list[i].append([1,1,1])
        elif color_type == 'black':
            for i in range(0, len(color_list)):
                color_list[i].append([0,0,0])
                color_list[i].append([0,0,0])
        elif color_type == 'red':
            for i in range(0, len(color_list)):
                color_list[i].append([1,0,0])
                color_list[i].append([1,0,0])
        elif color_type == 'green':
            for i in range(0, len(color_list)):
                color_list[i].append([0,1,0])
                color_list[i].append([0,1,0])
        elif color_type == 'blue':
            for i in range(0, len(color_list)):
                color_list[i].append([0,1,0])
                color_list[i].append([0,1,0])
        elif color_type == 'rainbow':
            color_list[0].append([0,0,0])
            color_list[0].append([1,0,0])

            color_list[1].append([0,0,0])
            color_list[1].append([0,1,0])

            color_list[2].append([0,0,0])
            color_list[2].append([0,0,1])

            color_list[3].append([1,0,0])
            color_list[3].append([1,1,0])

            color_list[4].append([1,0,0])
            color_list[4].append([1,0,1])

            color_list[5].append([0,1,0])
            color_list[5].append([1,1,0])

            color_list[6].append([0,1,0])
            color_list[6].append([0,1,1])

            color_list[7].append([0,0,1])
            color_list[7].append([1,0,1])

            color_list[8].append([0,0,1])
            color_list[8].append([0,1,1])

            color_list[9].append([1,1,0])
            color_list[9].append([1,1,1])

            color_list[10].append([1,0,1])
            color_list[10].append([1,1,1])

            color_list[11].append([0,1,1])
            color_list[11].append([1,1,1])

        else:
            print 'ERROR: WRONG COLOR!!!'
            exit()	

        return color_list

    def perpvector(self, v1, v2, v3):
            J = cross(v2 - v3, v2 - v1)+ v2
            J = (J - v2) /(sqrt(dot(J - v2, J - v2))) + v2
            return J
    
    def rotate(self, V, J, T):
            x = V[0]
            y = V[1]
            z = V[2]
            u = J[0]
            v = J[1]
            w = J[2]
            a = (u*(u*x + v*y + w*z) + (x * (v*v + w*w) - u *(v*y + w*z))*cos(T) + sqrt(u*u + v*v + w*w)*(-w*y + v*z)*sin(T))/(u*u + v*v + w*w)
            b = (v*(u*x + v*y + w*z) + (y * (u*u + w*w) - v *(u*x + w*z))*cos(T) + sqrt(u*u + v*v + w*w)*(w*x - u*z)*sin(T))/(u*u + v*v + w*w)
            c = (w*(u*x + v*y + w*z) + (z * (u*u + v*v) - w *(u*x + v*y))*cos(T) + sqrt(u*u + v*v + w*w)*(-v*x + u*y)*sin(T))/(u*u + v*v + w*w)
            return array([a, b, c])
    
    def vector_angle(self, a, b, c):
    
        # In case numpy.dot() returns larger than 1
        # and we cannot take acos() to that number
        acos_out_of_bound = 1.0
        v1 = a - b
        v2 = c - b
        v1 = v1 / sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
        v2 = v2 / sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
        dot_product = dot(v1,v2)
    
        if dot_product > acos_out_of_bound:
            dot_product = acos_out_of_bound
        if dot_product < -1.0 * acos_out_of_bound:
            dot_product = -1.0 * acos_out_of_bound
    
        return arccos(dot_product)
        
    def norm_vector(self, vector):
    
        length = sqrt((vector[0]*vector[0]) + (vector[1]*vector[1]) + (vector[2]*vector[2]))
        norm_vector = array([vector[0]/length,vector[1]/length,vector[2]/length])

        return norm_vector

# wrapper functions for pycell
def pycell_start(obj=None,name="cell"):

    newcell._obj = obj
    newcell._name = name

def pycell_auto(obj=None, radius = 0.1, color_type = 'white', name= "cell"):

    """
    use as 
        load test.pdb
        pycell_auto test
	or
	pycell_auto test, 0.05, black, cell

    """
    newcell._obj = obj
    newcell._name = name
    newcell.find_cell()
    newcell.init_cell()
    newcell.draw(float(radius), color_type)

def pycell_reset():
    del newcell
    newcell = pycell()

def pycell_setcell(a,b,c,alpha,beta,gamma):
    newcell.set_cell(float(a),float(b),float(c),float(alpha),float(beta),float(gamma))

def pycell_draw( radius = 0.1, color_type = "white", state = 1):
    newcell.draw(float(radius), color_type, int(state))

def pycell_setaligment( cell_type ):
    newcell.set_cell_aligment( cell_type )

def pycell_trajectory( filepath, radius = 0.1, color_type = 'white', name = "cell" ):
    newcell.draw_trajectory( filepath, float(radius), color_type, name)

# pycell object
newcell = pycell()

# pymol commands 
cmd.extend("pycell_start",pycell_start)
cmd.extend("pycell_auto",pycell_auto)
cmd.extend("pycell_setcell", pycell_setcell)
cmd.extend("pycell_draw", pycell_draw)
cmd.extend("pycell_setaligment", pycell_setaligment )
cmd.extend("pycell_reset", pycell_reset)
cmd.extend("pycell_trajectory", pycell_trajectory) 
