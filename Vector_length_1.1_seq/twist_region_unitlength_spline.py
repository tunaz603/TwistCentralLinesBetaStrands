# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 23:01:08 2016

@author: tunaz
"""


import sys

import re
import numpy
import math

from scipy import interpolate

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self) :
        for f in self.files:
            f.flush()
            
            
if len(sys.argv) == 1 or len(sys.argv) == 2 or len(sys.argv) > 3  :
    print "Usage: python twist_region.py pdbfilename.pdb direction"
    sys.exit(0)


lines_File = [lines1.rstrip('\n') for lines1 in open(sys.argv[1], 'r')]   


try:
    direction = int(sys.argv[2]) # give direction as input
    '''
    direction= 0 means forward direction, 
    direction = 1 means backword direction 
    '''
    if direction != 0 and direction != 1 :
        print "please input either 0 or 1. 0 = forward, 1 = backword"
        sys.exit(0)
except ValueError:  
        print "please input either 0 or 1. 0 = forward, 1 = backword"
        sys.exit(0)          

#sys.stdout = open("output_all.pdb", "w")
out_filename = 'twist_angles.txt'
outf = open(out_filename, 'w')

sys.stdout = open('output_pdb.pdb','w')
#outf2 = open(out_filename2, 'w')
#lines_File = [lines1.rstrip('\n') for lines1 in open('1aky.pdb')]    
#lines_File = [lines1.rstrip('\n') for lines1 in open('1aop.pdb')]  
#lines_File = [lines1.rstrip('\n') for lines1 in open('2p8y.pdb')] #index problem
#lines_File = [lines1.rstrip('\n') for lines1 in open('1s04.pdb')] #complex structure (multple same chain duplicate index)
#lines_File = [lines1.rstrip('\n') for lines1 in open('1atz.pdb')] 
#lines_File = [lines1.rstrip('\n') for lines1 in open('2qtr.pdb')] 
#lines_File = [lines1.rstrip('\n') for lines1 in open('1d5t.pdb')] 
#lines_File = [lines1.rstrip('\n') for lines1 in open('1elu.pdb')]
#lines_File = [lines1.rstrip('\n') for lines1 in open('1a12.pdb')]
#lines_File = [lines1.rstrip('\n') for lines1 in open('1b3a.pdb')]

#direction = 0 # direction= 0 means forward direction, direction = 1 means backword direction so change here        
            
def draw_interpolation(input_array):
    #print "input  ", input_array
    #print "length  ", len(input_array)
    if  len(input_array) > 2 :
        data = numpy.array(input_array)
        data = data.astype(float)
        
       # print data #ok
           
        data = data[0:100].transpose()
        #print "hi", data #ok
        #now we get all the knots and info about the interpolated spline
		#spacing = int(len(input_array)/.3)
        tck, u= interpolate.splprep(data,k=2)   
        #new = interpolate.splev(numpy.linspace(0,1,50), tck)
        #new = interpolate.splev(numpy.linspace(0,1,int(len(input_array)/.8)), tck) #unit length given 0.1
        #new = interpolate.splev(numpy.linspace(0,1,int(len(input_array)/.5)), tck) #we used it for ACM_BCB
        new = interpolate.splev(numpy.linspace(0, 1, int(len(input_array) /.4)), tck) #vector length 1.1

        #now lets plot it!
        
        #fig = plt.figure()
        #ax = Axes3D(fig)
        #ax.plot(data[0], data[1], data[2], label='originalpoints', lw =2, c='Dodgerblue')
        #ax.plot(new[0], new[1], new[2], label='fit', lw =2, c='red')
        #ax.legend()
        #plt.savefig('junk.png')
        #plt.show() 
        return new[0], new[1], new[2]             


#print (lines_File1)

direction_opposite = 1 #hydrogen bond 

count_beta_sheet = 0
count_beta_strand = 0
counter = 0
counter_index = []
problem = 0
a=[]
beta_sheet = []
atom = []
splitted_line = []
count_atom = 0
for sheet_line in lines_File: #find each line in all lines
     if sheet_line.startswith('SHEET'):     #find those lines start with "SHEET"
         #count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "SHEET"        
         #beta_sheet.append ( re.sub("\s\s+" , " ",sheet_line).split()) #ok for all except 1
         beta_sheet.append([sheet_line[0:6], sheet_line[7:10], sheet_line[11:14], sheet_line[14:16],sheet_line[17:20], sheet_line[21], sheet_line[22:27], sheet_line[28:31], sheet_line[32], sheet_line[33:37], sheet_line[37], sheet_line[38:40],sheet_line[41:45],sheet_line[45:48], sheet_line[49], sheet_line[50:55], sheet_line[56:60], sheet_line[60:63], sheet_line[64], sheet_line[65:69], sheet_line[69]])
     if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         #atom.append( [sheet_line[0:6], sheet_line[6:11], sheet_line[12:16], sheet_line[16:17],sheet_line[17:20], sheet_line[21:22], sheet_line[22:26], sheet_line[26:27], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54], sheet_line[46:54], sheet_line[54:60], sheet_line[60:66], sheet_line[76:78], sheet_line[78:80]])
         #atom.append ( re.sub("\s\s+" , " ",sheet_line).split())
#print (int(atom[1836][5].strip()))
#print (atom[1836])

strand_vector=[]
center_vector = []
CA_vector_between_strand = []
amino_acid = []


#common hydrogen bond
strand_vector_opposite=[]
center_vector_opposite = []
CA_vector_between_strand_opposite = []
amino_acid_opposite = []


poi = 0
track = 1

poi_opposite = 0 #hydrogen bond


for sheet_index in range (0,len(beta_sheet)): # handle if there is less than 3 atoms in a strand
    if int(beta_sheet[sheet_index][1]) == int(beta_sheet[sheet_index][3]):
        counter = counter + 1

#print ("number of beta sheet ",counter) #ok
out_char = "number of beta sheet " + str(counter) + '\n\n'
outf.write(out_char)            
for sheet_index in range (0,len(beta_sheet)): # handle if there is less than 3 atoms in a strand
    if abs((int(beta_sheet[sheet_index][9]))  - int(beta_sheet[sheet_index][6])) >= 2 :
        if int(beta_sheet[sheet_index][1]) == int(beta_sheet[sheet_index][3]):
            #counter = counter + 1
            if track == 0 :
                counter_index.append(int(beta_sheet[sheet_index][1]) - problem)
                track = 1
            else:   
                counter_index.append(int(beta_sheet[sheet_index][1]))
                track = 1
    else:
        #print "hi " #ok
        if  (int(beta_sheet[sheet_index][3]) > 2):           
             problem = problem + 1
             #print problem  #ok
             track = 0
            
            
#for i in range (0,len(counter_index)):
    #print ("counter_index ",counter_index[i]) #ok  
     
for sheet_index in range (0,len(beta_sheet)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
    atom_start_index=(int(beta_sheet[sheet_index][6]))
    #print("start index ",atom_start_index) #ok
    atom_end_index=(int(beta_sheet[sheet_index][9])) 
    #print("end index ",atom_end_index) #ok
    sheet_chain = beta_sheet[sheet_index][5]
    #print sheet_chain

    mid_point = []
    array_strand_index =[]
    amino_acid_array = []
    aa_index = []
    #array_strand_index = [[0 for i in range(2)] for i in range(2)] #2 should be change according to length of one beta sheet
    array_x = []
    array_y = []
    array_z = []
    CA_x = []
    CA_y = []
    CA_z = []
    
    #print (count_atom)
    
    for i in range (0,count_atom) : 
        
        if (atom_start_index <= int(atom[i][5].strip()) <=  atom_end_index)  :
            
            if  atom[i][2].strip() == 'N' or atom[i][2].strip() == 'C': 
                #print "duplicate ", int(atom[i][5].strip())
                
                if sheet_chain == atom[i][4].strip(): #this is needed to differentiate same startindex & endindex but different chains
                    array_strand_index.append(atom[i][5].strip())                  
                    array_x.append(atom[i][6].strip())
                    array_y.append(atom[i][7].strip())
                    array_z.append(atom[i][8].strip())
                    length=len(array_strand_index)  
            if atom[i][2].strip() == 'CA':
                if sheet_chain == atom[i][4].strip():
                    aa_index.append( i )             
                    amino_acid_array.append(atom[i][3].strip())        
                    CA_x.append(atom[i][6].strip())
                    CA_y.append(atom[i][7].strip())
                    CA_z.append(atom[i][8].strip())
                    length_CA = len (CA_x)
    #print ("strand index  ",array_strand_index) #ok
   # print (aa_index)
    #print ("Amino acid  ",amino_acid_array)  #ok        
   
   
   
    
    if direction == 0: #forward direction #ok
        midpoint_x_bef = []
        midpoint_y_bef = []
        midpoint_z_bef = []  
        
      
        
        midpoint_x = []
        midpoint_y = []
        midpoint_z = []
        
        new_x = []
        new_y = []
        new_z = []  
        
       
        #coef = [1,1,1,1,1,1,1,1]
        #coef = []
        
        c_index = 1
        n_index = 2
        while c_index < length-1 and n_index <  length-1 :  #mid point calculation
            #print ("c_index ",c_index, "n_index ", n_index)
            midpoint_x_bef.append((float(array_x[c_index]) + float(array_x[n_index]))/float(2.0))    #making midpoint       
            midpoint_y_bef.append((float(array_y[c_index]) + float(array_y[n_index]))/float(2.0))   #making midpoint 
            midpoint_z_bef.append((float(array_z[c_index]) + float(array_z[n_index]))/float(2.0))  #making midpoint 
           
            c_index = c_index + 2
            n_index = n_index +2  
        
        #print ("midpoint_x_bef  ", midpoint_x_bef ) 
        #print ("midpoint_y_bef  ", midpoint_y_bef  )
        #print ("midpoint_z_bef  ", midpoint_z_bef  )
            
        midpoint_draw = []
        draw_index = 0
        while draw_index < len (midpoint_x_bef) :
            midpoint_point = []
            mid_x = midpoint_x_bef [draw_index]
            mid_y = midpoint_y_bef [draw_index]
            mid_z = midpoint_z_bef [draw_index]           
            midpoint_point.append (mid_x)
            midpoint_point.append (mid_y)
            midpoint_point.append (mid_z)
            midpoint_draw.append (midpoint_point) 
            draw_index = draw_index + 1
        #print   (midpoint_draw)  #ok
        pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s6s10s2s3s'
        #shape = 'ATOM 1 N LEU A 81 25.865 39.612 19.376 1.00 53.28 N '
        shape = ['ATOM', '1','N', 'LEU', 'A', '81', '25.865', '39.612', '19.376']
        if len(midpoint_draw) <= 2:
             midpoint_x=midpoint_x_bef
             midpoint_y=midpoint_y_bef
             midpoint_z=midpoint_z_bef
             
             for i in range (len( midpoint_x)):
                poi = poi +1
                                           
                #print "%-6s%5s %4s %3s %s%4s    %8s%8s%8s\n"%tuple(shape) ok
                print "%-6s%5s %4s %3s %s%4s    %8s%8s%8s"%tuple(['ATOM',poi,'S','LEU', 'A', '81',round(midpoint_x[i],3),round(midpoint_y[i],3),round(midpoint_z[i],3)]) #print new atom position in PDB file
                
             
        else:    
            value = draw_interpolation (midpoint_draw)   #making 3D cubic spline interpolation for midpoints
            #print value[0] #ok
            midpoint_x=value[0]  #new midpoint after cubic spline
            midpoint_y=value[1]  #new midpoint after cubic spline
            midpoint_z=value[2]   #new midpoint after cubic spline
        
            
           
            pdb_length = len( midpoint_x)
          
            for i in range (pdb_length):
                poi = poi +1          
                print "%-6s%5s %4s %3s %s%4s    %8s%8s%8s"%tuple(['ATOM',poi,'S','LEU', 'A', '81',round(midpoint_x[i],3),round(midpoint_y[i],3),round(midpoint_z[i],3)]) #print new atom position in PDB file
                
       
       
        
        
       
        #print ("midpoint_x ", midpoint_x ) 
         
        
     
     
        length_midpoint_array = len(midpoint_x)
        vector_array =[]
        cen_array = []
        midpoint_index = 1
        while midpoint_index < length_midpoint_array:
            point = []
            cen_point = []
            temp_mid_index = midpoint_index - 1
            vec_x = midpoint_x [temp_mid_index] -  midpoint_x [temp_mid_index+1]  #making midpoint vector
            vec_y = midpoint_y [temp_mid_index] -  midpoint_y [temp_mid_index+1]  #making midpoint vector
            vec_z = midpoint_z [temp_mid_index] -  midpoint_z [temp_mid_index+1]  #making midpoint vector
            
            cen_x = (float(midpoint_x [temp_mid_index]) + float( midpoint_x [temp_mid_index+1])) / float(2.0)  #making center of midpoint vector
            cen_y =(float( midpoint_y [temp_mid_index]) + float(midpoint_y [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
            cen_z = (float(midpoint_z [temp_mid_index]) + float( midpoint_z [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
            
            #print ("vector x ",vec_x) 
           # print ("vector y ",vec_y) 
           # print ("vector z ",vec_z) 
            #make unit vector
            magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
            
           
            print("vector length ",magnitude)
            
            point.append(vec_x/magnitude) # unit vector
            point.append(vec_y/magnitude)  # unit vector
            point.append(vec_z/magnitude)  # unit vector
            
      
            #point.append(vec_x) #ok
            #point.append(vec_y)
            #point.append(vec_z)
      
         
            cen_point.append(cen_x) #ok
            cen_point.append(cen_y)
            cen_point.append(cen_z)
            
            
            
            
            vector_array.append(point)
            
            cen_array.append(cen_point)
            
            midpoint_index = midpoint_index +1
        #print "vector array ",vector_array
        
        
        aa = []
        aa_index = 1
        while aa_index < len(amino_acid_array) - 1:   # to get amino acid value
             aa.append(amino_acid_array[aa_index]) 
             aa_index = aa_index +1
        amino_acid.append(aa)
        
        #print ("x cordinate of CA in each strand" CA_x )
    
        #print("length ", length)
        #print ("CA length ",length_CA)
        CA_x_new = []
        CA_y_new = []
        CA_z_new = []
        CA_vector = []
        CA_index = 1
        while CA_index < length_CA - 1:
            CA_point = []
            CA_x_new=(float(CA_x[CA_index]))
            CA_y_new=(float(CA_y[CA_index]))
            CA_z_new=(float(CA_z[CA_index]))
            CA_point.append (CA_x_new)
            CA_point.append (CA_y_new)
            CA_point.append (CA_z_new)
            CA_vector.append (CA_point)   
            CA_index = CA_index + 1
        #draw_interpolation (CA_vector) #draw spline according to carbon alpha 

        CA_vector_between_strand.append(CA_vector)
        #draw_interpolation (CA_vector_between_strand)

         
    if direction_opposite == 1:   #backword direction
    
        midpoint_x_bef_opposite = []
        midpoint_y_bef_opposite = []
        midpoint_z_bef_opposite = []  

        midpoint_x_opposite = []
        midpoint_y_opposite = []
        midpoint_z_opposite = [] 
        
       
        
        c_index_opposite = length-3
        n_index_opposite = length-2
        while c_index_opposite > 0  and n_index_opposite > 0 : 
            #print ("c_index ",c_index, "n_index ", n_index)
            midpoint_x_bef_opposite.append((float(array_x[c_index_opposite]) + float(array_x[n_index_opposite]))/float(2.0))
            midpoint_y_bef_opposite.append((float(array_y[c_index_opposite]) + float(array_y[n_index_opposite]))/float(2.0))
            midpoint_z_bef_opposite.append((float(array_z[c_index_opposite]) + float(array_z[n_index_opposite]))/float(2.0)) 
            c_index_opposite = c_index_opposite - 2
            n_index_opposite = n_index_opposite - 2   
       
       
       
        midpoint_draw_opposite = []
        draw_index_opposite = len (midpoint_x_bef_opposite) -1
        while draw_index_opposite >= 0 :
            midpoint_point_opposite = []
            mid_x_opposite = midpoint_x_bef_opposite [draw_index_opposite]
            mid_y_opposite = midpoint_y_bef_opposite [draw_index_opposite]
            mid_z_opposite = midpoint_z_bef_opposite [draw_index_opposite]           
            midpoint_point_opposite.append (mid_x_opposite)
            midpoint_point_opposite.append (mid_y_opposite)
            midpoint_point_opposite.append (mid_z_opposite)
            midpoint_draw_opposite.append (midpoint_point_opposite) 
            draw_index_opposite = draw_index_opposite - 1
            
        #draw_interpolation (midpoint_draw)
       
        if len(midpoint_draw_opposite) <= 2:
             midpoint_x_opposite = midpoint_x_bef_opposite 
             midpoint_y_opposite = midpoint_y_bef_opposite
             midpoint_z_opposite = midpoint_z_bef_opposite
             
             for i in range (len( midpoint_x_opposite)):
                poi_opposite = poi_opposite +1
                                           
                #print "%-6s%5s %4s %3s %s%4s    %8s%8s%8s\n"%tuple(shape) ok
               # print "%-6s%5s %4s %3s %s%4s    %8s%8s%8s"%tuple(['ATOM',poi_opposite,'S','LEU', 'A', '81',round(midpoint_x_opposite[i],3),round(midpoint_y_opposite[i],3),round(midpoint_z_opposite[i],3)]) #print new atom position in PDB file
                
             
        else:    
            value_opposite = draw_interpolation (midpoint_draw_opposite)   #making 3D cubic spline interpolation for midpoints
            #print value_opposite[0] #ok
            midpoint_x_opposite = value_opposite[0]  #new midpoint after cubic spline
            midpoint_y_opposite = value_opposite[1]  #new midpoint after cubic spline
            midpoint_z_opposite = value_opposite[2]   #new midpoint after cubic spline
        
            
           
            pdb_length_opposite = len( midpoint_x_opposite)
          
            for i in range (pdb_length_opposite):
                poi_opposite = poi_opposite +1          
                #print "%-6s%5s %4s %3s %s%4s    %8s%8s%8s"%tuple(['ATOM',poi_opposite,'S','LEU', 'A', '81',round(midpoint_x_opposite[i],3),round(midpoint_y_opposite[i],3),round(midpoint_z_opposite[i],3)]) #print new atom position in PDB file
                
       
       
        midpoint_x_opposite = midpoint_x_bef_opposite
        midpoint_y_opposite = midpoint_y_bef_opposite
        midpoint_z_opposite = midpoint_z_bef_opposite
        
    
        #print ("midpoint_x ", midpoint_x )      
        length_midpoint_array_opposite = len (midpoint_x_opposite)   
        vector_array_opposite =[]
        cen_array_opposite = []
        
        #midpoint_index_opposite = length_midpoint_array
        #while midpoint_index > 1:
        
        midpoint_index_opposite = 1 #big change in this loop
        while midpoint_index_opposite < length_midpoint_array_opposite:
            point_opposite = []
            cen_point_opposite = []
            
            temp_mid_index_opposite = midpoint_index_opposite - 1
            #vec_x = midpoint_x [temp_mid_index] -  midpoint_x [temp_mid_index-1]
            #vec_y = midpoint_y [temp_mid_index] -  midpoint_y [temp_mid_index-1]
            #vec_z = midpoint_z [temp_mid_index] -  midpoint_z [temp_mid_index-1]
            
            vec_x_opposite = midpoint_x_opposite [temp_mid_index_opposite] -  midpoint_x_opposite [temp_mid_index_opposite+1]
            vec_y_opposite = midpoint_y_opposite [temp_mid_index_opposite] -  midpoint_y_opposite [temp_mid_index_opposite+1]
            vec_z_opposite = midpoint_z_opposite [temp_mid_index_opposite] -  midpoint_z_opposite [temp_mid_index_opposite+1]
            
            #print ("vector ",vec_x,vec_y,vec_z) 
            
            #cen_x = (float(midpoint_x [temp_mid_index]) + float( midpoint_x [temp_mid_index-1])) / float(2.0)  #making center of midpoint vector
            #cen_y =(float( midpoint_y [temp_mid_index]) + float(midpoint_y [temp_mid_index-1]))/float(2.0)   #making center of midpoint vector
            #cen_z = (float(midpoint_z [temp_mid_index]) + float( midpoint_z [temp_mid_index-1]))/float(2.0)   #making center of midpoint vector
            
            cen_x_opposite = (float(midpoint_x_opposite [temp_mid_index_opposite]) + float( midpoint_x_opposite [temp_mid_index_opposite+1])) / float(2.0)  #making center of midpoint vector
            cen_y_opposite =(float( midpoint_y_opposite [temp_mid_index_opposite]) + float(midpoint_y_opposite [temp_mid_index_opposite+1]))/float(2.0)   #making center of midpoint vector
            cen_z_opposite = (float(midpoint_z_opposite [temp_mid_index_opposite]) + float( midpoint_z_opposite [temp_mid_index_opposite+1]))/float(2.0)   #making center of midpoint vector
            
              #make unit vector
            magnitude_opposite = math.sqrt(pow(vec_x_opposite, 2) + pow(vec_y_opposite, 2) + pow(vec_z_opposite, 2))
            #print("vector length opposite ",magnitude_opposite)
            #point.append(vec_x/magnitude) #not ok
            #point.append(vec_y/magnitude)
            #point.append(vec_z/magnitude)
            
            '''
            point_opposite.append(vec_x_opposite) #ok
            point_opposite.append(vec_y_opposite)
            point_opposite.append(vec_z_opposite)
            '''
            point_opposite.append(vec_x_opposite/magnitude_opposite) #ok
            point_opposite.append(vec_y_opposite/magnitude_opposite)
            point_opposite.append(vec_z_opposite/magnitude_opposite)
            
            
            
            cen_point_opposite.append(cen_x_opposite) #ok
            cen_point_opposite.append(cen_y_opposite)
            cen_point_opposite.append(cen_z_opposite)
     
            vector_array_opposite.append(point_opposite)
            
            cen_array_opposite.append(cen_point_opposite)
            
            midpoint_index_opposite = midpoint_index_opposite + 1
            
            #midpoint_index = midpoint_index - 1 
        
        aa_opposite = []
        aa_index_opposite = len(amino_acid_array) - 2
        while aa_index_opposite > 0:   # to get amino acid value
             aa_opposite.append(amino_acid_array[aa_index_opposite]) 
             aa_index_opposite = aa_index_opposite - 1 
        amino_acid_opposite.append(aa_opposite)
   
   
        #print ("x cordinate of CA in each strand", CA_x )
    
        #print("length ", length)
        #print ("CA length ",length_CA)
        CA_x_new_opposite = []
        CA_y_new_opposite = []
        CA_z_new_opposite = []
        CA_vector_opposite = []
        CA_index_opposite = length_CA - 2
        while CA_index_opposite > 0 :
            CA_point_opposite = []
            CA_x_new_opposite=(float(CA_x[CA_index_opposite]))
            CA_y_new_opposite=(float(CA_y[CA_index_opposite]))
            CA_z_new_opposite=(float(CA_z[CA_index_opposite]))
            CA_point_opposite.append (CA_x_new_opposite)
            CA_point_opposite.append (CA_y_new_opposite)
            CA_point_opposite.append (CA_z_new_opposite)
            CA_vector_opposite.append (CA_point_opposite)   
            CA_index_opposite = CA_index_opposite - 1
        CA_vector_between_strand_opposite.append(CA_vector_opposite)
        
        
    strand_vector.append(vector_array)  
    
    center_vector.append(cen_array)
    
    #for hydrogen bond calculation
    strand_vector_opposite.append(vector_array_opposite)  
    
    center_vector_opposite.append(cen_array_opposite)
    
   #print(beta_sheet[sheet_index][1])  
    
    

 
   
    
    '''
    print ("vector ",vector_array)
    print (atom_start_index)
    print (array_strand_index) 
    print(length)
    
    print (array_x) 
    print (array_y) 
    print (array_z) 
    print (midpoint_x)
    print (midpoint_y)
    print (midpoint_z)
    print (length_midpoint_array)
    '''
print "\n"
length_strand_vector = len(strand_vector)
#print (strand_vector)
#print(amino_acid) # ok
#print ("vector strand length ",length_strand_vector)
#print ("CA_vector_between_strand ",len(CA_vector_between_strand)) # ok.  vector_strand length and CA_vector same
#print ("CA_vector_between_strand ",CA_vector_between_strand) 
#draw_interpolation (CA_vector_between_strand)

def distance(x,y):   
    return numpy.sqrt(numpy.sum((x-y)**2))


def find_nearest_distance(a,B, pair_aa1, pair_aa2):
    
    
    #print ("a ", a)
    min_point = B[0]
    #min_dist = abs(numpy.linalg.norm(a- numpy.array(B[0])))
    min_dist = math.sqrt ( pow((a[0]-B[0][0]),2) + pow((a[1]-B[0][1]),2) + pow((a[2]-B[0][2]),2)  )
    #min_dist = 12000
   # print ("mindistance ", min_dist)
    #print("len of B",len(B)) #ok
    aa_min_index = 0
    for i in range (0, len(B)): 
        #dist = []
                       
        b = numpy.array(B[i])        
        #print ("b ", b) #ok
        #dist=(numpy.linalg.norm(a-b))
        try_dist = math.sqrt ( pow((a[0]-b[0]),2) + pow((a[1]-b[1]),2) + pow((a[2]-b[2]),2)  )
        #print ("distance ", try_dist)
        #dist=abs(distance(a,b))
        if try_dist < min_dist:
            min_dist = try_dist
            min_point = b
            aa_min_index = i
 
    
    return (min_point, aa_min_index)


def calculate_twist_angle(a,minimum_dist_point,j,aa_min_index2):
        
    #print ("mindist point ", minimum_dist_point)
    dot_product = numpy.dot(a,minimum_dist_point) # A.B
    #print (dot_product)
       
    a_sum = pow(a[0],2)+pow(a[1],2)+pow(a[2],2)
    a_root=math.sqrt(a_sum)
    b_sum = pow(minimum_dist_point[0],2)+pow(minimum_dist_point[1],2)+pow(minimum_dist_point[2],2)
    b_root=math.sqrt(b_sum)
    a_b = a_root*b_root
    
    angle_dot = math.acos(dot_product/a_b)
    #print (angle_dot)
    twist_angle = math.degrees(angle_dot)
    
    '''
    ### cross product ###
    cross_product = numpy.cross(a,minimum_dist_point) # A*B
    #print (cross_product)
    sum_square_cross_product = pow(cross_product[0],2) +  pow(cross_product[1],2)+ pow(cross_product[2],2)
    #print (sum_square_cross_product)
    root_sum_square_cross_product = math.sqrt(sum_square_cross_product)
    #print (root_sum_square_cross_product)
   
   
   
    angle_cross = math.asin(root_sum_square_cross_product/a_b)
    #print (angle_cross)
    print(math.degrees(angle_cross))
    ######## end of cross product ###########
    '''
    if twist_angle > 90.:
        twist_angle = 180.0 - twist_angle
   # sys.stdout = original
   # print(" Amino Acid",pair_aa1[j], " and ", pair_aa2[aa_min_index2],":", "Twist angle::" , twist_angle) #necessary
    #print("Twist angle::" , round(twist_angle,3)) #necessary
    out_char = "Twist angle:: " + str(round(twist_angle,3)) + '\n'
    outf.write(out_char)
    return twist_angle
       # print ( "amino acid ", pair_aa1[j], " ", pair_aa2[aa_min_index2])
minimum_dist_point = []  
minimum_dist_point_dup = []
flag = 0
all_avg_twist_angle = []  
avg_count = 0    
avg_track = 1
pair_array = []
len_array = []
i = 0
while i < length_strand_vector : # pairwise calculation 
    #print ("beta sheet index", beta_sheet[i][1])
    #print ("beta strand", beta_sheet[i])
    
    #print(int(beta_sheet[i][1]), " : ", int(beta_sheet[i][3]))
    #print ("i ", i)
    #print ("flag ",flag)
    avg_twist_count = []
            
    if(int(beta_sheet[i][1])>=int(beta_sheet[i][3])):
        flag = int(beta_sheet[i][1])
        avg_count = 0
       
        #print ("avg count ",avg_count)
    else:
        flag = flag + 1
        avg_count = avg_count + 1
        A_t = []
        B_t= []
        A = []
        B = []
        
        #for hydrogen bond
        A_t_opposite = []
        B_t_opposite  = []
        A_opposite  = []
        B_opposite  = []
        
        alpha_A_t = []
        alpha_B_t = []
        alpha_A = []
        alpha_B = []
        
          #for hydrogen bond
        alpha_A_t_opposite = []
        alpha_B_t_opposite = []
        alpha_A_opposite = []
        alpha_B_opposite = []
        
        pair_aa1_t = []
        pair_aa2_t = []
        pair_aa1 = []
        pair_aa2 = []
        
         
        #for hydrogen bond
        pair_aa1_t_opposite = []
        pair_aa2_t_opposite = []
        pair_aa1_opposite = []
        pair_aa2_opposite = []
        
        
        index_vec = i
        A_t = strand_vector[i] #copy strand to A array
        B_t = strand_vector[index_vec+1] #copy strand to B array 
        
        #for hydrogen bond
        A_t_opposite = strand_vector_opposite[i] #copy strand to A array
        B_t_opposite = strand_vector_opposite[index_vec+1] #copy strand to B array 
        
        
        #alpha_A_t = CA_vector_between_strand[i] #copy strand to A array
        #alpha_B_t = CA_vector_between_strand[index_vec+1] #copy strand to B array 
        
         #for hydrogen bond
        #alpha_A_t_opposite = CA_vector_between_strand_opposite[i] #copy strand to A array (alpha-carbon) backword direction
        #alpha_B_t_opposite = CA_vector_between_strand_opposite[index_vec+1] #copy strand to B array (alpha-carbon) backword direction
        
        
        alpha_A_t = center_vector[i] #copy strand to A array
        alpha_B_t = center_vector[index_vec+1] #copy strand to B array 
        
         #for hydrogen bond
        alpha_A_t_opposite = center_vector_opposite[i] #copy strand to A array (alpha-carbon) backword direction
        alpha_B_t_opposite = center_vector_opposite[index_vec+1] #copy strand to B array (alpha-carbon) backword direction
        
        
        pair_aa1_t = amino_acid[i] #copy strand to A array
        pair_aa2_t = amino_acid[index_vec+1] #copy strand to B array 
        
        #hydrogen bond
        pair_aa1_t_opposite = amino_acid_opposite[i] #copy strand to A array
        pair_aa2_t_opposite = amino_acid_opposite[index_vec+1] #copy strand to B array 
        
        
      #  print ("A_t ",A_t)
       
      #  print ("B_t ",B_t)
       # sys.stdout = original
        #print ("beta strand 1", beta_sheet[i]) #necessary
        #print ("beta strand 2", beta_sheet[i+1]) #necessary
        out_char = "beta strand 1" + str(beta_sheet[i]) + '\n'
        outf.write(out_char)
        out_char = "beta strand 2" + str(beta_sheet[i+1]) + '\n'
        outf.write(out_char)
        #print ("A_t ",len(A_t))
       
        #print ("B_t ",len(B_t))
        if len(A_t) > len(B_t):
            A = B_t
            B = A_t
            
            #hydrogen bond
            A_opposite = B_t_opposite
            B_opposite = A_t_opposite
            
            alpha_A = alpha_B_t
            alpha_B = alpha_A_t
            
             #hydrogen bond
            alpha_A_opposite = alpha_B_t_opposite
            alpha_B_opposite = alpha_A_t_opposite
            
            pair_aa1 = pair_aa2_t
            pair_aa2 = pair_aa1_t
            
            #hydrogen bond
            pair_aa1_opposite = pair_aa2_t_opposite
            pair_aa2_opposite = pair_aa1_t_opposite
        else:
            A = A_t
            B = B_t
            
             #hydrogen bond
            A_opposite = A_t_opposite
            B_opposite = B_t_opposite
            
            alpha_A = alpha_A_t
            alpha_B = alpha_B_t 
            
            #hydrogen bond
            alpha_A_opposite = alpha_A_t_opposite
            alpha_B_opposite = alpha_B_t_opposite
            
            pair_aa1 = pair_aa1_t
            pair_aa2 = pair_aa2_t
            
             #hydrogen bond
            pair_aa1_opposite = pair_aa1_t_opposite
            pair_aa2_opposite = pair_aa2_t_opposite
        
        #aa_min_index2 = 0
        k = 0
        twist_angle_return = 0.0
        add_1 = 0 
        add_n = 0
        add_hp = 0 #hydrogen bond
        pick_index = 0
        
        '''for calculation of avg twist angle between two beta strands'''
        twist_count = 0
        sum_twist_angle = 0.0
        avg_twist_angle = []
        
        #taking minimum of 4 consecutive  angles' sum
        cons_twist_angle = []
       
        
        for j in range (0, len(A_opposite)):
            a = numpy.array(A[j])  
            
            #hydrogen bond
            a_opposite = numpy.array(A_opposite[j])
            
            alpha_a = numpy.array(alpha_A[j])  #this is needed for calculating center of vector in 3D space
            
            #hydrogen bond
            alpha_a_opposite = numpy.array(alpha_A_opposite[j])
           
               
            #minimum_dist_point, aa_min_index2 = find_nearest_distance(a, B, pair_aa1[j], pair_aa2)
            #print ("k ", k)
            if k == 0  : # determining first pair
              
                minimum_dist_point, aa_min_index2 = find_nearest_distance(alpha_a, alpha_B, pair_aa1[j], pair_aa2)#nearest neighbour checking it will work for determining first pair
                #minimum_dist_point, aa_min_index2 = find_nearest_distance(a, B, pair_aa1[j], pair_aa2)
               #hydrogen bond
                minimum_dist_point_opposite, aa_min_index2_opposite = find_nearest_distance(alpha_a_opposite, alpha_B_opposite, pair_aa1_opposite[j], pair_aa2_opposite)#nearest neighbour checking it will work for determining first pair
                #minimum_dist_point_opposite, aa_min_index2_opposite = find_nearest_distance(a_opposite, B_opposite, pair_aa1_opposite[j], pair_aa2_opposite)
               
               # sys.stdout = original
                print(" MIN_INDEX_Forward", aa_min_index2)   
                print(" MIN_INDEX_Backword", aa_min_index2_opposite)    
                #print(" MIN_POINT ", minimum_dist_point)
                #print(" MIN_B ", B[aa_min_index2])
                
                '''
                twist_angle_return = calculate_twist_angle(a,B[aa_min_index2],j,aa_min_index2)
                sum_twist_angle = sum_twist_angle + twist_angle_return
                twist_count = twist_count + 1             
                pick_index = aa_min_index2
                '''
                
                if aa_min_index2 == 0 and len(A) <= aa_min_index2_opposite: #hydrogen bond parallel and anti-parallel case equal length strand or small strand is inside of big strand (complete overlap)
                    add_1 = 1
                    
                    twist_angle_return = calculate_twist_angle(a,B[aa_min_index2],j,aa_min_index2)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1             
                    pick_index = aa_min_index2
                    
                if aa_min_index2 == 0  and len(A) > aa_min_index2_opposite: #hydrogen bond parallel and anti-parallel case   small strand is outside (partial overlap) of big strand
                    add_hp = 1
                    
                    twist_angle_return = calculate_twist_angle(a_opposite,B_opposite[aa_min_index2_opposite],j,aa_min_index2)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1             
                    pick_index = aa_min_index2_opposite
                    
                if aa_min_index2 != 0:
                    add_n = 1
                    twist_angle_return = calculate_twist_angle(a,B[aa_min_index2],j,aa_min_index2)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1             
                    pick_index = aa_min_index2
                
            else: #next pairs
                print aa_min_index2," ",len(B), " ",aa_min_index2_opposite
                if add_1 == 1 and k < len(B):
                    minimum_dist_point_dup = B[k]
                    pick_index = k
                    twist_angle_return = calculate_twist_angle(a,minimum_dist_point_dup,j,pick_index)
                    cons_twist_angle.append(twist_angle_return)   
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1
                    
                if add_hp == 1 and  int(beta_sheet[i+1][11])!=-1  and  pick_index < len(B_opposite) - 1:		#parallel  small strand is outside (partial overlap) of big strand                         
                   # j = pick_index - 1
                    pick_index = pick_index + 1
                  
                    print "k ",k
                    print "pick_index",pick_index
                    minimum_dist_point_dup = B_opposite[pick_index]
                    #print "minimum_dist_point_dup ",minimum_dist_point_dup #ok
                   # print "a_opposite ",a_opposite
            
                    twist_angle_return = calculate_twist_angle(a_opposite,minimum_dist_point_dup,j,pick_index)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1
                    
                if add_hp == 1 and  int(beta_sheet[i+1][11])==-1  and  pick_index != 0: # anti-parallel  small strand is outside (partial overlap) of big strand                         
                   # j = pick_index - 1
                    pick_index = pick_index - 1 #may be need to change later 
                  
                    print "k ",k
                    print "pick_index",pick_index
                    minimum_dist_point_dup = B_opposite[pick_index]
                    #print "minimum_dist_point_dup ",minimum_dist_point_dup #ok
                   # print "a_opposite ",a_opposite
                    
                   
                    twist_angle_return = calculate_twist_angle(a_opposite,minimum_dist_point_dup,j,pick_index)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1
                
                
                if add_n == 1 and  int(beta_sheet[i+1][11])!=-1 and k < len(B) - aa_min_index2: # Parallel beta strands 	
                    pick_index = pick_index + 1
                    minimum_dist_point_dup = B[pick_index]
                    twist_angle_return = calculate_twist_angle(a,minimum_dist_point_dup,j,pick_index)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1
                    
                if add_n == 1 and   int(beta_sheet[i+1][11])==-1  and pick_index!=0: # anti-Parallel beta strands 		
                    pick_index = pick_index - 1           
                    minimum_dist_point_dup = B[pick_index]
                    twist_angle_return = calculate_twist_angle(a,minimum_dist_point_dup,j,pick_index)
                    cons_twist_angle.append(twist_angle_return)   
                    sum_twist_angle = sum_twist_angle + twist_angle_return
                    twist_count = twist_count + 1
            
                                                                   
            k = k + 1
                 
        if  twist_count != 0:  
            avg_twist_angle.append (float(sum_twist_angle) / float(twist_count))  
           
            cons_sum_array = []
            if len(cons_twist_angle) > 4:
                for c in range (0,len(cons_twist_angle) - 3):
                    cons_sum = 0.0
                    for cc in range (c,c+4):       
                         cons_sum = cons_twist_angle[cc] + cons_sum
                    cons_sum_array.append(cons_sum)
                pair_array.append((round(min(cons_sum_array)  / 4.0,3)))
                len_array.append(len(cons_sum_array))
                #print "sum of consecutive 4 angles ", cons_sum_array
                #print "min of consecutive 4 angles ",min(cons_sum_array)
                out_char = "min of sum of consecutive 4 angles :: " + str(round(min(cons_sum_array),3)) + '\n'
                outf.write(out_char)
                out_char = "Twist per angle :: " + str(round(min(cons_sum_array)  / 4.0 ,3)) + '\n'
                outf.write(out_char)
            
            elif len(cons_twist_angle) == 4:
                for c in range (0,len(cons_twist_angle) - 2):
                    cons_sum = 0.0
                    for cc in range (c,c+3):       
                         cons_sum = cons_twist_angle[cc] + cons_sum
                    cons_sum_array.append(cons_sum)
                pair_array.append((round(min(cons_sum_array)  / 3.0,3)))
                len_array.append(len(cons_sum_array))
                #print "sum of consecutive 4 angles ", cons_sum_array
                #print "min of consecutive 4 angles ",min(cons_sum_array)
                out_char = "min of sum of consecutive 3 angles :: " + str(round(min(cons_sum_array),3)) + '\n'
                outf.write(out_char)  
                out_char = "Twist per angle :: " + str(round(min(cons_sum_array) / 3.0 ,3) ) + '\n'
                outf.write(out_char)
                
            elif len(cons_twist_angle) == 2 or len(cons_twist_angle) == 3 or len(cons_twist_angle) == 1:
               #print "smallest twist ",min(cons_twist_angle)
                pair_array.append((round(min(cons_twist_angle) , 3 )))
                len_array.append(len(cons_sum_array))
                out_char = "smallest angles :: " + str(round(min(cons_twist_angle),3)) + '\n'
                outf.write(out_char)
            avg_track = 1
        else:
            avg_track = 0
        
        for t in range (len(avg_twist_angle)):   
            #print ("avg_twist_angle between strands ", round(avg_twist_angle[t],3)) #necessary
            #print ( round(avg_twist_angle[t],3)) #necessary
            out_char = "avg_twist_angle between strands  " + str( round(avg_twist_angle[t],3) ) + '\n'
            outf.write(out_char)
            all_avg_twist_angle.append(avg_twist_angle[t])            
     
        #print ("...........................")#necessary
        #print cons_twist_angle#ok
      
        outf.write("................................................"+"\n")
   
    i = i +1
#print ("avg_twist_angle between strands ", all_avg_twist_angle) #ok     


'''for calculation of avg twist angle in a beta sheet'''
outf.write("\n\n")

twist_angle_beta_sheet = []
avg_avg_twist_angle_beta_sheet = []
jc = 0   
for i in range (len(counter_index)):  
    sum_avg_twist_angle = 0.0  
    j = jc
    #print jc, " ", counter_index[i]-1, " ",j #ok
    while j < (jc + counter_index[i]-1):
        #print "j ",j
        sum_avg_twist_angle = sum_avg_twist_angle + all_avg_twist_angle[j]  
        sum_avg_twist_angle = round(sum_avg_twist_angle,3)       
        j = j + 1    
    jc = jc + counter_index[i]-1
    #avg_avg_twist_angle_beta_sheet.append(sum_avg_twist_angle)
    twist_angle_beta_sheet.append(float(sum_avg_twist_angle)/float(counter_index[i]-1))
    #print ("avg_twist angle in all beta sheets ", round(twist_angle_beta_sheet[i],3)) #necessary
  
    #print ( round(twist_angle_beta_sheet[i],3)) 
    out_char = "avg_twist angle in all beta sheets " + str( round(twist_angle_beta_sheet[i],3)) + '\n'
    outf.write(out_char)
   
#print ("avg_twist angle in all beta sheets ", float(avg_avg_twist_angle_beta_sheet[0])/float(counter_index[0]-1))
#print pair_array #ok

outf.write("................................................"+"\n")

#print pair_array # ok
#print len_array #ok
'''
#### According to minimum value pair #############

kc = 0
for i in range (len(counter_index)):  
    j = kc
    min_arr = []
    while j < (kc + counter_index[i]-1):
        min_arr.append (pair_array[j])
        
        j = j + 1  
    kc = kc + counter_index[i]-1
    arr = numpy.array (min_arr) 
    ind1 , ind2 = arr.argsort()[:2] 
    min_val1 =   arr[ind1]
    min_val2 =   arr[ind2]
    #print min_val1, min_val2
    out_char = "Min pair val1, val2 :: " + str(min_val1) + ",  "+ str(min_val2)+'\n'
    outf.write(out_char)
    
    min_avg_twist = float(min_val1 + min_val2)/2.0
    #print min_avg_twist
    out_char = "Min Twist :: " + str(min_avg_twist) + '\n'
    outf.write(out_char)

outf.write("................................................"+"\n")

#### According to longest pair #############
kl = 0
for il in range (len(counter_index)):  
    jl = kl
    min_arr_l = []
    while jl < (kl + counter_index[il]-1):
        min_arr_l.append (len_array[jl])
        
        jl = jl + 1  
    kl = kl + counter_index[il]-1
    arr_len = numpy.array (len_array) 
    ind1 , ind2 = arr_len.argsort()[::-1][:2]
    max_val1 =   pair_array[ind1]
    max_val2 =   pair_array[ind2]
   # print ind1, ind2 #ok
    out_char = "longest pair val1, val2 :: " + str(max_val1) + ",  "+ str(max_val2)+'\n'
    outf.write(out_char)
    
    min_avg_twist_long_pair = float(max_val1 + max_val2)/2.0
    #print min_avg_twist_long_pair
    out_char = "Avg MinTwist according to longest pair :: " + str(min_avg_twist_long_pair) + '\n'
    outf.write(out_char)

'''


#print (beta_sheet)
#print (beta_sheet[0][6])
#print (beta_sheet[0][9])
#print (atom[0])
#print (atom[0][5])
#print (a)
#f.close()

    