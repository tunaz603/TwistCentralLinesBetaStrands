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
            
'''            
if len(sys.argv) == 1 or len(sys.argv) == 2 or len(sys.argv) > 3  :
    print "Usage: python twist_region.py pdbfilename.pdb direction"
    sys.exit(0)


lines_File = [lines1.rstrip('\n') for lines1 in open(sys.argv[1], 'r')]   


try:
    direction = int(sys.argv[2]) # give direction as input
 
    if direction != 0 and direction != 1 :
        print "please input either 0 or 1. 0 = forward, 1 = backword"
        sys.exit(0)
except ValueError:  
        print "please input either 0 or 1. 0 = forward, 1 = backword"
        sys.exit(0)          
'''
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
lines_File0 = [lines1.rstrip('\n') for lines1 in open('s1.pdb_data_vl_1.1.pdb')]
lines_File1 = [lines1.rstrip('\n') for lines1 in open('s2.pdb_data_vl_1.1.pdb')]
lines_File2 = [lines1.rstrip('\n') for lines1 in open('s3.pdb_data_vl_1.1.pdb')]
lines_File3 = [lines1.rstrip('\n') for lines1 in open('s4.pdb_data_vl_1.1.pdb')]
lines_File4 = [lines1.rstrip('\n') for lines1 in open('s5.pdb_data_vl_1.1.pdb')]
lines_File5 = [lines1.rstrip('\n') for lines1 in open('s6.pdb_data_vl_1.1.pdb')]
lines_File6 = [lines1.rstrip('\n') for lines1 in open('s7.pdb_data_vl_1.1.pdb')]


direction = 0 # direction= 0 means forward direction, direction = 1 means backword direction so change here        
            

direction_opposite = 1 #hydrogen bond 

count_beta_sheet = 0
count_beta_strand = 0
counter = 0
counter_index = []
problem = 0
a=[]

atom0 = []
atom1 = []
atom2 = []
atom3 = []
atom4 = []
atom5 = []
atom6 = []

count_atom = 0
for sheet_line in lines_File0: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom0.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         
         
         #atom.append( [sheet_line[0:6], sheet_line[6:11], sheet_line[12:16], sheet_line[16:17],sheet_line[17:20], sheet_line[21:22], sheet_line[22:26], sheet_line[26:27], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54], sheet_line[46:54], sheet_line[54:60], sheet_line[60:66], sheet_line[76:78], sheet_line[78:80]])
         #atom.append ( re.sub("\s\s+" , " ",sheet_line).split())
#print (int(atom[1836][5].strip()))
#print (atom[1836])

for sheet_line in lines_File1: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom1.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         
    
for sheet_line in lines_File2: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom2.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         
    
for sheet_line in lines_File3: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom3.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         

for sheet_line in lines_File4: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom4.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         
    
for sheet_line in lines_File5: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom5.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         

for sheet_line in lines_File6: #find each line in all lines
   if sheet_line.startswith('ATOM'):  #find those lines start with "ATOM"
         count_atom = count_atom+1
         #print (sheet_line)       #print the lines started with "ATOM"
         atom6.append( [sheet_line[:6], sheet_line[6:11], sheet_line[12:16], sheet_line[17:20], sheet_line[21], sheet_line[22:26], sheet_line[30:38], sheet_line[38:46], sheet_line[46:54]])
         
    
    

strand_vector=[]
center_vector = []
CA_vector_between_strand = []
amino_acid = []




poi = 0
track = 1

######## 1st
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom0)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom0[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom0[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom0[i][8].strip()))    #making midpoint 
   
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef


length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
vector_array =[]
cen_array = []
midpoint_index = 1
while midpoint_index < length_midpoint_array:
    point = []
    cen_point = []
    temp_mid_index = midpoint_index - 1
    #print ("mid x ", midpoint_x [temp_mid_index] , " ", midpoint_x [temp_mid_index+1]) 
    #print ("mid y ",midpoint_y [temp_mid_index] , " ", midpoint_y [temp_mid_index+1]) 
   # print ("mid z ",midpoint_z [temp_mid_index] , " ", midpoint_z [temp_mid_index+1]) 
    vec_x =  midpoint_x [temp_mid_index] -  midpoint_x [temp_mid_index+1]  #making midpoint vector
    vec_y =  midpoint_y [temp_mid_index] -  midpoint_y [temp_mid_index+1]   #making midpoint vector
    vec_z =  midpoint_z [temp_mid_index] -  midpoint_z [temp_mid_index+1] #making midpoint vector
    
    cen_x = (float(midpoint_x [temp_mid_index]) + float( midpoint_x [temp_mid_index+1])) / float(2.0)  #making center of midpoint vector
    cen_y =(float( midpoint_y [temp_mid_index]) + float(midpoint_y [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
    cen_z = (float(midpoint_z [temp_mid_index]) + float( midpoint_z [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
    
   # print ("vector x ",vec_x) 
   # print ("vector y ",vec_y) 
  #  print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag1 ",magnitude)
   # if magnitude!= 0:
    point.append(vec_x/round(magnitude,13)) # unit vector
    point.append(vec_y/round(magnitude,13))  # unit vector
    point.append(vec_z/round(magnitude,13))  # unit vector
    
  
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
        
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)

#print "center",center_vector
######## 2nd
    
    
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom1)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom1[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom1[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom1[i][8].strip()))    #making midpoint 
   
     
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef
length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
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
    #print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag2 ",magnitude)
    
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

    
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)

######## 3rd
    
    
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom2)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom2[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom2[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom2[i][8].strip()))    #making midpoint 
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef
     


length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
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
    #print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag3 ",magnitude)
    
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

    
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)


######## 4th
    
    
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom3)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom3[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom3[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom3[i][8].strip()))    #making midpoint 
   
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef
  

length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
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
    #print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag4 ",magnitude)
    
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

    
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)


######## 5th
    
    
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom4)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom4[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom4[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom4[i][8].strip()))    #making midpoint 
   
     
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef

length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
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
    #print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag5 ",magnitude)
    
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

    
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)

######## 6th
    
    
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom5)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom5[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom5[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom5[i][8].strip()))    #making midpoint 
   
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef

length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
vector_array =[]
cen_array = []
midpoint_index = 1
while midpoint_index < length_midpoint_array:
    point = []
    cen_point = []
    temp_mid_index = midpoint_index - 1
    
  #  print ("mid x ", midpoint_x [temp_mid_index] , " ", midpoint_x [temp_mid_index+1]) 
  #  print ("mid y ",midpoint_y [temp_mid_index] , " ", midpoint_y [temp_mid_index+1]) 
  #  print ("mid z ",midpoint_z [temp_mid_index] , " ", midpoint_z [temp_mid_index+1]) 
    vec_x = midpoint_x [temp_mid_index] -  midpoint_x [temp_mid_index+1]  #making midpoint vector
    vec_y = midpoint_y [temp_mid_index] -  midpoint_y [temp_mid_index+1]  #making midpoint vector
    vec_z = midpoint_z [temp_mid_index] -  midpoint_z [temp_mid_index+1]  #making midpoint vector
    
    cen_x = (float(midpoint_x [temp_mid_index]) + float( midpoint_x [temp_mid_index+1])) / float(2.0)  #making center of midpoint vector
    cen_y =(float( midpoint_y [temp_mid_index]) + float(midpoint_y [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
    cen_z = (float(midpoint_z [temp_mid_index]) + float( midpoint_z [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
    
    #print ("vector x ",vec_x) 
   # print ("vector y ",vec_y) 
    #print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag6 ",magnitude)
    
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

    
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)

   #print(beta_sheet[sheet_index][1])  
    
    

 ######## 7th
    
    
midpoint_x_bef = []
midpoint_y_bef = []
midpoint_z_bef = []  

midpoint_x = []
midpoint_y = []
midpoint_z = []

new_x = []
new_y = []
new_z = []  
        
for i in range (0,len(atom6)): # 5 should be number of beta stand in one beta sheet
    #print("love ",beta_sheet[sheet_index][10])
   
   

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
    
    
    if direction == 0: #forward direction #ok
       
        
              
        midpoint_x_bef.append(float(atom6[i][6].strip()))    #making midpoint       
        midpoint_y_bef.append(float(atom6[i][7].strip()))     #making midpoint 
        midpoint_z_bef.append(float(atom6[i][8].strip()))    #making midpoint 
   
        midpoint_x=midpoint_x_bef
        midpoint_y=midpoint_y_bef
        midpoint_z=midpoint_z_bef

length_midpoint_array = len(midpoint_x)
#print "length_midpoint_array ",length_midpoint_array# ok
vector_array =[]
cen_array = []
midpoint_index = 1
while midpoint_index < length_midpoint_array:
    point = []
    cen_point = []
    temp_mid_index = midpoint_index - 1
    print ("mid x ", midpoint_x [temp_mid_index] , " ", midpoint_x [temp_mid_index+1]) 
    print ("mid y ",midpoint_y [temp_mid_index] , " ", midpoint_y [temp_mid_index+1]) 
    print ("mid z ",midpoint_z [temp_mid_index] , " ", midpoint_z [temp_mid_index+1]) 
    vec_x = midpoint_x [temp_mid_index] -  midpoint_x [temp_mid_index+1]  #making midpoint vector
    vec_y = midpoint_y [temp_mid_index] -  midpoint_y [temp_mid_index+1]  #making midpoint vector
    vec_z = midpoint_z [temp_mid_index] -  midpoint_z [temp_mid_index+1]  #making midpoint vector
    
    cen_x = (float(midpoint_x [temp_mid_index]) + float( midpoint_x [temp_mid_index+1])) / float(2.0)  #making center of midpoint vector
    cen_y =(float( midpoint_y [temp_mid_index]) + float(midpoint_y [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
    cen_z = (float(midpoint_z [temp_mid_index]) + float( midpoint_z [temp_mid_index+1]))/float(2.0)   #making center of midpoint vector
    
    #print ("vector x ",vec_x) 
   # print ("vector y ",vec_y) 
    #print ("vector z ",vec_z) 
    #make unit vector
    magnitude = math.sqrt(pow(vec_x, 2) + pow(vec_y, 2) + pow(vec_z, 2))
    
   
    print("mag7 ",magnitude)
    
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

    
    
strand_vector.append(vector_array)  

center_vector.append(cen_array)


    
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
#print (center_vector)
#print len(center_vector)
#print(amino_acid) # ok
#print ("vector strand length ",length_strand_vector)
#print ("CA_vector_between_strand ",len(CA_vector_between_strand)) # ok.  vector_strand length and CA_vector same
#print ("CA_vector_between_strand ",CA_vector_between_strand) 
#draw_interpolation (CA_vector_between_strand)

def distance(x,y):   
    return numpy.sqrt(numpy.sum((x-y)**2))


def find_nearest_distance(a,B):
    
    
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
    #print  ("dot ",dot_product)
       
    a_sum = pow(a[0],2)+pow(a[1],2)+pow(a[2],2)
    a_root=math.sqrt(a_sum)
    b_sum = pow(minimum_dist_point[0],2)+pow(minimum_dist_point[1],2)+pow(minimum_dist_point[2],2)
    b_root=math.sqrt(b_sum)
    a_b = a_root*b_root
    #print "ab ",a_b
    angle_dot = math.acos(round(dot_product/a_b , 13))
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
while i < length_strand_vector-1  : # pairwise calculation 

    A_t = []
    B_t= []
    A = []
    B = []
    
    alpha_A_t = []
    alpha_B_t = []
    alpha_A = []
    alpha_B = []
    
    pair_aa1_t = []
    pair_aa2_t = []
    pair_aa1 = []
    pair_aa2 = []
    index_vec = i
    A_t = strand_vector[i] #copy strand to A array
    B_t = strand_vector[index_vec+1] #copy strand to B array 
    
    #alpha_A_t = CA_vector_between_strand[i] #copy strand to A array
    #alpha_B_t = CA_vector_between_strand[index_vec+1] #copy strand to B array 
    
    alpha_A_t = center_vector[i] #copy strand to A array
    alpha_B_t = center_vector[index_vec+1] #copy strand to B array 
   
   # print ("A_t ",A_t)
   
   # print ("B_t ",B_t)
   # sys.stdout = original
    #print ("beta strand 1", beta_sheet[i]) #necessary
    #print ("beta strand 2", beta_sheet[i+1]) #necessary
 
    #print ("A_t ",len(A_t))
   
    #print ("B_t ",len(B_t))
    if len(A_t) > len(B_t):
        A = B_t
        B = A_t
        
        alpha_A = alpha_B_t
        alpha_B = alpha_A_t
        
     
    else:
        A = A_t
        B = B_t
        
        alpha_A = alpha_A_t
        alpha_B = alpha_B_t 
        
     
   
    twist_angle_return = 0.0
    add_1 = 0 
    add_n = 0
    pick_index = 0
    
    '''for calculation of avg twist angle between two beta strands'''
    twist_count = 0
    sum_twist_angle = 0.0
    avg_twist_angle = []
    
     #taking minimum of 4 consecutive  angles' sum
    cons_twist_angle = []
   
    #nearest neighbour way
    a_start = numpy.array(alpha_A[0])  
    a_end = numpy.array(alpha_A[len(A) - 1])
    p1, aa_min_start = find_nearest_distance(a_start, alpha_B)
    p2, aa_min_end = find_nearest_distance(a_end, alpha_B)
    print "start ","end ",aa_min_start,' ',aa_min_end,' ', len (alpha_B)
    if aa_min_start < aa_min_end :
        k = aa_min_start
    else:
        k = aa_min_end
    for j in range (0, len(A)):
          if (k < len (B)):
                a = numpy.array(A[j])  
                alpha_a = numpy.array(alpha_A[j])  #this is needed for calculating center of vector in 3D space
                minimum_dist_point, aa_min_index2 = find_nearest_distance(alpha_a, alpha_B)
                twist_angle_return = calculate_twist_angle(a,B[aa_min_index2],j,aa_min_index2)
                cons_twist_angle.append(twist_angle_return)     
                sum_twist_angle = sum_twist_angle + twist_angle_return
                twist_count = twist_count + 1  
                k = k+1
        
        
    #print cons_twist_angle#ok
    cons_sum_array = []
    if len(cons_twist_angle) > 4:
        for c in range (0,len(cons_twist_angle) - 3):
            cons_sum = 0.0
            for cc in range (c,c+4):       
                 cons_sum = cons_twist_angle[cc] + cons_sum
            cons_sum_array.append(cons_sum)
        pair_array.append((round(min(cons_sum_array)  / 4.0,3)))
        len_array.append(len(cons_twist_angle))
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
        len_array.append(len(cons_twist_angle))
        #print "sum of consecutive 4 angles ", cons_sum_array
        #print "min of consecutive 4 angles ",min(cons_sum_array)
        out_char = "min of sum of consecutive 3 angles :: " + str(round(min(cons_sum_array),3)) + '\n'
        outf.write(out_char)  
        out_char = "Twist per angle :: " + str(round(min(cons_sum_array) / 3.0 ,3) ) + '\n'
        outf.write(out_char)
    
        
    elif len(cons_twist_angle) == 2 or len(cons_twist_angle) == 3 :
       #print "smallest twist ",min(cons_twist_angle)
        pair_array.append((round(min(cons_twist_angle) , 3 )))
        len_array.append(len(cons_twist_angle))
        out_char = "smallest angles :: " + str(round(min(cons_twist_angle),3)) + '\n'
        outf.write(out_char)
    
            
                   
              
 
    #print ("...........................")#necessary
    outf.write("................................................"+"\n")
   
    i = i + 1
#print ("avg_twist_angle between strands ", all_avg_twist_angle) #ok     

print pair_array # ok
print len_array #ok

#### According to minimum value pair #############
arr = numpy.array (pair_array) 
ind1 , ind2 = arr.argsort()[:2] 
min_val1 =   arr[ind1]
min_val2 =   arr[ind2]
#print min_val1, min_val2
out_char = "Min pair val1, val2 :: " + str(min_val1) + ",  "+ str(min_val2)+'\n'
outf.write(out_char)

min_avg_twist = float(min_val1 + min_val2)/2.0
#print min_avg_twist
out_char = "Avg Min Twist :: " + str(min_avg_twist) + '\n'
outf.write(out_char)

#### According to longest pair #############
arr_len = numpy.array (len_array) 
ind1 , ind2 = arr_len.argsort()[::-1][:2]
max_val1 =   pair_array[ind1]
max_val2 =   pair_array[ind2]
print ind1, ind2 #ok
out_char = "longest pair val1, val2 :: " + str(max_val1) + ",  "+ str(max_val2)+'\n'
outf.write(out_char)

min_avg_twist_long_pair = float(max_val1 + max_val2)/2.0
#print min_avg_twist_long_pair
out_char = "Avg MinTwist according to longest pair :: " + str(min_avg_twist_long_pair) + '\n'
outf.write(out_char)



