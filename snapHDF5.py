# Python HDF5 snapshot reader/writer
# (requirest PyTables http://www.pytables.org)
#
# import snapHDF5 as snap 
# header = snap.snapshot_header("snap_063") 
# mass = snap.read_block("snap_063", "MASS", parttype=5) 
#
# Mark Vogelsberger (mvogelsb@cfa.harvard.edu)

import numpy as np
import os
import sys
import math
import tables

############ 
#DATABLOCKS#
############
#descriptions of all datablocks -> add new datablocks here!
#format -> TAG:[HDF5_NAME, DIMENSION]
datablocks = {"POS ":["Coordinates",3], 
	      "VEL ":["Velocities",3],
	      "ID  ":["ParticleIDs",1],
	      "MASS":["Masses",1],
              "U   ":["InternalEnergy",1],
              "RHO ":["Density",1],
              "VOL ":["Volume",1],
              "CMCE":["Center-of-Mass",3],
              "AREA":["Surface Area",1],
              "NFAC":["Number of faces of cell",1],
              "NE  ":["ElectronAbundance",1],
              "NH  ":["NeutralHydrogenAbundance",1],
	      "HSML":["SmoothingLength",1],
              "SFR ":["StarFormationRate",1],
              "AGE ":["StellarFormationTime",1],
              "Z   ":["Metallicity",1],
	      "ACCE":["Acceleration",3],
              "VEVE":["VertexVelocity",3],
              "FACA":["MaxFaceAngle",1],              
	      "COOR":["CoolingRate",1],
              "POT ":["Potential",1],
	      "MACH":["MachNumber",1],
	      "DMHS":["DM Hsml",1],
	      "DMDE":["DM Density",1],
              "PHKE":["PHKey",1],
              "GROU":["GroupNr",1],
#SIDM (GADGET)
              "PTSU":["PSum",1],
              "DMNB":["DMNumNgb",1],
	      "NTSC":["NumTotalScatter",1],
              "SHSM":["SIDMHsml",1],
              "SRHO":["SIDMRho",1],
              "SVEL":["SVelDisp",1],
#SIDM (AREPO)
              "PTSU":["SIDM_Psum",1],
              "DMNB":["SIDM_NumNgb",1],
              "NTSC":["SIDM_NumTotalScatter",1],
              "SHSM":["SIDM_Hsml",1],
              "SRHO":["SIDM_Density",1],
              "SVEL":["SIDM_VelDisp",1],
#TRACER
              "TRFQ":["FluidQuantities", 3],
	      "TRNT":["NumTracers", 1],
	      "TRCE":["TracerField", 1],
#GFM 
	      "GAGE":["GFM_StellarFormationTime",1],
              "GIMA":["GFM_InitialMass",1],
              "GZ  ":["GFM_Metallicity",1],
              "GMET":["GFM_Metals",9],
              "GWHV":["GFM_WindHostVal",1],
              "GCOL":["GFM_CoolingRate",1],
              "GSPH":["GFM_StellarPhotometrics",8], #band luminosities: U, B, V, K, g, r, i, z
              "AGNR":["GFM_AGNRadiation",1],

#FM_RadFeedback
              "STRM":["StromgrenRadius",1],
              "TAUD":["RadFeedTau",1],
              "RFNN":["RadFeed_NumNgb",1],
              "RADE":["RadiationMomentumReleased",1],
              "RFNS":["NormSphRadFeedback",1],
              "RFCD":["RadCoolShutoffTime",1],
              "RFCG":["GasRadCoolShutoffTime",1],
              "VSNZ":["AvgMetalNgb",1],             #average metallicity of neighboring gas cells used in variable SN model
              "NGBM":["TotNgbMass",1],              #total mass of neighboring gas cells 
#GDE
              "CACO":["CausticCounter",1],
              "STDE":["StreamDensity",1],
              "PSDE":["PhaseSpaceDensity",1],
              "FLDE":["FlowDeterminant",1],
	      "TIPS":["TidalTensorPS",9],
              "DIPS":["DistortionTensorPS",36],
              "SHOR":["SheetOrientation", 9],
              "INDE":["InitDensity",1],
#ISM
              "BRD ":["BlastRadius", 1],
              "CSO ":["CoolShutoffTime", 1],
              "RadP":["RadiationPressureMoment",1], 
#SUBFIND
              "SFDE":["SubfindDensity", 1],
              "SFHS":["SubfindHsml", 1],
              "SFVD":["SubfindVelDisp", 1],
#BHs
              "BHMD":["BH_Mdot", 1],
              "BHHS":["BH_Hsml", 1],
              "BHMA":["BH_Mass", 1],
              "REBH":["RefBHCounter", 1],
	      "BHHG":["BH_HaloGasMass", 1],
              "BHCL":["BH_CoolingLuminosity", 1],
              "BHMR":["BH_Mdot_Radio", 1], 
              "BHPR":["BH_Progs", 1]}


#####################################################################################################################
#                                                    READING ROUTINES			                            #
#####################################################################################################################


########################### 
#CLASS FOR SNAPSHOT HEADER#
###########################  
class snapshot_header:
	def __init__(self, *args, **kwargs):
		if (len(args) == 1):
			filename = args[0]
 			if os.path.exists(filename):
			        curfilename=filename
			elif os.path.exists(filename+".hdf5"):
				curfilename = filename+".hdf5"
			elif os.path.exists(filename+".0.hdf5"): 
				curfilename = filename+".0.hdf5"
			else:	
				print "[error] file not found : ", filename
				sys.stdout.flush()
			    	sys.exit()
    
			f=tables.openFile(curfilename)
			self.npart = f.root.Header._v_attrs.NumPart_ThisFile 
			self.nall = f.root.Header._v_attrs.NumPart_Total
			self.nall_highword = f.root.Header._v_attrs.NumPart_Total_HighWord
			self.massarr = f.root.Header._v_attrs.MassTable 
			self.time = f.root.Header._v_attrs.Time 
			self.redshift = f.root.Header._v_attrs.Redshift 
			self.boxsize = f.root.Header._v_attrs.BoxSize
			self.filenum = f.root.Header._v_attrs.NumFilesPerSnapshot
			self.omega0 = f.root.Header._v_attrs.Omega0
			self.omegaL = f.root.Header._v_attrs.OmegaLambda
			self.hubble = f.root.Header._v_attrs.HubbleParam
			self.sfr = f.root.Header._v_attrs.Flag_Sfr 
			self.cooling = f.root.Header._v_attrs.Flag_Cooling
			self.stellar_age = f.root.Header._v_attrs.Flag_StellarAge
			self.metals = f.root.Header._v_attrs.Flag_Metals
			self.feedback = f.root.Header._v_attrs.Flag_Feedback
			self.double = f.root.Header._v_attrs.Flag_DoublePrecision #GADGET-2 change
			f.close()
		else:
			#read arguments
			self.npart = kwargs.get("npart")
			self.nall = kwargs.get("nall")
			self.nall_highword = kwargs.get("nall_highword")
			self.massarr = kwargs.get("massarr")
			self.time = kwargs.get("time")
			self.redshift = kwargs.get("redshift")
			self.boxsize = kwargs.get("boxsize")
			self.filenum = kwargs.get("filenum")
			self.omega0 = kwargs.get("omega0")
			self.omegaL = kwargs.get("omegaL")
			self.hubble = kwargs.get("hubble")
			self.sfr = kwargs.get("sfr")
			self.cooling = kwargs.get("cooling")
			self.stellar_age = kwargs.get("stellar_age")
			self.metals = kwargs.get("metals")
			self.feedback = kwargs.get("feedback")
			self.double = kwargs.get("double")

			#set default values
                        if (self.npart == None):
				self.npart = np.array([0,0,0,0,0,0], dtype="int32")
                        if (self.nall == None):				
				self.nall  = np.array([0,0,0,0,0,0], dtype="uint32")
                        if (self.nall_highword == None):				
				self.nall_highword = np.array([0,0,0,0,0,0], dtype="uint32")
                        if (self.massarr == None):
				self.massarr = np.array([0,0,0,0,0,0], dtype="float64")
                        if (self.time == None):				
				self.time = np.array([0], dtype="float64")
                        if (self.redshift == None):				
				self.redshift = np.array([0], dtype="float64")
                        if (self.boxsize == None):				
				self.boxsize = np.array([0], dtype="float64")
                        if (self.filenum == None):
				self.filenum = np.array([1], dtype="int32")
                        if (self.omega0 == None):
				self.omega0 = np.array([0], dtype="float64")
                        if (self.omegaL == None):
				self.omegaL = np.array([0], dtype="float64")
                        if (self.hubble == None):
				self.hubble = np.array([0], dtype="float64")
                        if (self.sfr == None):	
				self.sfr = np.array([0], dtype="int32")            
                        if (self.cooling == None):	
				self.cooling = np.array([0], dtype="int32")
                        if (self.stellar_age == None):	
				self.stellar_age = np.array([0], dtype="int32")
                        if (self.metals == None):	
				self.metals = np.array([0], dtype="int32")
                        if (self.feedback == None):	
				self.feedback = np.array([0], dtype="int32")
                        if (self.double == None):
				self.double = np.array([0], dtype="int32")
		
			
		

##############################
#READ ROUTINE FOR SINGLE FILE#
############################## 
def read_block_single_file(filename, block_name, dim2, parttype=-1, no_mass_replicate=False, fill_block_name="", slab_start=-1, slab_len=-1, verbose=False):

  if (verbose):
	  print "[single] reading file           : ", filename   	
	  print "[single] reading                : ", block_name
     	  sys.stdout.flush()
 
  head = snapshot_header(filename)
  npart = head.npart
  massarr = head.massarr
  nall = head.nall
  filenum = head.filenum
  doubleflag = head.double #GADGET-2 change
  #doubleflag = 0          #GADGET-2 change


  if (parttype!=-1):
	if (head.npart[parttype]==0):
		return [0, False]
  else:
	if (head.npart.sum()==0):
		return [0, False]
  del head


  #construct indices for partial access
  if (slab_start!=-1) & (slab_len!=-1):
  	data_slice = slice(slab_start, (slab_start+slab_len))
  else:
	data_slice = slice(None, None, 1)

  if (verbose):
	print "[single] data slice: ", data_slice
	sys.stdout.flush() 

  f=tables.openFile(filename)


  #read specific particle type (parttype>=0, non-default)
  if parttype>=0:
        if (verbose):
        	print "[single] parttype               : ", parttype 
		sys.stdout.flush()
	if ((block_name=="Masses") & (npart[parttype]>0) & (massarr[parttype]>0)):
	        if (verbose):
			print "[single] replicate mass block"
			sys.stdout.flush()	
		ret_val=np.repeat(massarr[parttype], npart[parttype])[data_slice]
        else:		
	  	part_name='PartType'+str(parttype)
	  	ret_val = f.root._f_getChild(part_name)._f_getChild(block_name)[data_slice]
        if (verbose):
        	print "[single] read particles (total) : ", ret_val.shape[0]/dim2
		sys.stdout.flush()

  #read all particle types (parttype=-1, default)
  if parttype==-1:
	first=True
	dim1=long(0)
	for parttype in range(0,6):
		part_name='PartType'+str(parttype)
		if (f.root.__contains__(part_name)):
			if (verbose):
				print "[single] parttype               : ", parttype
				print "[single] massarr                : ", massarr
				print "[single] npart                  : ", npart
				sys.stdout.flush()

		        #replicate mass block per default (unless no_mass_replicate is set)
	        	if ((block_name=="Masses") & (npart[parttype]>0) & (massarr[parttype]>0) & (no_mass_replicate==False)):
                       		if (verbose):
                               		print "[single] replicate mass block"
					sys.stdout.flush()
				if (first):
			        	data=np.repeat(massarr[parttype], npart[parttype])
					dim1+=long(data.shape[0])
					ret_val=data
					first=False
				else:	
					data=np.repeat(massarr[parttype], npart[parttype])
                                        dim1+=long(data.shape[0])
                                        ret_val=np.append(ret_val, data)
        	                if (verbose):
                	        	print "[single] read particles (total) : ", ret_val.shape[0]/dim2
					sys.stdout.flush()
                                if (doubleflag==0):
					ret_val=ret_val.astype("float32")


			#fill fill_block_name with zeros if fill_block_name is set and particle type is present and fill_block_name not already stored in file for that particle type
                        if ((block_name==fill_block_name) & (block_name!="Masses") & (npart[parttype]>0) & (f.root._f_getChild(part_name).__contains__(block_name)==False)):
                                if (verbose):
                                        print "[single] replicate block name", fill_block_name
					sys.stdout.flush()
                                if (first):
                                        data=np.repeat(0.0, npart[parttype]*dim2)
                                        dim1+=long(data.shape[0])
                                        ret_val=data
                                        first=False
                                else:
                                        data=np.repeat(0.0, npart[parttype]*dim2)
                                        dim1+=long(data.shape[0])
                                        ret_val=np.append(ret_val, data)
                                if (verbose):
                                        print "[single] read particles (total) : ", ret_val.shape[0]/dim2
					sys.stdout.flush()
                                if (doubleflag==0):
                                        ret_val=ret_val.astype("float32")

			#default: just read the block
			if (f.root._f_getChild(part_name).__contains__(block_name)):
				if (first):
					data=f.root._f_getChild(part_name)._f_getChild(block_name)[:]
					dim1+=long(data.shape[0])
					ret_val=data
					first=False
				else:
					data=f.root._f_getChild(part_name)._f_getChild(block_name)[:]
					dim1+=long(data.shape[0])
					ret_val=np.append(ret_val, data)
                		if (verbose):
                        		print "[single] read particles (total) : ", ret_val.shape[0]/dim2
					sys.stdout.flush()

	if ((dim1>0) & (dim2>1)):
		ret_val=ret_val.reshape(dim1,dim2)
                if (verbose):
                        print "[single] reshape done: ", ret_val.shape
			sys.stdout.flush()


  f.close()


  return [ret_val, True]

##############
#READ ROUTINE#
##############
def read_block(filename, block, parttype=-1, no_mass_replicate=False, fill_block="", slab_start=-1, slab_len=-1, verbose=False):
  if (verbose):
          print "reading block          : ", block
	  sys.stdout.flush()	

  if parttype not in [-1,0,1,2,3,4,5]:
    print "[error] wrong parttype given"
    sys.stdout.flush()
    sys.exit()

  curfilename=filename+".hdf5"

  if os.path.exists(curfilename):
    multiple_files=False
  elif os.path.exists(filename+".0"+".hdf5"):
    curfilename = filename+".0"+".hdf5"
    multiple_files=True
  else:
    print "[error] file not found : ", filename
    sys.stdout.flush()
    sys.exit()

  slabflag=False
  if ((slab_start!=-1) | (slab_len!=-1)):
	slabflag=True
	if (parttype==-1):
		print "[error] slabs only supported for specific parttype"
		sys.stdout.flush()
		sys.exit()

  head = snapshot_header(curfilename)
  filenum = head.filenum
  del head

 
  if (datablocks.has_key(block)):
        block_name=datablocks[block][0]
        dim2=datablocks[block][1]
        first=True
        if (verbose):
                print "Reading HDF5           : ", block_name
                print "Data dimension         : ", dim2
		print "Multiple file          : ", multiple_files
		print "Slab data              : ", slabflag
		sys.stdout.flush()
  else:
        print "[error] Block type ", block, "not known!"
	sys.stdout.flush()
        sys.exit()

  fill_block_name=""
  if (fill_block!=""):
        if (datablocks.has_key(fill_block)):
                fill_block_name=datablocks[fill_block][0]
                dim2=datablocks[fill_block][1]
                if (verbose):
                        print "Block filling active   : ", fill_block_name
			sys.stdout.flush()

  if (multiple_files):	
	if (slabflag==False):
		first=True
		dim1=long(0)
		for num in range(0,filenum):
			curfilename=filename+"."+str(num)+".hdf5"
			if (verbose):
				print "Reading file           : ", num, curfilename
				sys.stdout.flush() 
			if (first):
				data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start, slab_len, verbose)
				if (succ == True):
					dim1+=long(data.shape[0])
					ret_val = data
					first = False 
			else:	 
        	                data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start, slab_len, verbose)
				if (succ == True):
	                	        dim1+=long(data.shape[0])
					ret_val=np.append(ret_val, data)
        	        if (verbose):
				if (succ):
	                	        print "Read particles (total) : ", ret_val.shape[0]/dim2
					sys.stdout.flush()
				else:
					print "Read particles (total) : none"
					sys.stdout.flush()
			

	if (slabflag==True):
                off=slab_start
		left=slab_len
                first=True
                dim1=long(0)
                for num in range(0,filenum):
                        curfilename=filename+"."+str(num)+".hdf5"
			head = snapshot_header(curfilename)
			nloc=head.npart[parttype]
	                if (nloc > off):
        	                start = off
                	        if (nloc - off > left):
                        	        count = left
	                        else:
        	                        count = nloc - off
	                        if (verbose):
        	                        print "Reading file           : ", num, curfilename, start, count
					sys.stdout.flush()
                	        if (first):
                        	        data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start=start, slab_len=count, verbose=verbose)
					if (succ == True):
		                                dim1+=data.shape[0]
        		                        ret_val = data
                		                first = False
	                        else:
        	                        data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start=start, slab_len=count, verbose=verbose)
					if (succ == True):
	                	                dim1+=data.shape[0]
        	                	        ret_val=np.append(ret_val, data)
	                        if (verbose):
					if (succ):
	        	                        print "Read particles (total) : ", ret_val.shape[0]/dim2
						sys.stdout.flush()
					else:
						print "Read particles (total) : none"
						sys.stdout.flush()

                	        left -= count
	                        off += count
       	        	if (left==0):
               	        	break
	                off -= nloc



	if (verbose):
		print "all partial files read in"
		sys.stdout.flush()

	if ((dim1>0) & (dim2>1)):
		ret_val=ret_val.reshape(dim1,dim2)
                if (verbose):
                        print "Reshape done (total): ", ret_val.shape
			sys.stdout.flush()


  else:	
	ret_val, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start, slab_len, verbose)

  return ret_val


#############
#LIST BLOCKS#
#############
def list_blocks(filename, parttype=-1, verbose=False):
  
  f=tables.openFile(filename)
  for parttype in range(0,6):
  	part_name='PartType'+str(parttype)
        if (f.root.__contains__(part_name)):
        	print "Parttype contains : ", parttype
		print "-------------------"
		sys.stdout.flush()
		iter = it=datablocks.__iter__()
		next = iter.next()
		while (1):
			if (verbose):
				print "check ", next, datablocks[next][0]
				sys.stdout.flush()
			if (f.root._f_getChild(part_name).__contains__(datablocks[next][0])):
  				print next, datablocks[next][0]
				sys.stdout.flush()
			try:
				next=iter.next()
			except StopIteration:
				break	
  f.close() 

#################
#CONTAINS BLOCKS#
#################
def contains_block(filename, tag, parttype=-1, verbose=False):
  
  contains_flag=False
  f=tables.openFile(filename)
  for parttype in range(0,6):
        part_name='PartType'+str(parttype)
        if (f.root.__contains__(part_name)):
                iter = it=datablocks.__iter__()
                next = iter.next()
                while (1):
                        if (verbose):
                                print "check ", next, datablocks[next][0]
				sys.stdout.flush()
                        if (f.root._f_getChild(part_name).__contains__(datablocks[next][0])):
                                if (next.find(tag)>-1):
					contains_flag=True	
                        try:
                                next=iter.next()
                        except StopIteration:
                                break
  f.close() 
  return contains_flag

############
#CHECK FILE#
############
def check_file(filename):
  f=tables.openFile(filename)
  f.close()
                                                                                                                                                  






#####################################################################################################################
#                                                    WRITING ROUTINES    		                            #
#####################################################################################################################



#######################
#OPEN FILE FOR WRITING#
#######################
def openfile(filename, mode="w"):
	f=tables.openFile(filename, mode = mode)	 
	return f

############
#CLOSE FILE#
############
def closefile(f):
	f.close()

##############################
#WRITE SNAPSHOT HEADER OBJECT#
##############################
def writeheader(f, header):	
    	group_header=f.createGroup(f.root, "Header")
	group_header._v_attrs.NumPart_ThisFile=header.npart
	group_header._v_attrs.NumPart_Total=header.nall
	group_header._v_attrs.NumPart_Total_HighWord=header.nall_highword
	group_header._v_attrs.MassTable=header.massarr
	group_header._v_attrs.Time=header.time
	group_header._v_attrs.Redshift=header.redshift
	group_header._v_attrs.BoxSize=header.boxsize
	group_header._v_attrs.NumFilesPerSnapshot=header.filenum						
	group_header._v_attrs.Omega0=header.omega0							
	group_header._v_attrs.OmegaLambda=header.omegaL							
	group_header._v_attrs.HubbleParam=header.hubble	
	group_header._v_attrs.Flag_Sfr=header.sfr	
	group_header._v_attrs.Flag_Cooling=header.cooling
	group_header._v_attrs.Flag_StellarAge=header.stellar_age		
	group_header._v_attrs.Flag_Metals=header.metals		
	group_header._v_attrs.Flag_Feedback=header.feedback				
	group_header._v_attrs.Flag_DoublePrecision=header.double			

###############
#WRITE ROUTINE#
###############
def write_block(f, block, parttype, data):
	part_name="PartType"+str(parttype)
	if (f.root.__contains__(part_name)==False):
	    	group=f.createGroup(f.root, part_name)
	else:
		group=f.root._f_getChild(part_name)	
	
	if (datablocks.has_key(block)):
        	block_name=datablocks[block][0]
	        dim2=datablocks[block][1]		
		if (group.__contains__(block_name)==False):
			table=f.createArray(group, block_name, data)
		else:
			print "I/O block already written"
			sys.stdout.flush()
	else:
		print "Unknown I/O block"
		sys.stdout.flush()		

	



