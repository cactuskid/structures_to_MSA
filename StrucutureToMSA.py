#run fatcat with structures in a folder
#condense alignment into MSA
#it requires Fatcat and clustalo to work
#fatcat is the java based version included with the protein comparison tool
#from the PDB website.
#the phylogeny is derived from a distance matrix
#I used phylip kitsch to generate the newick tree. but any newick tree should work.
from Bio.PDB import *
import glob
import subprocess
import itertools
import shlex
from Bio import AlignIO , SeqIO
import pickle
import re
import numpy as np
import pandas as pd
import math
import pylab
import scipy.cluster.hierarchy as sch
from itertools import combinations, permutations
import scipy.spatial.distance as distutils
import multiprocessing as MP
import tempfile as tmp


#steps in the process
#compute all pairwise alignments with fatcat
grabPDBs = True
align = False
Fatcatalign = False
TMalign = False
#realign flexed structures with TM align
TMalign_from_Fatcat = True

#turn the output into fasta
parse = False
#use clustalo to create one sequence alignment of the strucutres
merge = False
#output a distance matrix based on SDM that you can use for phylogeny building
distmat = True
#use phylip kitsch after this step on the output : phylipmat.txt

#generate a figure of the distmat and the phylogenetic tree
show_distmat = False

#this is just a dumb hack to create non redundant fastas to analyze the merging of the intermediate alignments
resolve_names = False



#where a your structures? ... should all be single chain with decent homology. 
#filedir = '../domainIII/'
filedir = '../effHapViralNew/'
outpath = filedir
#where is fatcat?
FatcatPath =  './runFATCAT.sh' 
TMalignPath = './TMalign/TMalign'
#where is clustalo?
clustaloPath = 'clustalo'



def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'r') as f:
        return pickle.load(f)

def runFatcat(pargs):
	FatcatPath , struct1, struct2, outfile , outPDB = pargs
	
	args =  FatcatPath + ' -file1 ' + struct1 + ' -file2 ' + struct2 + '  -printFatCat -flexible True -outputPDB  -outFile '+ outPDB
	print args
	args = shlex.split(args)
	fatcatfile = open(outfile, 'w')
	p = subprocess.call(args  , stdout= fatcatfile )
	fatcatfile.close()
	return outfile
	

def runTMalign(pargs):
	FatcatPath , struct1, struct2, outfile , outPDB = pargs
	args =  TMalignPath + ' ' + struct1 + ' ' + struct2 
	print args
	args = shlex.split(args)
	TMfile = open(outfile, 'w')
	p = subprocess.call(args , stdout= TMfile )
	TMfile.close()
	return outfile
	

def runTMalign_FatcatOut(pargs):
	TMalignPath , Fatcatstruct, outfile  = pargs
	TMfile = open(outfile, 'w')
	"""
	TITLE  jFatCat_flexible 1.0 domain3_4hj1.pdb vs. domain3_4adi.pdb
	EXPDTA    NMR, 2 STRUCTURES
	MODEL      1
	"""
	coords = {}
	header = ''
	i =0
	start = False
	with open( Fatcatstruct , 'r' ) as PDBfile:
		for line in PDBfile:
			if start == False:
				if 'EXPDTA' in line:
					header += line.replace('2','1')
				else:	
					header += line
			if 'TITLE' in line:
				pdb2 = line.replace('TITLE  jFatCat_flexible 1.0' , '' ).split('vs.')[0].strip()
				pdb1 = line.replace('TITLE  jFatCat_flexible 1.0' , '' ).split('vs.')[1].strip()
			if 'MODEL' in line:
				start = True
				i+=1
				coords[i] = ''
				
			if start == True:
				coords [i] += line
	filehandle1 = tmp.NamedTemporaryFile( 'w' , dir = './')
	filehandle2 = tmp.NamedTemporaryFile( 'w' , dir = './')
	handles = [filehandle1, filehandle2]
	pdbs = [pdb1, pdb2]
	for i,struct in enumerate(coords):
		handles[i].write(header) 
		handles[i].write(coords[struct]) 
		handles[i].write('ENDMDL')	
	args =  TMalignPath + ' ' + handles[0].name + ' ' + handles[1].name 
	print args
	
	args = shlex.split(args)
	TMfile = open(outfile, 'w')
	p = subprocess.call(args , stdout= TMfile )
	
	for handle in handles:
		handle.close()
	TMfile.close()
	return outfile , pdbs

def testTMalign_fromFatacat():
	TMalignPath = './TMalign/TMalign'
	Fatcatstruct = '../effHapViralNew/1SVB_4ADI.factcatPDB'
	outfile = 'tmaligntest.txt'
	pargs =TMalignPath , Fatcatstruct, outfile  
	outfile, pdbs =runTMalign_FatcatOut(pargs)
	print pdbs 
	with open(outfile , 'r') as output:
		for line in output:
			print line
	
def testTMalign():
	pd1 = '../structural_alignment/effHapViralcontrol/1SVB.pdb'
	pd2 = '../structural_alignment/effHapViralcontrol/2ala.pdb'
	TMalignPath = './TMalign/TMalign'
	outpath = '../structural_alignment/effHapViralcontrol/'
	outfile = '../structural_alignment/effHapViralcontrol/TMaligntest.txt'

	runTMalign([TMalignPath,pd1, pd2, outfile , outpath])
	

def mergeAlign(clustalopath , fasta1, fasta2, outfile):
	args =  clustalopath + ' --force --p1 ' + fasta1 + ' --p2 ' + fasta2 + ' -o ' + outfile
	print args
	args = shlex.split(args)
	p = subprocess.call(args )

def FatcatToFasta(fatout, filename):
	print fatout 
	print filename
	alns = ['','']
	files = []
	BlockNext = False
	Blocknum = 0
	with open( fatout , 'r' ) as fatfile:
		for line in fatfile:
			try:
				print line
				if 'file from local' in line:
					files.append(line.split('/')[-1])	
				#format  = Chain 1:   16 PHC---------------SKTPIVRAQTSQNAMS-------RGMQMQFSIGLHT---------AVC----
				if 'Chain 1:' in line:
					alns[0] += line.split(':')[1].strip().split(' ')[1]
				if 'Chain 2:' in line:
					alns[1] += line.split(':')[1].strip().split(' ')[1]
			except:
				print 'line error:'
				print line
		alnstr1 = '>'+files[0][:-1] + 'Fatcat_block' +'\n' +alns[0] + '\n'
		alnstr2 = '>'+files[1][:-1] + 'Fatcat_block' +'\n' +alns[1] + '\n'
		print alnstr1
		print alnstr2
		handle = open(filename  , 'w')
		handle.write(alnstr1 + '\n')
		handle.write(alnstr2)
		handle.close()
		return filename


def FatcatToDistances(fatout):
	alns = ['','']
	files = []
	BlockNext = False
	Blocknum = 0
	pval = 0
	distances = {}
	counts = {}
	numbers = re.compile('[-+]?[0-9]*\.?[0-9]+.')
	integers = re.compile('\s[0-9]+\s')
	with open( fatout , 'r' ) as fatfile:
		for i,line in enumerate(fatfile):
			if 'P-value' in line:
				pval = float(line.split(' ')[1])
				#P-value 1.44e-03 Afp-num 65162 Identity 4.96% Similarity 16.08%
			if 'Align' in line:
				#Align 1SVB.pdb 395 with ememmodelFinalarabidop_chopped.pdb.pdb 467
				files = line.split(' ')
				pdbs = []
				for filename in files:
					if '.pdb' in filename.lower():
						pdbs.append(filename)	
				lengths = integers.findall(line)
				cleanlen = []
				for val in lengths:
					cleanlen.append(int(val.strip()))
			if 'Block' in line:
				#Block  0 afp 17 score 201.70 rmsd  3.92 gap 128 (0.48%)
				words = line.split('rmsd')[1]
				rmsd  = float(numbers.findall(words)[0].strip())
				distances[Blocknum] = rmsd
				Blocknum+=1
			if 'Block' not in line and 'Chain' not in line and i > 6:
				for j in range(len(distances)):
					if j not in counts:
						counts[j] = 0
					counts[j]+= line.count(str(j+1))
	return pdbs, cleanlen, distances, counts , pval
	
def TMaligntoDistance(TMalignfile):
	#parse TM align file and output parameters
	print TMalignfile
	"""
	format:
		 **************************************************************************
	 *                        TM-align (Version 20160521)                     *
	 * An algorithm for protein structure alignment and comparison            *
	 * Based on statistics:                                                   *
	 *       0.0 < TM-score < 0.30, random structural similarity              *
	 *       0.5 < TM-score < 1.00, in about the same fold                    *
	 * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)    *
	 * Please email your comments and suggestions to: zhng@umich.edu          *
	 **************************************************************************

	Name of Chain_1: ../structural_alignment/effHapViralcontrol/1SVB.pd
	Name of Chain_2: ../structural_alignment/effHapViralcontrol/2ala.pd
	Length of Chain_1:  395 residues
	Length of Chain_2:  384 residues

	Aligned length=  325, RMSD=   5.63, Seq_ID=n_identical/n_aligned= 0.062
	TM-score= 0.56731 (if normalized by length of Chain_1)
	TM-score= 0.57996 (if normalized by length of Chain_2)
	(You should use TM-score normalized by length of the reference protein)

	(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
	SRCTHLENRDFVTGTQG-TTRVTLVLELGGCVTITAEG-----K-PSMDVW-LDAIYQENPAKTREYCLHAKLSDTKVAARCPTMGPATLAEEHQGGTVCKRDQ-SDRGW-G-N-HC-GLFGKGSIVACVKAAC--EAKKKATGHVYDANKIVYTVKVEPHTGDYVAANET-HSGRKTASFTISSEKTILTMGEYGDVSLLCRVASGVDLAQTVILELDKTVEHLPTAWQVHRDWF-ND----LALPWKHE-GAQNWNNAERLVEFGAPHA--VKMDV-YNLGDQTGVLLKALAGVPV-------AHIEGTKYHLK-SGHVTCEVGLEK---LKMKG--LTYT-MCDKTKFTWKRAPTDSGHDTVVMEVTFS-GT-----KPCRIPVRAVAHGSPDVNVAMLITPNPTIENN-----GGGFIEMQLPPGDNIIY-V-------GELSHQWFQK-
					  .::::. .......::::::     : :::::: .::::.. :::: ::::::::: ::::::..::::.. . .::::.::::.: .:::: : : :: ::...:::::::.:::  ::::::::::...  .::::::::...         .:::::::..  .:::::: :  .:.::.::...:.....::::::   :  :: :::::.... ..    :::::::: .:::::: ::......:::  ::::. ..::.::.::.::.:....       ::::  ::::. .....::::..:   .....  .:.. ..... ...  .:.:.    .::....  .      ........         ...::::..::...      ...::::.:   ::... .       .....::..  
	-----------------YEHSTVM-PNVVGFPYKAHIERPGYSPLTLQMQVVETSLEPT-LNLE-YITCEYKTV-VPSPYVKCCGASEC-S-TKEKPDYQCKVYTGVYPFMWGGAYCFCDSENTQLSEAYVDRSDVCRHDHASAYKAHT--ASLKAKVRVMYG--------NVNQTVDVYVN--GDHAVTI-G--GTQFIFGPLSSAWTPFDNKIVVY---K--DE-VFNQDFPPYGSGQPGRFGDIQSRTVESNDLYA-NTALKLARPSPGMVHVPYTQTPSGFKYWLKEKGTALNTKAPFGCQIKTN--PVRAMNCAVGNIPVSMNLPDSAFTRIVEAPTIIDLTCT-VAT--CTHSS----DFGGVLT-LT-YKTNKNGDCSVHS---------HSNVATLQEATAKV-KTAGKVTLHFSTAS---ASPSFVVSLCSARATCSASCEPP-K



	"""
	#maybe incorporate more stuff into this later
	
	TMscores = []
	pdbs = []
	
	
		
	with open( TMalignfile, 'r') as output:
		for line in output:
			if 'Length of Chain_1:' in line:
				chain1len = int(line.split(':')[1].split()[0] )
			if 'Length of Chain_2:' in line:
				chain2len = int(line.split(':')[1].split()[0]	)
			if 'Name of Chain_' in line and '*' not in line :
				pdbs.append(line.split(':')[1].strip())
			if 'TM-score' in line and '*' not in line and ' normalized by length of the reference protein' not in line:
				score = float(line.split('=')[1].split( )[0])
				TMscores.append(score)
	if chain1len>chain2len:
		return TMscores[1], pdbs
	else :
		return TMscores[0], pdbs

def FatcatToDF(fatout ):
	print fatout
	alns = ['','']
	files = []
	pdbfile =fatout+ 'PDB'
	BlockNext = False
	Blocknum = 0
	pval = 0
	df1 = None
	df2 = None
	chain1_read = False
	chain2_read = False

	distances = {}
	counts = {}
	df1dict = {}
	df2dict = {}
	numbers = re.compile('[-+]?[0-9]*\.?[0-9]+.')
	integers = re.compile('\s[0-9]+\s')
	number = re.compile('[0-9]+')
	
	sequence = re.compile('[A-Z-]+')
	
	with open( fatout , 'r' ) as fatfile:
		for i,line in enumerate(fatfile):
			if 'Align' in line:
				#Align 1SVB.pdb 395 with ememmodelFinalarabidop_chopped.pdb.pdb 467
				files = line.split(' ')
				pdbs = []
				for filename in files:
					if '.pdb' in filename.lower():
						pdbs.append(filename)	
				lengths = integers.findall(line)
				cleanlen = []
				for val in lengths:
					cleanlen.append(int(val.strip()))

			
			if 'Block' in line:
				#Block  0 afp 17 score 201.70 rmsd  3.92 gap 128 (0.48%)
				words = line.split('rmsd')[1]
				rmsd  = float(numbers.findall(words)[0].strip())
				distances[Blocknum] = rmsd
				Blocknum+=1

			#Chain 1:   19 TWVDLVLEGDSCVTIMSK----DKPTIDVKMMNMEAANLAEVRSYCYLATVSDLSTKAACPTMGEAHNDK
	            #  			   111111111111111111    111111111111111 111111111111111111         11111
			#Chain 2:    3 HVTVIPNTVGVPYKTLVNRPGYSPMVLEMELLSVTLEPTLSLDYITCEYKTVIPSPYVKCCGTAECKDKN
			
			if chain1_read == True:
				#read block asiignment
				chain1_read = False
				blocks = line[14:]
				
			if 'Chain 1' in line:
				chain1_read = True
				cleanline = line.split(':')[1].strip()
				start1 = int(number.findall(cleanline)[0])
				alnseq1 = sequence.findall(cleanline)[0]

			if 'Chain 2' in line:
				chain2_read = True
				cleanline = line.split(':')[1].strip()
				start2 = int(number.findall(cleanline)[0])
				alnseq2 = sequence.findall(cleanline)[0]

			if chain2_read == True:
				chain2_read = False
				#parse the two chains and map residues to block, gap or aln positions 
				#residue count
				count1 = 0
				count2 = 0
				for i,char1 in enumerate(alnseq1):
					if i < len(alnseq2):
						Blocknum = blocks[i]
						char2 = alnseq2[i]
						skip = False
						if char1 == '-':
							# assign values for the residue that opened the gap
							row2 = [char2,'NONE', 'gap' , 'NONE' , 'NONE']	
							df2dict[start2+count2] = row2
							count2 +=1
							skip = True
						if char2 == '-':
							row1 = [char1,'NONE', 'gap' , 'NONE' ,  'NONE']	
							df1dict[start1+count1] = row1
							count1 +=1
							skip = True
						if skip == False:
							#spaces will be counted as aligned but with no blocknum
							row1 = [ char1, Blocknum , 'aln' , start2+count2 , char2 ]
							df1dict[start1+count1] = row1
							row2 = [ char2, Blocknum , 'aln' , start1+count1 , char1 ]
							df2dict[start2+count2] = row2
							count2 +=1
							count1 +=1
	
	df1 = pd.DataFrame.from_dict( df2dict , orient= 'index' )
	df2 = pd.DataFrame.from_dict( df1dict , orient= 'index' )

	df1.columns = ['res','block','status','match','matchres']
	df2.columns = ['res','block','status','match','matchres']

	
	return df1, df2 , pdbs, pdbfile

def norm2coords( args ):
	coords1 , coords2 , i , j = args
	NORM = np.linalg.norm(coords1[i,:]-coords2[j,:])
	return i , j , NORM

def qalnscore( args):
	i, j , index1a, index2a , index1b , index2b , Adistmat , Bdistmat = args 
	qaln = np.exp(- math.pow(Adistmat[index1a, index2a] - Bdistmat[index1b, index2b] ,2 )  / (2*math.pow(math.fabs(i-j),1.5 )  ) )
	return qaln

	
def calculateQh(pargs):
	pairwisePDB,df1,df2,pdbs = pargs
	#pairwise pdb is the aligned fatcat output
	#df1 should be the first structures dataframe showing where it aligned etc
	parser = PDBParser(PERMISSIVE=1)
	structure = parser.get_structure('Model 1' ,pairwisePDB)
	model1 = structure[0]
	model2 = structure[1]
	residues1 = model1.get_residues()
	residues2 = model2.get_residues()


	resdict1 = {}
	for res in residues1:
		if 'CA' in res:
			resdict1[res.id[1]] = res['CA'].coord
	
	resdict2 = {}
	for res in residues2:
		if 'CA' in res:
			resdict2[res.id[1]] = res['CA'].coord

	
	resdf1 = pd.DataFrame.from_dict( resdict2 , orient= 'index' )
	resdf2 = pd.DataFrame.from_dict( resdict1 , orient= 'index' )
	
	resdf1.columns = ['x','y','z']
	resdf2.columns = ['x','y','z']
	
	
	df1 =pd.concat([df1, resdf1], axis=1)
	df2 = pd.concat([df2, resdf2], axis =1)
	
	Acoords = resdf1[['x','y','z']].as_matrix()
	Bcoords = resdf2[['x','y','z']].as_matrix()
	
	
	pool = MP.Pool()
	
	Adistmat = np.zeros((len(resdf1.index), len(resdf1.index) ) )
	rundata = []
	for i , res1 in enumerate(resdf1.index):
		for j, res2 in enumerate(resdf1.index):
			if i < j :
				rundata.append([Acoords , Acoords , i , j])
	results = pool.map_async( norm2coords , rundata).get()
	for result in results:
		i,j,norm = result
		Adistmat[i,j] = norm
	Adistmat += Adistmat.T

	Bdistmat = np.zeros((len(resdf2.index), len(resdf2.index) ) )
	rundata = []
	for i , res1 in enumerate(resdf2.index):
		for j, res2 in enumerate(resdf2.index):
			if i < j :
				rundata.append([Bcoords , Bcoords , i , j])
	results = pool.map_async( norm2coords , rundata).get()
	for result in results:
		i,j,norm = result
		Bdistmat[i,j] = norm
	Bdistmat += Bdistmat.T

	print 'done distmats'

	alnDF2 = df2.query( 'status == "aln"')
	alnDF1 = df1.query( 'status == "aln"')
	
	resnums1 = resdf1.index.tolist()
	resnums2 = resdf2.index.tolist()
	qaln = 0
	rundata = []
	#some discontinuties in PDBs... check if the residue exists
	#residue start is rest at each fatcat line so the error shouldnt propagate much
	
	for i in alnDF1.index:
		for j in alnDF1.index:
			if i <= j-2:
				if i in resnums1 and j in resnums1:
					index1a = resnums1.index(i)
					index2a = resnums1.index(j)
					iprime = alnDF1['match'][i]
					jprime = alnDF1['match'][j]
					if iprime in resnums2 and jprime in resnums2:
						index1b = resnums2.index(iprime)
						index2b = resnums2.index(jprime)
						rundata.append([i, j , index1a, index2a , index1b , index2b , Adistmat , Bdistmat])
	results = pool.map_async( qalnscore , rundata).get()
	qaln = sum(results)
	print qaln
	print 'done qaln'
	gapDF2 = df2.query( 'status == "gap"')
	gapDF1 = df1.query( 'status == "gap"')
	qgap = 0
	g0 = 0
	g1 = 0 
	g2 = 0
	
	
	coords = []
	rundata = []
	#calculate the two sets of qgap scores
	for gapDF, alnDF, resn1, resn2 , Adist, Bdist in [(gapDF1,alnDF1,resnums1,resnums2,Adistmat , Bdistmat), (gapDF2,alnDF2,resnums2,resnums1, Bdistmat, Adistmat)]:
			
		print len(gapDF.index)
		print len(alnDF.index)
		for gapRes in gapDF.index:
			try:
				#closest aligned residue is edge1
				indexArray = np.asarray(alnDF.index)
				gapedge1 = indexArray[np.argmin(np.absolute(indexArray - gapRes )) ]
				#find opposite gap edge
				if gapedge1 > gapRes:
					otherside = indexArray[np.where(indexArray< gapRes)] 
					gapedge2 =  otherside[np.argmin(np.absolute( otherside - gapRes ) ) ]
				else:
					otherside = indexArray[np.where(indexArray>gapRes)] 
					gapedge2 =  otherside[np.argmin(np.absolute( otherside - gapRes ) ) ]
				gapedges = [gapedge1, gapedge2]
				if np.amin(np.absolute(np.asarray(gapedges)-gapRes) )> 2:
					g2+=1
				if np.amin(np.absolute(np.asarray(gapedges)-gapRes) ) == 2:
					g1+= 1
				if np.amin(np.absolute(np.asarray(gapedges)-gapRes) ) < 2:
					g0+= 1
				
				for alnRes in alnDF.index:
					index1a = resn1.index(gapRes)
					index2a = resn1.index(alnRes)
					#match j' and edges from prot A
					index1b = resn2.index(alnDF['match'][alnRes])
					for edge in gapedges:											
						index2b =  resn2.index(alnDF['match'][edge])
						rundata.append([resn1.index(gapRes), resn1.index(alnRes) , index1a, index2a , index1b , index2b , Adist , Bdist ])
						coords.append([gapRes,edge, alnRes])
			except:
				pass
				
		#calculate all and then use coords to find maxima before calculating the sum 
		results = pool.map_async(qalnscore, rundata).get()
		print len(results)
		scores = {}
		for i , score in enumerate(results):
			gapRes,edge, alnRes = coords[i]
			if (gapRes, edge, alnRes) not in scores:
				scores[(gapRes, edge, alnRes)]= score
			else:
				if scores[(gapRes, edge, alnRes)] < score:
					scores[(gapRes, edge, alnRes)] = score
		print sum(scores.values())
		qgap += sum(scores.values())
	pool.close()
	Naln = len(alnDF1.index)
	Ngap = len(gapDF1.index) + len(gapDF2.index)
	N = .5* (Naln -1)*(Naln -2) + g0*(Naln) + g1 *( Naln-1) + g2 * (Naln -2) + Ngap
	print {'N': N , 'Naln':Naln, 'Ngap' : Ngap, 'qaln': qaln, 'qgap':qgap , 'g0': g0 , 'g1':g1 , 'g2':g2}
	finalscore =  (qaln + qgap) / N
	print finalscore
	print 'DONE'
	return finalscore, N , qaln, qgap , [Naln, g0, g1, g2] ,pdbs


def testQHscore():
	df1, df2 , pdbs, pdbfile = FatcatToDF(filedir+'1SVB_4ADI.factcat')
	Qh = calculateQh([pdbfile, df1, df2])

if grabPDBs == True:
	structures = glob.glob(filedir + '*.pdb')
	print structures
	save_obj(structures, filedir + 'pdblist')
if align == True:
	print 'pairwise aligning structures'
	structures = load_obj(filedir +'pdblist')
	print structures
	pool = MP.Pool()
	if Fatcatalign == True :
		fatcatfiles = []
		runData = []
		for struct1,struct2 in itertools.combinations(structures, 2):
			outfile = outpath + struct1.split('/')[-1].replace('.pdb','')+ '_'+struct2.split('/')[-1].replace('.pdb','') + '.factcat'
			outPDB = outfile + 'PDB'
			runData.append([FatcatPath ,struct1,struct2, outfile, outPDB])
		
		print len(runData)
		print 'alignments to run'
		print 'fatcat align allvall'
		results = pool.map_async( runFatcat , runData, MP.cpu_count() ).get()
		save_obj( results, filedir +'fatcatfiles')

	if TMalign == True:
		TMalignfiles = []
		runData = []
		for struct1,struct2 in itertools.combinations(structures, 2):
			outfile = outpath + struct1.split('/')[-1].replace('.pdb','')+ '_'+struct2.split('/')[-1].replace('.pdb','') + '.TMalign'
			outPDB = outfile + 'PDB'
			runData.append([FatcatPath ,struct1,struct2, outfile, outPDB])
		print 'TMalign align allvall'
		results = pool.map_async( runTMalign , runData, MP.cpu_count() ).get()
		save_obj( results, filedir +'TMalignfiles')
	if TMalign_from_Fatcat == True:
		#run this after fatcat!
		runData = []
		alignments = load_obj(filedir + 'fatcatfiles')
		for aln in alignments:
			runData.append([TMalignPath , aln+'PDB' ,  aln+'TMRealign'])
		print runData
		results = pool.map_async( runTMalign_FatcatOut, runData, MP.cpu_count() ).get()
		save_obj( results, filedir +'TMRealignfiles')
	pool.close()
	
	print 'done aligning'
	

if parse == True:
	print filedir
	alignments = load_obj(filedir + 'fatcatfiles')
	fastas = []
	distances = {}
	print 'converting fatcat output to fasta'
	for aln in alignments:
		fastas.append(FatcatToFasta(aln,aln+'.fasta'))
		distances[aln] = FatcatToDistances(aln)
		
	print len(fastas)
	print 'alignments by fatcat'
	save_obj( fastas, filedir+'fastas')
	save_obj(  distances, filedir +'distances')
	print 'done writing fastas'

def distmat_to_txt( pdblist , distmat, filedir , name):
		#write out distmat in phylip compatible format
	outstr =' ' + str(len(pdblist)) + '\n'
	for i,pdb in enumerate(pdblist):
		if len(pdb)>10:
			namestr= pdb[0:10]
		if len(pdb)<10:
			namestr = pdb
			for pad in range(10 -len(pdb)):
				namestr += ' '
		outstr += namestr+ ' ' + np.array2string( distmat[i,:], formatter={'float_kind':lambda x: "%.2f" % x}).replace('[', '').replace(']', '')  + ' \n'
	print outstr
	handle = open(filedir + name + 'phylipmat.txt' , 'w')
	handle.write(outstr)
	handle.close()

	outstr = str(len(pdblist)) + '\n'
	for i,pdb in enumerate(pdblist):
		namestr = pdb.replace('.','').replace('_','')[0:20]
		outstr += namestr+ ' ' + np.array2string( distmat[i,:], formatter={'float_kind':lambda x: "%.2f" % x}).replace('[', '').replace(']', '').replace('\n', '')  + '\n'
	
	print outstr
	handle = open(filedir + name + 'fastmemat.txt' , 'w')
	handle.write(outstr)
	handle.close()
	
	
if distmat == True:
	
	distances = load_obj(filedir + 'distances')
	maxdist = 0
	pdblist = []
	for out in distances.values():
		
		pdblist += out[0]
		for dist in out[2]:
			if dist > maxdist:
				maxdist = dist
	pdblist = list(set(pdblist))
	print 'creating SDM distmat for Fatcat output'
	
	distmat = np.zeros((len(pdblist),len(pdblist)))
	pvalmat = np.zeros((len(pdblist),len(pdblist)))

	#calculate SDM based distmat
	for out in distances.values():
		lensmallprot = min(out[1])
		SRMSs = []
		count = 0
		pdbs = out[0]
		distances = out[2]
		counts = out[3]
		pval = out[4]
		for block in distances:
			count += counts[block]
			SRMSs.append(counts[block]*(1-distances[block]/maxdist))
		PFTE = count / lensmallprot
		avgSRMS = sum(SRMSs)/count
		w1 = (1 -PFTE + 1 -avgSRMS) /2 
		w2 = (PFTE +avgSRMS)/2
		sdm = -100* math.log(w1*PFTE + w2 *avgSRMS)
		distmat[ pdblist.index(pdbs[0]),pdblist.index(pdbs[1])] = sdm
		#low pval -> high score
		pvalmat[ pdblist.index(pdbs[0]),pdblist.index(pdbs[1])] = pval
	distmat += distmat.T
	#distmat /= np.amax(distmat)
	distmat_to_txt( pdblist , distmat, filedir , 'SDMmat')
	distmat_to_txt( pdblist , pvalmat, filedir , 'pvalmat')
	
	
	
	"""
	print 'creating QH and Qaln distmats for Fatcat output'
	alignments = load_obj(filedir + 'fatcatfiles')

	#qh based distmats
	QHmat = np.zeros((len(pdblist),len(pdblist)))
	Qalnmat = np.zeros((len(pdblist),len(pdblist)))
	
	runData = []
	for filename in alignments:
		df1,df2,pdbs,pdbfile = FatcatToDF(filename)
		runData.append( [pdbfile ,df1,df2 , pdbs])
	
	
	
	#pool = MP.Pool()
	#results = pool.map_async( calculateQh , runData ).get()
	for args in runData:
		results = calculateQh(args)
		finalscore, N , qaln, qgap , factors, pdbs = results
		Naln, g0, g1, g2 = factors
		QHmat[ pdblist.index(pdbs[0]),pdblist.index(pdbs[1])] = finalscore
		Qalnmat[ pdblist.index(pdbs[0]),pdblist.index(pdbs[1])] = qaln / N
	QHmat += QHmat.T
	Qalnmat += Qalnmat.T
	distmat_to_txt( pdblist , QHmat, filedir , 'QHmat')
	distmat_to_txt( pdblist , Qalnmat, filedir , 'Qalnmat')
	
	invQH = np.fill_diagonal(1/QHmat, 0 )
	invQaln = np.fill_diagonal(1/Qalnmat, 0 )
	distmat_to_txt( pdblist , invQH 	, filedir , 'invQHmat')
	distmat_to_txt( pdblist , invQaln , filedir , 'invQalnmat')
	
	save_obj( QHmat , filedir + 'QHmat')
	save_obj(   Qalnmat , filedir + 'Qalnmat')
	
	pdblist = load_obj(filedir + 'pdblist')
	print pdblist
	TMdistmat = np.zeros((len(pdblist),len(pdblist)))
	print 'creating TMalign based distmat'
	TMalignresults = load_obj(filedir + 'TMalignfiles')
	tmlist =[]
	namelist = []
	
	
	for filename in TMalignresults:
		distance, pdbs =TMaligntoDistance(filename)
		for filename in pdbs:
			if filename not in tmlist:
				tmlist.append(filename)
				namelist.append(filename.split('/')[-1])
	for filename in TMalignresults:
		distance, pdbs =TMaligntoDistance(filename)
		TMdistmat[ tmlist.index(pdbs[0]),tmlist.index(pdbs[1])] = distance
	TMdistmat += TMdistmat.T
	TMdistmat = (1-TMdistmat)
	np.fill_diagonal(TMdistmat, 0 )
	print namelist
	distmat_to_txt( namelist , TMdistmat, filedir , 'TMdistmat')
	save_obj(   TMdistmat , filedir + 'TMdistmat')
	
	"""
	
	tmlist = []
	namelist = []
	
	TMRealignresults = load_obj( filedir +'TMRealignfiles')
	TMdistmat = np.zeros((len(pdblist),len(pdblist)))
	print len(TMRealignresults)
	for filename,pdbs in TMRealignresults:
		for filename in pdbs:
			if filename not in tmlist:
				tmlist.append(filename)
				namelist.append(filename.split('/')[-1])

	TMdistmat = np.zeros((len(tmlist),len(tmlist)))
	
	for filename , pdbs in TMRealignresults:
		distance, pdbs2 =TMaligntoDistance(filename)
		TMdistmat[ tmlist.index(pdbs[0]),tmlist.index(pdbs[1])] = 1-distance
	TMdistmat += TMdistmat.T
	print TMdistmat
	distmat_to_txt( namelist , TMdistmat, filedir , 'TMRealigndistmat')
	save_obj(   TMdistmat , filedir + 'TMRealigndistmat')
	
	
	fig = pylab.figure(figsize=(8,8))
	axmatrix = fig.add_axes([0.1,0.1,.7,.7] )
	im = axmatrix.matshow(TMdistmat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu, interpolation = 'None')
	axmatrix.yaxis.tick_right()
	axmatrix.set_yticks( range(len(tmlist)) )
	axmatrix.set_xticks( range(len(tmlist)) )
	axmatrix.set_yticklabels( namelist , fontsize = 10 )
	axmatrix.set_xticklabels( namelist , fontsize = 10 , rotation = 90)
	
	pylab.show()


if show_distmat:
	pdblist = load_obj(filedir + 'pdblist')
	D= load_obj(filedir + 'distmat')
	# Compute and plot first dendrogram.

	#you need to calculate the phyllip kitsch outtree with the distmat before running this part.
	#it should be in the same folder. as the PDBS etc
	#this will output a graph with the distmat and the grouping given by phylip on the same graph


	fig = pylab.figure(figsize=(8,8))
	tree = PhyloTree( filedir +'outtree', sp_naming_function=None) 
	leaves = tree.get_leaf_names()
	idx_dict = {}


	for i,pdb in enumerate(leaves):
		idx_dict[pdb] = i
	dmat = np.zeros((len(leaves),len(leaves)))
	idx_labels = [idx_dict.keys()[idx_dict.values().index(i)] for i in range(0, len(idx_dict))]
	for l1,l2 in combinations(leaves,2):
		d = tree.get_distance(l1,l2)
		dmat[idx_dict[l1],idx_dict[l2]] = dmat[idx_dict[l2],idx_dict[l1]] = d
	schlink = sch.linkage(distutils.squareform(dmat),method='average',metric='euclidean')
	ax2 = fig.add_axes([0.09,0.1,0.2,0.6])
	
	Z2 = sch.dendrogram(schlink, orientation = 'right')#, labels = pdblist)#,  leaf_rotation=90 )
	ax2.set_yticks([])

	# Plot distance matrix.
	axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
	idx1 = Z2['leaves']
	D = D[idx1,:]
	D = D[:,idx1]
	pdblist = np.asarray(leaves)[idx1]
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu, interpolation = 'None')
	axmatrix.yaxis.tick_right()
	axmatrix.set_yticks( range(len(pdblist)) )
	axmatrix.set_xticks( range(len(pdblist)) )
	axmatrix.set_yticklabels( pdblist , fontsize = 10 )
	axmatrix.set_xticklabels( pdblist , fontsize = 10 , rotation = 90)
	
	# Plot colorbar.
	axcolor = fig.add_axes([0.3,0.07,0.6,0.02] )
	pylab.colorbar(im, cax=axcolor, orientation = 'horizontal')
	fig.subplots_adjust(right=5)
	pylab.show()

if merge == True:
	print alignments
	print 'merging all alignments'
	i = 0
	alignments = load_obj(filedir + 'fastas')
	alignMerge ={}
	alignMerge[i] = alignments
	while len(alignMerge[i]) > 2:
		#halve the alignment list with each iteration and grab the remainder for the next round
		alignMerge[i + 1] = []
		nextmerge = []
		print 'round ' + str(i)
		print alignMerge[i]
		for j,fasta in enumerate(alignMerge[i]):
			if len(nextmerge) == 2:
				mergeAlign(clustaloPath  ,nextmerge[0],nextmerge[1] , nextmerge[0] + str(i)+ '.fasta' )
				alignMerge[i+1].append(nextmerge[0] + str(i)+ '.fasta')
				nextmerge = []
				nextmerge.append(fasta)
			else:
				nextmerge.append(fasta)
		else:	
			alignMerge[i+1] += nextmerge
		i+=1
	finalFasta = alignMerge[alignMerge.keys()[-1]][0]
	#pare down the alignment so that it is non-redundant
	record_iterator = AlignIO.read(finalFasta, "fasta")
	writeout=[]
	titles = []
	print 'final fasta '
	print finalFasta

	for record in record_iterator:
		if record.id not in titles:
			titles.append(record.id)
			print record
			writeout.append(record)

	print len(writeout)

if resolve_names == True:
	fastas = glob.glob(filedir + '*fasta')
	for fasta in fastas:
		if 'nr' not in fasta:
			record_iterator = AlignIO.read(fasta, "fasta")
			writeout=[]
			redundant = []
			for i,record in enumerate(record_iterator):
					record.id = str(i)+ '_'+record.id
					record.description = ''
					writeout.append(record)
			handle = open(fasta + 'nr.fasta' , 'w')
			SeqIO.write( writeout, handle, 'fasta')
			handle.close()
