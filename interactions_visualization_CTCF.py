# coding=utf-8
import	sys
from PIL	import	Image, ImageDraw, ImageFont
#import	Image, ImageDraw, ImageFont
import	argparse
import	math
import numpy as np
from	scipy import spatial
import	random
## asdadcod ing=utf-8

def	read_data_CTCF(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	# The link coordinates and its stength
	Link = {}
	for i in x:
		xx = i.split('\t')
		#if xx[0] == Ch:
		Link[int(xx[1])] = [[int(xx[1]), int(xx[4]), 200]]
	return Link

def	read_data_CTCF_all(FF, Ch, ChrAccum):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	# The link coordinates and its stength
	Link = {}
	Acc = 0
	for i in x:
		xx = i.split('\t')
		prev = ChrAccum[xx[0]]
		if len(xx) == 6:
			Link[int(xx[1])] = [[prev + int(xx[1]), prev + int(xx[4]), int(xx[5])]]
		else:
			Link[int(xx[1])] = [[prev + int(xx[1]), prev + int(xx[4]), 200]]
	return Link

def	read_RNAseq(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	RNA = {}
	for i in x:
		xx = i.split('\t')
		if xx[0] == Ch:
			v = float(xx[5])
			if v > -1001.0:
				RNA[int(xx[1])] = v
	return RNA

def	read_ovRNAseq(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	ovRNA = {}
	GE = {}
	Sign_ovRNA = {}
	Sign_GE = {}
	for i in x:
		xx = i.split('\t')
		if xx[0] == Ch:
			if xx[3] == "yes":
				#print "yes", xx[5]
				ovRNA[int(xx[1])] = xx[3]
				vv = float(xx[4])
				if vv > 20.0:
					vv = 20.0
				GE[int(xx[1])] = vv
	return [ovRNA, GE]


# Read all the tss of all genes
def	read_tss(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	TSS = {}
	Gene = {}
	for i in x:
		xx = i.split('\t')
		#if xx[0] == Ch:
		if xx[0] == Ch or Ch == "all":
			# An artificial value of 1.0 in order to visualize
			# the TSS in the y = 1 line
			TSS[int(xx[1])] = 1.0
			Gene[int(xx[1])] = xx[4]
	return [TSS, Gene]

def	read_tss_all(FF, ChrAccum):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	TSS = {}
	Gene = {}
	for i in x:
		xx = i.split('\t')
		prev = ChrAccum[xx[0]]
		TSS[int(xx[1])+prev] = 1.0
		Gene[int(xx[1])+prev] = xx[4]
	return [TSS, Gene]

def	read_ovTSS(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	ovTSS = {}
	for i in x:
		xx = i.split('\t')
		if xx[0] == Ch:
			# An artificial value of 1.0 in order to visualize
			# the TSS in the y = 1 line
			if xx[3] == "yes\n":
				v = 1
			else:
				v = 0
			ovTSS[int(xx[1])] = v
	return ovTSS

"""
read_data_val.
FF is the file to be read
Ch is the chromosome to be inspected
Returns three dictionaries:
D: The value (strength or any other quantity).
L: The length of the region / peak.
E: The end location of the region / peak.
The indexes of the dictionary us the starting region

"""
def	read_data_val(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	#D: the strength of the peak
	D = {}
	# L: the length of the peak
	L = {}
	# E: end point of the peak
	E = {}
	P = {}
	for i in x:
		xx = i.split('\t')
		try:
			#prev = ChrAccum[xx[0]]
			prev = 0
			if xx[0] == Ch:
				#print "xx = ", xx
				# For now, consider only the peak starting point as the key
				D[int(xx[1])] = float(xx[3])
				L[int(xx[1])] = float(xx[2]) - float(xx[1])
				E[int(xx[1])] = float(xx[2])

				P[int(xx[1]) + prev] = [int(xx[1]), int(xx[2])]
		except:
			fl = 1

	#return [D, L, E]
	return [D, L, E, P]

"""
read_data_val_all
FF is the file to be read
Reads data from all chromosomes
Returns three dictionaries:
D: The value (strength or any other quantity).
L: The length of the region / peak.
E: The end location of the region / peak.
The index of the dictionary is the starting region

"""
def	read_data_val_all(FF, ChrAccum):
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	#D: the strength of the peak
	D = {}
	# L: the length of the peak
	L = {}
	# E: end point of the peak
	E = {}
	# P: The position of the peak
	P = {}
	for i in x:
		xx = i.split('\t')
		if xx[0] != "chrY":
			"""
			if xx[0] == "chr1":
				prev = 0
			else:
				prev = ChrAccum[xx[0]]
			"""
			prev = ChrAccum[xx[0]]
			#print "xx = ", xx
			# For now, consider only the peak starting point as the key
			D[int(xx[1]) + prev] = float(xx[3])
			L[int(xx[1]) + prev] = float(xx[2]) - float(xx[1])
			# This was the problem on figs. generated before 19.05.2015:
			# The ending of a peak was not calculated considering
			# the new coordinates
			E[int(xx[1]) + prev] = float(xx[2]) + prev
			# The peak location in the genome
			P[int(xx[1]) + prev] = [int(xx[1]), int(xx[2])]

	return [D, L, E, P]


def	read_tracks(FF):
	f = open(FF, "r")
	x = f.readlines()
	f.close()

	Tracks = []
	for i in range(len(x)):
		if x[i][0] != '#':
			xx = x[i].split('\t')
			yy = xx[3].split(',')
			col1 = []
			for j in yy:
				col1.append(int(j))
			yy = xx[8].split(',')
			col2 = []
			for j in yy:
				col2.append(int(j))
			tmp = [xx[0], xx[1], float(xx[2]), tuple(col1), float(xx[4]), xx[5], xx[6], float(xx[7]), tuple(col2), float(xx[9])]
			Tracks.append(tmp)
	return Tracks
		
	

def	read_num_motifs(FF, Ch):
	f = open(FF, "r")
	x = f.readlines()
	f.close()

	Motif = {}
	for i in x:
		xx = i.split('\t')
		if xx[0] == Ch:
			num = len(xx[3].split(','))
			Motif[int(xx[1])] = num
	#K = Motif.keys()
	return Motif




def	chromosome_length(FF):
	Size = {}
	f = open(FF, "r")
	x = f.readlines()
	f.close()
	Acc = {}
	T = 0
	for i in x:
		xx = i.split('\t')
		if xx[0] != "chrY" and xx[0] != "chrM":
			Size[xx[0]] = int(xx[1])
			Acc[xx[0]] = T
			T = T + int(xx[1])
	return [Size, Acc]

def	hashx(P, Coord, W):
	for i in range(W):
		if P <= Coord[i]:
			return i-1
			#return i
	return W-1

"""
Obtains the coordinates of the W points in the circle
"""
def	coordinates_circle(r, W):
	pi = 3.141592
	Loc = {}
	circum = 2.0 * pi
	inc_ang = circum / float(W)
	ang = 0.0
	k = 0
	while ang <= circum:
		x = math.cos(ang) * r + 1.0
		y = math.sin(ang) * r + 1.0
		Loc[k] = [x, y]
		k = k + 1
		ang = ang + inc_ang
	return Loc

def	transform(X, Xfac, W, Xstart, Xend, mn, mx):
	nX = W - (X - mn) * Xfac - Xstart
	return nX

def	transformx(X, Xfac, W, Xstart, Xend, mn, mx):
	nX = (X - mn) * Xfac + Xstart
	return nX

def	transformy(V, Yfac, mn, mx, Ystart, ps):
	#nV = ps - ((V - mn)/(mx - mn)) * Yfac + Ystart
	nV = ps - (V - mn) * Yfac + Ystart
	return nV

"""
def	rotate(ang, x, y):
	x1 = x * math.cos(ang) - y * math.sin(ang)
	y1 = x * math.sin(ang) + y * math.cos(ang)
	return [x1, y1]
"""

# Rotates point P around center point CP
def	rotate(ang, CP, P):
	rP = P[0] - CP[0], P[1] - CP[1]
	rP = (rP[0] * math.cos(ang) - rP[1] * math.sin(ang), rP[0] * math.sin(ang) + rP[1] * math.cos(ang))
	rP = rP[0] + CP[0], rP[1] + CP[1]
	return rP

def	draw_parabola(draw, Parab):
	for i in Parab:
		draw.point((i[0], i[1]), fill = (30, 90, 159))
		draw.point((i[0]+1, i[1]+1), fill = (30, 90, 159))


# This function obtains the parabola that has a as extreme points (x1, y1) and (x2, y2)
def	obtain_parabola(x1, y1, x2, y2, p, orient, LP):
	Pts = []
	x_ini = min([x1, x2])
	x_end = max([x1, x2])
	y_ini = min([y1, y2])
	y_end = max([y1, y2])
	xm = x_ini + abs(x1-x2)/2.0
	ym = y_ini
	x = x_ini
	ss = ((x1 - x2)**2.0 + (y1 - y2)**2)**0.5
	###p = (LP - (x_end - x_ini))/5.0
	###p = (0.2 - (x_end - x_ini))/30.0
	##p = ((x_end - x_ini))/10.0
	##p = (0.1 - abs(x_end - x_ini))/100.0
	p = ((3.0*ss))/50.0
	# Almost:
	#p = 0.0056 - (ss**0.5)/50.0
	##p = (0.1 - abs(x_end - x_ini))/100.0
	###p = (x_end - x_ini)/20.0
	###p = (1.01 - (ss / LP) + 0.00001) / 1200.0
	print "pp = ", p, x_end - x_ini, x_ini, x_end, ss, LP
	# yf: the size of the parabola in terms of y
	if y1 > 1.0:
		fac = -1.0
	else:
		fac = 1.0
	if orient == 'v':
		if y_ini < 1.0:
			fv = -1.0
			yf = abs(ym +  ((x_ini - xm) ** 2 / (4.0 * p ) + ym) )
		else:
			yf = abs( ((x_ini - xm) ** 2 / (4.0 * p ) + ym) - (0.0 + ym))
			fv = 1.0
	while x <= x_end + 0.00001:
	#while x <= x_end + 0.0011:
		if orient == 'h':
			if y_ini < 1.0:
				fv = 1.0
				yf = abs( ((x_ini - xm) ** 2 / (4.0 * p ) + ym) - (0.0 + ym))
			else:
				fv = -1.0
				yf = abs(ym +  ((x_ini - xm) ** 2 / (4.0 * p ) + ym) )
			#y = fac * math.sqrt(4.0 * p * (x - x_ini)) + ym
			y = fv * ( (x - xm) ** 2 / (4.0 * p)  + ym - yf)
		else:
			y = fv * ( (x - xm) ** 2 / (4.0 * p)  + ym - yf)
		Pts.append([x, y])
		#x = x + 0.0001
		x = x + 0.0001
	"""
	cp = x_ini + abs(x1-x2)/2.0, y_ini + abs(y1-y2)/2.0
	Pts.append([x1, y1])
	Pts.append([cp[0], y_max-0.1])
	Pts.append([x2, y2])
	"""
	return Pts

def	funccol(V):
	"""
	if V > 800:
		red = 150 + (V - 800)
		green = 0
		blue = 0
	else:
		if V > 600:
			red = 0
			green = 150 + (V - 600)
			blue = 0
		else:
			red = 0
			green = 0
			blue = 150 + (200 - V)
	"""

	"""
	red = 226
	green = 166
	blue = 25
	"""

	if V == 0:
		red = 254
		green = 0
		blue = 0
	else:
		if V == 1:
			red = 0
			green = 0
			blue = 250
		else:
			if V == 2:
				red = 154
				green = 48
				blue = 202
			else:
				if V == 3:
					red = 216
					green = 191
					blue = 216

	return (red, green, blue)

def	distance(x1, y1, x2, y2):
	dst = abs(x1-x2) + abs(y1-y2)
	return dst

def	color_chromosome(c):
	"""
	if c != 'chrX':
		n = int(c[3:5])
		if n < 6:
			col = (42*n, 2*n, n)
		else:
			if n < 12:
				col = (n, 42*n, 2*n)
			else:
				col = (2*n, 10 + n, 42*n)
	else:
		col = (100, 120, 49)
	"""
	r = int(random.random()*255)
	g = int(random.random()*255)
	b = int(random.random()*255)
	col = (r, g, b)
	return col
	
"""
Plots in a circular fashion some variable associated to each window
Dst is the distance to the circle
PT: plot type. 0 circular, 1 linear
WS: Window size for echa track
plAvg: Plot the averages for each chromosome and whole genome?
NormV: Normalization values in the form of [min, max]
"""
def	plot_data(Xstart, Xend, Ystart, Yend, imp, draw, nW, Val, Dst, col, r, yfac, PT, WS, plAvg, Info, ChrBorders, NormV, Label):
	W = nW
	H = nW
	#print "yfac = ", yfac
	#print "Dst = ", Dst
	#print "r = ", r
	fontv = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 46)
	fontu = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 26)
	fontx = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 78)
	fontc = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 46)
	if PT == 0:
		# Circular plot
		Diam = 2.0
		Xfac = float(W - Xstart - Xend) / Diam
		Yfac = float(H - Ystart - Yend) / Diam
		#print "Xfac = ", Xfac
		#print "Yfac = ", Yfac
		#Xfac = 0.31
		#Yfac = 0.31
		circ = 2.0 * 3.141592
		rx = r + Dst
		#print "rx = ", rx
		angfac = circ / nW
		for i in range(nW-1):
			ang1 = float(i) * angfac
			#ang2 = float(i+1) * angfac
			ang2 = ang1
			v1 = Val[i]/ yfac
			x1 = math.cos(ang1) * (rx + v1) + 1.0
			y1 = math.sin(ang1) * (rx + v1) + 1.0

			#v2 = Val[i+1] / yfac
			v2 = 0.0
			x2 = math.cos(ang2) * (rx + v2) + 1.0
			y2 = math.sin(ang2) * (rx + v2) + 1.0
			x1r = transformx(x1, Xfac, W, Xstart, Xend, 0.0, 0.0)
			y1r = transform(y1, Yfac, H, Ystart, Yend, 0.0, 0.0)
			x2r = transformx(x2, Xfac, W, Xstart, Xend, 0.0, 0.0)
			y2r = transform(y2, Yfac, H, Ystart, Yend, 0.0, 0.0)
			if i < 0:
				print "pt ", i, " Vali = ", Val[i], " v1 = ", v1, " x1 = ", x1, " y1 = ", y1, " v2 = ", v2, " x2 = ", x2, " y2 = ", y2, " x1r = ", x1r, " y1r = ", y1r, " x2r = ", x2r, " y2r = ", y2r
			if v1 != 0.0:
				draw.line([x1r, y1r, x2r, y2r], fill = col)
				draw.line([x1r, y1r-1, x2r, y2r-1], fill = col)
				draw.line([x1r, y1r-2, x2r, y2r-2], fill = col)

				draw.line([x1r+1, y1r, x2r+1, y2r], fill = col)
				draw.line([x1r-1, y1r, x2r-1, y2r], fill = col)
				draw.line([x1r+1, y1r+1, x2r+1, y2r+1], fill = col)
				draw.line([x1r-1, y1r-1, x2r-1, y2r-1], fill = col)
			else:
				draw.line([x1r, y1r, x2r, y2r], fill = (100,100,100))

			"""
			aa = (x1r + x2r ) / 2.0
			bb = (y1r + y2r ) / 2.0
			draw.line([aa, y1r, aa, y2r], fill = col)
			draw.line([aa, y1r-1, aa, y2r-1], fill = col)
			"""

		# The label. For the meantime, it is better to print it outside the circus
		"""
		y = math.sin(3.141592/2.0) * (rx + 0.02) + 1.0
		yr = transform(y, Yfac, H, Ystart, Yend, 0.0, 0.0)
		x = W / 2.0 - len(Label)/2.0 - Xstart
		draw.text((10, yr), Label, fill = col, font = fontu)
		"""
	else:
		# Linear plot
		if NormV == []:
			mx = -100000.0
			mn = 100000.0
			for k in range(nW):
				if Val[k] > mx:
					mx = Val[k]/(yfac/1.0)
				if Val[k] < mn:
					mn = Val[k]/(yfac/1.0)
		else:
			mn = NormV[0]
			mx = NormV[1]
		# WS is the y-size of the track
		#WS = 320
		Yfac = float(WS) / (mx-mn)
		#Yfac = float(W - Xstart - Xend) / (mx-mn)
		#rx: the starting y-location
		#rx = Ystart + r*100.0 + Dst*2200.0
		#rx = Ystart + r*1.0 + Dst*WS*6
		rx = Dst
		#print "Yfac = ", Yfac, mn, mx, rx
		Pp = int(H/2)
		for i in range(nW-1):
			# Plot the data
			v1 = Val[i]/ (yfac/1.0)
			v2 = Val[i+1]/ (yfac/1.0)
			#v1 = Val[i]
			#v2 = Val[i+1]
			y1r = transformy(v1, Yfac, mn, mx, Ystart, rx)
			y2r = transformy(v2, Yfac, mn, mx, Ystart, rx)
			#print v1, v2, y1r, y2r
			draw.line([Xstart + i, y1r, Xstart + i + 1, y2r], fill = col)
			draw.line([Xstart + i, y1r-1, Xstart + i + 1, y2r-1], fill = col)
		# Min and max
		#mm = mn + (mn + mx)/2.0
		mm = mn + (mx - mn)/2.0
		ymin = transformy(mn, Yfac, mn, mx, Ystart, rx)
		ymax = transformy(mx, Yfac, mn, mx, Ystart, rx) + 55
		# For log data, ymm = 0.0
		logD = 1
		if logD == 1:
			mm = 0.0
			# plot lines at y = 1 and y = -1
			y1 = transformy(1.0, Yfac, mn, mx, Ystart, rx)
			#draw.line([Xstart + 0, y1, W, y1], fill = 'black')
			#draw.line([Xstart + 0, y1-1, W, y1-1], fill = 'black')
			#draw.text((Xstart - 75, y1-20), "1", fill = 'black', font = fontv)
			# y = -1
			y_1 = transformy(-1.0, Yfac, mn, mx, Ystart, rx)
			draw.line([Xstart + 0, y_1, W, y_1], fill = 'black')
			draw.line([Xstart + 0, y_1-1, W, y_1-1], fill = 'black')
			#draw.text((Xstart - 75, y_1-20), "-1", fill = 'black', font = fontv)
		# The min 
		ymm = transformy(mm, Yfac, mn, mx, Ystart, rx) + 25
		# Draw the min.
		draw.text((Xstart - 85, ymin), str(mn)[0:4], fill = 'black', font = fontv)
		draw.line([Xstart + 0, ymin, Xstart + 9, ymin], fill = 'black')
		draw.line([Xstart + 0, ymin-1, Xstart + 9, ymin-1], fill = 'black')
		draw.line([Xstart + 0, ymin-2, Xstart + 9, ymin-2], fill = 'black')

		# The max
		draw.text((Xstart - 85, ymax), str(mx)[0:4], fill = 'black', font = fontv)
		draw.line([Xstart + 0, ymax, Xstart + 9, ymax], fill = 'black')
		draw.line([Xstart + 0, ymax+1, Xstart + 9, ymax+1], fill = 'black')
		draw.line([Xstart + 0, ymax+2, Xstart + 9, ymax+2], fill = 'black')

		# The mean
		#draw.text((Xstart - 75, ymm-45), str(mm)[0:2], fill = 'black', font = fontv)
		#draw.line([Xstart + 1, ymm, Xstart + 9, ymm], fill = 'black')
		#draw.line([Xstart + 1, ymm-1, Xstart + 9, ymm-1], fill = 'black')
		#draw.line([Xstart + 1, ymm+1, Xstart + 9, ymm+1], fill = 'black')

		# The y-axis
		#draw.line([Xstart - 2, ymin, Xstart - 2, ymax], fill = 'black')
		#draw.line([Xstart - 3, ymin, Xstart - 3, ymax], fill = 'black')
		#draw.line([Xstart - 4, ymin, Xstart - 4, ymax], fill = 'black')
		# The label
		#draw.text((nW/2, ymax+55), Label, fill = 'black', font = fontc)
		#draw.text((nW/2, ymax+55), Label, fill = 'black', font = fontx)
		if plAvg:
			# plot averages for each chromosome
			K = ChrBorders.keys()
			for i in K:
				y = Info[i]
				yr = transformy(y, Yfac, mn, mx, Ystart, rx) + 25
				ccol = color_chromosome(i)
				draw.line([Xstart + ChrBorders[i][0], yr, Xstart + ChrBorders[i][1], yr], fill = ccol)
				draw.line([Xstart + ChrBorders[i][0], yr+1, Xstart + ChrBorders[i][1], yr+1], fill = ccol)
				draw.line([Xstart + ChrBorders[i][0], yr+2, Xstart + ChrBorders[i][1], yr+2], fill = ccol)
		# plot the borders
		K = ChrBorders.keys()
		for i in K:
			draw.line([Xstart + ChrBorders[i][0], ymin, Xstart + ChrBorders[i][0], ymin-25], fill = 'black')
			draw.line([Xstart + ChrBorders[i][0] + 1, ymin, Xstart + ChrBorders[i][0] + 1, ymin-25], fill = 'black')
			txt = i[3:5]
			xxp = ChrBorders[i][0]  + (ChrBorders[i][1] - ChrBorders[i][0])/2.0 - 3
			#draw.text((Xstart + xxp, ymin - 25), txt, fill = 'black', font=fontc)

		#draw.rectangle([Xstart - 1, rx, W - 1, ymin], outline = 'black')


def	plot_labels(draw, NameTracks, Colors):
	y0 = 500
	x0 = 30
	fontv = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 18)
	for i in range(len(NameTracks)):
		draw.text((x0, y0), NameTracks[i], fill = Colors[i], font = fontv)
		draw.line([x0-3, y0+6, x0-20, y0+6], fill = Colors[i])
		draw.line([x0-3, y0+7, x0-20, y0+7], fill = Colors[i])
		y0 = y0 - 23

"""
plot_HS. This function plots a box surrounding a certain region (window). That region
can be a hot spot, for example.
"""
def	plot_HS(Xstart, Xend, Ystart, Yend, draw, nW, HS, Loc, Coord, Pos, Pos_i, Diam, r, Dst, col):
	#Diam = 2.0
	W = nW
	H = nW
	Xfac = float(W - Xstart - Xend) / Diam
	Yfac = float(H - Ystart - Yend) / Diam
	circ = 2.0 * 3.141592
	rx = r + Dst
	angfac = circ / nW
	#print "angfac = ", angfac
	fontv = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 22)
	for i in range(len(HS)):
		#"""
		ang1 = float(HS[i]-1) * angfac
		#print "i = ", i, HS[i], ang1
		v1 = 1.0/300.0
		x1 = math.cos(ang1) * (rx + v1) + 1.0
		y1 = math.sin(ang1) * (rx + v1) + 1.0
		x2 = math.cos(ang1) * (r + v1) + 1.0
		y2 = math.sin(ang1) * (r + v1) + 1.0
		x1r = transformx(x1, Xfac, W, Xstart, Xend, 0.0, 0.0)
		y1r = transform(y1, Yfac, H, Ystart, Yend, 0.0, 0.0)
		x2r = transformx(x2, Xfac, W, Xstart, Xend, 0.0, 0.0)
		y2r = transform(y2, Yfac, H, Ystart, Yend, 0.0, 0.0)
		draw.line([x1r, y1r-1, x2r, y2r-1], fill = col)
		#"""
		"""
		x1r = Pos[HS[i]-1][0]
		y1r = Pos[HS[i]-1][1]
		x2r = Pos[HS[i]+1][0]
		y2r = Pos[HS[i]+1][1]
		draw.line([x1r, y1r-1, x2r, y2r-1], fill = col)
		draw.ellipse([x1r-10, y1r-10, x2r+10, y2r+10], fill = col)
		print "xy = ", x1r, y1r, x2r, y2r
		"""
		#"""
		ang2 = float(HS[i]+1) * angfac
		v2 = 1.0/300.0
		a1 = math.cos(ang2) * (rx + v1) + 1.0
		b1 = math.sin(ang2) * (rx + v1) + 1.0
		a2 = math.cos(ang2) * (r + v1) + 1.0
		b2 = math.sin(ang2) * (r + v1) + 1.0
		a1r = transformx(a1, Xfac, W, Xstart, Xend, 0.0, 0.0)
		b1r = transform(b1, Yfac, H, Ystart, Yend, 0.0, 0.0)
		a2r = transformx(a2, Xfac, W, Xstart, Xend, 0.0, 0.0)
		b2r = transform(b2, Yfac, H, Ystart, Yend, 0.0, 0.0)
		draw.line([a1r, b1r-1, a2r, b2r-1], fill = col)
		draw.line([x1r, y1r-1, a1r, b1r-1], fill = col)
		draw.line([x2r, y2r-1, a2r, b2r-1], fill = col)
		#"""
		"""
		a1r = Pos[HS[i]-1][0]
		b1r = Pos[HS[i]-1][1]
		a2r = Pos[HS[i]+1][0]
		b2r = Pos[HS[i]+1][1]
		draw.line([a1r, b1r-1, a2r, b2r-1], fill = col)
		#draw.line([x1r, y1r-1, a1r, b1r-1], fill = col)
		#draw.line([x2r, y2r-1, a2r, b2r-1], fill = col)
		"""

def	too_close(a, L):
	for i in L:
		if abs(a-i) < 4:
			return 1
	return -1

"""
plot_gene
Display the name of the genes in the specified list HS
"""
def	plot_gene(Xstart, Ystart, Xend, Yend, Diam, r, Dst, nW, draw, Gene, HS):
	W = nW
	H = nW
	Xfac = float(W - Xstart - Xend) / Diam
	Yfac = float(H - Ystart - Yend) / Diam
	circ = 2.0 * 3.141592
	angfac = circ / nW
	#print "angfac = ", angfac
	rx = r + Dst
	fontv = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 22)
	Proc = []
	for i in HS:
		if too_close(i, Proc) != 1:
			# To avoid crowded names, just print one gene from the first window,
			# and ignore genes from closeby  windows
			ang1 = float(i) * angfac
			TSS = Gene[i]
			mn = min([1, len(TSS)])
			x1 = math.cos(ang1) * rx + 1.0
			y1 = math.sin(ang1) * rx + 1.0
			x1r = transformx(x1, Xfac, W, Xstart, Xend, 0.0, 0.0)
			y1r = transform(y1, Yfac, H, Ystart, Yend, 0.0, 0.0)
			for j in range(mn):
				draw.text((x1r, y1r), TSS[j], fill = 'black', font = fontv)
				if ang1 > 3.1415:
					y1r = y1r - 12
				else:
					y1r = y1r + 12
		Proc.append(i)

"""
plot_circle plots the circle representing the chromosome(s)
"""
def	plot_circle(Xstart, Xend, Ystart, Yend, draw, nW, radio, Ty, kfac, Label):
	W = nW
	H = nW
	# The diameter of the circle
	Diam = 2.0
	Xfac = float(W - Xstart - Xend) / Diam
	Yfac = float(H - Ystart - Yend) / Diam
	# The position in the display space
	Pos = {}
	# The position in the (0,0), (2,2) space
	Pos_i = {}
	# Plot the points representing chromosomal regions
	for i in range(nW):
		x = transformx(Loc[i][0], Xfac, W, Xstart, Xend, 0.0, 0.0)
		y = transform(Loc[i][1], Yfac, H, Ystart, Yend, 0.0, 0.0)
		Pos[i] = [x, y]
		Pos_i[i] = [Loc[i][0], Loc[i][1]]
		# Draw the position corresponding to the window
		#draw.ellipse((x-radio, y-radio, x+radio, y+radio), fill = (30, 130, 109))

	font3 = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 40)
	# The title
	draw.text((Xstart, 1), Label, fill = 'black', font=font3)

	font3 = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 22)
	font = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 22)

	if Ty == 1:
		# Plot the starting and ending locations
		[xs, ys] = Pos[0]
		#draw.text((xs-100, ys), '1st. window', fill = 'black', font=font)
		draw.text((xs-100, ys), '', fill = 'black', font=font)

		[xe, ye] = Pos[int(nW/2.0)]
		pm = int(kfac * int(nW/2.0))
		#draw.text((xe, ye), 'Center window: ' + str(pm), fill = 'black', font=font)
		draw.text((xe, ye), '', fill = 'black', font=font)
		#Region = [42501114, 44460904]
		"""
		Region = [46813776, 54647096]
		for rR in Region:
			winR = hashx(rR, Coord, W+1)
			print "winR = ", rR, winR
			[xe, ye] = Pos[winR]
			draw.text((xe, ye), str(rR), fill = 'black', font=font)
		"""

		#draw.line([xe, ye-3, xe-100, ye-3], fill = 'black')
		#draw.line([xe, ye-2, xe-100, ye-2], fill = 'black')

	return [Pos, Pos_i]

"""
plot_CTCF plots links between relevant locations
Xstart, Xend, Ystart, Yend: the coordinates of start and end of the canvas
nW: number of points
Links: Coordinates of the relevant positions, in the form:
	Link[loc] = [loci, endi, str]
	loci may be different to loc, endi is the ending location, and str strength is
	the strength of the interaction
Loc: the coordinates of the location in the circle space
Coord:
Ty: the type of the link (parabola, line)
"""
def	plot_CTCF(Xstart, Xend, Ystart, Yend, draw, nW, Links, Loc, Coord, Pos, Pos_i, Diam, Ty, Ch):
	W = nW
	H = nW

	Xfac = float(W - Xstart - Xend) / Diam
	Yfac = float(H - Ystart - Yend) / Diam
	#Xfac = float(W - Xstart - Xend) / 3.5
	#Yfac = float(H - Ystart - Yend) / 3.5
	font1 = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 22)
	font2 = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 18)
	font3 = ImageFont.truetype('/Library/Fonts/Tahoma.ttf', 40)

	#print "drawing TADs..."
	# Draw the links...
	K = Links.keys()

	# The largest parabola...
	LP = -10000.0
	for i in K:
		for j in range(len(Links[i])):
			Start = Links[i][j][0]
			End = Links[i][j][1]
			winS = hashx(float(Start), Coord, W+1)
			winE = hashx(float(End), Coord, W+1)
			[x1, y1] = Pos_i[winS]
			[x2, y2] = Pos_i[winE]
			Ss = ((x1 - x2)**2.0 + (y1 - y2)**2.0)
			if Ss > LP:
				LP = Ss

	# The windows in which the CTCF domains are located
	winCTCF = {}
	for i in K:
		#print "Li = ", Links[i]
		for j in range(len(Links[i])):
			Start = Links[i][j][0]
			End = Links[i][j][1]
			winS = hashx(float(Start), Coord, W+1)
			winE = hashx(float(End), Coord, W+1)
			winCTCF[i] = [winS, winE]
			#print "ws = ", winS, winE, Start, End
			#print "anchor 1 = ", Links[i][j][0]
			#print "anchor 2 = ", Links[i][j][1]
			[xS, yS] = Pos[winS]
			[xE, yE] = Pos[winE]
			# Draw the line linking the start and end locations of the CTCF border
			if Ty == 1:
				draw.line([xS, yS, xE, yE], fill = 'red')
			#draw.arc((int(xS), int(yS), int(xE), int(yE)), 0, 360, fill = 'red')
			else:
				[x1, y1] = Pos_i[winS]
				[x2, y2] = Pos_i[winE]
				#print "xxyy = ", x1, y1, x2, y2
				#cc = sys.stdin.read(1)
				dst = distance(x1, y1, x2, y2)
				"""
				print "p1 = ", x1, y1
				print "p2 = ", x2, y2
				print "Pos = ", Pos_i[winS], Pos_i[winE]
				print "dst = ", dst, winS, winE
				"""
				if dst > 0.003:
				#if dst > 0.01:
				##if dst > 0.03:
					# Fix the first point, rotate the second one
					if x2 > x1:
						# Fix the first point
						m = (y2 - y1) / (x2 - x1)
						ang = -math.atan(m)
						[x2a, y2a]= rotate(ang, (x1, y1), (x2, y2))
						x1a = x1
						y1a = y1
						rpx = x2a
						rpy = y2a
						FixP = (x1a, y1a)
					else:
						# Fix the second point
						m = (y1 - y2) / (x1 - x2)
						ang = -math.atan(m)
						[x1a, y1a] = rotate(ang, (x2, y2), (x1, y1))
						x2a = x2
						y2a = y2
						rpx = x1a
						rpy = y1a
						FixP = (x2a, y2a)
					#p = 0.1
					if dst > 1.5:
						p = 0.5
						print "dst = ", dst
					else:
						if Ch == "all" or Ch == "All":
							p = 0.0005
						else:
							p = 0.003
					#print "Obtaining parabola ", p
					# The slope the orientation
					if abs(m) < 8.2:
					#if abs(m) < 5.2:
						orient = 'v'
					else:
						orient = 'h'
					orient = 'v'
					Parab = obtain_parabola(x1a, y1a, x2a, y2a, p, orient, LP)
					#print "Parab = ", Parab
					# Rotate again the parabola
					rotParab = []
					rotParabPlot = []
					for k in Parab:
						nP = rotate(-ang, FixP, k)
						#nP = k
						rotParab.append(nP)
						nx = transformx(nP[0], Xfac, W, Xstart, Xend, 0.0, 0.0)
						ny = transform(nP[1], Yfac, W, Ystart, Yend, 0.0, 0.0)
						rotParabPlot.append([nx, ny])
					#draw_parabola(draw, rotParabPlot)
					#for k in range(len(rotParabPlot)-1):
					for k in range(3,len(rotParabPlot)-3):
						#draw.point((int(j[0]), int(j[1])), fill = (30, 90, 159))
						col = funccol(Links[i][j][2])
						#print "pr... = ", rotParabPlot[k], col
						draw.line([rotParabPlot[k][0], rotParabPlot[k][1], rotParabPlot[k+1][0], rotParabPlot[k+1][1]], fill = col)
						draw.line([rotParabPlot[k][0]+1, rotParabPlot[k][1], rotParabPlot[k+1][0]+1, rotParabPlot[k+1][1]], fill = col)
						draw.line([rotParabPlot[k][0], rotParabPlot[k][1]+1, rotParabPlot[k+1][0], rotParabPlot[k+1][1]+1], fill = col)
						draw.line([rotParabPlot[k][0], rotParabPlot[k][1]+2, rotParabPlot[k+1][0], rotParabPlot[k+1][1]+2], fill = col)

						draw.line([rotParabPlot[k][0]-1, rotParabPlot[k][1], rotParabPlot[k+1][0]-1, rotParabPlot[k+1][1]], fill = col)
						draw.line([rotParabPlot[k][0], rotParabPlot[k][1]-1, rotParabPlot[k+1][0], rotParabPlot[k+1][1]-1], fill = col)
						draw.line([rotParabPlot[k][0], rotParabPlot[k][1]+1, rotParabPlot[k+1][0], rotParabPlot[k+1][1]+1], fill = col)
						draw.line([rotParabPlot[k][0]+1, rotParabPlot[k][1]+1, rotParabPlot[k+1][0]+1, rotParabPlot[k+1][1]+1], fill = col)
						draw.line([rotParabPlot[k][0]-1, rotParabPlot[k][1]-1, rotParabPlot[k+1][0]-1, rotParabPlot[k+1][1]-1], fill = col)

	return winCTCF

	#return Pos

def	isNaN(v):
	return v != v

def	rearrange(n, DM2):
	Dt = [None] * n
	for i in range(n):
		Dt[i] = [0.0] * n
	s = 0
	for i in range(n):
		for j in range(i + 1, n):
			if i != j:
				if isNaN(DM2[s]):
					DM2[s] = 0.0
				Dt[i][j] = DM2[s]
				Dt[j][i] = DM2[s]
				#print "s = ", s, i, j, DM2[s]
				s = s + 1
				dd = int(Dt[i][j])
	return Dt


def	most_similar(W, Vars, Atts, nA, MS):
	# nA: number of attributes
	nA = len(Atts)
	# The variables
	Vars = zip(*Atts)
	# Distance matrix
	#print "computing distance"
	DMt = spatial.distance.pdist(np.asarray(Vars), 'euclidean')
	DM = rearrange(W, DMt)
	# Detect the regions with the highest absolute values...
	Sel = []
	for i in range(W):
		for j in range(i + 1, W):
			# Verify that the window has values different than zero in at least
			# a subset of the variables
			n1 = 0
			for k in range(nA):
				if Atts[k][i] != 0.0:
					n1 = n1 + 1
			if n1 > 6:
				# At least five variables have a value different than zero
				Sel.append([DM[i][j], i, j, Vars[i], Vars[j]])

	Sel_sort = sorted(Sel)
	R = Sel_sort[0:MS]
	WinMS = []
	for i in range(len(R)):
		WinMS.append(R[i][1])
		WinMS.append(R[i][2])

	return [R, WinMS]

def	convert(MS, kfac):
	L = {}
	for i in range(len(MS)):
		st = int(MS[i][1] * kfac)
		end = int(MS[i][2] * kfac)
		strength = int(MS[i][0]*100)
		L[st] = [[st, end, strength]]
	return L

def	dist(X, Y):
	DMt = spatial.distance.pdist(np.asarray([X, Y]), 'euclidean')
	return DMt[0]

def	num_zeros(X):
	nz = 0
	for i in X:
		if float(i) == 0.0:
			nz = nz + 1
	return nz
	

"""
hot_spots
# Att: rows are variables (number of tracks), columns are windows (W).
# Win: rows are windows (W), columns are variables.
This function identifies relevant windows (hot spots).
First, the average for each variable is found. Those windows with a significant different
value for that average are marked. Those windows are hotspots for that particular variable.
A global hotspot is a window marked as a hotspots for (almost) all variables.
"""
def	hot_spots(Att, Win, AD):
	# Obtain the average for each variable
	nA = len(Att)
	Avg = [0.0] * nA
	for i in range(nA):
		# nS: number of non-zero windows
		nS = 0.0
		for j in range(len(Att[i])):
			v = abs(Att[i][j])
			Avg[i] = Avg[i] + v
			if v > 0.0:
				nS = nS + 1.0
		if nS > 0.0:
			Avg[i] = Avg[i] / nS
		#Avg[i] = Avg[i] / len(Att[i])

	nw = len(Win)
	# At the begining, all regions are unremarkable for all variables
	Remarkable = [None] * nw
	for i in range(nw):
		Remarkable[i] = [0] * nA
	#tR:the number of variables for which region i has been marked as remarkable
	tR = [0] * nw

	HS = []
	if AD == 0:
		# A first way to detect outliers (hotspots) is to verify attribute by attribute
		# and keep the account of the number of those that are (significativelly) different
		# from the mean. Then if the number of different attributes is higher than
		# a threshold, mark the window as outlier.
		Theta = 1.5
		for i in range(nw):
			for j in range(nA):
				R = abs(Win[i][j]) / abs(Avg[j])
				if R > Theta:
					Remarkable[i][j] = 1
					tR[i] = tR[i] +  1

		# Gamma: the minimum number of variables for which a region has to be marked
		# as remarkable to be selected as a hotspot.
		#Gamma = 5
		Gamma = int(nA/2)
		for i in range(nw):
			if tR[i] > Gamma:
				HS.append(i)
	else:
		Rho = 12.0
		# If region i has a 0 on Omega or more attributes, then discard it.
		Omega = int(nA/2.0) - 7
		#Omega = int(nA/2.0) - 3
		for i in range(nw):
			# nz: the number of attributes equal to zero for region i.
			nz = num_zeros(Win[i])
			if nz < Omega:
				ds = dist(Win[i], Avg)
				if ds > Rho:
					HS.append(i)

	return [Avg, HS]


def	save_data(FF, Vars, Name):
	f = open(FF, "w")
	f.write("#")
	for n in Name:
		f.write(str(n) + "\t")
	f.write("\n")
	for i in Vars:
		for j in i:
			f.write(str(j) + "\t")
		f.write("\n")
	f.close()

def	save_stats(FF, AvgAll, TotAll, NpAll):
	K = AvgAll[0].keys()
	f = open(FF + "_avg.csv", "w")
	for i in K:
		f.write(i + "\t")
		for j in range(len(AvgAll)):
			f.write(str(AvgAll[j][i]) + "\t")
		f.write("\n")
	f.close()

	K = NpAll[0].keys()
	f = open(FF + "_dens.csv", "w")
	for i in K:
		f.write(i + "\t")
		for j in range(len(NpAll)):
			f.write(str(NpAll[j][i]) + "\t")
		f.write("\n")
	f.close()


def	previous_chromosome(a):
	if a == 'chr1':
		return "None"
	else:
		if a == 'chrX':
			return 'chr22'
		else:
			n = int(a[3:5])
			prev = 'chr' + str(n-1)
			return prev

def	next_chromosome(a):
	if a == 'chrX':
		return "none"
	else:
		if a == 'chr22':
			return 'chrX'
		else:
			n = int(a[3:5])
			prev = 'chr' + str(n+1)
			return prev

def	chromosome_borders(ChrAccum, Coord, W):
	ChrBorders = {}
	ChrStart = {}
	K = ChrAccum.keys()
	#print "K = ", K, len(K)
	#print "Coord = ", Coord
	for i in K:
		#print "i = ", i, ChrAccum[i]
		win = hashx(ChrAccum[i], Coord, W)
		#win = hashx(ChrAccum[i], Coord, W+1)
		#print "win = ", win
		ChrStart[i] = win
	for i in K:
		if i != 'chrX':
			nc = next_chromosome(i)
			ChrBorders[i] = [ChrStart[i], ChrStart[nc]-1]
		else:
			ChrBorders[i] = [ChrStart[i], W]
	return ChrBorders

		

def	which_chromosome(a, K, ChrBorders, W):
	for i in K:
		if a >= ChrBorders[i][0] and a <= ChrBorders[i][1]:
			return i
	return 'chrX'

"""
chromosome_stats obtains some statistics for each chromosome
"""
def	chromosome_stats(ChrBorders, ValX, NumPeaksX, W):
	# The number of attributes (tracks)
	N = len(ValX)
	NAll = []
	AvgAll = []
	TotAll = []
	K = ChrBorders.keys()
	for i in range(N):
		StChr = {}
		Np = {}
		Avg = {}
		Tot = {}
		for j in K:
			Np[j] = 0
			Avg[j] = 0
			Tot[j] = 0
		for j in range(W):
			ch = which_chromosome(j, K, ChrBorders, W)
			Np[ch] = Np[ch] + NumPeaksX[i][j]
			Tot[ch] = Tot[ch] + ValX[i][j]
		for j in K:
			if Np[j] > 0:
				Avg[j] = Tot[j] / Np[j]
		AvgAll.append(Avg)
		TotAll.append(Tot)
		NAll.append(Np)

	return [AvgAll, TotAll, NAll]

def	min_max(AvgAll):
	mn = 100000.0
	mx = -100000.0
	#print "ln = ", len(AvgAll)
	for i in range(len(AvgAll)):
		K = AvgAll[i].keys()
		#print "K = ", K
		for j in K:
			if AvgAll[i][j] < mn:
				mn = AvgAll[i][j]
			if AvgAll[i][j] > mx:
				mx = AvgAll[i][j]
	#return [mn, mx]
	return [mn, mx+0.01]
	
def	detect_outliers(avg_FC_V, ValV, NumPeaksV, ChrBorders, kfac, W, WhatPeaks, Theta):
	Outliers = {}
	s = 0
	K = ChrBorders.keys()
	for i in range(W):
		if avg_FC_V[i] > Theta or avg_FC_V[i] < -Theta:
			# The value in the window exceeds the threshold,
			# it is an outlier candidate
			for j in K:
				if i >= ChrBorders[j][0] and i <= ChrBorders[j][1]:
					# The window is in chromosome j
					cont = 1
					t = 0
					while cont == 1:
						Delta = (i - t) - ChrBorders[j][0]
						Loc = Delta * kfac
						if Loc <= WhatPeaks[i][0][4][1]:
							cont = 0
						else:
							t = t + 1
					Outliers[s] = [j, Delta, int(Loc), int(Loc + kfac), avg_FC_V[i], ValV[i], NumPeaksV[i], i, WhatPeaks[i]]
			s = s + 1
	return Outliers

def	save_outliers(FF, Outliers):
	K = Outliers.keys()
	f = open(FF, "w")
	for i in K:
		#print "i = ", i, Outliers[i]
		f.write(str(i) + "\t")
		L = Outliers[i]
		for j in range(len(L)-1):
			f.write(str(L[j]) + "\t")
		f.write(str(L[len(L)-1]) + "\n")
	f.close()

def	save_CTCF(Links, winCTCF, FF):
	K = Links.keys()
	f = open(FF, "w")
	for i in K:
		f.write(str(i) + "\t" + str(Links[i][0][1]) + "\t" + str(winCTCF[i][0]) + "\t" + str(winCTCF[i][1]) + "\n")
	f.close()

"""
python/Users/antonioneme/analysis/FAIRE/progs/interactions_visualization.py -i all_CTCF_sign_2h.csv -chr chr1 -RNA TS_all_2h.csv ... -out all_CTCF_sign_2h_chr1.png
The input file is of the form:
chr1	805136	805426	-1	-1	-1
chr1	875413	876476	839813	999431	200
...
chr1	1057502	1058838	1057205	1168373	200	1057403	1227872	300	913753

The first three columns indicate the location of the peak (or region of
interest). For the next data, each group of three columns indicates the
interacting regions and the strength of it.

-RNA is teh file containing RNA-seq data

Much more tracks can be specified... later.
"""

parser = argparse.ArgumentParser()
parser.add_argument('-i', action="store", dest = "d", help = "Input data")
parser.add_argument('-chr', action="store", dest = "ch", help = "The chromosome to be inspected")
parser.add_argument('-out', action="store", dest = "out", help = "The output file")
parser.add_argument('-stats', action="store", dest = "stats", help = "The output file")
parser.add_argument('-img', action="store", dest = "img", help = "The image")
parser.add_argument('-FAIRE', action="store", dest = "FAIRE", help = "Plot the average FC for FAIRE-seq peaks?")
parser.add_argument('-VDR', action="store", dest = "VDR", help = "Plot the average FC for VDR ChIP-seq peaks?")
parser.add_argument('-motifs', action="store", dest = "motifs", help = "The file containing the tracks for motifs")
parser.add_argument('-tss', action="store_true", dest = "tss", help = "Plot the number TSS in the window")
parser.add_argument('-ovtss', action="store", dest = "ovtss", help = "Plot the number of FAIRE-seq overlapping with TSS?")
parser.add_argument('-indep', action="store", dest = "indep", help = "The file containing the tracks of independent datasets")
parser.add_argument('-res', action="store", dest = "res", help = "The number of windows (pixels)")
parser.add_argument('-outliers', action="store", dest = "outliers", help = "The output file containing information about the outliers")
parser.add_argument('-CTCFloc', action="store", dest = "CTCFloc", help = "The output file containing the windows for the CTCF anchor sites")


args = parser.parse_args()



[Chrm, ChrAccum] = chromosome_length("/Users/antne/data/GENOME/chr_size.csv")
#print "ChrAccum = ", ChrAccum
#cc = sys.stdin.read(1)

if args.ch == "all" or args.ch == "All":
	Links = read_data_CTCF_all(args.d, args.ch, ChrAccum)
else:
	Links = read_data_CTCF(args.d, args.ch)

#print "Links = ", Links
#cc = sys.stdin.read(1)

# S: Chromosome size
if args.ch == "all" or args.ch == "All":
	S = 0
	K = Chrm.keys()
	for i in K:
		S = S + Chrm[i]
	All = True
else:
	S = float(Chrm[args.ch])
	All = False
#print "S = ", S

# mW: the amximum number fo windows
#mW = 7000
#mW = 5500
#mW = 2500
#mW = 1500
mW = int(args.res)
#mW = int(args.res) + 90 + 50 + 140

# For the same resolution:
#W = int(float(mW) * (S/float(249250621)) )

# For the same number of pixels:
W = mW
Wx = mW
Hx = mW

# for the same resolution:
#kfac = 249250621.0/float(mW)
# for the same number of pixels:
#kfac = S/float(mW)
# Used until 28.05.2015
kfac = S/float(Wx - 100 - 100)

# After the inconsistencies found in the CTCF data, kfac was changed
kfac = S/float(Wx)
print "kfac = ", kfac

#print "kfac = ", kfac
#print "W = ", W
#cc = sys.stdin.read(1)

Coord = {}

start = 0
for i in range(W+1):
	Coord[i] = start
	start = start + int(kfac)
#print "Coord = ", Coord[0], Coord[1], Coord[2], Coord[3], Coord[4], Coord[W-2], Coord[W-1]
#cc = sys.stdin.read(1)

# The index (starting points) of the peaks (regions of interest)
K = Links.keys()

#print "K = ", K
#print "Lk = ", Links[243507268]
# Loc: the coordinates in the circle of the W points
radius = 1.0
Loc = coordinates_circle(radius, W+1)
Rr = Loc.keys()
#print "coord = ", Rr[11], Loc[Rr[11]]
#print "num. points in circle = ", len(Loc)
#cc = sys.stdin.read(1)


imp = Image.new("RGB", (Wx, Hx), "white")
draw = ImageDraw.Draw(imp)
# PT is the plot type. 0 is circular, 1 is linear
PT = 0
#PT = 1
# This seems to work for nW = 2500
Diam = 2.0
print "MM = ", mW
if PT == 0:
	if mW <= 3500:
		"""
		Xstart = 550
		Xend = 550
		Ystart = 550
		Yend = 550
		"""
		Xstart = 710
		Xend = 710
		Ystart = 710
		Yend = 710
	# This seems to work for nW = 5500
	else:
		if mW == 5001:
			Xstart = 1000
			Xend = 1000
			Ystart = 1000
			Yend = 1000
		else:
			if mW == 1500:
				Xstart = 600
				Ystart = 600
				Xend = 600
				Yend = 600
			else:
				if mW >= 1000 and mW <= 7005:
				#if mW >= 10000:
					Xstart = 840
					Ystart = 840
					Xend = 840
					Yend = 840
				else:
					if mW > 7005:
						Xstart = 3200
						Ystart = 3200
						Xend = 3200
						Yend = 3200
	#[Pos, Pos_i] = plot_circle(Xstart, Xend, Ystart, Yend, draw, W+1, 2.0, 1, kfac, 'CTCF domains ' + args.ch)
	[Pos, Pos_i] = plot_circle(Xstart, Xend, Ystart, Yend, draw, W+1, 2.0, 1, kfac, '')
	#print "Coord = ", Coord
	"""
	print "Pos = ", Pos[10]
	print "Posi = ", Pos_i[10]
	print "Pos = ", Pos[1000]
	print "Posi = ", Pos_i[1000]
	Kkk = Links.keys()
	print "LK = ", sorted(Kkk)
	#print "Lk = ", Links
	#print "Links = ", Links[Kkk[0]], Links[Kkk[300]]
	cc = sys.stdin.read(1)
	"""
	winCTCF = plot_CTCF(Xstart, Xend, Ystart, Yend, draw, W+1, Links, Loc, Coord, Pos, Pos_i, Diam, 0, args.ch)
else:
	Xstart = 90
	Xend = 20
	Ystart = 50
	Yend = 50

save_CTCF(Links, winCTCF, args.CTCFloc)

ChrBorders = chromosome_borders(ChrAccum, Coord, W+1)

# Atts: The attributes for each window
Atts = []
# NameTracks: the name of the tracks
NameTracks = []
Colors = []
Indep = read_tracks(args.indep)

b = 2
# The number of independent tracks
NumTracks = len(Indep)
for track in Indep:
	if track[6] == "yes":
		NumTracks = NumTracks + 1
		b = 1
#NumTracks = NumTracks + b
# Window size for each track
WS = float(W) / NumTracks - 70
#print "NT / WS = ", NumTracks, WS
ValX = []
NumPeaksX = []
Avg_FC_X = []
for track in Indep:
	ValV = {}
	NumPeaksV = {}
	whatPeaks = {}
	# Can also the coverage be consireded?
	if not All:
		# Process only one chromosome
		[D, L, E, Ps] = read_data_val(track[1], args.ch)
		#[D, L, E] = read_data_val(track[1], args.ch)
		#print "HHH"
	else:
		# Process all chromosomes
		print "tr = ", track[1]
		[D, L, E, Ps] = read_data_val_all(track[1], ChrAccum)
	for i in range(W+1):
		ValV[i] = 0.0
		NumPeaksV[i] = 0
		whatPeaks[i] = []
	K = D.keys()
	PW = 0
	TEE = 0
	NS = 0
	print "tr = ", track[1]
	for i in K:
		win = hashx(float(i), Coord, W+1)
		win_end = hashx(float(E[i]), Coord, W+1)
		#if i == 1875102 or i == 2316401 or Ps[i][0] == 141721807:
		#	print "i = ", i, win, win_end, E[i], Ps[i]
		if win != win_end:
			# Since the peak is contained in two windows, it has to be
			# considered in both of them
			NumPeaksV[win] = NumPeaksV[win] + 1
			NumPeaksV[win_end] = NumPeaksV[win_end] + 1
			if All:
				whatPeaks[win].append([D[i], L[i], E[i], i, Ps[i]])
			else:
				#whatPeaks[win].append([D[i], L[i], E[i], i])
				whatPeaks[win].append([D[i], L[i], E[i], i, Ps[i]])
			# However, the strenghth contirbution....?
			# perhaps is to be considered equally in both windows
			ValV[win] = ValV[win] + D[i]
			ValV[win_end] = ValV[win_end] + D[i]
			if All:
				whatPeaks[win_end].append([D[i], L[i], E[i], i, Ps[i]])
			else:
				whatPeaks[win_end].append([D[i], L[i], E[i], i, Ps[i]])
				#whatPeaks[win_end].append([D[i], L[i], E[i], i])
			NS = NS + 1
		else:
			ValV[win] = ValV[win] + D[i]
			NumPeaksV[win] = NumPeaksV[win] + 1
			if All:
				whatPeaks[win].append([D[i], L[i], E[i], i, Ps[i]])
			else:
				#print "Ps = ", Ps[i]
				whatPeaks[win].append([D[i], L[i], E[i], i, Ps[i]])
				#whatPeaks[win].append([D[i], L[i], E[i], i])
		#print "i = ", i, win, win_end, NumPeaksV[win], NumPeaksV[win_end]
	#cc = sys.stdin.read(1)
	avg_FC_V = {}
	for i in range(W+1):
		if NumPeaksV[i] > 0:
			avg_FC_V[i] = ValV[i] / NumPeaksV[i]
		else:
			avg_FC_V[i] = 0.0

	#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, avg_FC_V, track[2], track[3], radius, track[4], PT, WS, plAvg, [], 'FC ' + track[5])
	Atts.append(avg_FC_V)
	NameTracks.append(track[5])
	#NameTracks.append('FC ' + track[5])
	#NameTracks.append(track[0])
	Colors.append(track[3])
	ValX.append(ValV)
	Avg_FC_X.append(avg_FC_V)

	#if track[6] == "yes":
		#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, NumPeaksV, track[7], track[8], radius, track[9], PT, WS, plAvg, [], 'Dens ' + track[5])
	Atts.append(NumPeaksV)
	NameTracks.append('Dens ' + track[5])
	Colors.append(track[8])
	NumPeaksX.append(NumPeaksV)
	if track[6] == "yes":
		# More than three peaks?
		ThetaT = 5.0
		Outliers1 = detect_outliers(NumPeaksV, ValV, NumPeaksV, ChrBorders, kfac, W+1, whatPeaks, ThetaT)
		#print "t = ", track[1]
		save_outliers(args.outliers + "_dens.csv", Outliers1)
		#save_outliers(args.outliers + "_dens_" + track[1], Outliers1)
	else:
		#ThetaT = 2.0
		ThetaT = 1.05
		Outliers1 = detect_outliers(avg_FC_V, ValV, NumPeaksV, ChrBorders, kfac, W+1, whatPeaks, ThetaT)
		save_outliers(args.outliers + "_FC.csv", Outliers1)
		#save_outliers(args.outliers + "_" + track[1], Outliers1)

# plAvg: display te average value for each chromosome?
#plAvg = True
plAvg = False
[AvgAll, TotAll, NpAll]  = chromosome_stats(ChrBorders, ValX, NumPeaksX, W+1)
#NV: the range for all tracks
NVval = min_max(Avg_FC_X)
#print "Range = ", NVval
NVDens = min_max(NumPeaksX)
s = 0
ss = 0
Y0 = W / NumTracks - Ystart - 15
#print "Y0 ", Y0, NumTracks, Ystart
for track in Indep:
	#print "plotting tr = ", track
	Start_pos = Y0 + ((WS + 8) * ss)
	#print "Avg = ", Avg_FC_X[s]
	#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, Avg_FC_X[s], Start_pos, track[3], radius, track[4], PT, WS, plAvg, AvgAll[s], ChrBorders, [], track[5])
	#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, Avg_FC_X[s], Start_pos, track[3], radius, track[4], PT, WS, plAvg, AvgAll[s], ChrBorders, NVval, track[5])
	# For circular plot:
	plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, Avg_FC_X[s], track[2], track[3], radius, track[4], PT, WS, plAvg, AvgAll[s], ChrBorders, [], track[5])
	ss = ss + 1
	if track[6] == "yes":
		Start_pos = Y0 + ((WS + 8) * ss)
		#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, NumPeaksX[s], Start_pos, track[8], radius, track[9], PT, WS, plAvg, NpAll[s], ChrBorders, [], 'Dens ' + track[5])
		#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, NumPeaksX[s], Start_pos, track[8], radius, track[9], PT, WS, plAvg, NpAll[s], ChrBorders, NVDens, '' + track[5])
		# For circular plot
		plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, NumPeaksX[s], track[7], track[8], radius, track[9], PT, WS, plAvg, NpAll[s], ChrBorders, [], 'Dens ' + track[5])
		ss = ss + 1
	s = s + 1

#"""
Motifs = read_tracks(args.motifs)
#print "Motifs = ", Motifs
Motif_num = []
for track in Motifs:
	Mtf = read_num_motifs(track[1], args.ch)
	Mtf_num = {}
	for i in range(W+1):
		Mtf_num[i] = 0.0
	K = Mtf.keys()
	for i in K:
		win = hashx(float(i), Coord, W+1)
		Mtf_num[win] = Mtf_num[win] + Mtf[i]
	#plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, Mtf_num, track[2], track[3], radius, track[4], PT, WS, plAvg, [], ChrBorders, [], 'Dens. ' + track[5])
	Atts.append(Mtf_num)
	#NameTracks.append(track[1])
	NameTracks.append(track[5])
	Colors.append(track[3])
	Motif_num.append(Mtf_num)
#"""
#"""
if args.tss == True:
	#[TSS, Gene] = read_tss('/Users/antonioneme/data/GENOME/xTSS.csv', args.ch)
	if args.ch != "all" and args.ch != "All":
		[TSS, Gene] = read_tss('/Users/antne/data/GENOME/xTSS_new_filtered.csv', args.ch)
	else:
		[TSS, Gene] = read_tss_all('/Users/antne/data/GENOME/xTSS_new_filtered.csv', ChrAccum)
	# There can be several TSS per window
	# TSS_num: number of TSS in the window
	TSS_num = {}
	Genes_in_Win = {}
	for i in range(W+1):
		TSS_num[i] = 0
		Genes_in_Win[i] = []
	for i in TSS:
		win = hashx(float(i), Coord, W+1)
		TSS_num[win] = TSS_num[win] + 1
		#print "i = ", i, win, TSS_num[win]
		#print "g = ", Gene[i]
		Genes_in_Win[win].append(Gene[i])

	plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, TSS_num, -0.13, (99, 100, 101), radius, 950.0, PT, WS, plAvg, [], ChrBorders, [], 'Num TSS')
	Atts.append(TSS_num)
	NameTracks.append('Num. TSS')
	Colors.append((99, 100, 101))

# Mark the chromosome borders...
CB = ChrAccum.keys()
print "CB = ", CB
print "ChrAccum = ", ChrAccum
ChrBorder = {}
for i in range(W+1):
	ChrBorder[i] = 0
for ch in CB:
	win = hashx(float(ChrAccum[ch]), Coord, W+1)
	print "winB = ", win
	if win < 0:
		win = 0
	ChrBorder[win] = ChrBorder[win] + 1
plot_data(Xstart, Xend, Ystart, Yend, imp, draw, W+1, ChrBorder, -0.001, (30, 25, 19), radius, 1.4, PT, WS, plAvg, [], ChrBorders, [], 'Chr borders')
	


##plot_labels(draw, NameTracks, Colors)

# nA: number of attributes
nA = len(Atts)
print "nA = ", nA
print "Tracks = ", NameTracks
At = []
for i in range(len(Atts)):
	tmp = []
	for j in range(W+1):
		tmp.append(Atts[i][j])
	At.append(tmp)
# At: rows are variables (number of tracks), columns  are windows (W+1).
# Vars: rows are windows (W+1), columns are the variables (number of tracks).
Vars = zip(*At)

# Identify hot spots ...
print "Number of attributes = ", nA
print "Number of windows (regions) = ", W+1

# Hot spots are anomalies. How to detect them?
AD = 1
[AvgVar, HS] = hot_spots(At, Vars, AD)
#print "HS = ", HS, len(HS), AvgVar
"""
for i in range(len(HS)):
	print Vars[HS[i]]
"""
mxR = 0.76 + 0.05
#plot_HS(Xstart, Xend, Ystart, Yend, draw, W+1, [2, 100, 200, 400, 900, 1200, 2400], Loc, Coord, Pos, Pos_i, Diam, radius-0.02, mxR, (0, 250, 10))
#plot_HS(Xstart, Xend, Ystart, Yend, draw, W+1, [278, 286, 2102, 2108, 900], Loc, Coord, Pos, Pos_i, Diam, radius-0.02, mxR, (0, 250, 10))
#plot_HS(Xstart, Xend, Ystart, Yend, draw, W+1, [1, 4000], Loc, Coord, Pos, Pos_i, Diam, radius-0.02, mxR, (0, 250, 10))

# Plot hot spots
if PT == 0:
	plot_HS(Xstart, Xend, Ystart, Yend, draw, W+1, HS, Loc, Coord, Pos, Pos_i, Diam, radius-0.02, mxR, (250, 0, 0))
	plot_gene(Xstart, Ystart, Xend, Yend, Diam, radius, -0.04, W+1, draw, Genes_in_Win, HS)

# Plot the genes in the hot spots HS
#plot_gene(Xstart, Ystart, Xend, Yend, Diam, radius, -0.04, W+1, draw, Genes_in_Win, HS)

# MS: the most similar windows, containing:
# [[dist, wini, winj, attsi, attsj], ...]
# where dist is the distance between attsi and attsj of wini and winj
# WinMS is the list containig all selected windows

"""
[MS, WinMS] = most_similar(W+1, Vars, At, nA, 5)
Links2 = convert(MS, kfac)
K = Links.keys()

# Plot the links between the similar windows
plot_CTCF(Xstart, Xend, Ystart, Yend, draw, W+1, Links2, Loc, Coord, Pos, Pos_i, Diam, 1)
"""

# Mark the similar windows
# Plot similar windows
#plot_HS(Xstart, Xend, Ystart, Yend, draw, W+1, WinMS, Loc, Coord, Pos, Pos_i, Diam, radius-0.02, mxR, (0, 0, 250))
#plot_gene(Xstart, Ystart, Xend, Yend, Diam, radius, -0.04, W+1, draw, Genes_in_Win, WinMS)

if All:
	#ChrBorders = chromosome_borders(ChrAccum, Coord, W+1)
	#print "Acc = ", ChrAccum
	#print "BRD = ", ChrBorders
	[AvgAll, TotAll, NpAll]  = chromosome_stats(ChrBorders, ValX, NumPeaksX, W+1)
	K = AvgAll[0].keys()
	"""
	T = 0
	for i in K:
		print i, AvgAll[0][i], NpAll[0][i]
		T = T + NpAll[0][i]
	print T
	"""
	if PT == 1:
		"""
		fontc = ImageFont.truetype('/Applications/iPhoto.app/Contents/Resources/Fonts/arial.ttf', 32)
		# Plot chromosome borders
		K = ChrAccum.keys()
		for i in K:
			win = hashx(ChrAccum[i], Coord, W+1)
			#draw.line([win, 1600, win, 1612], fill = 'black')
			#draw.line([Xstart + win, Ystart, Xstart + win, Hx - Yend], fill = 'black')
			draw.line([Xstart + win, 3, Xstart + win, Ystart-10], fill = 'black')
			txt = i[3:5]
			draw.text((Xstart + win + 3, 3), txt, fill = 'black', font=fontc)
		"""
		nothing = 0


# Vars includes all tracks plus tss if checked
save_data(args.out, Vars, NameTracks)
save_stats(args.stats, AvgAll, TotAll, NpAll)
#imp.show()
#imp.save(args.img)
imp.save(args.img, dpi=(5000,5000))
