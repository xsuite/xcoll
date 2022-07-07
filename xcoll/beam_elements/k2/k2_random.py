import numpy as np
from pathlib import Path
from xcoll import _pkg_root

__pseudo_rands__ = np.load(Path(_pkg_root,'beam_elements','k2','randoms.npy'))
__index__ = -1

cgen = np.zeros(200, dtype=np.float64)

def get_random():

    global __index__
    __index__ += 1
    thisran = __pseudo_rands__[__index__]
    return thisran


def get_random_ruth():
    global cgen
    ran = get_random()
    gap=0.010101010091602802 #1./99.
    gapinv=99.
    tleft= 2./99
    bright=97./99.
    gaps=tleft/49.
    gapins=1./gaps
    j = int(ran*gapinv) + 1
    if j < 3:
        j1 = int(ran *gapins)
        j = j1 + 101
        j = max(j,102)
        j = min(j,148)
        p = (ran -gaps*(j1-1)) * gapins
        a = (p+1.0) * cgen[j+1] - (p-2.0)*cgen[j-2]
        b = (p-1.0) * cgen[j-1] - p * cgen[j]
       
    elif j > 97:
        j1 = int((ran-bright)*gapins)
        j = j1 + 151
        j = max(j,152)
        j = min(j,198)
        p = ((ran-bright) - gaps*(j1-1)) * gapins
        a = (p+1.0) * cgen[j+1] - (p-2.0)*cgen[j-2]
        b = (p-1.0) * cgen[j-1] - p * cgen[j]
       
    else:
        p = (ran -gap*(j-1)) * gapinv
        a = (p+1.) * cgen[j+1] - (p-2.)*cgen[j-2]
        b = (p-1.) * cgen[j-1] - p * cgen[j]
    ran = ((a*p)*(p-1.))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5

    return ran




def initialise_random_ruth(material):
    global cgen
    cgen = np.load(Path(_pkg_root,'beam_elements','k2','cgen_' + material + '.npy'))


#func,ifunc,xlow,xhigh,nlo,nbins,tftot



# def rgs56p(f,a,b):

#     r1 = 1. 
#     hf = r1/2.
                                        
#     x5 = [4.6910077030668004e-02,
#     2.3076534494715846e-01,
#     5.0000000000000000e-01,
#     7.6923465505284154e-01,
#     9.5308992296933200e-01]

#     w5 = [1.1846344252809454e-01,
#     2.3931433524968324e-01,
#     2.8444444444444444e-01,
#     2.3931433524968324e-01,
#     1.1846344252809454e-0]
                                    
#     x6 = [3.3765242898423989e-02,
#     1.6939530676686775e-01,
#     3.8069040695840155e-01,
#     6.1930959304159845e-01,
#     8.3060469323313225e-01,
#     9.6623475710157601e-01]

#     w6 = [8.5662246189585178e-02,              
#     1.8038078652406930e-01,         
#     2.3395696728634552e-01,         
#     2.3395696728634552e-01,         
#     1.8038078652406930e-01,         
#     8.5662246189585178e-0]

#     rang = b-a
#     e5 = 0
#     e6 = 0

#     for i in range(1,6):
#         e5 = e5+np.double(w5[i-1]*f(a+rang*x5[i-1]))
#         e6 = e6+np.double(w6[i-1]*f(a+rang*x6[i-1]))

#     e6 = e6+np.double(w6[5]*f(a+rang*x6[5]))
#     res = np.real((np.double(hf)*(e6+e5))*np.double(rang))
#     err = np.real(abs((e6-e5)*np.double(rang)))

#     return res, err




# def funpct(func,ifunc,xlow,xhigh,nlo,nbins,tftot):

# # Array XFCUM is filled from NLO to NLO+NBINS, which makes
# # the number of values NBINS+1, or the number of bins NBINS
#     rteps=0.005
#     nz=10
#     maxz=20
#     nitmax=6
#     precis=1e-6
# #      DOUBLE PRECISION TPCTIL, TZ, TCUM, XINCR, DTABS,
# #     &  TINCR, TZMAX, XBEST, DTBEST, DTPAR2
# #
#     ierr = 0
#     if (tftot <= 0.): 
#         #go to 900
#         ierr = 1
#         return

#     xfcum = np.zeros(200, dtype=np.float64)

#     tpctil = tftot/np.real(nbins)
#     tz = tpctil/np.real(nz)
#     tzmax = tz * 2.
#     xfcum[nlo-1] = xlow
#     xfcum[nlo+nbins-1] = xhigh
#     x = xlow
#     f = func(x)

#     if (f < 0.):
#         #go to 900
#         ierr = 1
#         return ierr
# #   Loop over percentile bins 
   
#     # do 600 
#     for ibin in range(nlo, nlo+nbins-1):
#         tcum = 0.
#         x1 = x
#         f1 = f
#         dxmax = (xhigh -x) / nz
#         fmin = tz/dxmax
#         fminz = fmin
# #         Loop over trapezoids within a supposed percentil

#         # do 500 
#         for iz in range(1, maxz+1):

#             xincr = tz/max(f1,fmin,fminz)
#         #350
#             while True:
#                 x = x1 + xincr
#                 f = func(x)
#                 if (f < 0.):
#                     #go to 900
#                     ierr = 1
#                     return ierr
                
#                 tincr = ((x-x1) * 0.5) * (f+f1)
#                 if (tincr < tzmax):
#                     break
#                 xincr = xincr * 0.5
#         #370 continue
#             tcum = tcum + tincr
#             if (tcum >= tpctil*0.99):
#                 #go to 520
#                 break
#             fminz = (tz*f)/ (tpctil-tcum)
#             f1 = f
#             x1 = x
#         # 500 continue
#         # write(lout,*) ' FUNLUX:  WARNING. FUNPCT fails trapezoid.'
#         #         END OF TRAPEZOID LOOP
#         #         Adjust interval using Gaussian integration with
#         #             Newton corrections since F is the derivative
#     # 520 continue
#         x1 = xfcum[ibin-1]
#         xbest = x
#         dtbest = tpctil
#         tpart = tpctil
#     # Allow for maximum NITMAX more iterations on RADAPT
#         # do 550 
#         for ihome in range(1, nitmax+1):
#         #535
#             while True:
#                 xincr = (tpctil-tpart) / max(f,fmin)
#                 x = xbest + xincr
#                 x2 = x
#                 if (ihome>1 & x2==xbest):
#                     #write(lout,'(A,G12.3)') ' FUNLUX: WARNING from FUNPCT: insufficient precision at X=',x
#                     #go to 580
#                     break
#                 #endif
#                 refx = abs(x)+precis
#                 tpart2, uncert = radapt(func,x1,x2,1,rteps,0)
#                 dtpar2 = tpart2-tpctil
#                 dtabs = abs(dtpar2)
#                 if(abs(xincr)/refx < precis): 
#                     break
#                 if(dtabs < dtbest): 
#                     break
#                 xincr = xincr * 0.5
#                 #goto 535
#         #545 
#             dtbest = dtabs
#             xbest = x
#             tpart = tpart2
#             f = func(x)
#             if(f < 0.): 
#                 #goto 900
#                 ierr = 1
#                 return ierr
#             if(dtabs < rteps*tpctil):
#                 #goto 580
#                 break
#         # 550 continue
#         #write(lout,'(A,I4)') ' FUNLUX: WARNING from FUNPCT: cannot converge, bin',ibin

#     #580 continue
#         xincr = (tpctil-tpart) / max(f,fmin)
#         x = xbest + xincr
#         xfcum[ibin] = x
#         f = func(x)
#         if(f < 0.): 
#             # goto 900
#             ierr = 1
#             return ierr

# #600 continue
#     # END OF LOOP OVER BINS
#     x1 = xfcum[(nlo+nbins)-2]
#     x2 = xhigh
#     tpart, uncert = radapt(func,x1,x2,1,rteps,0)
#     aberr = abs(tpart-tpctil)/tftot
#     # WRITE(6,1001) IFUNC,XLOW,XHIGH
#     # if(aberr > rteps): 
#     #     write(lout,1002) aberr
#     return xfcum

# # #900 write(lout,1000) x,f
# #     ierr = 1
# #     return




# def radapt(f,a,b,nseg,reltol,abstol):

#     ndim=100
#     r1 = 1.
#     hf = r1/2.
#     nseg = 1
#     nsegd = 1
#     te = 0.0 
#     xlo = np.zeros(ndim)
#     xhi = np.zeros(ndim)
#     tval = np.zeros(ndim)
#     ters = np.zeros(ndim)
#     ibig=1
#     bige=0
#     res=0
#     err=0
#     tvals=0
#     root=0

#     #   dimension xlo(ndim),xhi(ndim),tval(ndim),ters(ndim)
#     nter = 0  # data nter /0/
#     flag = 0

#     if nseg <= 0:
#         if nter == 0:
#             nsegd = 1
#             flag = 2 #goto 2
#         if flag == 0:
#             tvals = 0
#             terss = 0
#             #1
#             for i in range(1,nter+1):
#                 tval[i-1], te = rgs56p(f, xlo[i-1], xhi[i-1])
#                 ters[i-1] = te**2
#                 tvals = tvals + tval[i-1]
#                 terss = terss + ters[i-1]
#                 #1  continue
#             root = np.sqrt(2*terss)
#             flag = 9 #goto 9
#     if flag == 0:
#         nsegd = min(nseg,ndim)
    
#     #2
#     if flag == 2:
#         xhib = a
#         bin = (b-a)/np.real(nsegd)
#         #3
#         for i in range(1,nsegd+1):
#             xlo[i-1] = xhib
#             xlob = xlo[i-1]
#             xhi[i-1] = xhib + bin
#             if i == nsegd:
#                 xhi[i-1] = b
#             xhib = xhi[i-1]
#             tval[i-1], te = rgs56p(f, xlob, xhib)
#             ters[i-1] = te**2
#             #3 continue
#         nter = nsegd
#         #4
#         for iter in range(1, ndim+1):
#             tvals = tval[0]
#             terss = ters[0]
#             #5
#             for i in range(2, nter+1):
#                 tvals = tvals + tval[i-1]
#                 terss = terss + ters[i-1]
#                 #5 continue
#             root = np.sqrt(2*terss)
#             if root <= abstol or root <= reltol*abs(tvals):
#                 break
#             if nter == ndim:
#                 break
#             bige = ters[0]
#             ibig = 1
#             #6
#             for i in range(2, nter+1):
#                 if ters[i-1] > bige:
#                     bige = ters[i-1]
#                     ibig = i
#                 #6 continue
#             nter = nter+1
#             xhi[nter-1] = xhi[ibig-1]
#             xnew = hf*(xlo[ibig-1] + xhi[ibig-1])
#             xhi[ibig-1] = xnew
#             xlo[nter-1] = xnew
#             tval[ibig-1], te = rgs56p(f, xlo[ibig-1], xhi[ibig-1])
#             ters[ibig-1] = te**2
#             tval[nter-1], te = rgs56p(f, xlo[nter-1], xhi[nter-1])
#             ters[nter-1] = te**2
#             #4 continue
#     #9
#     res=tvals
#     err=root
        
#     return res, err





# def funlz(func,x2low,x2high):
#     xlow = x2low
#     xhigh = x2high
#     xmid = xlow
#     found = 0
#     if func(xlow) > 0:
#         found = 120
#     else:
#         xmid = xhigh
#         if func(xhigh) > 0:
#             found = 50
#     #30
#     if found == 0:
#         for logn in range(1,8):
#             nslice = 2**logn
#             #20
#             for i in range(1, nslice+1, 2):
#                 xmid = xlow + (np.float64(i) * (xhigh-xlow)) / np.float64(nslice)
#                 if func(xmid) > 0:
#                     found = 50
#                     break
#             if found > 0:
#                 break
#             # 20 continue
#             # 30 continue

#     if found == 50:
#         xh = xmid
#         xl = xlow
#         #70
#         for k in range(1,21):
#             xnew = 0.5*(xh+xl)
#             if func(xnew) == 0:
#                 xl = xnew
#             else:
#                 xh = xnew
#             #70 continue
#         xlow = xl
#         found = 120
#     if found == 120:
#         if func(xhigh) > 0:
#             return xlow, xhigh
#         xl = xmid
#         xh = xhigh
#         #170
#         for k in range(1, 21):
#             xnew = 0.5*(xh+xl)
#             if func(xnew) == 0:
#                 xh = xnew
#             else:
#                 xl = xnew
#             #170 continue
#         xhigh = xh
#     if found == 0:
#         xlow = 0.
#         xhigh = 0.
#     return xlow, xhigh
        