#############################################################################
# SRWLib for Python: Undulator Utilities v 0.1
#############################################################################

from srwlib import *
#from copy import *

#****************************************************************************
def srwl_und_cor_fld_int(_mag3d, _dist_bw_kicks, _rms_len_kicks=0.05, _zc=0, _zcMesh=0, _zRange=0, _dupl=False):
    """
    Compensates 1st and 2nd integrals of undulator field by adding "kicks" before and after the undulator.
    :param _mag3d: 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _dist_bw_kicks: distance between kicks (in longitudinal direction) [m]
    :param _rms_len_kicks: RMS kick length [m]
    :param _zc: center position for the kick set (~center of indulator) [m]
    :param _zcMesh: longitudinal center position for the magnetic field mesh [m]
    :param _zRange: range of magnetic field to be taken into account [m]
    :param _dupl: duplicate the magnetic field object or not
    :returns: undated magnetic 3d field structure
    """

    if(_mag3d.nz <= 1): return _mag3d

    nz = _mag3d.nz
    zStep = _mag3d.rz/(nz - 1)

    #print(0.5*_dist_bw_kicks, _rms_len_kicks, _zc)

    halfDistBwKicks = 0.5*_dist_bw_kicks
    zIn = _zc - halfDistBwKicks #Position of Input Kick
    zOut = _zc + halfDistBwKicks #Position of Output Kick

    #zE = 0.5*_mag3d.rz #? End position of magnetic field
    #zB = -zE #? Start position of magnetic field
    halfFullRangeZ = 0.5*_mag3d.rz
    zB = _zcMesh - halfFullRangeZ #? End position of magnetic field
    zE = _zcMesh + halfFullRangeZ #? Start position of magnetic field

    print('zIn=', zIn, 'zOut=', zOut, 'zc=', _zc)
    print('zB=', zB, 'zE=', zE)

    zBr = zB
    zEr = zE
    if(_zRange > 0.):
        halfRange = 0.5*_zRange
        zBr = _zcMesh - halfRange
        zEr = _zcMesh + halfRange
    
    invDistBwKicks = 1./_dist_bw_kicks

    arBxOrig = _mag3d.arBx
    arByOrig = _mag3d.arBy
    manyTrPos = False
    if((_mag3d.nx > 1) or (_mag3d.ny > 1)):
        manyTrPos = True
        arBxOrig = array('d', [0]*nz) 
        arByOrig = array('d', [0]*nz) 

    halfRangeZ = 4.*_rms_len_kicks
    izInB = round((zIn - halfRangeZ - zB)/zStep)
    if(izInB < 0): izInB = 0
    izInE = round((zIn + halfRangeZ - zB)/zStep)
    if(izInE >= nz): izInE = nz - 1
    izOutB = round((zOut - halfRangeZ - zB)/zStep)
    if(izOutB < 0): izOutB = 0
    izOutE = round((zOut + halfRangeZ - zB)/zStep)
    if(izOutE >= nz): izOutE = nz - 1

    izBr = round((zBr - zB)/zStep)
    if(izBr < 0): izBr = 0
    izEr = round((zEr - zB)/zStep)
    if(izEr >= nz): izEr = nz - 1

    resMag3d = _mag3d
    if(_dupl == True): resMag3d = deepcopy(_mag3d)
    arBxRes = resMag3d.arBx
    arByRes = resMag3d.arBy

    #print(izInB, izInE, izOutB, izOutE)
    #print(izInB, izInE, izOutB, izOutE)
    print(izBr, izEr)
    
    perY = _mag3d.nx
    perZ = perY*_mag3d.ny
    for iy in range(_mag3d.ny):
        iy_perY = iy*perY
        for ix in range(_mag3d.nx):
            if(manyTrPos == True):
                arBxTot = _mag3d.arBx
                arByTot = _mag3d.arBy
                for iz in range(nz):
                    ofst = ix + iy_perY + iz*perZ
                    if(arBxTot != None):
                        arBxOrig = arBxTot[ofst]
                    if(arByTot != None):
                        arByOrig = arByTot[ofst]

            if(arBxOrig != None):
                arBxInt = uti_math.integ_array(arBxOrig, zStep, True)
                
                #auxI1X = -arBxInt[nz - 1] #[T.m]
                auxI1X_Er = arBxInt[izEr]
                auxI1X_Br = arBxInt[izBr]
                auxI1X = -(auxI1X_Er - auxI1X_Br) #[T.m]
                
                arBxInt = uti_math.integ_array(arBxInt, zStep, True)
                #auxI2X = -arBxInt[nz - 1] #[T.m^2]
                auxI2X = -(arBxInt[izEr] - arBxInt[izBr] - auxI1X_Br*(izEr - izBr)*zStep) #[T.m^2]
                
                kickI1Xin = (auxI2X - auxI1X*(zE - zOut))*invDistBwKicks
                kickI1Xout = (auxI1X - kickI1Xin)

                #print(auxI1X, 'T.m')
                #print(kickI1Xin*(zE - zIn), kickI1Xout*(zE - zOut))
                #I2test = kickI1Xin*(zE - zIn) + kickI1Xout*(zE - zOut)
                #print(auxI2X, I2test)

                B0kickX = kickI1Xin*0.3989422804/_rms_len_kicks
                z = zB + zStep*izInB
                for iz in range(izInB, izInE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zIn)/_rms_len_kicks
                    arBxRes[ofst] += B0kickX*exp(-0.5*t*t)
                    z += zStep
                B0kickX = kickI1Xout*0.3989422804/_rms_len_kicks
                z = zB + zStep*izOutB
                for iz in range(izOutB, izOutE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zOut)/_rms_len_kicks
                    arBxRes[ofst] += B0kickX*exp(-0.5*t*t)
                    z += zStep
                
            if(arByOrig != None):
                arByInt = uti_math.integ_array(arByOrig, zStep, True)
                
                #auxI1Y = -arByInt[nz - 1] #[T.m]
                auxI1Y_Er = arByInt[izEr]
                auxI1Y_Br = arByInt[izBr]
                auxI1Y = -(auxI1Y_Er - auxI1Y_Br) #[T.m]
                
                arByInt = uti_math.integ_array(arByInt, zStep, True)
                #auxI2Y = -arByInt[nz - 1] #[T.m^2]
                auxI2Y = -(arByInt[izEr] - arByInt[izBr] - auxI1Y_Br*(izEr - izBr)*zStep) #[T.m^2]

                kickI1Yin = (auxI2Y - auxI1Y*(zE - zOut))*invDistBwKicks
                kickI1Yout = (auxI1Y - kickI1Yin)

                print(auxI1Y, 'T.m')
                #I2test = kickI1Yin*(zE - zIn) + kickI1Yout*(zE - zOut)
                #print(auxI2Y, I2test)

                B0kickY = kickI1Yin*0.3989422804/_rms_len_kicks
                z = zB + zStep*izInB
                for iz in range(izInB, izInE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zIn)/_rms_len_kicks
                    arByRes[ofst] += B0kickY*exp(-0.5*t*t)
                    z += zStep
                B0kickY = kickI1Yout*0.3989422804/_rms_len_kicks
                z = zB + zStep*izOutB
                for iz in range(izOutB, izOutE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zOut)/_rms_len_kicks
                    arByRes[ofst] += B0kickY*exp(-0.5*t*t)
                    z += zStep
    return resMag3d

#****************************************************************************
def srwl_und_fld_add_const(_mag3d, _zcMesh=0, _zc=0, _zRange=0, _bx=0, _by=0, _bz=0, _dupl=False):
    """
    Adds constant magnetic field within given longitudinal range.
    :param _mag3d: 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _zcMesh: longitudinal center position for the magnetic field mesh [m]
    :param _zc: center position for the dipole field range [m]
    :param _zRange: range of magnetic field over which to add constant field [m]
    :param _bx: horizontal field component to add [T]
    :param _by: vertical field component to add [T]
    :param _bz: longitudinal field component to add [T]
    :param _dupl: duplicate the magnetic field object or not
    :returns: undated magnetic 3d field structure
    """
    
    if(_mag3d.nz <= 1): return _mag3d

    nz = _mag3d.nz
    zStep = _mag3d.rz/(nz - 1)
    inv_zStep = 0. if(zStep == 0.) else 1./zStep

    izAddB = 0
    izAddE = nz - 1

    if(_zRange > 0):
        halfFullRangeZ = 0.5*_mag3d.rz
        zB = _zcMesh - halfFullRangeZ #? End position of magnetic field
        zE = _zcMesh + halfFullRangeZ #? Start position of magnetic field

        halfRange = 0.5*_zRange
        zAddB = _zc - halfRange
        zAddE = _zc + halfRange

        izAddTestB = int(round((zAddB - zB)*inv_zStep))
        izAddTestE = int(round((zAddE - zB)*inv_zStep))

        if(izAddTestB <= izAddTestE):
            if((izAddB < izAddTestB) and (izAddTestB <= izAddE)): izAddB = izAddTestB
            if((izAddB <= izAddTestE) and (izAddTestE < izAddE)): izAddE = izAddTestE

    resMag3d = _mag3d
    if(_dupl == True): resMag3d = deepcopy(_mag3d)

    arBX = resMag3d.arBx
    arBY = resMag3d.arBy
    arBZ = resMag3d.arBz

    bxIsSet = False
    if(arBX != None):
        if(len(arBX) > 0): bxIsSet = True 
    byIsSet = False
    if(arBY != None):
        if(len(arBY) > 0): byIsSet = True
    bzIsSet = False
    if(arBZ != None):
        if(len(arBZ) > 0): bzIsSet = True
    
    for iz in range(izAddB, izAddE + 1):
        if(bxIsSet == True): arBX[iz] += _bx
        if(byIsSet == True): arBY[iz] += _by
        if(bzIsSet == True): arBZ[iz] += _bz

    return resMag3d

#****************************************************************************
def srwl_und_cut_fld(_mag3d, _res_rz, _zc=None, _dupl=False):
    """
    Cuts (truncates) undulator magnetic field. Assumes equidistant mesh.
    :param _mag3d: 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _res_rz: range vs longitudinal position to produce [m]
    :param _zc: center position for the kick set (~center of indulator) [m]
    :param _dupl: duplicate the magnetic field object or not
    :returns: updated 3d magnetic field structure
    """

    if(_mag3d.nz <= 1): return _mag3d

    nx = _mag3d.nx
    ny = _mag3d.ny
    nz = _mag3d.nz
    zStep = _mag3d.rz/(nz - 1)
    #zRange = zStep*(nz - 1)
    zRange = _mag3d.rz

    if((_res_rz <= 0) or (_res_rz >= zRange)):
        if(_dupl == False): return _mag3d
        else: return deepcopy(_mag3d)

    zBegOrig = -0.5*zRange
    zEndOrig = -zBegOrig

    perZ = nx*ny

    arBxOrig = _mag3d.arBx
    arByOrig = _mag3d.arBy
    arBzOrig = _mag3d.arBz

    if(_zc == None):
        #Find Max. field value and position
        Be2max = 0.
        i = 0
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    Bx = 0 if arBxOrig == None else arBxOrig[i]
                    By = 0 if arByOrig == None else arByOrig[i]
                    Bz = 0 if arBzOrig == None else arBzOrig[i]
                    curBe2 = Bx*Bx + By*By + Bz*Bz
                    if(Be2max < curBe2): Be2max = curBe2
                    i += 1

        #Find Longitudinal Center position of the field 
        Bthresh = sqrt(Be2max)*0.1 #to steer
        #print('Bthresh=', Bthresh)

        i = 0; izStart = 0
        wasBreak = False
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    Bx = 0 if arBxOrig == None else arBxOrig[i]
                    By = 0 if arByOrig == None else arByOrig[i]
                    Bz = 0 if arBzOrig == None else arBzOrig[i]
                    curB = sqrt(Bx*Bx + By*By + Bz*Bz)
                    #print('i=', i, 'curB=', curB)
                    
                    if(curB >= Bthresh):
                        izStart = iz
                        wasBreak = True
                        break
                    i += 1
                if(wasBreak): break
            if(wasBreak): break

        izEnd = nz - 1
        wasBreak = False
        for iz in range(nz - 1, -1, -1):
            iz_perZ = iz*perZ
            for iy in range(ny):
                iy_nx_p_iz_perZ = iy*nx + iz_perZ
                for ix in range(nx):
                    i = ix + iy_nx_p_iz_perZ
                    Bx = 0 if arBxOrig == None else arBxOrig[i]
                    By = 0 if arByOrig == None else arByOrig[i]
                    Bz = 0 if arBzOrig == None else arBzOrig[i]
                    curB = sqrt(Bx*Bx + By*By + Bz*Bz)
                    if(curB >= Bthresh):
                        izEnd = iz
                        wasBreak = True
                        break
                if(wasBreak): break
            if(wasBreak): break

        #print('izStart=', izStart, 'izEnd=', izEnd)

        zStart = zBegOrig + izStart*zStep
        zEnd = zBegOrig + izEnd*zStep
        _zc = 0.5*(zStart + zEnd)

    halfResRangeZ = 0.5*_res_rz
    zStartRes = _zc - halfResRangeZ
    if(zStartRes < zBegOrig): zStartRes = zBegOrig
    zEndRes = _zc + halfResRangeZ
    if(zEndRes > zEndOrig): zEndRes = zEndOrig

    #print('zc=', zc, 'zStartRes=', zStartRes, 'zEndRes=', zEndRes)

    zRangeRes = zEndRes - zStartRes

    izStartRes = int(round((zStartRes - zBegOrig)/zStep))
    if(izStartRes < 0): izStartRes = 0

    izEndRes = int(round((zEndRes - zBegOrig)/zStep))
    if(izEndRes >= nz): izEndRes = nz - 1

    if(izEndRes < izStartRes): izEndRes = izStartRes

    zRangeRes = (izEndRes - izStartRes)*zStep #OC04082016
    
    nzRes = izEndRes - izStartRes + 1
    nTot = perZ*nzRes
    auxList = [0]*nTot
    arBxRes = None if(arBxOrig == None) else array('d', auxList)
    arByRes = None if(arByOrig == None) else array('d', auxList)
    arBzRes = None if(arBzOrig == None) else array('d', auxList)

    for iz in range(izStartRes, izEndRes + 1):
        iz_perZ = iz*perZ
        iz0_perZ = (iz - izStartRes)*perZ
        for iy in range(ny):
            iy_nx_p_iz_perZ = iy*nx + iz_perZ
            iy_nx_p_iz0_perZ = iy*nx + iz0_perZ
            
            for ix in range(nx):
                i = ix + iy_nx_p_iz_perZ
                i0 = ix + iy_nx_p_iz0_perZ
                if(arBxOrig != None): arBxRes[i0] = arBxOrig[i]
                if(arByOrig != None): arByRes[i0] = arByOrig[i]
                if(arBzOrig != None): arBzRes[i0] = arBzOrig[i]

    if(_dupl==True):
        return SRWLMagFld3D(_arBx=arBxRes, _arBy=arByRes, _arBz=arBzRes, _nx=nx, _ny=ny, _nz=nzRes,
                            _rx=_mag3d.rx, _ry=_mag3d.ry, _rz=zRangeRes, _nRep=_mag3d.nRep, _interp=_mag3d.interp)
    else:
        _mag3d.arBx = arBxRes; _mag3d.arBy = arByRes; _mag3d.arBz = arBzRes; _mag3d.nz = nzRes; _mag3d.rz = zRangeRes

    return _mag3d

#****************************************************************************
def srwl_und_find_cen_len(_mag3d, relThrB=0.8):
    """
    Finds longitudinal center position and length of undulator (dominating field component) by analyzing its tabulated field.
    :param _mag3d: tabulated 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param relThrB: relative threshold to be used with respect to peak field for determining the center
    """

    if((_mag3d == None) or ((_mag3d.arBx == None) and (_mag3d.arBy == None))):
        raise Exception("Incorrect definition of the input tabulated undulator magnetic field")

    if((relThrB <= 0.) or (relThrB >= 1.)):
        raise Exception("relative threshold to be used with respect to peak field should be larger than 0 and less than 1")
    
    nx = _mag3d.nx
    ny = _mag3d.ny
    nz = _mag3d.nz
    if(nz <= 1): raise Exception("1D magnetic field with more than one grid point vs longitudinal position is expected")
    if((nx > 1) or (ny > 1)): raise Exception("1D magnetic field with more than one grid point vs longitudinal position and one point vs horizontal and vertical position is expected")

    BxIsDefined = False
    if(_mag3d.arBx != None):
        if(len(_mag3d.arBx) == nz): BxIsDefined = True

    ByIsDefined = False
    if(_mag3d.arBy != None):
        if(len(_mag3d.arBy) == nz): ByIsDefined = True

    if((BxIsDefined == False) and (ByIsDefined == False)):
        raise Exception("1D magnetic field data (vertical or horizontal component) are not defined")

    absBxMax = 0.
    absByMax = 0.
    for iz in range(nz):
        if(BxIsDefined):
            curAbsB = abs(_mag3d.arBx[iz])
            if(absBxMax < curAbsB): absBxMax = curAbsB
        if(ByIsDefined):
            curAbsB = abs(_mag3d.arBy[iz])
            if(absByMax < curAbsB): absByMax = curAbsB

    if((absBxMax <= 0.) and (absByMax <= 0.)):
        raise Exception("Non-zero 1D magnetic field data (vertical or horizontal component) are not defined")

    arB = None
    absBmax = 0.
    if(absByMax >= absBxMax):
        arB = _mag3d.arBy
        absBmax = absByMax
    else:
        arB = _mag3d.arBx
        absBmax = absBxMax

    absThreshB = relThrB*absBmax
    zHalfRange = 0.5*_mag3d.rz
    zThreshLeft = -zHalfRange
    zThreshRight = zHalfRange
    zStep = _mag3d.rz/(nz - 1)

    z = zThreshLeft
    for iz in range(nz):
        curAbsB = abs(arB[iz])
        if(curAbsB >= absThreshB):
            zThreshLeft = z
            break
        z += zStep

    nz_mi_1 = nz - 1
    z = zThreshRight
    for iz in range(nz):
        curAbsB = abs(arB[nz_mi_1 - iz])
        if(curAbsB >= absThreshB):
            zThreshRight = z
            break
        z -= zStep

    #print(zThreshLeft, zThreshRight)
    if(zThreshRight <= zThreshLeft):
        return 0., _mag3d.rz
    else:
        return 0.5*(zThreshRight + zThreshLeft), zThreshRight - zThreshLeft

#****************************************************************************
def srwl_und_fld_1d_mis(_mag3d, _per, _dg_by_len, _c1, _c2=0, _g0=0, _a=0, _y0=0, _dydz=0, _z0=None, _dupl=False):
    """
    Simulates effects of gap taper and electron "mis-steering" on magnetic field "seen" by the electron, using the formula:

    B(z) = B0(z)*exp(-(dg/(len*per))*(c1 - 2*c2*g0/per)*(z-z0) + c2*(dg/(len*per))^2*(z-z0)^2)*cosh((2*pi*a/per)*(y0+dydz*(z-z0)))
    
    :param _mag3d: tabulated 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _per: undulator period [m]
    :param _dg_by_len: ratio of gap variation betweein undulator exit and entrance to the undulator length (i.e. dg/len)
    :param _c1: coefficient before linear term in the argument of exponent describing the gap dependence (i.e. c1 in b0*exp(-c1*gap/per + c2*(gap/per)^2))
    :param _c2: coefficient before quadratic term in the argument of exponent describing the gap dependence (i.e. c2 in b0*exp(-c1*gap/per + c2*(gap/per)^2))
    :param _g0: approximate undulator gap [m]
    :param _a: coefficient in the argument of cosh defining dependence of the magnetic field on the vertical position (i.e. a in cosh(2*pi*a*y/per))
    :param _y0: initial vertical position of electron (i.e. at z0 or in the center of the mag3d range) [m]
    :param _dydz:  vertical angle of electron [rad]
    :param _z0: longitudinal position where y0 is defined (if z0 == None, then the center of the mag3d range is assumed) [m]
    :param _dupl: switch specifying whether the magnetic field object has to be duplicated or modified in place
    :returns: updated 3d magnetic field structure
    """

    if((_mag3d == None) or ((_mag3d.arBx == None) and (_mag3d.arBy == None)) or (_per <= 0.)):
        raise Exception("Incorrect definition of the input tabulated undulator magnetic field")
    
    resMag3D = deepcopy(_mag3d) if(_dupl == True) else _mag3d

    nx = resMag3D.nx
    ny = resMag3D.ny
    nz = resMag3D.nz

    if(nz <= 1): raise Exception("1D magnetic field with more than one grid point vs longitudinal position is expected")
    if((nx > 1) or (ny > 1)): raise Exception("1D magnetic field with more than one grid point vs longitudinal position and one point vs horizontal and vertical position is expected")

    arB = resMag3D.arBy
    BorigIsDefined = True
    if(arB == None):
        BorigIsDefined = False
    elif(len(arB) < nz):
        BorigIsDefined = False

    if(BorigIsDefined == False):
        arB = resMag3D.arBx
        if(arB != None):
            if(len(arB) == nz): BorigIsDefined = True

    if(BorigIsDefined == False): raise Exception("1D magnetic field data (vertical or horizontal component) is not defined")

    zStep = resMag3D.rz/(nz - 1)
    zRange = resMag3D.rz
    if(zRange <= 0): return resMag3D

    #z0 = _z0 if(_z0 != None) else 0.
    z0 = _z0 if(_z0 != None) else srwl_und_find_cen_len(resMag3D)[0]
    #print(z0)
    
    z = -0.5*zRange
    for iz in range(nz):
        dz = z - z0
        #print(dz)
       
        dz_dg_d_len_per = dz*_dg_by_len/_per
        argExp = -dz_dg_d_len_per*(_c1 + _c2*(dz_dg_d_len_per - 2*_g0/_per))
        argCosh = (6.283185307*_a/_per)*(_y0 + _dydz*dz)
        mult = exp(argExp)*cosh(argCosh)
        #print(mult)

        arB[iz] *= mult
        z += zStep

    return resMag3D

#****************************************************************************
#def srwl_fld_extrap_grad_off_mid_plane(_mag_mid, _ry, _ny, _grad_mult=1):
#    """
#    Extrapolates magnetic field off horizontal mid-plane based on field gradient in the mid-plane (e.g. for dipole with gradient or for a guad).
#    :param _mag_mid: input tabulated 3D magnetic field on 2D mesh vs x and z in the horizontal mid-plane in a container (object of SRWLMagFldC type)
#    :param _ry: vertical position range of the final 3D mesh [m]
#    :param _ny: number of points vs vertical position in the final 3D mesh
#    :param _grad_mult: a number the gradient has to be multiplied by
#    :returns: resulting exrapolated 3D magnetic field structure in a container (object of SRWLMagFldC type)
#    """

#    if((_mag_mid is None) or (_mag_mid.arMagFld is None) or (len(_mag_mid.arMagFld) != 1) or (_ry < 0.) or (_ny <= 1)):
#        raise Exception("Incorrect input parameters for magnetic field to be extrapolated")

#    fld3d_mid = _mag_mid.arMagFld[0]
#    #arBx_mid = fld3d_mid.arBx
#    arBy_mid = fld3d_mid.arBy
#    arBz_mid = fld3d_mid.arBz

#    xc = _mag_mid.arXc[0]
#    xStart = xc - 0.5*fld3d_mid.rx
#    nx = fld3d_mid.nx
#    xStep = fld3d_mid.rx/(nx - 1)

#    yc = 0
#    yStart = -0.5*_ry
#    ny = int(_ny)
#    yStep = _ry/(ny - 1)

#    zc = _mag_mid.arZc[0]
#    zStart = zc - 0.5*fld3d_mid.rz
#    nz = fld3d_mid.nz
#    zStep = fld3d_mid.rz/(nz - 1)

#    nTotRes = int(nx*ny*nz)
#    arBxRes = array('d', [0]*nTotRes)
#    arByRes = array('d', [0]*nTotRes)
#    arBzRes = None if(arBz_mid is None) else array('d', [0]*nTotRes)
    
#    nx_mi_1 = nx - 1
#    #print(nx, xStart, xStep, ny, yStart, yStep, nz, zStart, zStep)
    
#    for iz in range(nz):
#        iz_nx = iz*nx
#        iz_nx_ny = iz*nx*ny
#        for ix in range(nx):
#            b1y = 0; b2y = 0; b0y = 0; dx = xStep
#            if(ix == 0):
#                b1y = arBy_mid[iz_nx]
#                b2y = arBy_mid[iz_nx + 1]
#                b0y = b1y
#            elif(ix == nx_mi_1):
#                b1y = arBy_mid[iz_nx + nx_mi_1 - 1]
#                b2y = arBy_mid[iz_nx + nx_mi_1]
#                b0y = b2y
#            else:
#                b1y = arBy_mid[iz_nx + ix - 1]
#                b0y = arBy_mid[iz_nx + ix]
#                b2y = arBy_mid[iz_nx + ix + 1]
#                dx = 2*xStep
                
#            curGrad = _grad_mult*(b2y - b1y)/dx
#            y = yStart
#            for iy in range(ny):
#                ofst = iz_nx_ny + iy*nx + ix
#                arBxRes[ofst] = y*curGrad
#                arByRes[ofst] = b0y
#                y += yStep

#    fld3dRes = SRWLMagFld3D(arBxRes, arByRes, arBzRes, nx, ny, nz, fld3d_mid.rx, _ry, fld3d_mid.rz)
#    #print(fld3dRes.nx, fld3dRes.rx, fld3dRes.ny, fld3dRes.ry, fld3dRes.nz, fld3dRes.rz)
#    return SRWLMagFldC(fld3dRes, xc, yc, zc)

#****************************************************************************
##def srwl_fld_extrap_grad_curv(_mag_curv_mid, _cen_trj_data, _xi, _xf, _nx, _yi, _yf, _ny, _zi, _zf, _nz, _grad_mult=1):
##    """
##    Extrapolates magnetic field on rectangulat 3D mesh in and off horizontal mid-plane based on field gradient in the mid-plane (e.g. for dipole with gradient or for a guad).
##    :param _mag_curv_mid: input tabulated 3D magnetic field on 2D mesh vs x and z in Natural frame in the horizontal mid-plane in a container (object of SRWLMagFldC type)
##    :param _cen_trj_data: central trajectory data defining the Natural frame of the input magnetic field (_cen_trj_data[0] is ct, _cen_trj_data[1] is x. _cen_trj_data[2] is z)
##    :param _xi: horizontal initial position of the final 3D mesh [m]
##    :param _xf: horizontal final position of the final 3D mesh [m]
##    :param _nx: number of points vs horizontal position in the final 3D mesh
##    :param _yi: vertical initial position of the final 3D mesh [m]
##    :param _yf: vertical initial position of the final 3D mesh [m]
##    :param _ny: number of points vs vertical position in the final 3D mesh
##    :param _zi: longitudinal initial position of the final 3D mesh [m]
##    :param _zf: longitudinal final position of the final 3D mesh [m]
##    :param _nz: number of points vs longitudinal position in the final 3D mesh
##    :param _grad_mult: a number the gradient has to be multiplied by
##    :returns: exrapolated 3D magnetic field structure in a container (object of SRWLMagFldC type)
##    """
##
##    if((_mag_curv_mid is None) or (_mag_curv_mid.arMagFld is None) or (len(_mag_curv_mid.arMagFld) < 1) or (_cen_trj_data is None) or (len(_cen_trj_data) < 3)
##       or (_nx < 1) or (_ny < 1) or (_nz < 1)):
##        raise Exception("Incorrect input parameters for magnetic field to be extrapolated")
##    #print('In srwl_fld_extrap_grad_curv')
##
##    fld3d_curv_mid = _mag_curv_mid.arMagFld[0]
##    #arBx_curv_mid = fld3d_curv_mid.arBx
##    arBy_curv_mid = fld3d_curv_mid.arBy
##    arBz_curv_mid = fld3d_curv_mid.arBz

#*********************************Entry point (for command-prompt calls)
if __name__=="__main__":
    import optparse
    p = optparse.OptionParser()
    p.add_option('--cor', dest='do_cor_fld_int', action='store_true', default=False, help="correct 1st and 2nd field integrals of tabulated undulator magnetic field")
    p.add_option('--conv', dest='do_conv_fld', action='store_true', default=False, help="convert magnetic measurement data file to SRW ASCII format")
    #p.add_option('--trun', dest='do_trunc_fld', action='store_true', default=False, help="truncate magnetic field (i.e. reduce range) vs longitudinal position")
    #p.add_option('--efg', dest='do_ext_fld_grd', action='store_true', default=False, help="extrapolate 3D magnetic measurement data out of horizontal mid-plane using gradient")
    #p.add_option('--efgc', dest='do_ext_fld_grd_curv', action='store_true', default=False, help="extrapolate 3D magnetic measurement data out of horizontal mid-plane using gradient, assuming input field data is gived for curved path along central trajectory")

    p.add_option('--ifn', dest='ifn', type='string', default='', help="input field file name")
    p.add_option('--ofn', dest='ofn', type='string', default='', help="output field file name")
    p.add_option('--itfn', dest='itfn', type='string', default='', help="input trajectory file name")
    p.add_option('--otfn', dest='otfn', type='string', default='', help="output trajectory file name")

    p.add_option('--ozc', dest='ozc', type="float", default=1.e+23, metavar="NUMBER", help="center longitudinal coordinate of the output magnetic field [m]") #OC22062017
    
    p.add_option('--dbwk', dest='dbwk', type="float", default=2.8, metavar="NUMBER", help="distance between correcting kicks [m]")
    p.add_option('--rmsk', dest='rmsk', type="float", default=0.05, metavar="NUMBER", help="RMS length of the correcting kicks [m]")

    p.add_option('--rz', dest='rz', type="float", default=0., metavar="NUMBER", help="new longitudinal range of magnetic field after truncation [m] (effective if > 0; the center field will be taken from zk)")

    p.add_option('--zk', dest='zk', type="float", default=1.e+23, metavar="NUMBER", help="longitudinal position of the correcting kicks [m]")
    p.add_option('--z0t', dest='z0t', type="float", default=1.e+23, metavar="NUMBER", help="initial longitudinal position for trajectory calculation [m]")
    p.add_option('--x0t', dest='x0t', type="float", default=0., metavar="NUMBER", help="initial horizontal position for trajectory calculation [m]")
    p.add_option('--y0t', dest='y0t', type="float", default=0., metavar="NUMBER", help="initial vertical position for trajectory calculation [m]")
    p.add_option('--xp0t', dest='xp0t', type="float", default=0., metavar="NUMBER", help="initial horizontal angle for trajectory calculation [m]")
    p.add_option('--yp0t', dest='yp0t', type="float", default=0., metavar="NUMBER", help="initial vertical angle for trajectory calculation [m]")
    p.add_option('--elen', dest='elen', type="float", default=3., metavar="NUMBER", help="electron energy [GeV]")
    p.add_option('--dbx', dest='dbx', type="float", default=0., metavar="NUMBER", help="horizontal magnetic field component to add [T]")
    p.add_option('--dby', dest='dby', type="float", default=0., metavar="NUMBER", help="vertical magnetic field component to add [T]")
    p.add_option('--dbz', dest='dbz', type="float", default=0., metavar="NUMBER", help="longitudinal magnetic field component to add [T]")
    p.add_option('--dzc', dest='dzc', type="float", default=0., metavar="NUMBER", help="center longitudinal position for adding constant magnetic field [m]")
    p.add_option('--dzr', dest='dzr', type="float", default=0., metavar="NUMBER", help="range of longitudinal position within which to add constant magnetic field [m]")

    #p.add_option('--xi', dest='xi', type="float", default=0., metavar="NUMBER", help="horizontal initial position of the final 3D mesh [m]")
    #p.add_option('--xf', dest='xf', type="float", default=0., metavar="NUMBER", help="horizontal final position of the final 3D mesh [m]")
    #p.add_option('--nx', dest='nx', type="float", default=0., metavar="NUMBER", help="number of points vs horizontal position in the final 3D mesh")
    #p.add_option('--yi', dest='yi', type="float", default=0., metavar="NUMBER", help="vertical initial position of the final 3D mesh [m]")
    #p.add_option('--yf', dest='yf', type="float", default=0., metavar="NUMBER", help="vertical final position of the final 3D mesh [m]")
    #p.add_option('--ny', dest='ny', type="float", default=0., metavar="NUMBER", help="number of points vs vertical position in the final 3D mesh")
    #p.add_option('--zi', dest='zi', type="float", default=0., metavar="NUMBER", help="longitudinal initial position of the final 3D mesh [m]")
    #p.add_option('--zf', dest='zf', type="float", default=0., metavar="NUMBER", help="longitudinal final position of the final 3D mesh [m]")    
    #p.add_option('--nz', dest='nz', type="float", default=0., metavar="NUMBER", help="number of points vs longitudinal position in the final 3D mesh")
 
    opt, args = p.parse_args()

    #if(opt.do_ext_fld_grd): #Extrapolate 3D magnetic measurement data out of horizontal mid-plane using gradient
    #    fldMidPl = srwl_uti_read_mag_fld_3d(opt.ifn)
    #    fldRes = srwl_fld_extrap_grad_off_mid_plane(fldMidPl, opt.ry, opt.ny)
    #    fld3dRes = fldRes.arMagFld[0]
    #    fld3dRes.save_ascii(opt.ofn, fldRes.arXc[0], fldRes.arYc[0], fldRes.arZc[0])
    #    sys.exit()

    #if(opt.do_ext_fld_grd_curv): #Extrapolate 3D magnetic measurement data out of horizontal mid-plane using gradient, assuming input field data is gived for curved path along central trajectory"
    #    fldCurvMidPl = srwl_uti_read_mag_fld_3d(opt.ifn)
    #    cenTrajDataAll = srwl_uti_read_data_cols(opt.itfn, '\t', _i_col_start=0, _i_col_end=6, _n_line_skip=1)
    #    #print(len(cenTrajData))
    #    cenTrajData_ct_x_z = [cenTrajDataAll[0], cenTrajDataAll[1], cenTrajDataAll[5]]
    #    fldRes = srwl_fld_extrap_grad_curv(fldCurvMidPl, cenTrajData_ct_x_z, opt.xi, opt.xf, opt.nx, opt.yi, opt.yf, opt.ny, opt.zi, opt.zf, opt.nz)
        
    #    #fld3dRes = fldRes.arMagFld[0]
    #    #fld3dRes.save_ascii(opt.ofn, fldRes.arXc[0], fldRes.arYc[0], fldRes.arZc[0])
    #    sys.exit()

    #print(opt.ifn)
    fldCnt = None

    if(opt.do_conv_fld): #Convertion from magnetic measurements data format
        fldDataCols = srwl_uti_read_data_cols(opt.ifn, '\t')
        
        arZ = array('d', fldDataCols[0])
        for i in range(len(arZ)):
            arZ[i] *= 0.001
        
        arBX = array('d', fldDataCols[1])
        arBY = array('d', fldDataCols[2])
        arBZ = array('d', fldDataCols[3])

        zNp = len(arZ)
        zSt = arZ[0] #*0.001
        zFi = arZ[zNp - 1] #*0.001

        #Field tabulated on Irregular Mesh:
        #fldCntIrreg = SRWLMagFldC(SRWLMagFld3D(arBX, arBY, arBZ, 1, 1, zNp, 0., 0., zFi - zSt, _interp=2, _arX=None, _arY=None, _arZ=arZ), 0., 0., 0.5*(zSt + zFi))
        fldCntIrreg = SRWLMagFldC(SRWLMagFld3D(arBX, arBY, arBZ, 1, 1, zNp, 0., 0., zFi - zSt, _interp=2, _arX=None, _arY=None, _arZ=arZ), 0., 0., 0.)

        #Re-calculating the Field on Regular Mesh:
        #fldCnt = SRWLMagFldC(SRWLMagFld3D(arBX, arBY, arBZ, 1, 1, zNp, 0., 0., zFi - zSt), 0., 0., zcRes) #OC22062017
        fldCnt = SRWLMagFldC(SRWLMagFld3D(arBX, arBY, arBZ, 1, 1, zNp, 0., 0., zFi - zSt), 0., 0., 0.5*(zSt + zFi))
        #fldCnt = SRWLMagFldC(SRWLMagFld3D(arBX, arBY, arBZ, 1, 1, zNp, 0., 0., zFi - zSt), 0., 0., 0.)
        srwl.CalcMagnField(fldCnt, fldCntIrreg, [0])
        
        #fldCnt = SRWLMagFldC(SRWLMagFld3D(array('d', fldDataCols[1]), array('d', fldDataCols[2]), array('d', fldDataCols[3]), 1, 1, zNp, 0., 0., zFi - zSt), 0., 0., 0.5*(zSt + zFi))
        #print(1)
        
    if(fldCnt == None): fldCnt = srwl_uti_read_mag_fld_3d(opt.ifn)
    
    fld3d = fldCnt.arMagFld[0]
    zc = fldCnt.arZc[0]
    if(opt.ozc < 1.e+23):
        zc = opt.ozc
        fldCnt.arZc[0] = opt.ozc
    
    print('   zc=', zc)

    if((opt.dbx != 0.) or (opt.dby != 0.) or (opt.dbz != 0.)):
        srwl_und_fld_add_const(fld3d, fldCnt.arZc[0], opt.dzc, opt.dzr, opt.dbx, opt.dby, opt.dbz)
        #print(2)

    zk = opt.zk
    
    if(opt.do_cor_fld_int):
        #zc = fldCnt.arZc[0]
        if(zk >= 1.e+22):
            zk, auxUndLen = srwl_und_find_cen_len(fld3d)
            #print('zk=', zk)
            
        srwl_und_cor_fld_int(fld3d, opt.dbwk, opt.rmsk, zk, zc)
        #print(3)

    if(opt.rz > 0.): #Cutting the field
        if(zk >= 1.e+22):
            zk, auxUndLen = srwl_und_find_cen_len(fld3d)
            
        srwl_und_cut_fld(fld3d, opt.rz, zk)
        zc += zk #Because after the cutting zk becomes 0!
        #print(4)

    if(len(opt.ofn) > 0):
        
        fld3d.save_ascii(opt.ofn, 0., 0., zc)
        print('zc=', zc)
        #print(5)

    if(len(opt.otfn) > 0):
        #z0 = opt.z0t if(opt.z0t != 1.e+23) else fldCnt.arZc[0] - 0.5*fld3d.rz
        z0 = opt.z0t if(opt.z0t != 1.e+23) else zc - 0.5*fld3d.rz
        
        elec = SRWLParticle(_x=opt.x0t, _y=opt.y0t, _z=z0, _xp=opt.xp0t, _yp=opt.yp0t, _gamma=opt.elen/0.51099890221e-03)
        traj = SRWLPrtTrj(_np=fld3d.nz, _ctStart=0, _ctEnd=fld3d.rz, _partInitCond=elec)
        #print(fld3d.arBy)
        srwl.CalcPartTraj(traj, fldCnt, [1])
        traj.save_ascii(opt.otfn)
        #print(6)
