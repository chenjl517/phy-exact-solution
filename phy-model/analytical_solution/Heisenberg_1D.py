"""
High Temperature Expansion for XXX Heisenberg chain (external field h =0)
Integral Equation Generates High-Temperature Expansion of the Heisenberg Chain
DOI: 10.1103/PhysRevLett.89.117201
"""

import numpy as np

def calc_free_energy(J, T, n=10):
    x = J / T

    f_div_T = 0
    for i in range(n):
        f_div_T += calc_coeff(i) * x ** i
    f = f_div_T * T

    ham_bias = -1.0 / 4

    f = f + ham_bias

    return f

#-------------------------------------
import numpy as np

def calc_coeff(n):
    if(n==0):
        return -np.log(2*np.cosh(0))
    elif(n==1):
        return 1/4
    else:
        alpha_n = [0, 0, 6, -36, -360, 7200, 15120, -1848672, 11426688, 594846720, -11558004480, 
                   -199812856320, 10106191180800, 19376365252608, -9289795522775040, 121944211136778240,
                    8791781390116945920, -310402124957945954304, -7225535925744106143744, 643407197363813620776960,
                    96147483542540314214400, -1279121513829538179364945920, 27962069861743501862336200704, 2398518627113966015427501883392,
                    -129834725539335848980192847462400, -3493877000064415911285457158144000, 484455914465376683487755420408217600,
                    -1592964364128699671723658807556964352, -1659222341377723674454893065936371187712, 53827694891210973745020673240061454581760,
                    5090517962961447184851808942927438864711680, -388446833192725659973817494776649147157053440, -11028658525378274359384407389662010654796546048,
                    2255854109806569120380670028755308714386167693312, -18878702580622989070078793482363993425701807063040, -11721570087037734701356860480896473609272542720163840,
                    527183038642469386328859769728396518893185382125404160, 52548749252010967993948480669499712309299856992951599104,
                    -5382365237582398925074954773487741035075672601159589167104, -150281021589219619860159284209265140804955107276364175114240,
                    44482678475307391762958213932359681737961800852665438128046080, -670778300712303276022754187872671936481343744506621812675706880,
                    -323311185126253530334911142992092649497388429937499549362387156992, 19271391500067613736198673193545354611765664770995927250862568636416,
                    1963797073102024140530884388201619017857423297479642613447074848440320, -261757449501391383349154989821467694962901907072496780634929173522022400,
                    -6715036186134671522475926929150627328836680429076585020863949295740518400, 2897640509780835688484069216581936412870887902144153768250804439297575354368,
                    -67884583842448252705729493380589916284543590592089243545625816663055684075520, -27839667354545349510646207322267839343449407617219100508752156300912706663219200,
                    2130333568970965233678580974509426707048585535358286474483865352948915543579033600,  216827879657653769500650387534438914339017251205042192428944424756474567924803174400]
        if(n > len(alpha_n)):
            return 0
        else:
            return alpha_n[n] / (4**n * np.math.factorial(n) * n *(1-n))