cez: 12/14/21: since I couldn't get NPP to stay constant on darwin3,
cez: back to redoing this here. New Yields from Barbara Bayer:
cez:    AOA: 0.096 mol NH4/mol C fixed (table 2): mean of 0.082 and 0.11 
cez:    NOB: 0.042 mol NO2/mol C fixed (table 2): for NOB 
cez: ALSO equal affinity for NOB and AOA (NOfac = tentothethird = 1), and 
cez: ALSO make nitrifiers (all microbes) temp sensitive 
cez: ALSO fix negligible bug: resptoNO2 = growth*(1. _d 0 / ynh4-R_NC(j)) 
cez: 1/3/21: no temp, but high yields (didn't change O2 yields --OK ) 
cez: 1/3/21: no temp, but low yields (didn't change O2 yields --OK ) 
cez: 1/3/21: add nitrifier temp back in (also changed reminTempfunc to
cez reminTempFunc but fortran is case insensitive so this doesn't matter)
cez: 6/22/17: from codeNew6BbothOM_pn30_f3_kd01_LOWERyno2_fourG_adjT,
cez: take out temp relationship for nitrification rate/N uptake: (Horak)
cez: 6/12/17: adjust pnmax for nitrifiers to match 50.8 at 30C
cez: 5/13/17: all one OM-consuming population: B eat both POM and DOM
cez:          taking out hydrolosis of POM to DOM for now (aka alpha=1)
cez: 12/21/16: finding the right vmax
cez: 11/29/16: now, use allometric theory for NOB bigger: 
cez:    -->10^(1/3)*kN and 10^(-1/3)*Vmax
cez: 11/16/16: here, make NOB bigger and lower affinity
cez: 11/16/16: back to f=3 for yields, pnmax=50
cez: 11/12/16: new yields
cez: pnmax is 10, so really 50, (equiv to 50 in 1D model)
cez: 11/10/16: realized y*pnmax was off by a factor of 5 (so was 100!
cez: (so max growth rate was 5x higher than i thought - put this 5 in)
cez: 11/2/16: changed pnmax to 80 (from 20) for 1Z model
cez: 10/6/16: fixed prodO2 bug (was zeroing out by bacteria!!!)
cez: 10/3/16: track the uptake of NO3 by pp somehow: use "Denit" diagn
cez: 9/25/16: CNred (so use CNold) and f=0.02. also ynd=0.13 (from 0.11)
cez: 9/19/16: update nitrification stoich. f=0.03
cez: 5/19/16: get 6B going again with params from May work for 4B
cez: 12/29/15: and add new O2 prod by pp
cez: 12/26/15: add nitrifiers as j=11 and 12
cez: 12/15/15: modify to match updates in Add2B (POM,growth diag)
cez: 10/27/15: add two denitrifiers in, so four total! 
cez: this is for adding 2 types, so B are j=7 and 8
cez: just POM and DOM eating aerobes
cez: gud_model_4 types: has DOM and POM aerobes and nitrifiers:
cez, cont"  DOM denitrifers (to N2) commented out here. 
cez: gud_model_ALL: make this the code to start from for all types (can
change j to turn on/off).
cez: possibilities here:
cez: PA heterotrophic (eats POM), FL heterotrophic (eats DOM),
cez: Nitrification step 1 and 2, Denitrification full (to N2).
cez: codeNitr: from code3, make c02 and c03 nitrifiers
cez: code 3: from code2, replace DOM remin, too.
C cez: this was modified to include bacteria starting 8/15/14 
C $Header$
C $Name$
C TODO:
C - carbon
C - cdom
C - conservation of Chl?
C - PIC: grazers must have R_PICPOC = 0!
C   need to set all R_PICPOC to zero to emulate monod without carbon
C   (maybe not)
C   with carbon???

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_MODEL
C !INTERFACE: ==========================================================
      SUBROUTINE GUD_MODEL(
     I     Ptr,
     U     gTr,
     O     chlout, diags,
     I     PAR, photoTempFunc, reminTempFunc, uptakeTempFunc,
#ifdef GUD_DEBUG
     I     iG, jG, k, dT,
#endif
     I     myTime,myIter,myThid)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_DIAGS.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  Ptr    :: gud model tracers
C  PAR    :: PAR in uEin/s/m2
C         :: either non-spectral (tlam=1) or waveband total
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      _RL Ptr(nGud)
      _RL PAR(nlam)
      _RL photoTempFunc(nplank)
      _RL reminTempFunc
      _RL uptakeTempFunc
      INTEGER myThid, myIter
      _RL myTime
#ifdef GUD_DEBUG
      INTEGER iG, jG, k
      _RL dT
#endif

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: accumulates computed tendencies
      _RL gTr(nGud)

C !OUTPUT PARAMETERS: ==================================================
C  chlout :: computed acclimated chlorophyll if not dynamic
      _RL chlout(nPhoto)
      _RL diags(gud_nDiag)
CEOP

#ifdef ALLOW_GUD

c !LOCAL VARIABLES: ====================================================
      INTEGER j, l

      _RL DIC
      _RL NH4
      _RL NO2
      _RL NO3
      _RL PO4
      _RL SiO2
      _RL FeT
      _RL DOC
      _RL DON
      _RL DOP
      _RL DOFe
      _RL POC
      _RL PON
      _RL POP
      _RL POSi
      _RL POFe
      _RL PIC
      _RL O2

      _RL X
      _RL Qc
      _RL Qn
#ifdef GUD_ALLOW_PQUOTA
      _RL Qp
#endif
#ifdef GUD_ALLOW_SIQUOTA
      _RL Qsi
#endif
#ifdef GUD_ALLOW_FEQUOTA
      _RL Qfe
#endif
#ifdef GUD_ALLOW_CHLQUOTA
      _RL QChl
#endif

      _RL regQ

      _RL limitpCO2
      _RL limitNH4
      _RL limitNO2
      _RL limitNO3
      _RL fracNH4
      _RL fracNO2
      _RL fracNO3
      _RL limitn
      _RL limitp
      _RL limitsi
      _RL limitfe
      _RL limitnut
      _RL limitI
      _RL ngrow
cez add POM and O2 limits
      _RL limitPON
      _RL limitPOC
      _RL limitPOP
      _RL limitPOFe
      _RL limitDON
      _RL limitDOC
      _RL limitDOP
      _RL limitDOFe
      _RL limitO2
cez end add

      _RL uptakeDIC
      _RL uptakeNH4
      _RL uptakeNO2
      _RL uptakeNO3
      _RL uptakeN
      _RL uptakePO4
      _RL uptakeSiO2
      _RL uptakeFeT
      _RL consumDIC
      _RL consumDIC_PIC
      _RL consumNH4
      _RL consumNO2
      _RL consumNO3
      _RL consumPO4
      _RL consumSiO2
      _RL consumFeT

cez add POM uptake and consumption: 
      _RL uptakePON
      _RL uptakePOP
      _RL uptakePOC
      _RL uptakePOFe
      _RL uptakeDON
      _RL uptakeDOP
      _RL uptakeDOC
      _RL uptakeDOFe
      _RL uptakeO2
      
      _RL respPON
      _RL respPOP
      _RL respPOC
      _RL respPOFe
      _RL respDON
      _RL respDOP
      _RL respDOC
      _RL respDOFe
      _RL resptoNO2
      _RL resptoNO3
      
      _RL hydrolPON
      _RL hydrolPOP
      _RL hydrolPOC
      _RL hydrolPOFe
      _RL solubilPON
      _RL solubilPOP
      _RL solubilPOC
      _RL solubilPOFe

      _RL consumPON
      _RL consumPOP
      _RL consumPOC
      _RL consumPOFe
      _RL consumDON
      _RL consumDOP
      _RL consumDOC
      _RL consumDOFe
      _RL consumO2
cez add O2 production by phytopl:
      _RL prodO2
cez end add
      
      _RL inhibNH4

      _RL alpha_I
      _RL alpha_I_growth
      _RL PCm
      _RL PC
      _RL acclim
      _RL chl2c
      _RL growth
      _RL rhochl
      _RL Ek
      _RL EkoverE

      _RL synthChl

      _RL reminDOC
      _RL reminDON
      _RL reminDOP
      _RL reminDOFe
      _RL reminPOC
      _RL reminPON
      _RL reminPOP
      _RL reminPOSi
      _RL reminPOFe
      _RL disscPIC

      _RL prodNO2
      _RL prodNO3

      _RL PARtot

      _RL tmp

#ifdef GUD_ALLOW_CDOM
      _RL CDOM
      _RL reminPOC_CDOM
      _RL reminPON_CDOM
      _RL reminPOP_CDOM
      _RL reminPOFe_CDOM
      _RL degrCDOM_DOC
      _RL degrCDOM_DON
      _RL degrCDOM_DOP
      _RL degrCDOM_DOFe
#endif

#ifdef GUD_ALLOW_DENIT
      _RL denit, denitNH4
#endif

      _RL EPS

cez   add parameters
      _RL ksatPON
cez: uptake maximum is for all OM here:
      _RL pmaxOM
      _RL ksatDON
      _RL pcoefO2
      _RL pmaxDINden
      _RL pmaxDIN
      _RL convpmaxDIN
      _RL tentothethird
      _RL ksatDIN
      _RL d_p
      _RL yod
      _RL yoe
      _RL ynd
      _RL yne1
      _RL yne2
      _RL ynh4
      _RL yonh4
      _RL yno2
      _RL yono2
      _RL ynh4anmx
      _RL yno2anmx
      _RL alpha
      _RL R_O2prod
cez-end
      PARAMETER (EPS=1D-38)

C==== precompute a few things ==========================================

cez   set parameter values
      ksatPON=0.1 _d 0
      ksatDON=0.1 _d 0
      pmaxOM=1. _d 0/86400 _d 0
      pcoefO2=2329.100 _d 0/86400 _d 0
      pmaxDINden=20. _d 0/86400 _d 0
cez: this value is so that at 30C, pmaxDIN is 50.8. (at 20C, it's 32.4)
cez: see gud_tempfunc: tempcoeff is 0.8
cez: put in rate at 30C:
      pmaxDIN=50.8 _d 0/86400 _d 0
cez: use this to convert pmax*y since y's are in N/N
      convpmaxDIN=0.2 _d 0
cez: size 
cez      tentothethird=2.154434690031884 _d 0
cez: 1.414 gives a 2-fold decrease in affinity = sqrt(2)       
cez      tentothethird=1.414 _d 0
cez: Equal affinity:
      tentothethird=1. _d 0
      ksatDIN=0.133 _d 0
      alpha=2.0 _d 0
cez: 7/28/16: change B to C5O2H7N (R+McC)=20*16 (before:467)
cez: 7/28/16: change nitrifier stoich f to data-based 0.02
cez: 9/19/16: add d_p, the denominator for biomass synth norm. to P      
cez: should consistent with R_PC (320 for C:P=80, 467 for C:P of 106) 
      d_p=467. _d 0 
      yod=0.14 _d 0
cez: 9/19/16: here, R of C:P is of Bo (C:N of 5, Zimmerman 2014, R*McC 2014)
cez: SHOULD put this in below as this R_PC(j).
cez: checked - this is right, but need to keep d_p consistent with R_PC-
!should link or keep assigned here.
      yoe=yod/d_p*4. _d 0/(1. _d 0-yod)/R_PC(7)*1. _d 0 !
      R_O2prod=1. _d 0/yoe/(1. _d 0/yod-1)*1. _d 0
      ynd=0.13 _d 0
!NO3 --> NO2:
      yne1=ynd/d_p*2. _d 0/(1. _d 0-ynd)/R_PC(7)*1. _d 0
!NO2 --> N2:
      yne2=ynd/d_p*3. _d 0/(1. _d 0-ynd)/R_PC(7)*1. _d 0
!cez      fnh4=0.02 _d 0
cez 1/4/21:      ynh4=1. _d 0 / 112. _d 0 * 5. _d 0 !1/(1+d_p/6/fnh4)/R_PC(7)*1. _d 0 !units: mol C biomass/mol NH4
      ynh4 = 0.096 _d 0 !mol NH4/mol C fixed BBayer mean of 0.082 and 0.110 AOA
      yonh4=1. _d 0 / 162. _d 0 * 5. _d 0 !fnh4/d_p*4/(1-fnh4)/R_PC(7)*1. _d 0 !units: mol C biom/mol O2
cez      fno2=0.02 _d 0
cez 1/4/21:      yno2=1. _d 0 / 334. _d 0 * 5. _d 0 !1/(1+d_p/2/fno2)/R_PC(7)*1 _d 0
      yno2 = 0.042 _d 0 !mol NO2/mol C fixed BBayer for NOB
cez      yno2=1. _d 0 / 223. _d 0 * 5. _d 0 !1/(1+d_p/2/fno2)/R_PC(7)*1 _d 0
cez: next, try a quarter of it: (NOTE: now not e balanced. mock mixotr)
cez      yno2=1. _d 0 / 168. _d 0 * 5. _d 0 !1/(1+d_p/2/fno2)/R_PC(7)*1 _d 0
      yono2=1. _d 0 / 162. _d 0 * 5. _d 0 !fno2/d_p*4/(1-fno2)/R_PC(7)*1. _d 0
!cez: anammox:
      ynh4anmx=1. _d 0 / 223. _d 0 * 5. _d 0 !1/(1+d_p/2/fno2)/R_PC(7)*1 _d 0
      yno2anmx=1. _d 0 / 216. _d 0 * 5. _d 0 !1/(1+d_p/2/fno2)/R_PC(7)*1 _d 0
cez end param add

      PARtot = SUM(PAR)

C==== make all bio fields non-negative and compute quotas ==============

      DIC  = MAX(0., Ptr(iDIC))
      NH4  = MAX(0., Ptr(iNH4))
      NO2  = MAX(0., Ptr(iNO2))
      NO3  = MAX(0., Ptr(iNO3))
      PO4  = MAX(EPS, Ptr(iPO4))
      SiO2 = MAX(EPS, Ptr(iSiO2))
      FeT  = MAX(EPS, Ptr(iFeT))

      DOC  = MAX(0., Ptr(iDOC))
      DON  = MAX(0., Ptr(iDON))
      DOP  = MAX(0., Ptr(iDOP))
      DOFe = MAX(0., Ptr(iDOFe))

      PIC  = MAX(0., Ptr(iPIC))
      POC  = MAX(0., Ptr(iPOC))
      PON  = MAX(0., Ptr(iPON))
      POP  = MAX(0., Ptr(iPOP))
      POSi = MAX(0., Ptr(iPOSi))
      POFe = MAX(0., Ptr(iPOFe))
#ifdef GUD_ALLOW_CARBON
      O2   = MAX(0., Ptr(iO2))
#endif
#ifdef GUD_ALLOW_CDOM
      CDOM = MAX(0., Ptr(iCDOM))
#endif

      chlout = 0.0 _d 0

      consumDIC  = 0.0 _d 0
      consumDIC_PIC = 0.0 _d 0
      consumNH4  = 0.0 _d 0
      consumNO2  = 0.0 _d 0
      consumNO3  = 0.0 _d 0
      consumPO4  = 0.0 _d 0
      consumSiO2 = 0.0 _d 0
      consumFeT  = 0.0 _d 0
cez add
      consumPON  = 0.0 _d 0
      consumPOP  = 0.0 _d 0
      consumPOC = 0.0 _d 0
      consumPOFe  = 0.0 _d 0
      consumDON  = 0.0 _d 0
      consumDOP  = 0.0 _d 0
      consumDOC = 0.0 _d 0
      consumDOFe  = 0.0 _d 0
      consumO2 = 0.0 _d 0
      reminPON  = 0.0 _d 0
      reminPOP  = 0.0 _d 0
      reminPOC = 0.0 _d 0
      reminPOFe  = 0.0 _d 0
      reminDON  = 0.0 _d 0
      reminDOP  = 0.0 _d 0
      reminDOC = 0.0 _d 0
      reminDOFe  = 0.0 _d 0
      solubilPON  = 0.0 _d 0
      solubilPOP  = 0.0 _d 0
      solubilPOC = 0.0 _d 0
      solubilPOFe  = 0.0 _d 0
      prodNO2 = 0.0 _d 0
      prodNO3 = 0.0 _d 0
      prodO2 = 0.0 _d 0
cez end add

      DO l=1,gud_nDiag
        diags(l) = 0.0
      ENDDO

C fixed carbon quota, for now 1.0 (may change later)
C other elements: get quota from corresponding ptracer or set to fixed
c ratio if not variable.
      DO j = iMinPhoto, iMaxPhoto
        X = MAX(0.0, Ptr(ic+j-1))
        Qc = 1.0

cez: ========PHYTOPLANKTON===================================

cez: start allowing this part to only run for the two real phytopl
cez        if (j.eq.1.or.j.eq.6.or.j.eq.2.or.j.eq.3) then
cez : phytopl are 1-6:
cez        if (j.ne.7.or.j.ne.8) then       
      if (j.eq.1.or.j.eq.2.or.j.eq.3.or.j.eq.4.or.j.eq.5.or.j.eq.6) then
C==== uptake and nutrient limitation ===================================
C       for quota elements, growth is limiteed by available quota,
C       for non-quota elements, by available nutrients in medium

C       to not use PO4, ..., set ksatPO4=0 and Vmax_PO4=0 (if GUD_ALLOW_PQUOTA)
C       or R_PC=0 (if not)
C       the result will be limitp = 1, uptakePO4 = 0

c phosphorus

         limitp = PO4/(PO4 + ksatPO4(j))
#ifdef GUD_ALLOW_PQUOTA
         Qp  = MAX(EPS*R_PC(j), Ptr(ip+j-1)) / MAX(EPS, X)
         regQ = MAX(0., MIN(1., (Qpmax(j)-Qp)/
     &                          (Qpmax(j)-Qpmin(j)) ))
         uptakePO4 = Vmax_PO4(j) * limitp * regQ *
     &                        uptakeTempFunc * X
c        normalized Droop limitation
         limitp = MAX(0., MIN(1., (1.0-Qpmin(j)/Qp)/
     &                            (1.0-Qpmin(j)/Qpmax(j))))
#endif

c silica
         limitsi = SiO2/(SiO2 + ksatSiO2(j))
#ifdef GUD_ALLOW_SIQUOTA
         Qsi = MAX(EPS*R_SiC(j), Ptr(isi+j-1)) / MAX(EPS, X)
         regQ = MAX(0., MIN(1., (Qsimax(j) - Qsi)/
     &                          (Qsimax(j) - Qsimin(j)) ))
         uptakeSiO2 = Vmax_SiO2(j) * limitsi * regQ *
     &                                  uptakeTempFunc * X

c        linear limitation
         limitsi = MAX(0., MIN(1., (Qsi - Qsimin(j))/
     &                          (Qsimax(j) - Qsimin(j)) ))
#endif

c iron
         limitfe = FeT/(FeT + ksatFeT(j))
#ifdef GUD_ALLOW_FEQUOTA
         Qfe = MAX(EPS*R_FeC(j), Ptr(ife+j-1)) / MAX(EPS, X)
         regQ = MAX(0., MIN(1., (Qfemax(j)-Qfe)/
     &                          (Qfemax(j)-Qfemin(j)) ))
         uptakeFeT = Vmax_FeT(j) * limitfe * regQ *
     &                        uptakeTempFunc * X

c        normalized Droop limitation
         limitfe = MAX(0., MIN(1., (1.0-Qfemin(j)/Qfe)/
     &                             (1.0-Qfemin(j)/Qfemax(j))))
#endif

c nitrogen
#ifdef GUD_ALLOW_NQUOTA
c       have nitrogen quota
        inhibNH4 = EXP(-amminhib(j)*NH4)
        limitNH4 = NH4/(NH4 + ksatNH4(j))
        limitNO2 = NO2/(NO2 + ksatNO2(j))*inhibNH4
        limitNO3 = NO3/(NO3 + ksatNO3(j))*inhibNH4
        Qn  = MAX(EPS*R_NC(j), Ptr(in+j-1)) / MAX(EPS, X)
        regQ = MAX(0., MIN(1., (Qnmax(j)-Qn)/
     &                         (Qnmax(j)-Qnmin(j)) ))
        uptakeNH4 = Vmax_NH4(j)*limitNH4*regQ*uptakeTempFunc*X
        uptakeNO2 = Vmax_NO2(j)*limitNO2*regQ*uptakeTempFunc*X
        uptakeNO3 = Vmax_NO3(j)*limitNO3*regQ*uptakeTempFunc*X
#ifdef GUD_ALLOW_FEQUOTA
        uptakeNO3 = uptakeNO3 * limitfe
#endif

        uptakeN = MAX(uptakeNH4 + uptakeNO2 + uptakeNO3,
     &                Vmax_N(j)*regQ*uptakeTempFunc*X*diazo(j))

c       linear limitation
        limitn = MAX(0., MIN(1., (Qn - Qnmin(j))/
     &                        (Qnmax(j) - Qnmin(j)) ))
#else /* not GUD_ALLOW_NQUOTA */
        Qn = R_NC(j)
        inhibNH4 = EXP(-amminhib(j)*NH4)
        limitNH4 = useNH4(j)*NH4/(NH4 + ksatNH4(j))
        limitNO2 = useNO2(j)*NO2/
     &   (NO2 + combNO(j)*(NO3 + ksatNO3(j) - ksatNO2(j)) + ksatNO2(j))*
     &   inhibNH4
        limitNO3 = useNO3(j)*NO3/
     &   (combNO(j)*NO2 + NO3 + ksatNO3(j))*inhibNH4
        limitn = limitNH4 + limitNO2 + limitNO3
C       normalize to sum (approx) 1
        fracNH4 = limitNH4/(limitn + EPS)
        fracNO2 = limitNO2/(limitn + EPS)
        fracNO3 = limitNO3/(limitn + EPS)
C if diazo, all fracN* == 0 but want no N limitation
        limitn = MIN(1.0, limitn + diazo(j))
#endif /* GUD_ALLOW_NQUOTA */
c        IF (limitn .GT. 0.0) THEN
c          ngrow = ((10*4+2)/(10*4 + 2*limitNH4/limitn +
c     &                       8*limitNO2/limitn + 10*limitNO3/limitn))
c        ELSE
        ngrow = 1.0
c        ENDIF
        
        limitnut = MIN(limitn, limitp, limitsi)
#ifndef GUD_ALLOW_FEQUOTA
        limitnut = MIN(limitnut, limitfe)
#endif

        limitpCO2 = 1.

C==== growth ===========================================================
#ifdef GUD_ALLOW_GEIDER

        alpha_I = 0 _d 0
        DO l = 1, nlam
          alpha_I = alpha_I + alphachl(j,l)*PAR(l)
#ifdef GUD_DEBUG
      if (ig.eq.iDEBUG.and.jg.eq.jDEBUG) then
       print*,'alphaI',myiter,k,j,l,alpha_I,alphachl(j,l),
     &       PAR(l)
      ENDIF
#endif
        ENDDO
C       NB: for quota, PCmax(j) = Vmax_c(j)
        PCm = PCmax(j)*limitnut*photoTempFunc(j)*limitpCO2
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG)
     & print*,'pcm',myiter,k,j,PCmax(j),limitnut,phototempfunc(j),
     &      limitpco2
#endif

        IF (PCm .GT. 0.0) THEN
          acclim = MAX(chl2cmin(j), MIN(chl2cmax(j),
     &             chl2cmax(j)/(1+(chl2cmax(j)*alpha_I)/(2*PCm)) ))
        ELSE
          acclim = chl2cmin(j)
        ENDIF

#ifdef GUD_ALLOW_CHLQUOTA
        QChl=MAX(EPS*R_ChlC(j),Ptr(ichl+j-iMinPhoto))/MAX(EPS,X)
C       quotas are already relative to carbon
        chl2c = QChl
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG)
     & print*,'chl2c',myiter,k,j,chl2c,EPS,R_ChlC(j),
     &     Ptr(ichl+j-iMinPhoto),X

#endif
#else
        chl2c = acclim
#endif

        alpha_I_growth = alpha_I
C a la quota
#ifdef GUD_ALLOW_FEQUOTA
        alpha_I_growth = alpha_I_growth*limitfe
#endif

C       carbon-specific growth rate
C       PC = PCm*(1-EXP(-alpha_I_growth*chl2c/MAX(EPS, PCm)))
        IF (PCm .GT. 0.0 .AND. PARtot .GT. PARmin) THEN
          PC = PCm*(1-EXP(-alpha_I_growth*chl2c/PCm))
        ELSE
          PC = 0.0 _d 0
        ENDIF
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG)
     & print*,'pc',myiter,k,j,PCm,alpha_I_growth,chl2c,PC
#endif

        IF (inhibcoef_geid(j) .GT. 0.0) THEN
C         "total" PAR:
          tmp = alpha_I/alpha_mean(j)
          Ek = PCm/(chl2c*alpha_mean(j))
          EkoverE = Ek / tmp
          IF (tmp .GE. Ek) THEN
            PC = PC*EkoverE*inhibcoef_geid(j)
          ENDIF
        ENDIF

#else /* not GUD_ALLOW_GEIDER */

        IF (PARtot .GT. PARmin) THEN
C         only 1 waveband without GUD_ALLOW_GEIDER
          limitI = (1.0 _d 0 - EXP(-PARtot*ksatPAR(j)))*
     &             EXP(-PARtot*kinhPAR(j)) * normI(j)
          PC = PCmax(j)*limitnut*limitI*photoTempFunc(j)*limitpCO2
        ELSE
          PC = 0.0 _d 0
        ENDIF
        synthChl = 0.0

#endif /* GUD_ALLOW_GEIDER */

        growth = PC*ngrow*X

        uptakeDIC = growth

C non-quota elements are taken up with fixed stoichiometry
#ifndef GUD_ALLOW_NQUOTA
          uptakeN = growth*R_NC(j)
          uptakeNH4 = uptakeN*fracNH4
          uptakeNO2 = uptakeN*fracNO2
          uptakeNO3 = uptakeN*fracNO3
#endif
#ifndef GUD_ALLOW_PQUOTA
         uptakePO4 = growth*R_PC(j)
#endif
#ifndef GUD_ALLOW_SIQUOTA
         uptakeSiO2 = growth*R_SiC(j)
#endif
#ifndef GUD_ALLOW_FEQUOTA
         uptakeFeT = growth*R_FeC(j)
#endif

C==== chlorophyll ======================================================
#ifdef GUD_ALLOW_GEIDER
#ifdef GUD_ALLOW_CHLQUOTA
#ifdef GUD_ALLOW_NQUOTA
C         Geider 1998
          IF (alpha_I*chl2c .GT. 0.0 _d 0) THEN
C            rhochl = Chl2Nmax/(alpha_I*chl2c)*ngrow ???
            rhochl = Chl2Nmax*PC*ngrow/(alpha_I*chl2c)
          ELSE
            rhochl = Chl2Nmax
          ENDIF
          uptakeDIC = uptakeDIC - synthcost*uptakeN
          synthChl = rhochl*uptakeN
#else
#ifdef GUD_GEIDER_RHO_SYNTH
          IF (alpha_I .GT. 0.0 _d 0 .AND. acclim .GT. 0.0 _d 0) THEN
            rhochl = Chl2Cmax(j)*PC*ngrow/(alpha_I*acclim)
          ELSE
            rhochl = 0.0 _d 0    ! should be Chl2Cmax(j) ?????
          ENDIF
          synthChl = rhochl*growth +
     &                        acclimtimescl(j)*(acclim-chl2c)*X
#else
          synthChl = acclim*growth +
     &                        acclimtimescl(j)*(acclim-chl2c)*X
#endif
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
       print*,'rhochl',myiter,k,j,rhochl,chl2cmax(j),pc,alpha_I,
     &      acclim,ngrow
       print*,'synth',myiter,k,j,synthchl,rhochl,growth,
     &       acclimtimescl(j),acclim,chl2c,X
      ENDIF
#endif
#endif /* GUD_ALLOW_NQUOTA */
#else /* GUD_ALLOW_CHLQUOTA */
         chlout(j-iMinPhoto+1) = X*Qc*chl2c
         synthChl = 0.0
#endif /* GUD_ALLOW_CHLQUOTA */
#endif /* GUD_ALLOW_GEIDER */
C=======================================================================
cez: here's the consumption part originally, but i'm moving consumption below the ifs (but w/in the j loop), allowing phytopl and bacteria cumulative effects.

        diags(iPP) = diags(iPP) + growth
        IF (j .LE. nPPplank) THEN
cez: take out the following and replace with specific growth rate        
cez: diags(iPPplank+j-1) = diags(iPPplank+j-1) + growth
        diags(iPPplank+j-1) = diags(iPPplank+j-1) + PC*ngrow
        
        ENDIF
        IF (diazo(j) .GT. 0.0 _d 0) THEN
         diags(iNfix)=diags(iNfix)+uptakeN-uptakeNH4-uptakeNO2-uptakeNO3
        ENDIF

#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG.and.k.eq.kDEBUG) THEN
       print*,'uptakeNO3',myiter,j,X,NO3,uptakeNO3*dT,uptakeN*dT,
     &       fracNO3,limitn,limitNH4,limitNO2,limitNO3,
     &  useNO3(j),combNO(j),NO2,ksatNO3(j),inhibNH4
      ENDIF
#endif
C=======================================================================

       gTr(ic+j-1)=gTr(ic+j-1)  + uptakeDIC
#ifdef GUD_ALLOW_NQUOTA
       gTr(in+j-1)=gTr(in+j-1)  + uptakeN
#endif
#ifdef GUD_ALLOW_PQUOTA
       gTr(ip+j-1)=gTr(ip+j-1)  + uptakePO4
#endif
#ifdef GUD_ALLOW_SIQUOTA
       gTr(isi+j-1)=gTr(isi+j-1) + uptakeSiO2
#endif
#ifdef GUD_ALLOW_FEQUOTA
       gTr(ife+j-1)=gTr(ife+j-1) + uptakeFeT
#endif
#ifdef GUD_ALLOW_CHLQUOTA
       gTr(iChl+j-iMinPhoto)=gTr(iChl+j-iMinPhoto) + synthChl
#endif


#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
       print*,'uptake',myiter,k,j,
     &     uptakeDIC,
     &     uptakeNH4,
     &     uptakeNO2,
     &     uptakeNO3,
     &     uptakeN,
     &     uptakePO4,
     &     uptakeSiO2,
     &     uptakeFeT
      ENDIF
#endif

cez: here is where the j loop ended originally:

cez zero out the new consum, etc. possibilities for bacteria/archaea:
            uptakePOC = 0. _d 0
            uptakePON = 0. _d 0
            uptakePOP = 0. _d 0
            uptakePOFe = 0. _d 0
            hydrolPOC = 0. _d 0
            hydrolPON = 0. _d 0
            hydrolPOP = 0. _d 0
            hydrolPOFe = 0. _d 0
            uptakeDOC = 0. _d 0
            uptakeDON = 0. _d 0
            uptakeDOP = 0. _d 0
            uptakeDOFe = 0. _d 0
            uptakeO2 = 0. _d 0
            respPOC = 0. _d 0
            respPON = 0. _d 0
            respPOP = 0. _d 0
            respPOFe = 0. _d 0
            respDOC = 0. _d 0
            respDON = 0. _d 0
            respDOP = 0. _d 0
            respDOFe = 0. _d 0
            resptoNO2 = 0. _d 0
            resptoNO3 = 0. _d 0
      
cez: add the oxygen production given the new ratio that balances B O2
consumption
            prodO2=prodO2 + growth*R_O2prod 

cez 10/3/16: coopt Denit diagn. to track Vno3P
        diags(iDenit) = diags(iDenit) + uptakeNO3

cez: end the phytoplankton part:
cez: ======== BACTERIA/ARCHAEA ====================================
cez: all bacteria types run through this portion of the code:
cez: if j IS 7, 8, 9, or 10: 
        else

cez: zero out everything first:

cez: original pp uptake:
            uptakeDIC = 0. _d 0 
            uptakeN = 0. _d 0
            uptakeNH4 = 0. _d 0 
            uptakeNO2 = 0. _d 0 
            uptakeNO3 = 0. _d 0 
            uptakePO4 = 0. _d 0 
            uptakeSiO2 = 0. _d 0 
            uptakeFeT = 0. _d 0 
cez: new pp O2 production:
cez: 10/6/16: need to NOT zero this out here!!
cez            prodO2 = 0. _d 0
cez: added substrates for bacteria/archaea:
            uptakePOC = 0. _d 0
            uptakePON = 0. _d 0
            uptakePOP = 0. _d 0
            uptakePOFe = 0. _d 0
            hydrolPOC = 0. _d 0
            hydrolPON = 0. _d 0
            hydrolPOP = 0. _d 0
            hydrolPOFe = 0. _d 0
            uptakeDOC = 0. _d 0
            uptakeDON = 0. _d 0
            uptakeDOP = 0. _d 0
            uptakeDOFe = 0. _d 0
            uptakeO2 = 0. _d 0
            respPOC = 0. _d 0
            respPON = 0. _d 0
            respPOP = 0. _d 0
            respPOFe = 0. _d 0
            respDOC = 0. _d 0
            respDON = 0. _d 0
            respDOP = 0. _d 0
            respDOFe = 0. _d 0
            resptoNO2 = 0. _d 0
            resptoNO3 = 0. _d 0
           
cez: now add specifics for each type:
cez using limitnut as a placeholder, but units of these limits are 1/day (=growth rates)

cez: ===B-TYPE:=== aerobic POM and DOM consuming:
        if (j.eq.7) then       
                limitPON = yod*pmaxOM*PON/(PON+ksatPON)
                limitPOC = yod*pmaxOM*POC/(POC+ksatPON/R_NC(j))
                limitPOP = yod*pmaxOM*POP/(POP+ksatPON/R_NC(j)*R_PC(j))
                limitPOFe = yod*pmaxOM*POFe/(POFe+
     &                       ksatPON/R_NC(j)*R_FeC(j))
                limitDON = yod*pmaxOM*DON/(DON+ksatDON)
                limitDOC = yod*pmaxOM*DOC/(DOC+ksatDON/R_NC(j))
                limitDOP = yod*pmaxOM*DOP/(DOP+ksatDON/R_NC(j)*R_PC(j))
                limitDOFe = yod*pmaxOM*DOFe/(DOFe+
     &                       ksatDON/R_NC(j)*R_FeC(j))
                limitO2 = yoe*pcoefO2*O2

cez: growth is limited collectively by both DOM and POM for each element
                limitnut = MIN(1.,limitPON+limitDON, limitPOC+limitDOC, 
     &                           limitPOP+limitDOP, limitPOFe+limitDOFe,
     &                           limitO2)

cez: growth = u*B, in terms of C (X is mmol/m3 C biomass) 
                growth = limitnut*X*reminTempFunc
                gTr(ic+j-1)=gTr(ic+j-1) + growth
cez: uptake for POM is fraction of POM of total OM pool for each element 
                uptakePOC = growth/yod *
     &                      limitPOC/(limitPOC+limitDOC+EPS)
                uptakePON = growth/yod*R_NC(j) * 
     &                      limitPON/(limitPON+limitDON+EPS)
                uptakePOP = growth/yod*R_PC(j) *
     &                      limitPOP/(limitPOP+limitDOP+EPS)
                uptakePOFe = growth/yod*R_FeC(j) *
     &                       limitPOFe/(limitPOFe+limitDOFe+EPS)
                uptakeDOC = growth/yod *
     &                      limitDOC/(limitPOC+limitDOC+EPS)
                uptakeDON = growth/yod*R_NC(j) * 
     &                      limitDON/(limitPON+limitDON+EPS)
                uptakeDOP = growth/yod*R_PC(j) *
     &                      limitDOP/(limitPOP+limitDOP+EPS)
                uptakeDOFe = growth/yod*R_FeC(j) *
     &                       limitDOFe/(limitPOFe+limitDOFe+EPS)
                           
cez -oxygen is only used for the part of POC that is metabolized:
                uptakeO2 = growth/yoe

cez: Take out hydrolysis of POM to DOM for this model
cez This is the part of POM that is hydrolized into DOM:
cez                hydrolPOC = (alpha- 1. _d 0)*growth/yod
cez                hydrolPON = hydrolPOC*R_NC(j)
cez                hydrolPOP = hydrolPOC*R_PC(j)
cez                hydrolPOFe = hydrolPOC*R_FeC(j)
        
cez These are the bacteria products for remineralization of POM:
                respPOC = growth*(1. _d 0 / yod - 1. _d 0) * 
     &                      limitPOC/(limitPOC+limitDOC+EPS)
                respPON = growth*(1. _d 0 / yod - 1. _d 0)*R_NC(j) * 
     &                      limitPON/(limitPON+limitDON+EPS)
                respPOP = growth*(1. _d 0 / yod - 1. _d 0)*R_PC(j) *
     &                      limitPOP/(limitPOP+limitDOP+EPS)
                respPOFe = growth*(1. _d 0 / yod - 1. _d 0)*R_FeC(j) *
     &                       limitPOFe/(limitPOFe+limitDOFe+EPS)
                respDOC = growth*(1. _d 0 / yod - 1. _d 0) *
     &                      limitDOC/(limitPOC+limitDOC+EPS)
                respDON = growth*(1. _d 0 / yod - 1. _d 0)*R_NC(j) * 
     &                      limitDON/(limitPON+limitDON+EPS)
                respDOP = growth*(1. _d 0 / yod - 1. _d 0)*R_PC(j) *
     &                      limitDOP/(limitPOP+limitDOP+EPS)
                respDOFe = growth*(1. _d 0 / yod - 1. _d 0)*R_FeC(j) *
     &                       limitDOFe/(limitPOFe+limitDOFe+EPS)


cez: ===B-TYPE:=== nitrate reducing (NO3->NO2) POM and DOM consuming:
        else if (j.eq.8) then
                limitPON = ynd*pmaxOM*PON/(PON+ksatPON)
                limitPOC = ynd*pmaxOM*POC/(POC+ksatPON/R_NC(j))
                limitPOP = ynd*pmaxOM*POP/(POP+ksatPON/R_NC(j)*R_PC(j))
                limitPOFe = ynd*pmaxOM*POFe/(POFe+
     &                       ksatPON/R_NC(j)*R_FeC(j))
                limitDON = ynd*pmaxOM*DON/(DON+ksatDON)
                limitDOC = ynd*pmaxOM*DOC/(DOC+ksatDON/R_NC(j))
                limitDOP = ynd*pmaxOM*DOP/(DOP+ksatDON/R_NC(j)*R_PC(j))
                limitDOFe = ynd*pmaxOM*DOFe/(DOFe+
     &                       ksatDON/R_NC(j)*R_FeC(j))
                limitNO3 = yne1*pmaxDINden*NO3/(NO3 + ksatDIN)

cez: growth is limited collectively by both DOM and POM for each element
                limitnut = MIN(1.,limitPON+limitDON, limitPOC+limitDOC, 
     &                           limitPOP+limitDOP, limitPOFe+limitDOFe,
     &                           limitNO3)

cez: growth = u*B, in terms of C (X is mmol/m3 C biomass) 
                growth = limitnut*X*reminTempFunc
                gTr(ic+j-1)=gTr(ic+j-1) + growth
cez: uptake for POM is fraction of POM of total OM pool for each element 
                uptakePOC = growth/ynd *
     &                      limitPOC/(limitPOC+limitDOC+EPS)
                uptakePON = growth/ynd*R_NC(j) * 
     &                      limitPON/(limitPON+limitDON+EPS)
                uptakePOP = growth/ynd*R_PC(j) *
     &                      limitPOP/(limitPOP+limitDOP+EPS)
                uptakePOFe = growth/ynd*R_FeC(j) *
     &                       limitPOFe/(limitPOFe+limitDOFe+EPS)
                uptakeDOC = growth/ynd *
     &                      limitDOC/(limitPOC+limitDOC+EPS)
                uptakeDON = growth/ynd*R_NC(j) * 
     &                      limitDON/(limitPON+limitDON+EPS)
                uptakeDOP = growth/ynd*R_PC(j) *
     &                      limitDOP/(limitPOP+limitDOP+EPS)
                uptakeDOFe = growth/ynd*R_FeC(j) *
     &                       limitDOFe/(limitPOFe+limitDOFe+EPS)

                uptakeNO3 = growth/yne1

cez These are the bacteria products for remineralization:
                respPOC = growth*(1. _d 0 / ynd - 1. _d 0) * 
     &                      limitPOC/(limitPOC+limitDOC+EPS)
                respPON = growth*(1. _d 0 / ynd - 1. _d 0)*R_NC(j) * 
     &                      limitPON/(limitPON+limitDON+EPS)
                respPOP = growth*(1. _d 0 / ynd - 1. _d 0)*R_PC(j) *
     &                      limitPOP/(limitPOP+limitDOP+EPS)
                respPOFe = growth*(1. _d 0 / ynd - 1. _d 0)*R_FeC(j) *
     &                       limitPOFe/(limitPOFe+limitDOFe+EPS)
                respDOC = growth*(1. _d 0 / ynd - 1. _d 0) *
     &                      limitDOC/(limitPOC+limitDOC+EPS)
                respDON = growth*(1. _d 0 / ynd - 1. _d 0)*R_NC(j) * 
     &                      limitDON/(limitPON+limitDON+EPS)
                respDOP = growth*(1. _d 0 / ynd - 1. _d 0)*R_PC(j) *
     &                      limitDOP/(limitPOP+limitDOP+EPS)
                respDOFe = growth*(1. _d 0 / ynd - 1. _d 0)*R_FeC(j) *
     &                       limitDOFe/(limitPOFe+limitDOFe+EPS)

cez: add NO2 production from nitrate reduction:
                resptoNO2 = growth/yne1

cez: ===B-TYPE:=== denitrifying (NO2->N2) POM and DOM consuming 
        else if (j.eq.9) then
                limitPON = ynd*pmaxOM*PON/(PON+ksatPON)
                limitPOC = ynd*pmaxOM*POC/(POC+ksatPON/R_NC(j))
                limitPOP = ynd*pmaxOM*POP/(POP+ksatPON/R_NC(j)*R_PC(j))
                limitPOFe = ynd*pmaxOM*POFe/(POFe+
     &                       ksatPON/R_NC(j)*R_FeC(j))
                limitDON = ynd*pmaxOM*DON/(DON+ksatDON)
                limitDOC = ynd*pmaxOM*DOC/(DOC+ksatDON/R_NC(j))
                limitDOP = ynd*pmaxOM*DOP/(DOP+ksatDON/R_NC(j)*R_PC(j))
                limitDOFe = ynd*pmaxOM*DOFe/(DOFe+
     &                       ksatDON/R_NC(j)*R_FeC(j))
                limitNO2 = yne2*pmaxDINden*NO2/(NO2 + ksatDIN)

cez: growth is limited collectively by both DOM and POM for each element
                limitnut = MIN(1.,limitPON+limitDON, limitPOC+limitDOC, 
     &                           limitPOP+limitDOP, limitPOFe+limitDOFe,
     &                           limitNO2)

cez: growth = u*B, in terms of C (X is mmol/m3 C biomass) 
                growth = limitnut*X*reminTempFunc
                gTr(ic+j-1)=gTr(ic+j-1) + growth
cez: uptake for POM is fraction of POM of total OM pool for each element 
                uptakePOC = growth/ynd *
     &                      limitPOC/(limitPOC+limitDOC+EPS)
                uptakePON = growth/ynd*R_NC(j) * 
     &                      limitPON/(limitPON+limitDON+EPS)
                uptakePOP = growth/ynd*R_PC(j) *
     &                      limitPOP/(limitPOP+limitDOP+EPS)
                uptakePOFe = growth/ynd*R_FeC(j) *
     &                       limitPOFe/(limitPOFe+limitDOFe+EPS)
                uptakeDOC = growth/ynd *
     &                      limitDOC/(limitPOC+limitDOC+EPS)
                uptakeDON = growth/ynd*R_NC(j) * 
     &                      limitDON/(limitPON+limitDON+EPS)
                uptakeDOP = growth/ynd*R_PC(j) *
     &                      limitDOP/(limitPOP+limitDOP+EPS)
                uptakeDOFe = growth/ynd*R_FeC(j) *
     &                       limitDOFe/(limitPOFe+limitDOFe+EPS)

                uptakeNO2 = growth/yne2

cez: no N2 production yet:  add N2 production from denitr:
cez:                resptoN2 = growth/yne2 (would be this)
 
cez These are the bacteria products for remineralization:
                respPOC = growth*(1. _d 0 / ynd - 1. _d 0) * 
     &                      limitPOC/(limitPOC+limitDOC+EPS)
                respPON = growth*(1. _d 0 / ynd - 1. _d 0)*R_NC(j) * 
     &                      limitPON/(limitPON+limitDON+EPS)
                respPOP = growth*(1. _d 0 / ynd - 1. _d 0)*R_PC(j) *
     &                      limitPOP/(limitPOP+limitDOP+EPS)
                respPOFe = growth*(1. _d 0 / ynd - 1. _d 0)*R_FeC(j) *
     &                       limitPOFe/(limitPOFe+limitDOFe+EPS)
                respDOC = growth*(1. _d 0 / ynd - 1. _d 0) *
     &                      limitDOC/(limitPOC+limitDOC+EPS)
                respDON = growth*(1. _d 0 / ynd - 1. _d 0)*R_NC(j) * 
     &                      limitDON/(limitPON+limitDON+EPS)
                respDOP = growth*(1. _d 0 / ynd - 1. _d 0)*R_PC(j) *
     &                      limitDOP/(limitPOP+limitDOP+EPS)
                respDOFe = growth*(1. _d 0 / ynd - 1. _d 0)*R_FeC(j) *
     &                       limitDOFe/(limitPOFe+limitDOFe+EPS)

cez: ===B-TYPE:=== anammox!! 
        else if (j.eq.10) then
cez                limitNH4 = ynh4anmx*convpmaxDIN*pmaxDIN/tentothethird*
cez     &                   NH4/(NH4 + ksatDIN*tentothethird)*reminTempFunc
cez                limitNO2 = yno2anmx*convpmaxDIN*pmaxDIN/tentothethird*
cez     &                   NO2/(NO2 + ksatDIN*tentothethird)*reminTempFunc
                limitNH4 = ynh4anmx*convpmaxDIN*pmaxDIN*
     &                   NH4/(NH4 + ksatDIN)*reminTempFunc
                limitNO2 = yno2anmx*convpmaxDIN*pmaxDIN*
     &                   NO2/(NO2 + ksatDIN)*reminTempFunc
                limitnut = MIN(1.,limitNH4, limitNO2)
                
                growth = limitnut*X
                gTr(ic+j-1)=gTr(ic+j-1) + growth
              
                uptakeNH4 = growth/ynh4anmx
                uptakeNO2 = growth/yno2anmx
                uptakeDIC = growth
                uptakePO4 = uptakeDIC*R_PC(j)
                uptakeFeT = uptakeDIC*R_FeC(j)
                
cez: tally up denitrification:
cez        diags(iDenit) = denit
cez                diags(iDenit) = diags(iDenit) + uptakeNO3

cez: ===B-TYPE:=== NH4 oxidizer:
        else if (j.eq.11) then
                limitNH4 = ynh4*convpmaxDIN*pmaxDIN*NH4/(NH4+ksatDIN)*
     &                      reminTempFunc 
                limitO2 = yonh4*pcoefO2*O2
                limitnut = MIN(1.,limitNH4, limitO2)
          
                growth = limitnut*X
                gTr(ic+j-1)=gTr(ic+j-1) + growth
          
                uptakeNH4 = growth/ynh4
                uptakeDIC = growth
                uptakePO4 = uptakeDIC*R_PC(j)
                uptakeFeT = uptakeDIC*R_FeC(j)
                uptakeO2 = growth/yonh4
cez:                resptoNO2 = growth*(1. _d 0 / ynh4 - 1. _d 0)
                resptoNO2 = growth*(1. _d 0 / ynh4 - R_NC(j))
        
cez: ===B-TYPE:=== NO2 oxidizer:
        else if (j.eq.12) then        
                limitNO2 = yno2*convpmaxDIN*pmaxDIN/tentothethird*
     &                     NO2/(NO2+ksatDIN*tentothethird)*
     &                      reminTempFunc 
cez: no size:               limitNO2 = yno2*convpmaxDIN*pmaxDIN*NO2/(NO2+ksatDIN)
                limitO2 = yono2*pcoefO2*O2
                limitnut = MIN(1.,limitNO2, limitO2)
          
                growth = limitnut*X
                gTr(ic+j-1)=gTr(ic+j-1) + growth
                
                uptakeNO2 = growth/yno2
                uptakeDIC = growth
                uptakePO4 = uptakeDIC*R_PC(j)
                uptakeFeT = uptakeDIC*R_FeC(j)
                uptakeO2 = growth/yono2
cez:                resptoNO3 = growth*(1. _d 0 / yno2 - 1. _d 0)  
                resptoNO3 = growth*(1. _d 0 / yno2 - R_NC(j))  
 
cez: end the specifics for each type:                
        end if

        IF (j .LE. nPPplank) THEN
cez counterintuitively, "growth" = u*B = mol/m3/day and limitnut = 1/day
cez          diags(iPPplank+j-1) = diags(iPPplank+j-1) + growth
          diags(iPPplank+j-1) = diags(iPPplank+j-1) + limitnut
        ENDIF

cez end the bacteria/archaea code:
        end if 


C========== CUMULATIVE CONSUM, REMIN, and PROD ================================
cez: now, within the j loop still, accumulate consumption by pp and bacter (which was originally here for pp), the remineralization and production (named historically to refer to OM and DIN, respectively) by bacteria/archaea.

cez: words: uptake, resp, solubil, and respto all are collected here

cez: this next bit (9 lines) was the original phytoplankton code:
cez: now, the uptake of NO2, NO3, and DIC, etc. by bact. are possible
C=======================================================================
        consumDIC_PIC = consumDIC_PIC + uptakeDIC*R_PICPOC(j)
        consumDIC  = consumDIC  + uptakeDIC 
        consumNH4  = consumNH4  + uptakeNH4 
        consumNO2  = consumNO2  + uptakeNO2 
        consumNO3  = consumNO3  + uptakeNO3 
        consumPO4  = consumPO4  + uptakePO4 
        consumSiO2 = consumSiO2 + uptakeSiO2
        consumFeT  = consumFeT  + uptakeFeT
 
cez add B consum and accumulating remin, and prod:
        consumO2 = consumO2 + uptakeO2

cez this part is necessary, since denitrifiers also consume DOM!:
        consumDOC = consumDOC + uptakeDOC
        consumDON = consumDON + uptakeDON
        consumDOP = consumDOP + uptakeDOP
        consumDOFe = consumDOFe + uptakeDOFe

cez this part is for now unnecessary, since only this one type consumes POM:
        consumPOC = consumPOC + uptakePOC
        consumPON = consumPON + uptakePON
        consumPOP = consumPOP + uptakePOP
        consumPOFe = consumPOFe + uptakePOFe

        reminPOC = reminPOC + respPOC
        reminPON = reminPON + respPON
        reminPOP = reminPOP + respPOP
        reminPOFe = reminPOFe + respPOFe
        
cez        solubilPOC = solubilPOC + hydrolPOC
cez        solubilPON = solubilPON + hydrolPON
cez        solubilPOP = solubilPOP + hydrolPOP
cez        solubilPOFe = solubilPOFe + hydrolPOFe

        reminDOC = reminDOC + respDOC
        reminDON = reminDON + respDON
        reminDOP = reminDOP + respDOP
        reminDOFe = reminDOFe + respDOFe

cez:  nitrifying parameterization here
        prodNO2 = prodNO2 + resptoNO2
        prodNO3 = prodNO3 + resptoNO3

C     j loop end
      ENDDO

      gTr(iDIC )=gTr(iDIC ) - consumDIC - consumDIC_PIC
      gTr(iNH4 )=gTr(iNH4 ) - consumNH4
      gTr(iNO2 )=gTr(iNO2 ) - consumNO2
      gTr(iNO3 )=gTr(iNO3 ) - consumNO3
      gTr(iPO4 )=gTr(iPO4 ) - consumPO4
      gTr(iSiO2)=gTr(iSiO2) - consumSiO2
      gTr(iFeT )=gTr(iFeT ) - consumFeT

cez Take out the following, replace with bacteria reminP products defined above in j loop
cez      reminDOC  = reminTempFunc*KDOC *DOC
cez      reminDON  = reminTempFunc*KDON *DON
cez      reminDOP  = reminTempFunc*KDOP *DOP
cez      reminDOFe = reminTempFunc*KDOFe*DOFe
cez      reminPOC  = reminTempFunc*KPOC *POC
cez      reminPON  = reminTempFunc*KPON *PON
cez      reminPOP  = reminTempFunc*KPOP *POP
      reminPOSi = reminTempFunc*KPOSi*POSi
cez      reminPOFe = reminTempFunc*KPOFe*POFe
      disscPIC = Kdissc*PIC
      
c nitrogen chemistry
c NH4 -> NO2 -> NO3 by bacterial action
cez replace this by nitrifier prodNO2 and prodNO3 above
cez      prodNO2 = knita*NH4
cez      prodNO3 = knitb*NO2
cez      IF (PAR_oxi .NE. 0.0 _d 0) THEN
cezc        print*,'prod',prodNO2,PARtot,PAR_oxi
cez        prodNO2 = prodNO2*MAX(0.0, 1.0 - PARtot/PAR_oxi)
cez        prodNO3 = prodNO3*MAX(0.0, 1.0 - PARtot/PAR_oxi)
cez      ENDIF

#ifdef GUD_ALLOW_CDOM
      reminPOP_CDOM = fracCDOM*reminPOP
      reminPOC_CDOM = R_CP_CDOM*reminPOP_CDOM
      reminPON_CDOM = R_NP_CDOM*reminPOP_CDOM
      reminPOFe_CDOM = R_FeP_CDOM*reminPOP_CDOM
c degradation of  CDOM - high when bleached by light
      degrCDOM_DOP = reminTempFunc*CDOM*
     &     (CDOMdegrd+CDOMbleach*MIN(1.0 _d 0, PARtot/PARCDOM))
      degrCDOM_DOC  = R_CP_CDOM  * degrCDOM_DOP
      degrCDOM_DON  = R_NP_CDOM  * degrCDOM_DOP
      degrCDOM_DOFe = R_FeP_CDOM * degrCDOM_DOP
#endif

#ifdef GUD_ALLOW_DENIT
cez      IF (O2 .LT. O2crit .AND. NO3 .LT. NO3crit) THEN
cezc no remineralization for N, P, Fe (not Si?)
cez        reminDOC = 0. _d 0
cez        reminDOP  = 0. _d 0
cez        reminDON  = 0. _d 0
cez        reminDOFe = 0. _d 0
cez        reminPOC = 0. _d 0
cez        reminPOP  = 0. _d 0
cez        reminPON  = 0. _d 0
cez        reminPOFe = 0. _d 0
cez#ifdef GUD_ALLOW_CDOM
cez        reminPOC_cdom = 0. _d 0
cez        reminPOP_cdom  = 0. _d 0
cez        reminPON_cdom  = 0. _d 0
cez        reminPOFe_cdom = 0. _d 0
cez        degrCDOM_DOP = reminTempFunc*CDOM*
cez     &                 CDOMbleach*MIN(1.0 _d 0, PARtot/PARCDOM)
cez        degrCDOM_DON  = R_NP_CDOM  * degrCDOM_DOP
cez        degrCDOM_DOFe = R_FeP_CDOM * degrCDOM_DOP
cez        degrCDOM_DOC  = R_CP_CDOM  * degrCDOM_DOP
cez#endif
cez      ENDIF
#endif /* GUD_ALLOW_DENIT */

C==== tendencies =======================================================

#ifdef GUD_ALLOW_CARBON
c production of O2 by photosynthesis
cez take out old account of O2 production for new balanced:
cez      gTr(iO2  )=gTr(iO2  ) + R_OP*consumPO4
      gTr(iO2  )=gTr(iO2  ) + prodO2
c loss of O2 by remineralization
cez:      IF (O2 .GT. O2crit) THEN
cez        gTr(iO2)=gTr(iO2) - R_OP*reminDOP
        gTr(iO2)=gTr(iO2) - consumO2
#ifndef GUD_ALLOW_CDOM
cez modify but check this- consumO2 and R_OP still need to match:
cez        gTr(iO2)=gTr(iO2) - R_OP*reminPOP
cez        gTr(iO2)=gTr(iO2) - consumO2
#endif
cez      ENDIF

      gTr(iALK)=gTr(iALK) - (prodNO3 - consumNO3)
     &                    - 2.0 _d 0*(consumDIC_PIC - disscPIC)
#endif /* GUD_ALLOW_CARBON */

      gTr(iDIC )=gTr(iDIC ) + reminDOC + disscPIC
cez: consumption of NH4 and NO2 is already accounted for above
cez      gTr(iNH4 )=gTr(iNH4 ) + reminDON - prodNO2
cez      gTr(iNO2 )=gTr(iNO2 ) + prodNO2 - prodNO3
      gTr(iNH4 )=gTr(iNH4 ) + reminDON
      gTr(iNO2 )=gTr(iNO2 ) + prodNO2
      gTr(iNO3 )=gTr(iNO3 ) + prodNO3
cez 10/3/16: bc i don't know if DenitN is different than Denit, comment out
cez      diags(iDenitN) = 0.0
#ifdef GUD_ALLOW_DENIT
cez      IF (O2 .LT. O2crit) THEN
cez        denitNH4 = reminDON
cez        denit = denit_NP*reminDOP
cez#ifndef GUD_ALLOW_CDOM
cez        denitNH4 = denitNH4 + reminPON
cez        denit = denit + denit_NP*reminPOP
cez#endif
cez        diags(iDenit) = denit
cez        gTr(iNH4)=gTr(iNH4) - denitNH4
cez        gTr(iNO3)=gTr(iNO3) - denit_NO3/denit_np*denit
cez        gTr(iALK)=gTr(iALK) + denit_NO3/denit_np*denit
cez        diags(iDenitN) = denitNH4 + denit_NO3/denit_np*denit
cez      ENDIF
#endif /* GUD_ALLOW_DENIT */
#ifdef GUD_DEBUG
cez      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG.and.k.eq.kDEBUG) THEN
cez       print*,'consumNO3',myiter,k,NO3,consumNO3*10800
cez #ifdef GUD_ALLOW_DENIT
cez       IF (O2 .LT. O2crit)
cez     &       print*,'denitNO3',denit_NO3/denit_np*denit*10800
cez #endif
cez      ENDIF
#endif

      gTr(iPO4 )=gTr(iPO4 ) + reminDOP
      gTr(iFeT )=gTr(iFeT ) + reminDOFe
      gTr(iSiO2)=gTr(iSiO2)             + reminPOSi

cez- hack mods begin: replace remin uptake with B consumption:
cez      gTr(iDOC )=gTr(iDOC ) - reminDOC
cez      gTr(iDON )=gTr(iDON ) - reminDON
cez      gTr(iDOP )=gTr(iDOP ) - reminDOP
cez      gTr(iDOFe)=gTr(iDOFe) - reminDOFe

cez DOC is created by #4 PA-assoc solubilization and consumed by #5
cez      gTr(iDOC )=gTr(iDOC ) + solubilPOC - consumDOC
cez      gTr(iDON )=gTr(iDON ) + solubilPON - consumDON
cez      gTr(iDOP )=gTr(iDOP ) + solubilPOP - consumDOP
cez      gTr(iDOFe)=gTr(iDOFe) + solubilPOFe - consumDOFe
      gTr(iDOC )=gTr(iDOC ) - consumDOC
      gTr(iDON )=gTr(iDON ) - consumDON
      gTr(iDOP )=gTr(iDOP ) - consumDOP
      gTr(iDOFe)=gTr(iDOFe) - consumDOFe
       
       gTr(iPIC )=gTr(iPIC ) - disscPIC
cez      gTr(iPOC )=gTr(iPOC ) - reminPOC
cez      gTr(iPON )=gTr(iPON ) - reminPON       
cez      gTr(iPOP )=gTr(iPOP ) - reminPOP       
       gTr(iPOSi)=gTr(iPOSi) - reminPOSi
cez      gTr(iPOFe)=gTr(iPOFe) - reminPOFe
       
       gTr(iPOC )=gTr(iPOC ) - consumPOC
       gTr(iPON )=gTr(iPON ) - consumPON
       gTr(iPOP )=gTr(iPOP ) - consumPOP
       gTr(iPOFe)=gTr(iPOFe) - consumPOFe
       
#ifdef GUD_ALLOW_CDOM
      gTr(iDOC )=gTr(iDOC ) + reminPOC  - reminPOC_CDOM  + degrCDOM_DOC
      gTr(iDON )=gTr(iDON ) + reminPON  - reminPON_CDOM  + degrCDOM_DON
      gTr(iDOP )=gTr(iDOP ) + reminPOP  - reminPOP_CDOM  + degrCDOM_DOP
      gTr(iDOFe)=gTr(iDOFe) + reminPOFe - reminPOFe_CDOM + degrCDOM_DOFe

      gTr(iCDOM)=gTr(iCDOM) + reminPOP_CDOM - degrCDOM_DOP
#else
      gTr(iDIC )=gTr(iDIC ) + reminPOC
      gTr(iNH4 )=gTr(iNH4 ) + reminPON
      gTr(iPO4 )=gTr(iPO4 ) + reminPOP
      gTr(iFeT )=gTr(iFeT ) + reminPOFe
#endif /* GUD_ALLOW_CDOM */

#endif /* ALLOW_GUD */

      RETURN
      END SUBROUTINE
