$PARAM @annotated
KDT     : 0.0006583  : Duodenum to plasma transferrin rate constant (1/h)
KDC     : 0.00001396 : Duodenum to plasma citrate rate constant (1/h)
KPTPC   : 0.9191     : Plasma transferrin to plasma citrate rate constant (1/h)
KPCPT   : 14.442     : Plasma citrate to plasma transferrin rate constant (1/h)
KPCIC   : 20.00      : Plasma citrate to interstitial fluid citrate rate constant (1/h)
KICPC   : 5.000      : Interstitial fluid citrate to plasma citrate rate constant (1/h)
KICIT   : 18.03      : Interstitial fluid citrate to interstitial fluid transferrin rate constant (1/h)
KITIC   : 3.676      : Interstitial fluid transferrin to interstitial fluid citrate rate constant (1/h)
KLIVIN  : 1.084      : Plasma transferrin to liver rate constant (1/h)
KLIVOUT : 0.1521     : Liver to plasma transferrin rate constant (1/h)
KMUSIN  : 0.1074     : Interstitial fluid transferrin to muscle rate constant (1/h)
KMUSOUT : 0.001845   : Muscle to interstitial fluid transferrin rate constant (1/h)
KBONIN  : 0.2135     : Interstitial fluid citrate to bone rate constant (1/h)
KBONOUT : 0.00003583 : Bone to interstitial fluid citrate rate constant (1/h)
KU      : 2.752      : Plasma citrate to urine rate constant (1/h)
KTD     : 0.01444    : Plasma transferrin to duodenum rate constant (1/h)
KSD     : 4.000      : Stomach to duodenum rate constant (1/h)
KDR     : 0.2708     : Duodenum to residual intestinal track rate constant (1/h)
KRS     : 0.02778    : Residual intestinal track to faeces rate constant (1/h)

$CMT @annotated
STO : Stomach
DUO : Duodenum
RES : Residual intestinal track
PT  : Plasma transferrin
PC  : Plasma citrate
IT  : Interstitial fluid transferrin
IC  : Interstitial fluid citrate
LIV : Liver
MUS : Muscle
BON : Bone
U   : Urine
FAE : Faeces

$GLOBAL
#define PLASMA (PT + PC)
#define INTERSTITUALFLUID (IC + IT)
#define RETENTION (STO + DUO + RES + PT + PC + IT + IC + LIV + MUS + BON)

$ODE
dxdt_STO = -KSD*STO;
dxdt_DUO = KSD*STO - KDR *DUO - KDC*DUO - KDT*DUO + KTD*PT;
dxdt_RES = KDR*DUO - KRS*RES;
dxdt_PT  = KDT*DUO - KTD*PT - KPTPC*PT + KPCPT*PC - KLIVIN*PT + KLIVOUT*LIV;
dxdt_PC  = KDC*DUO - KPCPT*PC + KPTPC*PT - KPCIC*PC + KICPC*IC - KU*PC;
dxdt_IC  = -KICPC*IC + KPCIC*PC - KICIT*IC + KITIC*IT - KBONIN*IC + KBONOUT*BON;
dxdt_IT  = -KITIC*IT + KICIT*IC - KMUSIN*IT + KMUSOUT*MUS;
dxdt_LIV = KLIVIN*PT - KLIVOUT*LIV;
dxdt_MUS = KMUSIN*IT - KMUSOUT*MUS;
dxdt_BON = KBONIN*IC - KBONOUT*BON;
dxdt_U   = KU*PC;
dxdt_FAE = KRS*RES; 

$CAPTURE STO DUO RES FAE PT PC IT IC LIV MUS BON U PLASMA INTERSTITUALFLUID RETENTION