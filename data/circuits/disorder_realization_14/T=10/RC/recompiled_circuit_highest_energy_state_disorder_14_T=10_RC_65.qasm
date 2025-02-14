OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90647107) q[0];
sx q[0];
rz(4.8405092) q[0];
sx q[0];
rz(9.7066896) q[0];
rz(-2.6126722) q[1];
sx q[1];
rz(-1.6493874) q[1];
sx q[1];
rz(1.563501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603005) q[0];
sx q[0];
rz(-1.2055802) q[0];
sx q[0];
rz(-2.9486548) q[0];
rz(0.27484244) q[2];
sx q[2];
rz(-1.5656398) q[2];
sx q[2];
rz(2.8066563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7719745) q[1];
sx q[1];
rz(-1.6965515) q[1];
sx q[1];
rz(-2.5096276) q[1];
rz(2.7759477) q[3];
sx q[3];
rz(-2.4198101) q[3];
sx q[3];
rz(0.66058285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6866744) q[2];
sx q[2];
rz(-0.29759559) q[2];
sx q[2];
rz(2.5285524) q[2];
rz(-0.4736627) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(-1.4325498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365725) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(-0.75102425) q[0];
rz(0.48149064) q[1];
sx q[1];
rz(-2.057169) q[1];
sx q[1];
rz(2.1717333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3595122) q[0];
sx q[0];
rz(-0.94612288) q[0];
sx q[0];
rz(-1.1083353) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9578147) q[2];
sx q[2];
rz(-0.2193976) q[2];
sx q[2];
rz(0.32599005) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0326421) q[1];
sx q[1];
rz(-1.3130929) q[1];
sx q[1];
rz(2.2326681) q[1];
rz(-2.5447846) q[3];
sx q[3];
rz(-1.1973901) q[3];
sx q[3];
rz(1.663409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91886175) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(2.6027423) q[2];
rz(-0.017596267) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(-0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70055145) q[0];
sx q[0];
rz(-0.38527641) q[0];
sx q[0];
rz(-3.0803296) q[0];
rz(1.5047269) q[1];
sx q[1];
rz(-0.69925362) q[1];
sx q[1];
rz(1.1963074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6590779) q[0];
sx q[0];
rz(-1.3318303) q[0];
sx q[0];
rz(-0.84521291) q[0];
rz(2.8792686) q[2];
sx q[2];
rz(-1.50365) q[2];
sx q[2];
rz(-0.1131499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7972826) q[1];
sx q[1];
rz(-2.4786665) q[1];
sx q[1];
rz(-1.0271038) q[1];
x q[2];
rz(0.34516224) q[3];
sx q[3];
rz(-0.90800873) q[3];
sx q[3];
rz(-1.8927285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68941826) q[2];
sx q[2];
rz(-1.8345366) q[2];
sx q[2];
rz(0.43251953) q[2];
rz(-2.050926) q[3];
sx q[3];
rz(-0.64041036) q[3];
sx q[3];
rz(-1.0961016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5525621) q[0];
sx q[0];
rz(-0.93695372) q[0];
sx q[0];
rz(2.9402148) q[0];
rz(0.51678139) q[1];
sx q[1];
rz(-2.7613381) q[1];
sx q[1];
rz(2.1929599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059221642) q[0];
sx q[0];
rz(-1.2223402) q[0];
sx q[0];
rz(-2.2414464) q[0];
rz(-pi) q[1];
rz(-2.2791512) q[2];
sx q[2];
rz(-1.3459599) q[2];
sx q[2];
rz(2.1612957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3624787) q[1];
sx q[1];
rz(-2.6370905) q[1];
sx q[1];
rz(0.027536784) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6679156) q[3];
sx q[3];
rz(-1.0518952) q[3];
sx q[3];
rz(3.1176381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2739233) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(-2.5841827) q[2];
rz(-3.0607767) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(-1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4163365) q[0];
sx q[0];
rz(-1.6660322) q[0];
sx q[0];
rz(1.0303372) q[0];
rz(-0.15580767) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(0.20419289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9568528) q[0];
sx q[0];
rz(-2.889688) q[0];
sx q[0];
rz(-0.27916081) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7014808) q[2];
sx q[2];
rz(-1.3595264) q[2];
sx q[2];
rz(1.7344765) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0597282) q[1];
sx q[1];
rz(-1.7769741) q[1];
sx q[1];
rz(-2.6463406) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24498265) q[3];
sx q[3];
rz(-2.3534273) q[3];
sx q[3];
rz(0.80313166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2121409) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(0.41352752) q[2];
rz(-0.84189576) q[3];
sx q[3];
rz(-1.7588408) q[3];
sx q[3];
rz(0.97257417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464762) q[0];
sx q[0];
rz(-2.0501417) q[0];
sx q[0];
rz(0.0055775642) q[0];
rz(-1.7976409) q[1];
sx q[1];
rz(-2.9426212) q[1];
sx q[1];
rz(-0.70095789) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85421766) q[0];
sx q[0];
rz(-1.1297884) q[0];
sx q[0];
rz(-2.1483351) q[0];
rz(-pi) q[1];
rz(1.0089325) q[2];
sx q[2];
rz(-2.418926) q[2];
sx q[2];
rz(0.36107963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.720916) q[1];
sx q[1];
rz(-1.6117967) q[1];
sx q[1];
rz(-2.8662205) q[1];
rz(-pi) q[2];
rz(-0.72356059) q[3];
sx q[3];
rz(-0.71092194) q[3];
sx q[3];
rz(-1.5462745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31412101) q[2];
sx q[2];
rz(-2.8077329) q[2];
sx q[2];
rz(-1.1673048) q[2];
rz(2.2589034) q[3];
sx q[3];
rz(-0.76789951) q[3];
sx q[3];
rz(-2.816443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2348787) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(3.0562905) q[0];
rz(0.75434297) q[1];
sx q[1];
rz(-0.5032379) q[1];
sx q[1];
rz(1.0275966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9011544) q[0];
sx q[0];
rz(-1.580997) q[0];
sx q[0];
rz(-2.9290694) q[0];
x q[1];
rz(0.0079389056) q[2];
sx q[2];
rz(-2.5449341) q[2];
sx q[2];
rz(-2.3158429) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4484549) q[1];
sx q[1];
rz(-1.692728) q[1];
sx q[1];
rz(-1.8961402) q[1];
rz(-2.2619234) q[3];
sx q[3];
rz(-0.68118775) q[3];
sx q[3];
rz(-2.3613514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.1520285) q[2];
sx q[2];
rz(-2.1542408) q[2];
sx q[2];
rz(-2.7222471) q[2];
rz(0.63465214) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(1.1254719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3318876) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(0.69212717) q[0];
rz(2.5634815) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(-1.7587761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22631881) q[0];
sx q[0];
rz(-1.5902441) q[0];
sx q[0];
rz(1.5032728) q[0];
rz(-pi) q[1];
rz(-2.369129) q[2];
sx q[2];
rz(-2.303745) q[2];
sx q[2];
rz(-0.39847429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9005712) q[1];
sx q[1];
rz(-1.3374778) q[1];
sx q[1];
rz(1.0020578) q[1];
rz(-pi) q[2];
rz(2.2949831) q[3];
sx q[3];
rz(-1.203152) q[3];
sx q[3];
rz(1.6366307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3708923) q[2];
sx q[2];
rz(-0.53052491) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(0.51618451) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(-1.9019351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4952963) q[0];
sx q[0];
rz(-1.4343867) q[0];
sx q[0];
rz(0.37149757) q[0];
rz(2.3067572) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(3.1239948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3018349) q[0];
sx q[0];
rz(-2.5389414) q[0];
sx q[0];
rz(-2.0728803) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8210635) q[2];
sx q[2];
rz(-2.2907718) q[2];
sx q[2];
rz(-1.4158451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0549893) q[1];
sx q[1];
rz(-1.2003683) q[1];
sx q[1];
rz(-0.61062529) q[1];
x q[2];
rz(2.5334355) q[3];
sx q[3];
rz(-1.3038692) q[3];
sx q[3];
rz(-1.1692695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1060433) q[2];
sx q[2];
rz(-2.273166) q[2];
sx q[2];
rz(-0.10031984) q[2];
rz(2.912168) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(-2.6917698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68596524) q[0];
sx q[0];
rz(-0.28814155) q[0];
sx q[0];
rz(-0.57383865) q[0];
rz(-1.0539184) q[1];
sx q[1];
rz(-1.046448) q[1];
sx q[1];
rz(0.67646772) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94714245) q[0];
sx q[0];
rz(-1.8136171) q[0];
sx q[0];
rz(0.28698289) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0208548) q[2];
sx q[2];
rz(-0.46560198) q[2];
sx q[2];
rz(0.6731205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4852462) q[1];
sx q[1];
rz(-2.606483) q[1];
sx q[1];
rz(-2.9675304) q[1];
rz(-pi) q[2];
rz(2.9432137) q[3];
sx q[3];
rz(-1.0470445) q[3];
sx q[3];
rz(-2.3959121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68710589) q[2];
sx q[2];
rz(-0.7578308) q[2];
sx q[2];
rz(-1.6472316) q[2];
rz(-0.86862653) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(-2.6715265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070521991) q[0];
sx q[0];
rz(-1.949953) q[0];
sx q[0];
rz(-1.9198887) q[0];
rz(1.3337878) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(-3.0467544) q[2];
sx q[2];
rz(-0.76764501) q[2];
sx q[2];
rz(-2.1697247) q[2];
rz(-1.3258237) q[3];
sx q[3];
rz(-0.180937) q[3];
sx q[3];
rz(-2.141249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
