OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(4.1788221) q[0];
sx q[0];
rz(9.0691789) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(5.2390715) q[1];
sx q[1];
rz(10.704438) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0922001) q[0];
sx q[0];
rz(-2.85673) q[0];
sx q[0];
rz(-2.0869135) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6718164) q[2];
sx q[2];
rz(-2.6758286) q[2];
sx q[2];
rz(0.85675203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(2.5024662) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7198556) q[3];
sx q[3];
rz(-0.30656439) q[3];
sx q[3];
rz(-0.43678624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(2.546229) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-0.48746902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66747626) q[0];
sx q[0];
rz(-1.6352909) q[0];
sx q[0];
rz(1.4897896) q[0];
rz(-pi) q[1];
rz(2.3357453) q[2];
sx q[2];
rz(-2.5666551) q[2];
sx q[2];
rz(2.2573059) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6007874) q[1];
sx q[1];
rz(-1.4122737) q[1];
sx q[1];
rz(1.7608789) q[1];
rz(-1.3120193) q[3];
sx q[3];
rz(-0.32734713) q[3];
sx q[3];
rz(-2.5395288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-0.22274676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6762786) q[0];
sx q[0];
rz(-1.9440117) q[0];
sx q[0];
rz(2.7200384) q[0];
rz(1.1995302) q[2];
sx q[2];
rz(-1.3604593) q[2];
sx q[2];
rz(0.51031761) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45914868) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(-3.0241806) q[1];
rz(-pi) q[2];
rz(-2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(-2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(2.0866701) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.614914) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(-0.99292361) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22569457) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(0.29633488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11776609) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(-1.0002506) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6357189) q[3];
sx q[3];
rz(-0.40588356) q[3];
sx q[3];
rz(2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.2801923) q[2];
rz(-0.81930339) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(-0.76830307) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(2.6884902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74005175) q[0];
sx q[0];
rz(-1.4972367) q[0];
sx q[0];
rz(-3.1042276) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1066936) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(-0.81894433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0616152) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(-1.7916937) q[1];
rz(0.46053967) q[3];
sx q[3];
rz(-1.0214877) q[3];
sx q[3];
rz(-2.3159546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(-0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0244004) q[0];
sx q[0];
rz(-1.2478561) q[0];
sx q[0];
rz(-2.6431503) q[0];
rz(-0.061821826) q[2];
sx q[2];
rz(-2.0197868) q[2];
sx q[2];
rz(-2.2668554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7442419) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(-0.97126295) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7558277) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(-1.827084) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(1.1869173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531187) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(-0.32062809) q[0];
rz(1.600012) q[2];
sx q[2];
rz(-1.6299106) q[2];
sx q[2];
rz(-2.7616449) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.742332) q[1];
sx q[1];
rz(-2.2982344) q[1];
sx q[1];
rz(2.9620693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39077057) q[3];
sx q[3];
rz(-2.4739389) q[3];
sx q[3];
rz(2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(-1.5073744) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1925826) q[0];
sx q[0];
rz(-0.858925) q[0];
sx q[0];
rz(-2.683995) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.048642283) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(0.7359879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.028138782) q[1];
sx q[1];
rz(-0.33543643) q[1];
sx q[1];
rz(1.0949275) q[1];
x q[2];
rz(-0.52328531) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(0.39213359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-1.1351599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410256) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(-0.39649773) q[0];
rz(-pi) q[1];
rz(0.029149292) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(2.0707891) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4371722) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(1.0874332) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9477356) q[3];
sx q[3];
rz(-1.5514246) q[3];
sx q[3];
rz(1.3458348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-2.7897575) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-2.9454254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0862741) q[0];
sx q[0];
rz(-1.0090758) q[0];
sx q[0];
rz(-2.4713211) q[0];
rz(1.2162131) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(1.2186288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4539459) q[1];
sx q[1];
rz(-1.5861662) q[1];
sx q[1];
rz(0.83666283) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65532834) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(2.4650246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(-1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186196) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(0.20028533) q[2];
sx q[2];
rz(-2.9030048) q[2];
sx q[2];
rz(1.5422274) q[2];
rz(-3.0951981) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
