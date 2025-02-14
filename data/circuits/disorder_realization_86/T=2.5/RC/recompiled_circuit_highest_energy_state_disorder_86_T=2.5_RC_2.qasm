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
rz(0.64492172) q[0];
sx q[0];
rz(-0.46907297) q[0];
sx q[0];
rz(-0.99035779) q[0];
rz(-1.1733836) q[1];
sx q[1];
rz(-0.85964179) q[1];
sx q[1];
rz(0.8492066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17371236) q[0];
sx q[0];
rz(-3.1076508) q[0];
sx q[0];
rz(-1.4871661) q[0];
rz(-pi) q[1];
rz(2.2706896) q[2];
sx q[2];
rz(-1.3063626) q[2];
sx q[2];
rz(-1.3178133) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1942595) q[1];
sx q[1];
rz(-2.4521356) q[1];
sx q[1];
rz(1.1442776) q[1];
rz(-pi) q[2];
rz(-1.3732713) q[3];
sx q[3];
rz(-0.53086262) q[3];
sx q[3];
rz(1.498675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2363756) q[2];
sx q[2];
rz(-1.507501) q[2];
sx q[2];
rz(-2.6017453) q[2];
rz(1.5652462) q[3];
sx q[3];
rz(-2.7326475) q[3];
sx q[3];
rz(1.7319771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09963116) q[0];
sx q[0];
rz(-0.67622447) q[0];
sx q[0];
rz(-1.3270295) q[0];
rz(-2.3565893) q[1];
sx q[1];
rz(-1.4229341) q[1];
sx q[1];
rz(-0.50055093) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7206524) q[0];
sx q[0];
rz(-2.5897103) q[0];
sx q[0];
rz(-3.0938074) q[0];
rz(-2.9517216) q[2];
sx q[2];
rz(-1.9315757) q[2];
sx q[2];
rz(-1.9343513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5633531) q[1];
sx q[1];
rz(-2.1643442) q[1];
sx q[1];
rz(-2.7337573) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0359382) q[3];
sx q[3];
rz(-2.3243679) q[3];
sx q[3];
rz(1.5907839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0565722) q[2];
sx q[2];
rz(-2.4041921) q[2];
sx q[2];
rz(-2.6596587) q[2];
rz(2.7815172) q[3];
sx q[3];
rz(-1.1432546) q[3];
sx q[3];
rz(-2.9329407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8423186) q[0];
sx q[0];
rz(-2.0330918) q[0];
sx q[0];
rz(-2.6639248) q[0];
rz(-2.281669) q[1];
sx q[1];
rz(-1.0483024) q[1];
sx q[1];
rz(2.8544676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.053517) q[0];
sx q[0];
rz(-1.4447803) q[0];
sx q[0];
rz(-0.066670316) q[0];
rz(-pi) q[1];
rz(1.7874009) q[2];
sx q[2];
rz(-1.2239211) q[2];
sx q[2];
rz(1.6190478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17086731) q[1];
sx q[1];
rz(-1.5795465) q[1];
sx q[1];
rz(2.2544202) q[1];
x q[2];
rz(-1.7805598) q[3];
sx q[3];
rz(-0.38897369) q[3];
sx q[3];
rz(-0.28023187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3890248) q[2];
sx q[2];
rz(-1.7914881) q[2];
sx q[2];
rz(0.019850578) q[2];
rz(-0.94414532) q[3];
sx q[3];
rz(-0.70725924) q[3];
sx q[3];
rz(-0.39785644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95553628) q[0];
sx q[0];
rz(-0.63025403) q[0];
sx q[0];
rz(1.3007042) q[0];
rz(0.4862673) q[1];
sx q[1];
rz(-1.9673037) q[1];
sx q[1];
rz(-0.035331443) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5946884) q[0];
sx q[0];
rz(-0.42922089) q[0];
sx q[0];
rz(-1.5080601) q[0];
rz(-pi) q[1];
rz(-0.39725077) q[2];
sx q[2];
rz(-1.78671) q[2];
sx q[2];
rz(-1.4891032) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7481193) q[1];
sx q[1];
rz(-1.5474209) q[1];
sx q[1];
rz(-1.5651902) q[1];
rz(1.0159303) q[3];
sx q[3];
rz(-1.0175704) q[3];
sx q[3];
rz(-2.826626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6639634) q[2];
sx q[2];
rz(-2.5783381) q[2];
sx q[2];
rz(0.25838724) q[2];
rz(2.431331) q[3];
sx q[3];
rz(-2.5370772) q[3];
sx q[3];
rz(1.5822423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6579987) q[0];
sx q[0];
rz(-2.3778264) q[0];
sx q[0];
rz(0.64436954) q[0];
rz(-0.98980347) q[1];
sx q[1];
rz(-2.3376696) q[1];
sx q[1];
rz(-2.03233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98609867) q[0];
sx q[0];
rz(-1.7664096) q[0];
sx q[0];
rz(-1.0367111) q[0];
rz(-pi) q[1];
rz(1.057187) q[2];
sx q[2];
rz(-1.6053146) q[2];
sx q[2];
rz(2.2396954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68298816) q[1];
sx q[1];
rz(-0.79958497) q[1];
sx q[1];
rz(1.7728642) q[1];
rz(-pi) q[2];
rz(1.8568138) q[3];
sx q[3];
rz(-0.83854616) q[3];
sx q[3];
rz(-2.5570208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41205078) q[2];
sx q[2];
rz(-1.5770301) q[2];
sx q[2];
rz(0.61005074) q[2];
rz(-0.080502056) q[3];
sx q[3];
rz(-0.1463612) q[3];
sx q[3];
rz(3.1350873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7215111) q[0];
sx q[0];
rz(-0.93669909) q[0];
sx q[0];
rz(-1.6625846) q[0];
rz(1.9526019) q[1];
sx q[1];
rz(-0.34591302) q[1];
sx q[1];
rz(-1.4422653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1816671) q[0];
sx q[0];
rz(-0.64347351) q[0];
sx q[0];
rz(-1.305278) q[0];
rz(-pi) q[1];
rz(3.0411554) q[2];
sx q[2];
rz(-2.5578376) q[2];
sx q[2];
rz(0.69086087) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3644818) q[1];
sx q[1];
rz(-1.2419257) q[1];
sx q[1];
rz(-0.84192217) q[1];
rz(-pi) q[2];
rz(-2.6053794) q[3];
sx q[3];
rz(-1.271476) q[3];
sx q[3];
rz(-2.0423391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46001616) q[2];
sx q[2];
rz(-0.67395335) q[2];
sx q[2];
rz(0.67503929) q[2];
rz(-0.33440822) q[3];
sx q[3];
rz(-1.6655917) q[3];
sx q[3];
rz(2.7736751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48208958) q[0];
sx q[0];
rz(-0.98564321) q[0];
sx q[0];
rz(-3.0159045) q[0];
rz(2.9478759) q[1];
sx q[1];
rz(-0.88231641) q[1];
sx q[1];
rz(1.9901336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14979449) q[0];
sx q[0];
rz(-2.3201482) q[0];
sx q[0];
rz(-0.69355884) q[0];
rz(0.27982462) q[2];
sx q[2];
rz(-1.2132267) q[2];
sx q[2];
rz(1.2841061) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9509778) q[1];
sx q[1];
rz(-2.0101317) q[1];
sx q[1];
rz(-2.5530035) q[1];
rz(2.5732521) q[3];
sx q[3];
rz(-2.3821444) q[3];
sx q[3];
rz(-2.8041149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2842399) q[2];
sx q[2];
rz(-1.7809296) q[2];
sx q[2];
rz(2.8947158) q[2];
rz(-2.2829368) q[3];
sx q[3];
rz(-2.1571428) q[3];
sx q[3];
rz(0.27913678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7849279) q[0];
sx q[0];
rz(-2.6981638) q[0];
sx q[0];
rz(0.32522935) q[0];
rz(2.6204956) q[1];
sx q[1];
rz(-0.63322133) q[1];
sx q[1];
rz(-0.31347832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8756008) q[0];
sx q[0];
rz(-0.82339215) q[0];
sx q[0];
rz(-0.70988795) q[0];
rz(-1.3455079) q[2];
sx q[2];
rz(-1.5048001) q[2];
sx q[2];
rz(-0.57541945) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9814138) q[1];
sx q[1];
rz(-1.9367083) q[1];
sx q[1];
rz(-1.9311848) q[1];
rz(1.2282335) q[3];
sx q[3];
rz(-0.54751626) q[3];
sx q[3];
rz(-1.3207796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43361214) q[2];
sx q[2];
rz(-1.2182451) q[2];
sx q[2];
rz(1.8617967) q[2];
rz(2.4158939) q[3];
sx q[3];
rz(-1.1596707) q[3];
sx q[3];
rz(-0.80997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63876605) q[0];
sx q[0];
rz(-1.9483197) q[0];
sx q[0];
rz(2.9916812) q[0];
rz(1.8984849) q[1];
sx q[1];
rz(-1.5538235) q[1];
sx q[1];
rz(1.4814203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4325162) q[0];
sx q[0];
rz(-0.066532739) q[0];
sx q[0];
rz(2.3155022) q[0];
rz(-pi) q[1];
x q[1];
rz(0.051542087) q[2];
sx q[2];
rz(-1.8111501) q[2];
sx q[2];
rz(2.2759144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77566389) q[1];
sx q[1];
rz(-1.3559113) q[1];
sx q[1];
rz(-2.752423) q[1];
x q[2];
rz(1.3263014) q[3];
sx q[3];
rz(-1.6942548) q[3];
sx q[3];
rz(-2.3924912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82568613) q[2];
sx q[2];
rz(-1.76182) q[2];
sx q[2];
rz(-2.1659577) q[2];
rz(-1.1286831) q[3];
sx q[3];
rz(-2.5895139) q[3];
sx q[3];
rz(0.25477195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9752556) q[0];
sx q[0];
rz(-2.14125) q[0];
sx q[0];
rz(2.9794203) q[0];
rz(1.8099248) q[1];
sx q[1];
rz(-0.63260308) q[1];
sx q[1];
rz(3.0955637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9105658) q[0];
sx q[0];
rz(-1.5767158) q[0];
sx q[0];
rz(-1.5886515) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.677605) q[2];
sx q[2];
rz(-1.9205689) q[2];
sx q[2];
rz(-2.2217563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6190336) q[1];
sx q[1];
rz(-1.4231893) q[1];
sx q[1];
rz(1.1655432) q[1];
x q[2];
rz(-0.33933731) q[3];
sx q[3];
rz(-1.7320314) q[3];
sx q[3];
rz(1.1649023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33596805) q[2];
sx q[2];
rz(-0.38463548) q[2];
sx q[2];
rz(1.9920721) q[2];
rz(2.6286821) q[3];
sx q[3];
rz(-1.3866813) q[3];
sx q[3];
rz(0.30470595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.649986) q[0];
sx q[0];
rz(-2.2254324) q[0];
sx q[0];
rz(-1.9944763) q[0];
rz(3.0803549) q[1];
sx q[1];
rz(-1.1920659) q[1];
sx q[1];
rz(1.7658284) q[1];
rz(2.3394924) q[2];
sx q[2];
rz(-1.987793) q[2];
sx q[2];
rz(-1.1588617) q[2];
rz(-1.1652395) q[3];
sx q[3];
rz(-2.3074987) q[3];
sx q[3];
rz(-2.6354811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
