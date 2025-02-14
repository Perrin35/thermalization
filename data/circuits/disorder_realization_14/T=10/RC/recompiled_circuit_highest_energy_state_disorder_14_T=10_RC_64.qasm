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
rz(-1.4426761) q[0];
sx q[0];
rz(0.28191167) q[0];
rz(0.52892041) q[1];
sx q[1];
rz(-1.4922053) q[1];
sx q[1];
rz(1.5780916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96254565) q[0];
sx q[0];
rz(-1.3907278) q[0];
sx q[0];
rz(1.942304) q[0];
x q[1];
rz(0.27484244) q[2];
sx q[2];
rz(-1.5759528) q[2];
sx q[2];
rz(-2.8066563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7719745) q[1];
sx q[1];
rz(-1.4450412) q[1];
sx q[1];
rz(0.63196504) q[1];
rz(-pi) q[2];
rz(0.6880349) q[3];
sx q[3];
rz(-1.8092938) q[3];
sx q[3];
rz(1.1900589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6866744) q[2];
sx q[2];
rz(-2.8439971) q[2];
sx q[2];
rz(0.61304027) q[2];
rz(-2.6679299) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(-1.7090428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365725) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(-0.75102425) q[0];
rz(2.660102) q[1];
sx q[1];
rz(-2.057169) q[1];
sx q[1];
rz(0.96985936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.65362) q[0];
sx q[0];
rz(-0.75838415) q[0];
sx q[0];
rz(-2.5874596) q[0];
rz(-pi) q[1];
rz(-1.6115236) q[2];
sx q[2];
rz(-1.3551522) q[2];
sx q[2];
rz(-3.0037896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9958933) q[1];
sx q[1];
rz(-2.4384192) q[1];
sx q[1];
rz(1.165676) q[1];
x q[2];
rz(2.5447846) q[3];
sx q[3];
rz(-1.9442026) q[3];
sx q[3];
rz(-1.4781836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2227309) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(0.53885031) q[2];
rz(-0.017596267) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(2.2331451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(3.0803296) q[0];
rz(-1.6368658) q[1];
sx q[1];
rz(-2.442339) q[1];
sx q[1];
rz(1.9452852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4825148) q[0];
sx q[0];
rz(-1.8097623) q[0];
sx q[0];
rz(0.84521291) q[0];
x q[1];
rz(-2.8792686) q[2];
sx q[2];
rz(-1.6379426) q[2];
sx q[2];
rz(3.0284428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34431008) q[1];
sx q[1];
rz(-2.4786665) q[1];
sx q[1];
rz(-1.0271038) q[1];
rz(-1.9798093) q[3];
sx q[3];
rz(-0.73506415) q[3];
sx q[3];
rz(-2.4215557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68941826) q[2];
sx q[2];
rz(-1.8345366) q[2];
sx q[2];
rz(2.7090731) q[2];
rz(1.0906667) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(-2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5525621) q[0];
sx q[0];
rz(-2.2046389) q[0];
sx q[0];
rz(2.9402148) q[0];
rz(-2.6248113) q[1];
sx q[1];
rz(-2.7613381) q[1];
sx q[1];
rz(-0.94863272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059221642) q[0];
sx q[0];
rz(-1.9192524) q[0];
sx q[0];
rz(2.2414464) q[0];
x q[1];
rz(-0.29251137) q[2];
sx q[2];
rz(-0.88380948) q[2];
sx q[2];
rz(-0.77921898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81056726) q[1];
sx q[1];
rz(-1.0665035) q[1];
sx q[1];
rz(1.585998) q[1];
rz(-2.620655) q[3];
sx q[3];
rz(-1.4864941) q[3];
sx q[3];
rz(-1.4985633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8676694) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(2.5841827) q[2];
rz(3.0607767) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7252561) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-1.0303372) q[0];
rz(0.15580767) q[1];
sx q[1];
rz(-1.067433) q[1];
sx q[1];
rz(-2.9373998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4725114) q[0];
sx q[0];
rz(-1.3288427) q[0];
sx q[0];
rz(-1.4999963) q[0];
x q[1];
rz(1.7014808) q[2];
sx q[2];
rz(-1.7820662) q[2];
sx q[2];
rz(1.4071161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0818644) q[1];
sx q[1];
rz(-1.3646185) q[1];
sx q[1];
rz(0.49525201) q[1];
rz(-pi) q[2];
rz(-0.77301003) q[3];
sx q[3];
rz(-1.7436308) q[3];
sx q[3];
rz(-2.5484249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.92945176) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(2.7280651) q[2];
rz(-2.2996969) q[3];
sx q[3];
rz(-1.3827518) q[3];
sx q[3];
rz(-2.1690185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464762) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(-0.0055775642) q[0];
rz(-1.7976409) q[1];
sx q[1];
rz(-2.9426212) q[1];
sx q[1];
rz(2.4406348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537121) q[0];
sx q[0];
rz(-1.0544262) q[0];
sx q[0];
rz(-0.51306458) q[0];
x q[1];
rz(0.43918886) q[2];
sx q[2];
rz(-2.1648295) q[2];
sx q[2];
rz(2.8043945) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29418221) q[1];
sx q[1];
rz(-2.8632616) q[1];
sx q[1];
rz(2.9918482) q[1];
rz(0.57311975) q[3];
sx q[3];
rz(-1.1240715) q[3];
sx q[3];
rz(-0.61441159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8274716) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(1.9742879) q[2];
rz(-2.2589034) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(0.32514969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90671396) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(-3.0562905) q[0];
rz(-0.75434297) q[1];
sx q[1];
rz(-0.5032379) q[1];
sx q[1];
rz(2.1139961) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32815702) q[0];
sx q[0];
rz(-1.3582843) q[0];
sx q[0];
rz(-1.5603609) q[0];
rz(-0.59664388) q[2];
sx q[2];
rz(-1.5663356) q[2];
sx q[2];
rz(-0.73847929) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3049406) q[1];
sx q[1];
rz(-1.8936367) q[1];
sx q[1];
rz(0.12860997) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6646752) q[3];
sx q[3];
rz(-2.0774042) q[3];
sx q[3];
rz(-1.5445386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9895642) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(0.41934553) q[2];
rz(0.63465214) q[3];
sx q[3];
rz(-0.80497634) q[3];
sx q[3];
rz(2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80970508) q[0];
sx q[0];
rz(-0.87600791) q[0];
sx q[0];
rz(-0.69212717) q[0];
rz(-2.5634815) q[1];
sx q[1];
rz(-0.41936857) q[1];
sx q[1];
rz(-1.7587761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.064474) q[0];
sx q[0];
rz(-3.0713284) q[0];
sx q[0];
rz(-1.2901359) q[0];
rz(-0.77246364) q[2];
sx q[2];
rz(-2.303745) q[2];
sx q[2];
rz(-2.7431184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67685157) q[1];
sx q[1];
rz(-0.60982043) q[1];
sx q[1];
rz(-1.1552385) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6666497) q[3];
sx q[3];
rz(-0.90427665) q[3];
sx q[3];
rz(-0.24196821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3708923) q[2];
sx q[2];
rz(-0.53052491) q[2];
sx q[2];
rz(-1.716506) q[2];
rz(-2.6254081) q[3];
sx q[3];
rz(-2.1631212) q[3];
sx q[3];
rz(-1.2396575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6462964) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(0.37149757) q[0];
rz(2.3067572) q[1];
sx q[1];
rz(-0.8131665) q[1];
sx q[1];
rz(0.017597839) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71401063) q[0];
sx q[0];
rz(-1.0508063) q[0];
sx q[0];
rz(0.31975694) q[0];
rz(-0.73569466) q[2];
sx q[2];
rz(-1.7580877) q[2];
sx q[2];
rz(-3.1295619) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0549893) q[1];
sx q[1];
rz(-1.2003683) q[1];
sx q[1];
rz(0.61062529) q[1];
rz(-0.60815717) q[3];
sx q[3];
rz(-1.8377234) q[3];
sx q[3];
rz(-1.9723231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1060433) q[2];
sx q[2];
rz(-0.86842662) q[2];
sx q[2];
rz(-0.10031984) q[2];
rz(0.22942461) q[3];
sx q[3];
rz(-1.0474297) q[3];
sx q[3];
rz(0.44982287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4556274) q[0];
sx q[0];
rz(-0.28814155) q[0];
sx q[0];
rz(-0.57383865) q[0];
rz(-2.0876743) q[1];
sx q[1];
rz(-1.046448) q[1];
sx q[1];
rz(2.4651249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1944502) q[0];
sx q[0];
rz(-1.8136171) q[0];
sx q[0];
rz(-0.28698289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25679882) q[2];
sx q[2];
rz(-1.9635824) q[2];
sx q[2];
rz(-1.8671672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6868848) q[1];
sx q[1];
rz(-2.096972) q[1];
sx q[1];
rz(-1.4684907) q[1];
x q[2];
rz(2.1031688) q[3];
sx q[3];
rz(-1.7422973) q[3];
sx q[3];
rz(-2.4166783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4544868) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(-1.494361) q[2];
rz(0.86862653) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(-0.47006616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-3.0710707) q[0];
sx q[0];
rz(-1.949953) q[0];
sx q[0];
rz(-1.9198887) q[0];
rz(1.3337878) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(-1.6619353) q[2];
sx q[2];
rz(-0.80747928) q[2];
sx q[2];
rz(0.840431) q[2];
rz(-1.815769) q[3];
sx q[3];
rz(-2.9606557) q[3];
sx q[3];
rz(1.0003436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
