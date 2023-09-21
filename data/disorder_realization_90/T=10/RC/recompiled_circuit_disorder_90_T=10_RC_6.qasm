OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4092769) q[0];
sx q[0];
rz(-1.9415932) q[0];
sx q[0];
rz(2.2229574) q[0];
rz(-pi) q[1];
rz(-1.4423223) q[2];
sx q[2];
rz(-0.80768425) q[2];
sx q[2];
rz(-1.8345923) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7510208) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(-1.8608171) q[1];
rz(-pi) q[2];
rz(-1.617241) q[3];
sx q[3];
rz(-1.9825476) q[3];
sx q[3];
rz(2.9573033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13880754) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(2.7677317) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.194838) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9276792) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(2.7620865) q[0];
x q[1];
rz(1.3090918) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(1.2145834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5217168) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(0.39342777) q[1];
x q[2];
rz(-0.11081103) q[3];
sx q[3];
rz(-0.46289819) q[3];
sx q[3];
rz(-3.0447931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(0.17671281) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(1.8331029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.20298) q[0];
sx q[0];
rz(-0.70883195) q[0];
sx q[0];
rz(2.3908486) q[0];
rz(-pi) q[1];
rz(0.53738014) q[2];
sx q[2];
rz(-1.3668622) q[2];
sx q[2];
rz(-0.98762074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32654999) q[1];
sx q[1];
rz(-0.76154852) q[1];
sx q[1];
rz(1.9425068) q[1];
rz(-1.4806467) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999009) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(-2.2256057) q[0];
rz(-2.6782716) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(2.0844918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(1.5989499) q[0];
rz(2.1787203) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(2.059666) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4790198) q[1];
sx q[1];
rz(-0.23049196) q[1];
sx q[1];
rz(1.5107933) q[1];
x q[2];
rz(1.3435059) q[3];
sx q[3];
rz(-2.3152581) q[3];
sx q[3];
rz(-2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(-0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3559568) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(2.6729029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5993354) q[0];
sx q[0];
rz(-2.242803) q[0];
sx q[0];
rz(-2.6319648) q[0];
x q[1];
rz(-0.30829633) q[2];
sx q[2];
rz(-1.7436276) q[2];
sx q[2];
rz(-2.1413213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9710755) q[1];
sx q[1];
rz(-1.3689539) q[1];
sx q[1];
rz(1.438373) q[1];
rz(-pi) q[2];
rz(-1.6563985) q[3];
sx q[3];
rz(-0.75827956) q[3];
sx q[3];
rz(0.12320923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(0.86432499) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(-0.29278452) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(0.19700225) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-2.8053455) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42423781) q[0];
sx q[0];
rz(-2.3399118) q[0];
sx q[0];
rz(-2.5108811) q[0];
x q[1];
rz(2.5480812) q[2];
sx q[2];
rz(-2.0715908) q[2];
sx q[2];
rz(1.8261432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.044144883) q[1];
sx q[1];
rz(-2.1598585) q[1];
sx q[1];
rz(-1.5138022) q[1];
x q[2];
rz(-1.7429966) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3939312) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(-1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31539311) q[0];
sx q[0];
rz(-2.0059735) q[0];
sx q[0];
rz(2.7639593) q[0];
rz(-1.7933493) q[2];
sx q[2];
rz(-1.5587274) q[2];
sx q[2];
rz(-2.0633069) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0371857) q[1];
sx q[1];
rz(-1.4763325) q[1];
sx q[1];
rz(0.027992804) q[1];
rz(-3.0435211) q[3];
sx q[3];
rz(-1.5925928) q[3];
sx q[3];
rz(-1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(-2.705412) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(-0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748579) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(-1.0048237) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(2.337713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893387) q[0];
sx q[0];
rz(-2.5544871) q[0];
sx q[0];
rz(0.071285204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5875823) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(0.49288921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1357321) q[1];
sx q[1];
rz(-1.5809959) q[1];
sx q[1];
rz(-1.482153) q[1];
rz(-1.6611093) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(-0.46935287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-2.1311029) q[2];
rz(2.5381952) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(2.7340775) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-0.75072748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6347046) q[0];
sx q[0];
rz(-1.8543715) q[0];
sx q[0];
rz(-2.8371235) q[0];
x q[1];
rz(2.8620371) q[2];
sx q[2];
rz(-1.4923555) q[2];
sx q[2];
rz(-2.560175) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0093065) q[1];
sx q[1];
rz(-2.7704151) q[1];
sx q[1];
rz(2.1585141) q[1];
rz(-1.5446072) q[3];
sx q[3];
rz(-1.2695754) q[3];
sx q[3];
rz(-1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0687381) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.9753974) q[2];
rz(1.3646305) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(-1.3903842) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.6814544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4936192) q[0];
sx q[0];
rz(-2.4612392) q[0];
sx q[0];
rz(-0.86662678) q[0];
rz(-pi) q[1];
rz(-1.6129458) q[2];
sx q[2];
rz(-0.50439207) q[2];
sx q[2];
rz(-1.1500037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32669386) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(-1.7705998) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4246033) q[3];
sx q[3];
rz(-2.2319712) q[3];
sx q[3];
rz(2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.4617408) q[2];
rz(2.0215624) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(-0.35342758) q[2];
sx q[2];
rz(-0.8197439) q[2];
sx q[2];
rz(2.843429) q[2];
rz(0.47386668) q[3];
sx q[3];
rz(-2.4110473) q[3];
sx q[3];
rz(-1.9163781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];