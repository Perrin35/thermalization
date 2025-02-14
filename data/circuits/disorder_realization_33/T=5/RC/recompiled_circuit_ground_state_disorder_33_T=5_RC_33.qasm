OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44242087) q[0];
sx q[0];
rz(-2.3306263) q[0];
sx q[0];
rz(2.6851658) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(2.1906617) q[1];
sx q[1];
rz(9.6074109) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5125677) q[0];
sx q[0];
rz(-2.1898666) q[0];
sx q[0];
rz(2.7946212) q[0];
x q[1];
rz(-0.86948189) q[2];
sx q[2];
rz(-1.1966146) q[2];
sx q[2];
rz(2.2187198) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6485893) q[1];
sx q[1];
rz(-1.7436639) q[1];
sx q[1];
rz(-1.6088617) q[1];
x q[2];
rz(2.6218518) q[3];
sx q[3];
rz(-1.9155353) q[3];
sx q[3];
rz(-0.23101364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7543588) q[2];
sx q[2];
rz(-1.0789472) q[2];
sx q[2];
rz(2.8299502) q[2];
rz(1.257487) q[3];
sx q[3];
rz(-2.8897372) q[3];
sx q[3];
rz(-0.18865147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99609128) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(1.2392932) q[0];
rz(-2.7481825) q[1];
sx q[1];
rz(-2.0937803) q[1];
sx q[1];
rz(-1.0308824) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1686954) q[0];
sx q[0];
rz(-0.6873695) q[0];
sx q[0];
rz(2.1363791) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8003045) q[2];
sx q[2];
rz(-1.82845) q[2];
sx q[2];
rz(-2.9837556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0075025) q[1];
sx q[1];
rz(-2.1975027) q[1];
sx q[1];
rz(0.43372633) q[1];
rz(-pi) q[2];
rz(0.64149789) q[3];
sx q[3];
rz(-2.201295) q[3];
sx q[3];
rz(1.3288095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3139412) q[2];
sx q[2];
rz(-2.2590019) q[2];
sx q[2];
rz(2.7871056) q[2];
rz(0.83141023) q[3];
sx q[3];
rz(-2.3737213) q[3];
sx q[3];
rz(-1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671234) q[0];
sx q[0];
rz(-2.9614083) q[0];
sx q[0];
rz(-2.6572976) q[0];
rz(0.88227415) q[1];
sx q[1];
rz(-2.2703998) q[1];
sx q[1];
rz(2.5994515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061848693) q[0];
sx q[0];
rz(-1.925191) q[0];
sx q[0];
rz(2.8097879) q[0];
rz(-pi) q[1];
rz(-0.32140478) q[2];
sx q[2];
rz(-1.633145) q[2];
sx q[2];
rz(-0.37918175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0391991) q[1];
sx q[1];
rz(-2.3619283) q[1];
sx q[1];
rz(-1.8236266) q[1];
rz(-pi) q[2];
x q[2];
rz(0.063850689) q[3];
sx q[3];
rz(-1.333916) q[3];
sx q[3];
rz(0.76343918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.813628) q[2];
sx q[2];
rz(-0.70983228) q[2];
sx q[2];
rz(2.1072809) q[2];
rz(-1.3616925) q[3];
sx q[3];
rz(-1.9618278) q[3];
sx q[3];
rz(-0.55083197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4807602) q[0];
sx q[0];
rz(-1.7518504) q[0];
sx q[0];
rz(2.0939636) q[0];
rz(-1.3230336) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(2.7563162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.229363) q[0];
sx q[0];
rz(-1.1426539) q[0];
sx q[0];
rz(1.0083126) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6462506) q[2];
sx q[2];
rz(-2.2860043) q[2];
sx q[2];
rz(-2.070141) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4874607) q[1];
sx q[1];
rz(-1.7791554) q[1];
sx q[1];
rz(-1.4980493) q[1];
rz(0.88726081) q[3];
sx q[3];
rz(-1.6561015) q[3];
sx q[3];
rz(2.4312756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72383991) q[2];
sx q[2];
rz(-0.94261348) q[2];
sx q[2];
rz(2.8670132) q[2];
rz(-0.56139055) q[3];
sx q[3];
rz(-2.0761108) q[3];
sx q[3];
rz(-1.5076465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0059589) q[0];
sx q[0];
rz(-2.5649286) q[0];
sx q[0];
rz(2.2744001) q[0];
rz(2.7658956) q[1];
sx q[1];
rz(-0.58940327) q[1];
sx q[1];
rz(1.9120749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.00261) q[0];
sx q[0];
rz(-2.0073237) q[0];
sx q[0];
rz(-1.3912203) q[0];
rz(0.30699344) q[2];
sx q[2];
rz(-2.0684264) q[2];
sx q[2];
rz(2.1461329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41709677) q[1];
sx q[1];
rz(-2.1388106) q[1];
sx q[1];
rz(-0.54996164) q[1];
x q[2];
rz(-0.8939871) q[3];
sx q[3];
rz(-0.5115307) q[3];
sx q[3];
rz(0.77340305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.025297252) q[2];
sx q[2];
rz(-1.9090434) q[2];
sx q[2];
rz(3.0787025) q[2];
rz(-0.40999117) q[3];
sx q[3];
rz(-2.4153109) q[3];
sx q[3];
rz(-2.9819152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9922239) q[0];
sx q[0];
rz(-0.069644444) q[0];
sx q[0];
rz(2.3058291) q[0];
rz(-1.1831076) q[1];
sx q[1];
rz(-1.2703905) q[1];
sx q[1];
rz(-2.3979208) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2398674) q[0];
sx q[0];
rz(-2.5081303) q[0];
sx q[0];
rz(2.6153436) q[0];
x q[1];
rz(-1.9299149) q[2];
sx q[2];
rz(-0.15378498) q[2];
sx q[2];
rz(-1.6395456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6059123) q[1];
sx q[1];
rz(-2.4125068) q[1];
sx q[1];
rz(-1.8465183) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5178237) q[3];
sx q[3];
rz(-1.3154239) q[3];
sx q[3];
rz(-2.6290174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30367294) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(-0.79353235) q[2];
rz(-0.41905904) q[3];
sx q[3];
rz(-1.4895118) q[3];
sx q[3];
rz(0.52217531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.94988743) q[0];
sx q[0];
rz(-2.1623623) q[0];
sx q[0];
rz(-1.1631843) q[0];
rz(2.361182) q[1];
sx q[1];
rz(-2.8051832) q[1];
sx q[1];
rz(1.6132678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96161973) q[0];
sx q[0];
rz(-0.69937569) q[0];
sx q[0];
rz(-1.6716624) q[0];
rz(-pi) q[1];
rz(-2.0397908) q[2];
sx q[2];
rz(-2.7534979) q[2];
sx q[2];
rz(2.3562252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53832053) q[1];
sx q[1];
rz(-1.5495991) q[1];
sx q[1];
rz(1.7689598) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9933543) q[3];
sx q[3];
rz(-1.3179639) q[3];
sx q[3];
rz(-1.2444843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10716001) q[2];
sx q[2];
rz(-1.9623423) q[2];
sx q[2];
rz(-3.072928) q[2];
rz(0.56985235) q[3];
sx q[3];
rz(-2.6611501) q[3];
sx q[3];
rz(0.3705875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079085199) q[0];
sx q[0];
rz(-0.67254368) q[0];
sx q[0];
rz(1.8261209) q[0];
rz(1.2830118) q[1];
sx q[1];
rz(-0.43671572) q[1];
sx q[1];
rz(-0.049093094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0822711) q[0];
sx q[0];
rz(-1.647864) q[0];
sx q[0];
rz(-1.4801816) q[0];
x q[1];
rz(-2.6025199) q[2];
sx q[2];
rz(-2.1207715) q[2];
sx q[2];
rz(-2.0520738) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.093542592) q[1];
sx q[1];
rz(-0.88359264) q[1];
sx q[1];
rz(2.6395413) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2146739) q[3];
sx q[3];
rz(-1.4689494) q[3];
sx q[3];
rz(2.1702914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7580238) q[2];
sx q[2];
rz(-0.87791666) q[2];
sx q[2];
rz(-0.54086584) q[2];
rz(-1.0848378) q[3];
sx q[3];
rz(-2.4386051) q[3];
sx q[3];
rz(1.7613523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.908602) q[0];
sx q[0];
rz(-2.1451696) q[0];
sx q[0];
rz(-3.0365699) q[0];
rz(-0.54166334) q[1];
sx q[1];
rz(-0.88625208) q[1];
sx q[1];
rz(2.7856316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.134577) q[0];
sx q[0];
rz(-1.4316214) q[0];
sx q[0];
rz(3.1037168) q[0];
x q[1];
rz(-0.84635205) q[2];
sx q[2];
rz(-2.0931819) q[2];
sx q[2];
rz(-1.236793) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4527725) q[1];
sx q[1];
rz(-2.4311307) q[1];
sx q[1];
rz(-0.97150357) q[1];
rz(-pi) q[2];
rz(3.0186986) q[3];
sx q[3];
rz(-1.8884522) q[3];
sx q[3];
rz(-1.4113219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2124704) q[2];
sx q[2];
rz(-1.4996303) q[2];
sx q[2];
rz(-1.8222202) q[2];
rz(1.8170554) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(0.96778473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20973715) q[0];
sx q[0];
rz(-3.1071438) q[0];
sx q[0];
rz(1.4631648) q[0];
rz(-1.6819008) q[1];
sx q[1];
rz(-1.9821143) q[1];
sx q[1];
rz(0.34585888) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5992085) q[0];
sx q[0];
rz(-0.092723474) q[0];
sx q[0];
rz(-1.0371764) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9458617) q[2];
sx q[2];
rz(-0.32265857) q[2];
sx q[2];
rz(-3.0374073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.085850216) q[1];
sx q[1];
rz(-2.3338103) q[1];
sx q[1];
rz(2.4959688) q[1];
rz(-pi) q[2];
rz(-1.1283952) q[3];
sx q[3];
rz(-0.69303362) q[3];
sx q[3];
rz(-2.5062163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86089245) q[2];
sx q[2];
rz(-1.2714551) q[2];
sx q[2];
rz(-0.26091179) q[2];
rz(1.4704618) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(1.7117333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83229257) q[0];
sx q[0];
rz(-2.0335048) q[0];
sx q[0];
rz(-0.19620398) q[0];
rz(0.55083864) q[1];
sx q[1];
rz(-1.6549587) q[1];
sx q[1];
rz(0.5400198) q[1];
rz(2.9020799) q[2];
sx q[2];
rz(-0.85776599) q[2];
sx q[2];
rz(0.89606482) q[2];
rz(0.30487074) q[3];
sx q[3];
rz(-0.75406995) q[3];
sx q[3];
rz(-2.0537805) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
