OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.83752051) q[0];
sx q[0];
rz(-0.66523319) q[0];
sx q[0];
rz(1.2256149) q[0];
rz(-0.6839112) q[1];
sx q[1];
rz(-2.749935) q[1];
sx q[1];
rz(2.439523) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4944899) q[0];
sx q[0];
rz(-0.70870525) q[0];
sx q[0];
rz(-2.9604218) q[0];
rz(-pi) q[1];
x q[1];
rz(1.443154) q[2];
sx q[2];
rz(-2.9753471) q[2];
sx q[2];
rz(-3.018741) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7154123) q[1];
sx q[1];
rz(-2.3857834) q[1];
sx q[1];
rz(0.73985696) q[1];
rz(-pi) q[2];
rz(2.2625173) q[3];
sx q[3];
rz(-1.723515) q[3];
sx q[3];
rz(2.9756551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2012653) q[2];
sx q[2];
rz(-1.5387115) q[2];
sx q[2];
rz(-2.8931457) q[2];
rz(-2.0072319) q[3];
sx q[3];
rz(-0.13392197) q[3];
sx q[3];
rz(2.0094481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092875384) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(-1.5485113) q[0];
rz(-2.6379207) q[1];
sx q[1];
rz(-1.2232989) q[1];
sx q[1];
rz(2.7647387) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232613) q[0];
sx q[0];
rz(-2.0428162) q[0];
sx q[0];
rz(-1.5492155) q[0];
x q[1];
rz(-0.81266788) q[2];
sx q[2];
rz(-2.118131) q[2];
sx q[2];
rz(0.081204435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7347757) q[1];
sx q[1];
rz(-1.8207014) q[1];
sx q[1];
rz(-2.6491449) q[1];
rz(-0.21864076) q[3];
sx q[3];
rz(-0.34014103) q[3];
sx q[3];
rz(-2.4778333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5571931) q[2];
sx q[2];
rz(-0.69584766) q[2];
sx q[2];
rz(0.86563555) q[2];
rz(-1.4731167) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(0.98751155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.57403785) q[0];
sx q[0];
rz(-1.2873298) q[0];
sx q[0];
rz(-2.0903184) q[0];
rz(1.4569262) q[1];
sx q[1];
rz(-1.3070062) q[1];
sx q[1];
rz(1.5164703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74836377) q[0];
sx q[0];
rz(-0.81384515) q[0];
sx q[0];
rz(2.7500344) q[0];
rz(-pi) q[1];
rz(0.90041884) q[2];
sx q[2];
rz(-1.0068276) q[2];
sx q[2];
rz(-1.054686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.17371) q[1];
sx q[1];
rz(-0.62025242) q[1];
sx q[1];
rz(2.2642676) q[1];
x q[2];
rz(0.72318913) q[3];
sx q[3];
rz(-0.08130493) q[3];
sx q[3];
rz(1.7618084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39525825) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(-2.9193817) q[2];
rz(0.5018417) q[3];
sx q[3];
rz(-1.4072489) q[3];
sx q[3];
rz(-0.80880222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.2489081) q[0];
sx q[0];
rz(-2.6539256) q[0];
sx q[0];
rz(-1.9768313) q[0];
rz(-2.248863) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(-2.8597615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696071) q[0];
sx q[0];
rz(-1.0907249) q[0];
sx q[0];
rz(-1.7064554) q[0];
rz(-pi) q[1];
rz(0.99165062) q[2];
sx q[2];
rz(-2.1985612) q[2];
sx q[2];
rz(2.6038458) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52543312) q[1];
sx q[1];
rz(-1.5010009) q[1];
sx q[1];
rz(1.7341803) q[1];
rz(0.025679703) q[3];
sx q[3];
rz(-2.0168983) q[3];
sx q[3];
rz(3.1334973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.07936) q[2];
sx q[2];
rz(-2.5793109) q[2];
sx q[2];
rz(-2.6083561) q[2];
rz(2.0594635) q[3];
sx q[3];
rz(-2.5366668) q[3];
sx q[3];
rz(2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.767652) q[0];
sx q[0];
rz(-2.7426608) q[0];
sx q[0];
rz(1.3237413) q[0];
rz(-0.81870493) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(2.0735819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49379738) q[0];
sx q[0];
rz(-1.6500705) q[0];
sx q[0];
rz(-0.1597516) q[0];
rz(2.2862412) q[2];
sx q[2];
rz(-0.87424874) q[2];
sx q[2];
rz(2.0717422) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6883478) q[1];
sx q[1];
rz(-1.4679586) q[1];
sx q[1];
rz(2.4842841) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16509861) q[3];
sx q[3];
rz(-2.0378651) q[3];
sx q[3];
rz(-1.9593173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1800804) q[2];
sx q[2];
rz(-1.1137806) q[2];
sx q[2];
rz(-2.2799344) q[2];
rz(1.0263475) q[3];
sx q[3];
rz(-1.15871) q[3];
sx q[3];
rz(2.0238743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79477972) q[0];
sx q[0];
rz(-2.3220799) q[0];
sx q[0];
rz(2.2487707) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-2.4632958) q[1];
sx q[1];
rz(3.0111664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.923188) q[0];
sx q[0];
rz(-1.7291862) q[0];
sx q[0];
rz(0.10938258) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4685235) q[2];
sx q[2];
rz(-2.4166738) q[2];
sx q[2];
rz(2.6672305) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32864704) q[1];
sx q[1];
rz(-2.9025902) q[1];
sx q[1];
rz(1.6944597) q[1];
rz(-pi) q[2];
rz(-2.942286) q[3];
sx q[3];
rz(-2.350507) q[3];
sx q[3];
rz(-1.6419322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22455198) q[2];
sx q[2];
rz(-1.2265393) q[2];
sx q[2];
rz(-0.87635931) q[2];
rz(-1.8996436) q[3];
sx q[3];
rz(-2.2010937) q[3];
sx q[3];
rz(2.131264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7414311) q[0];
sx q[0];
rz(-0.09859666) q[0];
sx q[0];
rz(-2.7778991) q[0];
rz(0.74186507) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(-2.738764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14430732) q[0];
sx q[0];
rz(-1.8107521) q[0];
sx q[0];
rz(1.994094) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3260457) q[2];
sx q[2];
rz(-2.6167653) q[2];
sx q[2];
rz(1.2957525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5093358) q[1];
sx q[1];
rz(-1.6922573) q[1];
sx q[1];
rz(-3.0640825) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8419344) q[3];
sx q[3];
rz(-2.1880272) q[3];
sx q[3];
rz(2.8927876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3591298) q[2];
sx q[2];
rz(-0.368258) q[2];
sx q[2];
rz(-1.5943607) q[2];
rz(0.83526978) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(-0.48867759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35992026) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(-0.51505995) q[0];
rz(1.7022279) q[1];
sx q[1];
rz(-1.8481588) q[1];
sx q[1];
rz(0.39915592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24267749) q[0];
sx q[0];
rz(-1.5726046) q[0];
sx q[0];
rz(-0.13092069) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24243124) q[2];
sx q[2];
rz(-2.2848086) q[2];
sx q[2];
rz(-0.27513181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6281575) q[1];
sx q[1];
rz(-1.6353459) q[1];
sx q[1];
rz(1.8316395) q[1];
x q[2];
rz(-1.5695482) q[3];
sx q[3];
rz(-1.6364508) q[3];
sx q[3];
rz(-1.1923238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1601552) q[2];
sx q[2];
rz(-1.258054) q[2];
sx q[2];
rz(2.0951648) q[2];
rz(2.4753172) q[3];
sx q[3];
rz(-2.8847238) q[3];
sx q[3];
rz(2.090914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1995354) q[0];
sx q[0];
rz(-3.0362447) q[0];
sx q[0];
rz(-2.5355205) q[0];
rz(0.82707682) q[1];
sx q[1];
rz(-1.3304293) q[1];
sx q[1];
rz(-1.3333295) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925643) q[0];
sx q[0];
rz(-1.6743393) q[0];
sx q[0];
rz(-2.7825481) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6624062) q[2];
sx q[2];
rz(-2.3385417) q[2];
sx q[2];
rz(2.242964) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.044802) q[1];
sx q[1];
rz(-1.4540744) q[1];
sx q[1];
rz(-2.2404284) q[1];
rz(-pi) q[2];
rz(-2.6831818) q[3];
sx q[3];
rz(-1.6473149) q[3];
sx q[3];
rz(-1.0118226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1313608) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(-0.38491797) q[2];
rz(-0.63079232) q[3];
sx q[3];
rz(-1.8600978) q[3];
sx q[3];
rz(1.7990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.887562) q[0];
sx q[0];
rz(-1.3986724) q[0];
sx q[0];
rz(1.4400462) q[0];
rz(-0.74132672) q[1];
sx q[1];
rz(-1.4011551) q[1];
sx q[1];
rz(-0.25873605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3869303) q[0];
sx q[0];
rz(-1.5748236) q[0];
sx q[0];
rz(-1.8167102) q[0];
rz(-pi) q[1];
rz(2.0889241) q[2];
sx q[2];
rz(-1.8473704) q[2];
sx q[2];
rz(-0.58108789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9548134) q[1];
sx q[1];
rz(-2.3329321) q[1];
sx q[1];
rz(-1.0891799) q[1];
rz(-pi) q[2];
rz(0.40269884) q[3];
sx q[3];
rz(-2.3220563) q[3];
sx q[3];
rz(1.0454901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72186333) q[2];
sx q[2];
rz(-1.1121007) q[2];
sx q[2];
rz(-1.8224243) q[2];
rz(2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(0.54660249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7405613) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(-1.7026547) q[1];
sx q[1];
rz(-0.66822744) q[1];
sx q[1];
rz(-2.6034036) q[1];
rz(-1.5580675) q[2];
sx q[2];
rz(-0.71011484) q[2];
sx q[2];
rz(-1.9488123) q[2];
rz(-2.1323754) q[3];
sx q[3];
rz(-1.781096) q[3];
sx q[3];
rz(-2.5566035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
