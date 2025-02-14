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
rz(2.9958148) q[0];
sx q[0];
rz(4.4895953) q[0];
sx q[0];
rz(9.0564981) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(0.84288016) q[1];
sx q[1];
rz(9.1755484) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0195784) q[0];
sx q[0];
rz(-1.459157) q[0];
sx q[0];
rz(0.48099244) q[0];
x q[1];
rz(2.8299242) q[2];
sx q[2];
rz(-2.083792) q[2];
sx q[2];
rz(-2.5530961) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.492384) q[1];
sx q[1];
rz(-1.3351775) q[1];
sx q[1];
rz(-0.71028965) q[1];
rz(-0.35855583) q[3];
sx q[3];
rz(-1.4643781) q[3];
sx q[3];
rz(-2.2582683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.92980802) q[2];
sx q[2];
rz(-0.70694184) q[2];
sx q[2];
rz(-0.24016538) q[2];
rz(0.54720488) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(-2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7996247) q[0];
sx q[0];
rz(-2.7662321) q[0];
sx q[0];
rz(-0.64390916) q[0];
rz(-2.0945235) q[1];
sx q[1];
rz(-1.1771076) q[1];
sx q[1];
rz(-2.8783096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851259) q[0];
sx q[0];
rz(-1.5432285) q[0];
sx q[0];
rz(0.03027244) q[0];
rz(-pi) q[1];
rz(-1.6356757) q[2];
sx q[2];
rz(-1.5047964) q[2];
sx q[2];
rz(-0.22021107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66052478) q[1];
sx q[1];
rz(-0.98036375) q[1];
sx q[1];
rz(-2.6284638) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41384048) q[3];
sx q[3];
rz(-0.97025052) q[3];
sx q[3];
rz(-2.5418856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(-0.72570428) q[2];
rz(0.30722412) q[3];
sx q[3];
rz(-1.8351464) q[3];
sx q[3];
rz(-2.0871625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82045186) q[0];
sx q[0];
rz(-1.4844002) q[0];
sx q[0];
rz(-2.3725574) q[0];
rz(-3.0342195) q[1];
sx q[1];
rz(-1.4391876) q[1];
sx q[1];
rz(-0.69951397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25722593) q[0];
sx q[0];
rz(-2.6623137) q[0];
sx q[0];
rz(0.81247654) q[0];
rz(0.6147521) q[2];
sx q[2];
rz(-1.5745133) q[2];
sx q[2];
rz(-2.2509509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86068639) q[1];
sx q[1];
rz(-1.6526892) q[1];
sx q[1];
rz(-2.0211092) q[1];
x q[2];
rz(2.3885257) q[3];
sx q[3];
rz(-2.1478466) q[3];
sx q[3];
rz(-2.6073539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37381441) q[2];
sx q[2];
rz(-2.5641597) q[2];
sx q[2];
rz(-2.2334297) q[2];
rz(2.5670037) q[3];
sx q[3];
rz(-1.8879954) q[3];
sx q[3];
rz(-0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.813756) q[0];
sx q[0];
rz(-1.8151374) q[0];
sx q[0];
rz(2.6326219) q[0];
rz(0.058188997) q[1];
sx q[1];
rz(-1.6845082) q[1];
sx q[1];
rz(1.0675272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539538) q[0];
sx q[0];
rz(-1.937307) q[0];
sx q[0];
rz(-2.3838759) q[0];
x q[1];
rz(1.9827438) q[2];
sx q[2];
rz(-1.6399523) q[2];
sx q[2];
rz(1.032304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6696251) q[1];
sx q[1];
rz(-1.4635007) q[1];
sx q[1];
rz(1.4221684) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.068693585) q[3];
sx q[3];
rz(-1.8634081) q[3];
sx q[3];
rz(-3.0076487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1099781) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(2.1515089) q[2];
rz(1.7456985) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(1.8949738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2331959) q[0];
sx q[0];
rz(-1.7447423) q[0];
sx q[0];
rz(-1.0594077) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.6672986) q[1];
sx q[1];
rz(1.4310736) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9741668) q[0];
sx q[0];
rz(-1.7210809) q[0];
sx q[0];
rz(1.7420064) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0756948) q[2];
sx q[2];
rz(-1.3422728) q[2];
sx q[2];
rz(-1.3040207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5268577) q[1];
sx q[1];
rz(-0.6893553) q[1];
sx q[1];
rz(-2.7562642) q[1];
rz(0.21417136) q[3];
sx q[3];
rz(-1.7167544) q[3];
sx q[3];
rz(-1.0792866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3502256) q[2];
sx q[2];
rz(-0.75608772) q[2];
sx q[2];
rz(-1.6737326) q[2];
rz(2.098846) q[3];
sx q[3];
rz(-1.0576893) q[3];
sx q[3];
rz(-2.3500672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4333711) q[0];
sx q[0];
rz(-1.8480166) q[0];
sx q[0];
rz(-2.9300387) q[0];
rz(-0.64282203) q[1];
sx q[1];
rz(-1.1245518) q[1];
sx q[1];
rz(2.0932253) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38163588) q[0];
sx q[0];
rz(-0.71361976) q[0];
sx q[0];
rz(-2.9455393) q[0];
x q[1];
rz(1.6238113) q[2];
sx q[2];
rz(-1.3919733) q[2];
sx q[2];
rz(0.12763466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.590946) q[1];
sx q[1];
rz(-2.5027486) q[1];
sx q[1];
rz(0.039012564) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95825742) q[3];
sx q[3];
rz(-1.6102092) q[3];
sx q[3];
rz(-2.266021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50531498) q[2];
sx q[2];
rz(-1.2039528) q[2];
sx q[2];
rz(-2.8875276) q[2];
rz(0.17357477) q[3];
sx q[3];
rz(-0.55145276) q[3];
sx q[3];
rz(1.180163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5467065) q[0];
sx q[0];
rz(-0.43989023) q[0];
sx q[0];
rz(0.79676262) q[0];
rz(0.88675371) q[1];
sx q[1];
rz(-1.795307) q[1];
sx q[1];
rz(2.5409882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85854218) q[0];
sx q[0];
rz(-1.5158537) q[0];
sx q[0];
rz(0.29015975) q[0];
x q[1];
rz(1.7037292) q[2];
sx q[2];
rz(-0.80604751) q[2];
sx q[2];
rz(1.9423167) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1315508) q[1];
sx q[1];
rz(-0.89183346) q[1];
sx q[1];
rz(0.20874397) q[1];
x q[2];
rz(1.0375627) q[3];
sx q[3];
rz(-1.8775465) q[3];
sx q[3];
rz(0.83707929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0299006) q[2];
sx q[2];
rz(-1.7039958) q[2];
sx q[2];
rz(2.0491484) q[2];
rz(1.7732636) q[3];
sx q[3];
rz(-0.60951257) q[3];
sx q[3];
rz(2.7294559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71548897) q[0];
sx q[0];
rz(-1.1142718) q[0];
sx q[0];
rz(-0.45147595) q[0];
rz(0.077797912) q[1];
sx q[1];
rz(-1.0323689) q[1];
sx q[1];
rz(-0.66158867) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1790893) q[0];
sx q[0];
rz(-2.668619) q[0];
sx q[0];
rz(-2.0684558) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0691452) q[2];
sx q[2];
rz(-1.5094286) q[2];
sx q[2];
rz(-2.5453486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3319131) q[1];
sx q[1];
rz(-0.47985489) q[1];
sx q[1];
rz(-1.2820246) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39644661) q[3];
sx q[3];
rz(-1.3260734) q[3];
sx q[3];
rz(0.94371599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0852069) q[2];
sx q[2];
rz(-1.2678009) q[2];
sx q[2];
rz(-2.3351604) q[2];
rz(-1.8993529) q[3];
sx q[3];
rz(-0.16568383) q[3];
sx q[3];
rz(-0.14036673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30271444) q[0];
sx q[0];
rz(-2.0063722) q[0];
sx q[0];
rz(-2.7022527) q[0];
rz(1.2002523) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(2.1913948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.882975) q[0];
sx q[0];
rz(-1.1990917) q[0];
sx q[0];
rz(0.88754334) q[0];
x q[1];
rz(0.15301159) q[2];
sx q[2];
rz(-1.8540314) q[2];
sx q[2];
rz(-2.6161414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2216144) q[1];
sx q[1];
rz(-1.6755381) q[1];
sx q[1];
rz(-1.9280738) q[1];
x q[2];
rz(2.1113273) q[3];
sx q[3];
rz(-1.9196307) q[3];
sx q[3];
rz(0.23603786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83303678) q[2];
sx q[2];
rz(-0.50806442) q[2];
sx q[2];
rz(-1.1888095) q[2];
rz(3.0360119) q[3];
sx q[3];
rz(-0.9413541) q[3];
sx q[3];
rz(-2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078876615) q[0];
sx q[0];
rz(-2.833241) q[0];
sx q[0];
rz(-2.4421332) q[0];
rz(1.4106916) q[1];
sx q[1];
rz(-2.3532093) q[1];
sx q[1];
rz(-0.20111045) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5749436) q[0];
sx q[0];
rz(-2.0127949) q[0];
sx q[0];
rz(-0.7676564) q[0];
x q[1];
rz(-2.2314448) q[2];
sx q[2];
rz(-2.4322699) q[2];
sx q[2];
rz(0.34914474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0929151) q[1];
sx q[1];
rz(-0.72775048) q[1];
sx q[1];
rz(1.9074832) q[1];
rz(-pi) q[2];
rz(-1.6612382) q[3];
sx q[3];
rz(-1.3916257) q[3];
sx q[3];
rz(-2.1018525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19645277) q[2];
sx q[2];
rz(-1.8818734) q[2];
sx q[2];
rz(-0.067616612) q[2];
rz(-2.0130646) q[3];
sx q[3];
rz(-1.0849378) q[3];
sx q[3];
rz(-1.6245925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6832798) q[0];
sx q[0];
rz(-2.0742317) q[0];
sx q[0];
rz(-1.7924894) q[0];
rz(1.6364527) q[1];
sx q[1];
rz(-0.22947336) q[1];
sx q[1];
rz(-0.089692399) q[1];
rz(-1.7791053) q[2];
sx q[2];
rz(-2.1635593) q[2];
sx q[2];
rz(0.61232739) q[2];
rz(2.454461) q[3];
sx q[3];
rz(-0.41768585) q[3];
sx q[3];
rz(-1.0286898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
