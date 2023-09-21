OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5776423) q[0];
sx q[0];
rz(-0.9194153) q[0];
sx q[0];
rz(1.7786068) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3876786) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-2.4726601) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2540993) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(-1.9894132) q[1];
rz(-0.9382117) q[3];
sx q[3];
rz(-1.0217561) q[3];
sx q[3];
rz(-0.051018056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(-2.1263188) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.9706443) q[0];
sx q[0];
rz(-1.0785372) q[0];
sx q[0];
rz(1.2600391) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(-0.84567436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2226919) q[0];
sx q[0];
rz(-2.0881623) q[0];
sx q[0];
rz(-1.067724) q[0];
x q[1];
rz(1.4405865) q[2];
sx q[2];
rz(-1.7034334) q[2];
sx q[2];
rz(-0.33078937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3413275) q[1];
sx q[1];
rz(-0.53293537) q[1];
sx q[1];
rz(-2.76782) q[1];
x q[2];
rz(2.1917079) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(1.8267531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95124328) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(0.7730661) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(1.0864331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13585424) q[0];
sx q[0];
rz(-1.5691301) q[0];
sx q[0];
rz(1.1093596) q[0];
rz(-0.92347446) q[2];
sx q[2];
rz(-1.4787276) q[2];
sx q[2];
rz(1.1084686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.31049) q[1];
sx q[1];
rz(-1.5655787) q[1];
sx q[1];
rz(-2.8527841) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9757189) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(-2.3111642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(1.4455459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915582) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(0.66288373) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6608597) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(-2.4243674) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85019894) q[1];
sx q[1];
rz(-1.351601) q[1];
sx q[1];
rz(1.6367957) q[1];
rz(2.8533832) q[3];
sx q[3];
rz(-1.0681579) q[3];
sx q[3];
rz(-2.1484745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(-2.1172822) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720649) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(-3.1307427) q[0];
x q[1];
rz(-2.1335667) q[2];
sx q[2];
rz(-0.78283435) q[2];
sx q[2];
rz(0.40751878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8728719) q[1];
sx q[1];
rz(-2.1495021) q[1];
sx q[1];
rz(-2.242356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5913127) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(-2.9328336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-2.234941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(3.0153826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.135658) q[0];
sx q[0];
rz(-1.6023484) q[0];
sx q[0];
rz(0.21261442) q[0];
x q[1];
rz(0.38527617) q[2];
sx q[2];
rz(-0.81309536) q[2];
sx q[2];
rz(-2.3713881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0969442) q[1];
sx q[1];
rz(-1.3828039) q[1];
sx q[1];
rz(-2.2915927) q[1];
x q[2];
rz(-0.19212171) q[3];
sx q[3];
rz(-1.2479094) q[3];
sx q[3];
rz(-2.9411112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(2.3664756) q[2];
rz(2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-2.0194139) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37724272) q[0];
sx q[0];
rz(-0.61315216) q[0];
sx q[0];
rz(0.22293798) q[0];
rz(-pi) q[1];
rz(1.6227116) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.3494929) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0911078) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(2.1662103) q[1];
x q[2];
rz(2.8758994) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(-0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7155834) q[0];
sx q[0];
rz(-1.4953488) q[0];
sx q[0];
rz(-1.5619318) q[0];
rz(-pi) q[1];
rz(-2.7996054) q[2];
sx q[2];
rz(-2.7139398) q[2];
sx q[2];
rz(-0.74117408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65054446) q[1];
sx q[1];
rz(-0.77495134) q[1];
sx q[1];
rz(0.52853711) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7730764) q[3];
sx q[3];
rz(-2.5296122) q[3];
sx q[3];
rz(-2.4806541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(2.0020206) q[2];
rz(-1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9386439) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(-2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7302007) q[0];
sx q[0];
rz(-3.054266) q[0];
sx q[0];
rz(-1.9090396) q[0];
x q[1];
rz(-1.7622856) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(-2.7188403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.081454885) q[1];
sx q[1];
rz(-0.36768915) q[1];
sx q[1];
rz(2.0245488) q[1];
rz(-pi) q[2];
rz(-0.41096656) q[3];
sx q[3];
rz(-2.4604359) q[3];
sx q[3];
rz(0.70511234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-2.7900556) q[2];
rz(2.0848138) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(0.25340733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8966658) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(0.11163296) q[0];
rz(-pi) q[1];
rz(-2.9741653) q[2];
sx q[2];
rz(-1.9055467) q[2];
sx q[2];
rz(-1.4565005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3897755) q[1];
sx q[1];
rz(-0.42418617) q[1];
sx q[1];
rz(0.61702375) q[1];
rz(1.1425584) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(0.67728562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59166756) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-0.7111711) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-2.860643) q[2];
sx q[2];
rz(-1.8200257) q[2];
sx q[2];
rz(-0.65489468) q[2];
rz(-0.9234225) q[3];
sx q[3];
rz(-2.2150726) q[3];
sx q[3];
rz(1.1005145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];