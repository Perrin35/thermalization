OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(-1.0153271) q[0];
sx q[0];
rz(0.46749687) q[0];
rz(2.6214018) q[1];
sx q[1];
rz(-1.7953035) q[1];
sx q[1];
rz(-2.4612114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8898534) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(-2.2962909) q[0];
rz(-pi) q[1];
rz(-1.275185) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(2.0762028) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89741325) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(-1.0531055) q[1];
rz(-pi) q[2];
rz(0.47547961) q[3];
sx q[3];
rz(-2.7694422) q[3];
sx q[3];
rz(-1.1572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9900069) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(2.4123689) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4366348) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(1.0013642) q[0];
rz(2.9691866) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(-0.52406812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66602089) q[0];
sx q[0];
rz(-2.4040939) q[0];
sx q[0];
rz(-2.0185508) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42627724) q[2];
sx q[2];
rz(-2.2239416) q[2];
sx q[2];
rz(-0.041235812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.14428917) q[1];
sx q[1];
rz(-1.3376437) q[1];
sx q[1];
rz(0.38645978) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7254132) q[3];
sx q[3];
rz(-2.1025476) q[3];
sx q[3];
rz(1.0242517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(1.8015507) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(-2.952081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7127842) q[0];
sx q[0];
rz(-1.3002987) q[0];
sx q[0];
rz(0.74032797) q[0];
rz(-1.6763716) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(-2.2764652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1358007) q[1];
sx q[1];
rz(-1.6264919) q[1];
sx q[1];
rz(-2.9018351) q[1];
rz(-pi) q[2];
rz(2.965851) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(-0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(1.3872046) q[2];
rz(2.710279) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65961924) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(-1.5959928) q[0];
rz(-2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(-1.0901573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6369612) q[0];
sx q[0];
rz(-1.9177027) q[0];
sx q[0];
rz(2.7142801) q[0];
rz(0.81972576) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(2.0476598) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5174597) q[1];
sx q[1];
rz(-2.1412666) q[1];
sx q[1];
rz(2.0681417) q[1];
rz(-pi) q[2];
rz(-2.7146043) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1426992) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(2.2284171) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49044931) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8613375) q[0];
sx q[0];
rz(-2.4147408) q[0];
sx q[0];
rz(1.0759541) q[0];
rz(1.9636376) q[2];
sx q[2];
rz(-0.71937865) q[2];
sx q[2];
rz(0.17578416) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40844803) q[1];
sx q[1];
rz(-1.1322347) q[1];
sx q[1];
rz(1.1327729) q[1];
x q[2];
rz(2.6036161) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(0.81975421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(-1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(-0.67289105) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2900989) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(1.7593988) q[0];
x q[1];
rz(0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(1.6753472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91025483) q[1];
sx q[1];
rz(-0.64097039) q[1];
sx q[1];
rz(2.7892968) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17500413) q[3];
sx q[3];
rz(-2.724218) q[3];
sx q[3];
rz(-1.1328896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.7395082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3081449) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(1.863198) q[0];
rz(-pi) q[1];
rz(0.76485302) q[2];
sx q[2];
rz(-1.7424889) q[2];
sx q[2];
rz(-0.21966759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40560383) q[1];
sx q[1];
rz(-2.7466008) q[1];
sx q[1];
rz(2.7337381) q[1];
rz(-pi) q[2];
rz(0.98251179) q[3];
sx q[3];
rz(-1.4559828) q[3];
sx q[3];
rz(-1.7867775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4222251) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(-0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.592214) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-0.54135281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4810825) q[0];
sx q[0];
rz(-2.5773002) q[0];
sx q[0];
rz(0.020945992) q[0];
rz(-pi) q[1];
rz(-0.45221046) q[2];
sx q[2];
rz(-1.2621677) q[2];
sx q[2];
rz(1.2436777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7975446) q[1];
sx q[1];
rz(-0.86595264) q[1];
sx q[1];
rz(-2.123358) q[1];
rz(-pi) q[2];
rz(2.7612711) q[3];
sx q[3];
rz(-1.5690648) q[3];
sx q[3];
rz(-2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.4452971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8953067) q[0];
sx q[0];
rz(-0.9406957) q[0];
sx q[0];
rz(2.2340328) q[0];
x q[1];
rz(-2.2575349) q[2];
sx q[2];
rz(-2.3381655) q[2];
sx q[2];
rz(-0.49056177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5053133) q[1];
sx q[1];
rz(-1.5863824) q[1];
sx q[1];
rz(-0.63025766) q[1];
rz(0.85083665) q[3];
sx q[3];
rz(-0.81504226) q[3];
sx q[3];
rz(-1.9581865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(2.5342069) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(-0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8059175) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(0.6138531) q[0];
rz(-1.0461668) q[1];
sx q[1];
rz(-2.8764953) q[1];
sx q[1];
rz(2.7493431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989202) q[0];
sx q[0];
rz(-0.93085104) q[0];
sx q[0];
rz(-1.0660893) q[0];
rz(-pi) q[1];
rz(2.1456111) q[2];
sx q[2];
rz(-2.1157584) q[2];
sx q[2];
rz(1.121322) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2892462) q[1];
sx q[1];
rz(-0.95818555) q[1];
sx q[1];
rz(-0.5457408) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20499968) q[3];
sx q[3];
rz(-1.5424171) q[3];
sx q[3];
rz(1.8207707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-2.0742119) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(-2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(0.36841064) q[2];
sx q[2];
rz(-0.80130063) q[2];
sx q[2];
rz(-0.038566312) q[2];
rz(2.5946887) q[3];
sx q[3];
rz(-2.1756267) q[3];
sx q[3];
rz(-0.970943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
