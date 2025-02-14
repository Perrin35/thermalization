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
rz(1.3801112) q[0];
sx q[0];
rz(-1.5910281) q[0];
sx q[0];
rz(-0.38183364) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5474553) q[0];
sx q[0];
rz(-2.1004803) q[0];
sx q[0];
rz(2.269777) q[0];
rz(-pi) q[1];
rz(0.37021356) q[2];
sx q[2];
rz(-1.5523124) q[2];
sx q[2];
rz(3.0188675) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3081835) q[1];
sx q[1];
rz(-1.7004968) q[1];
sx q[1];
rz(2.7144578) q[1];
x q[2];
rz(2.1006868) q[3];
sx q[3];
rz(-1.9090416) q[3];
sx q[3];
rz(1.6774981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.8483346) q[2];
sx q[2];
rz(-1.5748242) q[2];
rz(-1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67398706) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(1.9655193) q[0];
rz(-3.0740956) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(2.8228021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8628937) q[0];
sx q[0];
rz(-0.79813939) q[0];
sx q[0];
rz(-0.44095914) q[0];
x q[1];
rz(2.6574357) q[2];
sx q[2];
rz(-0.47704298) q[2];
sx q[2];
rz(-0.31443982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76004878) q[1];
sx q[1];
rz(-2.4439619) q[1];
sx q[1];
rz(0.95894869) q[1];
x q[2];
rz(3.0158668) q[3];
sx q[3];
rz(-1.8549524) q[3];
sx q[3];
rz(-1.3181669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0789644) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(-2.6641565) q[2];
rz(1.7834974) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(-3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(1.1392449) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.258305) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8157522) q[0];
sx q[0];
rz(-1.6597972) q[0];
sx q[0];
rz(2.1584216) q[0];
rz(-pi) q[1];
rz(-2.0482721) q[2];
sx q[2];
rz(-2.4585063) q[2];
sx q[2];
rz(-2.8841599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.660608) q[1];
sx q[1];
rz(-2.0913731) q[1];
sx q[1];
rz(1.03517) q[1];
x q[2];
rz(2.0800679) q[3];
sx q[3];
rz(-1.2767999) q[3];
sx q[3];
rz(-1.7727838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4262126) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(2.0646084) q[2];
rz(-1.8814603) q[3];
sx q[3];
rz(-2.1144805) q[3];
sx q[3];
rz(0.25585678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.857321) q[0];
sx q[0];
rz(-1.9673328) q[0];
sx q[0];
rz(-0.76048365) q[0];
rz(0.35176945) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(-2.9913091) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12211299) q[0];
sx q[0];
rz(-1.5727336) q[0];
sx q[0];
rz(-1.5674855) q[0];
x q[1];
rz(1.1614805) q[2];
sx q[2];
rz(-2.7285353) q[2];
sx q[2];
rz(1.5395791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6071668) q[1];
sx q[1];
rz(-2.1517506) q[1];
sx q[1];
rz(2.7423285) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3129381) q[3];
sx q[3];
rz(-2.0153196) q[3];
sx q[3];
rz(0.049098102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.030423) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(-0.58491582) q[2];
rz(3.0758744) q[3];
sx q[3];
rz(-1.1394371) q[3];
sx q[3];
rz(-1.0544624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86972648) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(-2.4140893) q[0];
rz(-1.8456521) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(-2.4268699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489577) q[0];
sx q[0];
rz(-1.4997109) q[0];
sx q[0];
rz(-0.22146225) q[0];
rz(-1.6863912) q[2];
sx q[2];
rz(-1.8697479) q[2];
sx q[2];
rz(1.7345127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59698518) q[1];
sx q[1];
rz(-1.3860354) q[1];
sx q[1];
rz(-1.9562592) q[1];
x q[2];
rz(-1.0017702) q[3];
sx q[3];
rz(-2.4845893) q[3];
sx q[3];
rz(3.0083619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1456445) q[2];
sx q[2];
rz(-1.5278634) q[2];
sx q[2];
rz(1.145251) q[2];
rz(-2.9389985) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(0.32874671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458493) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(-0.16786815) q[0];
rz(-0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(1.3289183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58121366) q[0];
sx q[0];
rz(-2.3920076) q[0];
sx q[0];
rz(-2.4102274) q[0];
rz(-pi) q[1];
rz(1.6199048) q[2];
sx q[2];
rz(-1.430776) q[2];
sx q[2];
rz(2.5852709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5649453) q[1];
sx q[1];
rz(-2.7912346) q[1];
sx q[1];
rz(0.74585657) q[1];
rz(1.331948) q[3];
sx q[3];
rz(-1.3535366) q[3];
sx q[3];
rz(2.1715733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7557175) q[2];
sx q[2];
rz(-1.8942602) q[2];
sx q[2];
rz(2.7937549) q[2];
rz(0.31473413) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(-1.1381963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1494074) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(2.1904679) q[0];
rz(0.29257193) q[1];
sx q[1];
rz(-0.89076275) q[1];
sx q[1];
rz(-0.94400418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2601449) q[0];
sx q[0];
rz(-1.8788188) q[0];
sx q[0];
rz(-0.016168895) q[0];
x q[1];
rz(2.0538141) q[2];
sx q[2];
rz(-2.000914) q[2];
sx q[2];
rz(-1.5120235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71783644) q[1];
sx q[1];
rz(-1.593381) q[1];
sx q[1];
rz(-1.3976239) q[1];
rz(-pi) q[2];
rz(1.0349501) q[3];
sx q[3];
rz(-2.6831045) q[3];
sx q[3];
rz(1.5036086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6501179) q[2];
sx q[2];
rz(-2.398141) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(-0.024638351) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(0.98954454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80855075) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(0.1380052) q[0];
rz(2.0637312) q[1];
sx q[1];
rz(-1.4638008) q[1];
sx q[1];
rz(0.93625751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662921) q[0];
sx q[0];
rz(-1.202824) q[0];
sx q[0];
rz(-1.247768) q[0];
rz(-pi) q[1];
rz(-0.54674863) q[2];
sx q[2];
rz(-0.90675747) q[2];
sx q[2];
rz(-0.67686096) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.028658) q[1];
sx q[1];
rz(-1.3880265) q[1];
sx q[1];
rz(-3.1266646) q[1];
x q[2];
rz(2.8209723) q[3];
sx q[3];
rz(-2.4302539) q[3];
sx q[3];
rz(2.5742755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1428947) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(-3.0190492) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(-0.43012777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1634624) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(-2.7769856) q[0];
rz(1.3847903) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(0.91420954) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96959463) q[0];
sx q[0];
rz(-1.9978317) q[0];
sx q[0];
rz(-0.14694218) q[0];
rz(2.4070074) q[2];
sx q[2];
rz(-1.5365639) q[2];
sx q[2];
rz(1.6357259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89489854) q[1];
sx q[1];
rz(-1.7082038) q[1];
sx q[1];
rz(0.64260428) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8960885) q[3];
sx q[3];
rz(-2.5872018) q[3];
sx q[3];
rz(1.7911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.30283516) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(-1.6287104) q[2];
rz(0.58879876) q[3];
sx q[3];
rz(-1.9353119) q[3];
sx q[3];
rz(0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47422472) q[0];
sx q[0];
rz(-2.5741757) q[0];
sx q[0];
rz(-0.2844511) q[0];
rz(-2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(-2.7213352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7922554) q[0];
sx q[0];
rz(-0.47264034) q[0];
sx q[0];
rz(-1.0853173) q[0];
x q[1];
rz(1.1070232) q[2];
sx q[2];
rz(-0.83597413) q[2];
sx q[2];
rz(-0.42586621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0033231) q[1];
sx q[1];
rz(-1.4821293) q[1];
sx q[1];
rz(-1.6900464) q[1];
rz(-pi) q[2];
rz(-2.5758366) q[3];
sx q[3];
rz(-1.8744191) q[3];
sx q[3];
rz(-0.88241258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9922716) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(-2.1743656) q[2];
rz(2.1238756) q[3];
sx q[3];
rz(-1.8571721) q[3];
sx q[3];
rz(2.4093936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.9781072) q[0];
sx q[0];
rz(-2.0757984) q[0];
sx q[0];
rz(2.1708873) q[0];
rz(-0.12374395) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(-0.32044784) q[2];
sx q[2];
rz(-1.6978227) q[2];
sx q[2];
rz(2.6200079) q[2];
rz(0.10931482) q[3];
sx q[3];
rz(-2.8277665) q[3];
sx q[3];
rz(1.2841429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
