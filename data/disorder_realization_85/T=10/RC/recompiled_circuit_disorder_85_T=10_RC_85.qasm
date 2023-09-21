OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(-1.7655656) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-1.1321804) q[0];
sx q[0];
rz(2.2493275) q[0];
x q[1];
rz(-1.0797834) q[2];
sx q[2];
rz(-2.2248189) q[2];
sx q[2];
rz(-2.430928) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.260533) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(1.1313603) q[1];
rz(-pi) q[2];
rz(-0.81935482) q[3];
sx q[3];
rz(-1.6091533) q[3];
sx q[3];
rz(-2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(0.72431272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99740072) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(-0.86667592) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19947796) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(2.2897838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38072941) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(1.0357344) q[1];
rz(-pi) q[2];
rz(1.3734666) q[3];
sx q[3];
rz(-1.8036246) q[3];
sx q[3];
rz(-2.3372834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(0.36188564) q[2];
rz(-0.13606717) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(1.746159) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.9225072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1299382) q[0];
sx q[0];
rz(-1.219795) q[0];
sx q[0];
rz(0.28052335) q[0];
x q[1];
rz(-1.9169541) q[2];
sx q[2];
rz(-2.4097754) q[2];
sx q[2];
rz(3.0423622) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-0.59191275) q[1];
sx q[1];
rz(2.0134258) q[1];
rz(-pi) q[2];
rz(-0.015467042) q[3];
sx q[3];
rz(-1.3028212) q[3];
sx q[3];
rz(1.0166575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(2.3148361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087238833) q[0];
sx q[0];
rz(-2.0728489) q[0];
sx q[0];
rz(2.0421844) q[0];
rz(-pi) q[1];
rz(2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(1.3354288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0485059) q[1];
sx q[1];
rz(-0.48492453) q[1];
sx q[1];
rz(2.4946458) q[1];
x q[2];
rz(1.2590253) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(0.28516969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42442214) q[0];
sx q[0];
rz(-1.6132014) q[0];
sx q[0];
rz(-2.0457343) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6599777) q[2];
sx q[2];
rz(-0.72313213) q[2];
sx q[2];
rz(-0.30578223) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3045029) q[1];
sx q[1];
rz(-1.2577406) q[1];
sx q[1];
rz(0.00043991107) q[1];
x q[2];
rz(0.22484803) q[3];
sx q[3];
rz(-0.41089155) q[3];
sx q[3];
rz(1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(2.6873798) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(2.9845797) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(1.3670115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0317504) q[0];
sx q[0];
rz(-1.6037918) q[0];
sx q[0];
rz(-1.5026291) q[0];
rz(-0.579367) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(0.52057779) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4821266) q[1];
sx q[1];
rz(-0.30059338) q[1];
sx q[1];
rz(1.8468922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8664687) q[3];
sx q[3];
rz(-2.776006) q[3];
sx q[3];
rz(0.58055731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(-2.6223555) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0872333) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(-1.5646294) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.621775) q[2];
sx q[2];
rz(-1.4547252) q[2];
sx q[2];
rz(-2.799084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.891174) q[1];
sx q[1];
rz(-2.182057) q[1];
sx q[1];
rz(0.82959081) q[1];
rz(-2.6550754) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-2.2138514) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-2.4954605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883811) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(-1.7779011) q[0];
x q[1];
rz(-2.6693194) q[2];
sx q[2];
rz(-1.4204645) q[2];
sx q[2];
rz(-0.78778247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7039973) q[1];
sx q[1];
rz(-2.5181209) q[1];
sx q[1];
rz(0.61203476) q[1];
x q[2];
rz(-0.61026791) q[3];
sx q[3];
rz(-0.90859298) q[3];
sx q[3];
rz(2.3074647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(-0.10841766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3603044) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(-1.4772619) q[0];
x q[1];
rz(-0.046594521) q[2];
sx q[2];
rz(-1.7492883) q[2];
sx q[2];
rz(-0.5878693) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(-0.91067578) q[1];
x q[2];
rz(-2.2750862) q[3];
sx q[3];
rz(-1.4276854) q[3];
sx q[3];
rz(-2.0742311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(-2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(0.055158786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3653152) q[0];
sx q[0];
rz(-2.758983) q[0];
sx q[0];
rz(-0.17392735) q[0];
rz(-pi) q[1];
rz(2.5818985) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(0.44911227) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-2.2844237) q[1];
sx q[1];
rz(-2.39141) q[1];
rz(-pi) q[2];
rz(1.4961365) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(-0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72538439) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(2.9329119) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(-0.29253929) q[2];
sx q[2];
rz(-2.8164622) q[2];
sx q[2];
rz(-2.5278317) q[2];
rz(-2.901554) q[3];
sx q[3];
rz(-1.5272899) q[3];
sx q[3];
rz(1.8361113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
