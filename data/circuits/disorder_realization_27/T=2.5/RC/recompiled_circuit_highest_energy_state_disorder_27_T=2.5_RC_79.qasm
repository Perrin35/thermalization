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
rz(-1.6838411) q[0];
sx q[0];
rz(-0.10310752) q[0];
sx q[0];
rz(-2.4007894) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(3.8771602) q[1];
sx q[1];
rz(9.1280042) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73707092) q[0];
sx q[0];
rz(-1.9849791) q[0];
sx q[0];
rz(-2.7254654) q[0];
rz(-pi) q[1];
rz(-0.27539092) q[2];
sx q[2];
rz(-2.4045334) q[2];
sx q[2];
rz(3.0050957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1512013) q[1];
sx q[1];
rz(-1.5784653) q[1];
sx q[1];
rz(0.5504613) q[1];
rz(-pi) q[2];
rz(-2.0127679) q[3];
sx q[3];
rz(-1.4085341) q[3];
sx q[3];
rz(-1.6236537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2510117) q[2];
sx q[2];
rz(-0.0958395) q[2];
sx q[2];
rz(1.4061692) q[2];
rz(0.74875325) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(-0.21193084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9894079) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(0.33921355) q[0];
rz(0.57505125) q[1];
sx q[1];
rz(-1.1851858) q[1];
sx q[1];
rz(-2.7599879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0038892) q[0];
sx q[0];
rz(-0.48332542) q[0];
sx q[0];
rz(2.0893731) q[0];
rz(2.8782042) q[2];
sx q[2];
rz(-1.2424505) q[2];
sx q[2];
rz(1.245861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4424628) q[1];
sx q[1];
rz(-0.3990631) q[1];
sx q[1];
rz(1.5956053) q[1];
rz(-pi) q[2];
rz(-2.6032002) q[3];
sx q[3];
rz(-0.88148553) q[3];
sx q[3];
rz(-1.6834843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40994689) q[2];
sx q[2];
rz(-2.4531328) q[2];
sx q[2];
rz(-0.48639578) q[2];
rz(-2.5977123) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(0.16415183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.630702) q[0];
sx q[0];
rz(-2.7257901) q[0];
sx q[0];
rz(2.1422332) q[0];
rz(-0.73127812) q[1];
sx q[1];
rz(-2.6731698) q[1];
sx q[1];
rz(2.9670002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098487735) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(1.5652324) q[0];
rz(-pi) q[1];
x q[1];
rz(1.393494) q[2];
sx q[2];
rz(-1.7288704) q[2];
sx q[2];
rz(-0.47249913) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.021589702) q[1];
sx q[1];
rz(-1.6232511) q[1];
sx q[1];
rz(-1.4823556) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1483795) q[3];
sx q[3];
rz(-0.49634051) q[3];
sx q[3];
rz(-1.5143192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38873765) q[2];
sx q[2];
rz(-2.2896705) q[2];
sx q[2];
rz(-1.277415) q[2];
rz(-0.82469213) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(-2.9741014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.326062) q[0];
sx q[0];
rz(-0.46634665) q[0];
sx q[0];
rz(1.1482358) q[0];
rz(0.3321906) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(-1.0848328) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058131071) q[0];
sx q[0];
rz(-1.8634184) q[0];
sx q[0];
rz(-0.14614848) q[0];
rz(2.3246194) q[2];
sx q[2];
rz(-1.8228056) q[2];
sx q[2];
rz(-1.9793881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9124604) q[1];
sx q[1];
rz(-2.3082151) q[1];
sx q[1];
rz(1.6417437) q[1];
rz(-pi) q[2];
rz(1.1365436) q[3];
sx q[3];
rz(-1.1449185) q[3];
sx q[3];
rz(-1.1970465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3499902) q[2];
sx q[2];
rz(-0.76452667) q[2];
sx q[2];
rz(-0.0066268607) q[2];
rz(-2.918112) q[3];
sx q[3];
rz(-2.1036744) q[3];
sx q[3];
rz(0.82911432) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70581907) q[0];
sx q[0];
rz(-2.3888102) q[0];
sx q[0];
rz(-1.9594132) q[0];
rz(-1.8403107) q[1];
sx q[1];
rz(-0.92295206) q[1];
sx q[1];
rz(0.15929793) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2625339) q[0];
sx q[0];
rz(-2.8516483) q[0];
sx q[0];
rz(-0.13061173) q[0];
rz(-pi) q[1];
rz(1.2550687) q[2];
sx q[2];
rz(-1.3906533) q[2];
sx q[2];
rz(0.9100998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22896773) q[1];
sx q[1];
rz(-1.1054313) q[1];
sx q[1];
rz(1.0292589) q[1];
x q[2];
rz(-1.0491269) q[3];
sx q[3];
rz(-0.73724607) q[3];
sx q[3];
rz(0.76747074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84676877) q[2];
sx q[2];
rz(-0.20192768) q[2];
sx q[2];
rz(-1.3429886) q[2];
rz(0.35074562) q[3];
sx q[3];
rz(-1.1387768) q[3];
sx q[3];
rz(2.8158367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.89966929) q[0];
sx q[0];
rz(-0.82829183) q[0];
sx q[0];
rz(-2.9747466) q[0];
rz(0.22653656) q[1];
sx q[1];
rz(-1.8140565) q[1];
sx q[1];
rz(-2.66364) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83376965) q[0];
sx q[0];
rz(-1.2306899) q[0];
sx q[0];
rz(-0.97619636) q[0];
rz(-pi) q[1];
rz(-1.208838) q[2];
sx q[2];
rz(-0.88506341) q[2];
sx q[2];
rz(0.63206965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93343268) q[1];
sx q[1];
rz(-2.0590977) q[1];
sx q[1];
rz(0.6330755) q[1];
rz(-pi) q[2];
rz(-2.7109954) q[3];
sx q[3];
rz(-1.1599419) q[3];
sx q[3];
rz(1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53092521) q[2];
sx q[2];
rz(-1.764856) q[2];
sx q[2];
rz(1.2923856) q[2];
rz(0.90062201) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6677299) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(-1.5245755) q[0];
rz(1.698311) q[1];
sx q[1];
rz(-2.2459005) q[1];
sx q[1];
rz(-0.5009833) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058280073) q[0];
sx q[0];
rz(-2.4078363) q[0];
sx q[0];
rz(2.1408129) q[0];
x q[1];
rz(-1.8800182) q[2];
sx q[2];
rz(-1.0337325) q[2];
sx q[2];
rz(1.352965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0308548) q[1];
sx q[1];
rz(-2.1424865) q[1];
sx q[1];
rz(-0.60380377) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37338169) q[3];
sx q[3];
rz(-0.33003673) q[3];
sx q[3];
rz(-0.16597834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3524126) q[2];
sx q[2];
rz(-0.12426201) q[2];
sx q[2];
rz(-2.6782356) q[2];
rz(-3.0739259) q[3];
sx q[3];
rz(-1.8384408) q[3];
sx q[3];
rz(0.1629924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8891334) q[0];
sx q[0];
rz(-1.2754138) q[0];
sx q[0];
rz(-0.5109936) q[0];
rz(2.4976318) q[1];
sx q[1];
rz(-2.0168346) q[1];
sx q[1];
rz(0.28800979) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0181684) q[0];
sx q[0];
rz(-0.19561681) q[0];
sx q[0];
rz(0.19851144) q[0];
rz(-pi) q[1];
rz(-2.5019849) q[2];
sx q[2];
rz(-2.158463) q[2];
sx q[2];
rz(1.4195051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0102699) q[1];
sx q[1];
rz(-1.2966709) q[1];
sx q[1];
rz(2.5585973) q[1];
rz(-pi) q[2];
rz(-2.8657416) q[3];
sx q[3];
rz(-0.39182651) q[3];
sx q[3];
rz(-2.0656757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94706941) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(0.21491773) q[2];
rz(-1.8042709) q[3];
sx q[3];
rz(-0.53842068) q[3];
sx q[3];
rz(-0.18068331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27338481) q[0];
sx q[0];
rz(-2.2053563) q[0];
sx q[0];
rz(-1.0816164) q[0];
rz(-0.99010211) q[1];
sx q[1];
rz(-1.5568045) q[1];
sx q[1];
rz(2.8066011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31934798) q[0];
sx q[0];
rz(-1.2480422) q[0];
sx q[0];
rz(0.035495338) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1352409) q[2];
sx q[2];
rz(-1.1489023) q[2];
sx q[2];
rz(2.5308501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44030467) q[1];
sx q[1];
rz(-3.0106643) q[1];
sx q[1];
rz(-0.71016117) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1331691) q[3];
sx q[3];
rz(-1.1698186) q[3];
sx q[3];
rz(3.0110735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0830903) q[2];
sx q[2];
rz(-2.6164656) q[2];
sx q[2];
rz(2.7527909) q[2];
rz(0.36924103) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9487172) q[0];
sx q[0];
rz(-0.95272869) q[0];
sx q[0];
rz(0.70190758) q[0];
rz(-1.6096055) q[1];
sx q[1];
rz(-0.67370266) q[1];
sx q[1];
rz(-0.38175499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8030539) q[0];
sx q[0];
rz(-0.10333867) q[0];
sx q[0];
rz(-0.65729143) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7306546) q[2];
sx q[2];
rz(-1.4559846) q[2];
sx q[2];
rz(-1.4474335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.020479105) q[1];
sx q[1];
rz(-1.319863) q[1];
sx q[1];
rz(-0.24848715) q[1];
x q[2];
rz(-0.79910652) q[3];
sx q[3];
rz(-1.0630084) q[3];
sx q[3];
rz(0.72572177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4664885) q[2];
sx q[2];
rz(-2.1110822) q[2];
sx q[2];
rz(0.64811903) q[2];
rz(-1.0068007) q[3];
sx q[3];
rz(-1.1454134) q[3];
sx q[3];
rz(3.0692611) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61160144) q[0];
sx q[0];
rz(-1.7399104) q[0];
sx q[0];
rz(-2.4763784) q[0];
rz(0.29162677) q[1];
sx q[1];
rz(-1.3529774) q[1];
sx q[1];
rz(-1.1381961) q[1];
rz(1.9251346) q[2];
sx q[2];
rz(-1.0783429) q[2];
sx q[2];
rz(1.1781296) q[2];
rz(-2.9838647) q[3];
sx q[3];
rz(-1.1554416) q[3];
sx q[3];
rz(-0.7994061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
