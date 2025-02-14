OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8785414) q[0];
sx q[0];
rz(-0.6113373) q[0];
sx q[0];
rz(-1.6482469) q[0];
rz(0.85236323) q[1];
sx q[1];
rz(4.8893856) q[1];
sx q[1];
rz(8.1946876) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2821101) q[0];
sx q[0];
rz(-1.9766269) q[0];
sx q[0];
rz(0.52284436) q[0];
rz(-1.1793433) q[2];
sx q[2];
rz(-1.149893) q[2];
sx q[2];
rz(2.1564623) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7820471) q[1];
sx q[1];
rz(-0.57947928) q[1];
sx q[1];
rz(2.2234099) q[1];
x q[2];
rz(1.8147179) q[3];
sx q[3];
rz(-2.8496242) q[3];
sx q[3];
rz(-2.1413745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7386231) q[2];
sx q[2];
rz(-1.3961184) q[2];
sx q[2];
rz(1.9782861) q[2];
rz(0.93825424) q[3];
sx q[3];
rz(-2.3353751) q[3];
sx q[3];
rz(-1.4279385) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52629483) q[0];
sx q[0];
rz(-1.2328923) q[0];
sx q[0];
rz(1.1616608) q[0];
rz(1.42234) q[1];
sx q[1];
rz(-1.6484345) q[1];
sx q[1];
rz(3.0551522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0491701) q[0];
sx q[0];
rz(-1.0455275) q[0];
sx q[0];
rz(2.9591333) q[0];
rz(-pi) q[1];
rz(0.49098726) q[2];
sx q[2];
rz(-1.3730959) q[2];
sx q[2];
rz(-3.1065772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.045043) q[1];
sx q[1];
rz(-0.031686671) q[1];
sx q[1];
rz(-0.50848728) q[1];
rz(-pi) q[2];
rz(-1.2041041) q[3];
sx q[3];
rz(-2.9584998) q[3];
sx q[3];
rz(-2.2568373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5689759) q[2];
sx q[2];
rz(-2.6628351) q[2];
sx q[2];
rz(-3.0968481) q[2];
rz(0.71988002) q[3];
sx q[3];
rz(-1.6241112) q[3];
sx q[3];
rz(-1.6911136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2429263) q[0];
sx q[0];
rz(-1.9654322) q[0];
sx q[0];
rz(-1.9157008) q[0];
rz(-2.2236845) q[1];
sx q[1];
rz(-1.2914912) q[1];
sx q[1];
rz(1.8570522) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51581406) q[0];
sx q[0];
rz(-2.6373793) q[0];
sx q[0];
rz(1.4674241) q[0];
x q[1];
rz(1.5628556) q[2];
sx q[2];
rz(-0.94915945) q[2];
sx q[2];
rz(-1.7253184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0902024) q[1];
sx q[1];
rz(-1.0388684) q[1];
sx q[1];
rz(2.4064785) q[1];
rz(-1.8281047) q[3];
sx q[3];
rz(-0.9633102) q[3];
sx q[3];
rz(1.846254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7104177) q[2];
sx q[2];
rz(-1.147889) q[2];
sx q[2];
rz(0.61719027) q[2];
rz(0.91444531) q[3];
sx q[3];
rz(-0.72943288) q[3];
sx q[3];
rz(2.7845553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0964088) q[0];
sx q[0];
rz(-3.0953188) q[0];
sx q[0];
rz(0.98840493) q[0];
rz(1.4811966) q[1];
sx q[1];
rz(-0.55440569) q[1];
sx q[1];
rz(-0.60715094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6364215) q[0];
sx q[0];
rz(-1.2984411) q[0];
sx q[0];
rz(0.26557458) q[0];
rz(-pi) q[1];
rz(-0.49272196) q[2];
sx q[2];
rz(-2.4396076) q[2];
sx q[2];
rz(-1.6175601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1616832) q[1];
sx q[1];
rz(-1.9075127) q[1];
sx q[1];
rz(-2.0585795) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0801717) q[3];
sx q[3];
rz(-0.20397025) q[3];
sx q[3];
rz(-0.56817618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0854757) q[2];
sx q[2];
rz(-1.3885219) q[2];
sx q[2];
rz(-0.2307387) q[2];
rz(-2.3750316) q[3];
sx q[3];
rz(-2.9975588) q[3];
sx q[3];
rz(-1.8421596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615999) q[0];
sx q[0];
rz(-1.5138641) q[0];
sx q[0];
rz(3.0738714) q[0];
rz(-1.3041152) q[1];
sx q[1];
rz(-1.8505406) q[1];
sx q[1];
rz(-1.7052604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50698603) q[0];
sx q[0];
rz(-1.9546967) q[0];
sx q[0];
rz(-0.83457729) q[0];
rz(-pi) q[1];
rz(-3.0313644) q[2];
sx q[2];
rz(-2.0617635) q[2];
sx q[2];
rz(1.5256745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41714121) q[1];
sx q[1];
rz(-1.0994689) q[1];
sx q[1];
rz(0.36496867) q[1];
x q[2];
rz(2.1077638) q[3];
sx q[3];
rz(-2.7269533) q[3];
sx q[3];
rz(1.6920167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0962254) q[2];
sx q[2];
rz(-2.2787091) q[2];
sx q[2];
rz(0.53675845) q[2];
rz(1.9219575) q[3];
sx q[3];
rz(-2.2370971) q[3];
sx q[3];
rz(1.9513244) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9959975) q[0];
sx q[0];
rz(-0.97766101) q[0];
sx q[0];
rz(-1.4033432) q[0];
rz(0.504269) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(2.9020342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5546416) q[0];
sx q[0];
rz(-1.7514146) q[0];
sx q[0];
rz(0.15968542) q[0];
x q[1];
rz(-2.1920565) q[2];
sx q[2];
rz(-0.64709787) q[2];
sx q[2];
rz(2.4280649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9592683) q[1];
sx q[1];
rz(-0.9048008) q[1];
sx q[1];
rz(1.4592501) q[1];
x q[2];
rz(2.2277545) q[3];
sx q[3];
rz(-0.60937928) q[3];
sx q[3];
rz(2.1003124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2704894) q[2];
sx q[2];
rz(-1.5364545) q[2];
sx q[2];
rz(2.8883873) q[2];
rz(0.95868715) q[3];
sx q[3];
rz(-0.91512338) q[3];
sx q[3];
rz(1.2747214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36509982) q[0];
sx q[0];
rz(-2.4204142) q[0];
sx q[0];
rz(1.391063) q[0];
rz(0.59025383) q[1];
sx q[1];
rz(-1.9275459) q[1];
sx q[1];
rz(-0.50277695) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816526) q[0];
sx q[0];
rz(-1.5852196) q[0];
sx q[0];
rz(-1.676286) q[0];
rz(-pi) q[1];
rz(-2.6739357) q[2];
sx q[2];
rz(-2.4291647) q[2];
sx q[2];
rz(-1.4234655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.630809) q[1];
sx q[1];
rz(-1.8165339) q[1];
sx q[1];
rz(2.3151957) q[1];
rz(-pi) q[2];
rz(1.8842949) q[3];
sx q[3];
rz(-1.0491706) q[3];
sx q[3];
rz(-0.23594638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7949924) q[2];
sx q[2];
rz(-1.9575926) q[2];
sx q[2];
rz(0.45197519) q[2];
rz(-2.82708) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(-3.0934635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5101584) q[0];
sx q[0];
rz(-0.87961125) q[0];
sx q[0];
rz(1.6360224) q[0];
rz(0.10558852) q[1];
sx q[1];
rz(-2.0629203) q[1];
sx q[1];
rz(0.67965913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1733579) q[0];
sx q[0];
rz(-2.2799186) q[0];
sx q[0];
rz(-2.3308332) q[0];
rz(-pi) q[1];
rz(-2.1887052) q[2];
sx q[2];
rz(-1.2373239) q[2];
sx q[2];
rz(3.0798387) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42592827) q[1];
sx q[1];
rz(-1.1008465) q[1];
sx q[1];
rz(-3.1380395) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7962436) q[3];
sx q[3];
rz(-1.3890919) q[3];
sx q[3];
rz(-2.7971706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6203561) q[2];
sx q[2];
rz(-0.53400365) q[2];
sx q[2];
rz(0.83068577) q[2];
rz(-1.0384809) q[3];
sx q[3];
rz(-2.4887648) q[3];
sx q[3];
rz(1.7415107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251544) q[0];
sx q[0];
rz(-1.1517628) q[0];
sx q[0];
rz(1.0333767) q[0];
rz(1.0721463) q[1];
sx q[1];
rz(-1.9629581) q[1];
sx q[1];
rz(0.56195608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852405) q[0];
sx q[0];
rz(-2.4187208) q[0];
sx q[0];
rz(-1.3600574) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2267954) q[2];
sx q[2];
rz(-1.8567698) q[2];
sx q[2];
rz(-0.15483072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9431253) q[1];
sx q[1];
rz(-2.2930305) q[1];
sx q[1];
rz(-0.46015443) q[1];
rz(1.5132684) q[3];
sx q[3];
rz(-1.5794601) q[3];
sx q[3];
rz(-1.6144365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0434887) q[2];
sx q[2];
rz(-1.6489112) q[2];
sx q[2];
rz(1.4081504) q[2];
rz(-3.0812541) q[3];
sx q[3];
rz(-1.3121759) q[3];
sx q[3];
rz(-1.8797125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14335808) q[0];
sx q[0];
rz(-0.77333212) q[0];
sx q[0];
rz(-1.3704569) q[0];
rz(1.0073608) q[1];
sx q[1];
rz(-1.3887364) q[1];
sx q[1];
rz(-3.0000684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816982) q[0];
sx q[0];
rz(-2.0260677) q[0];
sx q[0];
rz(-2.0596402) q[0];
rz(-1.6312286) q[2];
sx q[2];
rz(-1.5055033) q[2];
sx q[2];
rz(-3.1217965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87301577) q[1];
sx q[1];
rz(-2.0076224) q[1];
sx q[1];
rz(-0.42307968) q[1];
rz(-0.022079682) q[3];
sx q[3];
rz(-1.0754657) q[3];
sx q[3];
rz(1.2091523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1141899) q[2];
sx q[2];
rz(-1.9582615) q[2];
sx q[2];
rz(1.2664504) q[2];
rz(3.0053084) q[3];
sx q[3];
rz(-2.2253021) q[3];
sx q[3];
rz(0.19792476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4895353) q[0];
sx q[0];
rz(-2.1283524) q[0];
sx q[0];
rz(1.5555489) q[0];
rz(2.1736705) q[1];
sx q[1];
rz(-0.5236917) q[1];
sx q[1];
rz(2.0655469) q[1];
rz(-2.3308771) q[2];
sx q[2];
rz(-1.6494807) q[2];
sx q[2];
rz(-2.6930566) q[2];
rz(-1.380873) q[3];
sx q[3];
rz(-1.6217762) q[3];
sx q[3];
rz(-3.1358596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
