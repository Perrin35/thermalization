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
rz(2.4616315) q[0];
sx q[0];
rz(-2.2357219) q[0];
sx q[0];
rz(-2.048197) q[0];
rz(1.5054585) q[1];
sx q[1];
rz(-1.9688164) q[1];
sx q[1];
rz(-2.7170031) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6494006) q[0];
sx q[0];
rz(-1.369242) q[0];
sx q[0];
rz(-3.1026476) q[0];
rz(-1.7248254) q[2];
sx q[2];
rz(-1.5206778) q[2];
sx q[2];
rz(0.94238867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72480768) q[1];
sx q[1];
rz(-0.25836408) q[1];
sx q[1];
rz(-1.2577673) q[1];
rz(1.7261476) q[3];
sx q[3];
rz(-2.4253824) q[3];
sx q[3];
rz(2.9368212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37497416) q[2];
sx q[2];
rz(-1.582229) q[2];
sx q[2];
rz(1.9826822) q[2];
rz(-0.22289395) q[3];
sx q[3];
rz(-1.736172) q[3];
sx q[3];
rz(0.29002732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463773) q[0];
sx q[0];
rz(-1.615849) q[0];
sx q[0];
rz(-2.0624397) q[0];
rz(1.7201299) q[1];
sx q[1];
rz(-2.454897) q[1];
sx q[1];
rz(3.0770643) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20708974) q[0];
sx q[0];
rz(-0.10182589) q[0];
sx q[0];
rz(-1.0290716) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05414836) q[2];
sx q[2];
rz(-2.6647419) q[2];
sx q[2];
rz(2.3091174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56386212) q[1];
sx q[1];
rz(-0.81474333) q[1];
sx q[1];
rz(-1.9063063) q[1];
x q[2];
rz(-2.9546176) q[3];
sx q[3];
rz(-0.83010736) q[3];
sx q[3];
rz(-2.376698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1530389) q[2];
sx q[2];
rz(-1.7861853) q[2];
sx q[2];
rz(1.1174196) q[2];
rz(-0.85876632) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(-0.43706885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3808463) q[0];
sx q[0];
rz(-2.8570211) q[0];
sx q[0];
rz(2.9685156) q[0];
rz(2.1848047) q[1];
sx q[1];
rz(-1.0280949) q[1];
sx q[1];
rz(2.079336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05325584) q[0];
sx q[0];
rz(-0.8475786) q[0];
sx q[0];
rz(3.1331314) q[0];
x q[1];
rz(-0.71346475) q[2];
sx q[2];
rz(-2.3629945) q[2];
sx q[2];
rz(0.32003357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.838321) q[1];
sx q[1];
rz(-1.0365937) q[1];
sx q[1];
rz(-0.81373416) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5161602) q[3];
sx q[3];
rz(-2.1533826) q[3];
sx q[3];
rz(-2.7435477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.411285) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(0.74472204) q[2];
rz(-2.5005285) q[3];
sx q[3];
rz(-1.791626) q[3];
sx q[3];
rz(0.49066576) q[3];
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
rz(-1.4500126) q[0];
sx q[0];
rz(-2.5847961) q[0];
sx q[0];
rz(2.2963754) q[0];
rz(-0.40329626) q[1];
sx q[1];
rz(-1.5396298) q[1];
sx q[1];
rz(0.32803112) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16763891) q[0];
sx q[0];
rz(-0.22952794) q[0];
sx q[0];
rz(0.91201241) q[0];
rz(-0.26706605) q[2];
sx q[2];
rz(-0.24644463) q[2];
sx q[2];
rz(1.3738969) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58807948) q[1];
sx q[1];
rz(-2.5633286) q[1];
sx q[1];
rz(2.9248613) q[1];
rz(-pi) q[2];
rz(2.8680236) q[3];
sx q[3];
rz(-2.0988587) q[3];
sx q[3];
rz(-0.4116962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8959117) q[2];
sx q[2];
rz(-2.3978105) q[2];
sx q[2];
rz(0.15920676) q[2];
rz(-0.016247449) q[3];
sx q[3];
rz(-2.1135606) q[3];
sx q[3];
rz(0.73133674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943587) q[0];
sx q[0];
rz(-1.1927698) q[0];
sx q[0];
rz(-1.0900981) q[0];
rz(-3.0116426) q[1];
sx q[1];
rz(-0.81472412) q[1];
sx q[1];
rz(1.5257588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9989706) q[0];
sx q[0];
rz(-2.4283613) q[0];
sx q[0];
rz(-0.66340982) q[0];
x q[1];
rz(1.9984869) q[2];
sx q[2];
rz(-2.9614355) q[2];
sx q[2];
rz(-3.0017972) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8061749) q[1];
sx q[1];
rz(-1.5433784) q[1];
sx q[1];
rz(-2.5078678) q[1];
x q[2];
rz(-3.106254) q[3];
sx q[3];
rz(-1.3900847) q[3];
sx q[3];
rz(-0.23739761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7502363) q[2];
sx q[2];
rz(-2.3199234) q[2];
sx q[2];
rz(2.6643122) q[2];
rz(-0.20398772) q[3];
sx q[3];
rz(-2.9450649) q[3];
sx q[3];
rz(0.6240713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9452962) q[0];
sx q[0];
rz(-2.3260703) q[0];
sx q[0];
rz(-0.77769172) q[0];
rz(2.5241191) q[1];
sx q[1];
rz(-1.6505046) q[1];
sx q[1];
rz(-2.2106574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6361178) q[0];
sx q[0];
rz(-1.9808928) q[0];
sx q[0];
rz(2.0652886) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7226376) q[2];
sx q[2];
rz(-2.0134996) q[2];
sx q[2];
rz(2.1328164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1738064) q[1];
sx q[1];
rz(-1.5266982) q[1];
sx q[1];
rz(-1.2323061) q[1];
rz(-1.0369426) q[3];
sx q[3];
rz(-2.4483238) q[3];
sx q[3];
rz(1.1806837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65333873) q[2];
sx q[2];
rz(-1.3335756) q[2];
sx q[2];
rz(-0.2674357) q[2];
rz(0.93287647) q[3];
sx q[3];
rz(-1.8578015) q[3];
sx q[3];
rz(0.4944087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701508) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(-0.66584051) q[0];
rz(2.937607) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(-1.1873672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68722938) q[0];
sx q[0];
rz(-2.6626514) q[0];
sx q[0];
rz(2.9402551) q[0];
rz(-pi) q[1];
rz(-0.52296488) q[2];
sx q[2];
rz(-2.209216) q[2];
sx q[2];
rz(2.0948727) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0625035) q[1];
sx q[1];
rz(-2.3983208) q[1];
sx q[1];
rz(-0.39689245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7660308) q[3];
sx q[3];
rz(-0.25849202) q[3];
sx q[3];
rz(-0.067276567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54619706) q[2];
sx q[2];
rz(-0.51837817) q[2];
sx q[2];
rz(2.4714244) q[2];
rz(-2.8403122) q[3];
sx q[3];
rz(-1.2944841) q[3];
sx q[3];
rz(0.41845751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8381074) q[0];
sx q[0];
rz(-1.0024339) q[0];
sx q[0];
rz(-2.0759034) q[0];
rz(2.4844555) q[1];
sx q[1];
rz(-0.70825759) q[1];
sx q[1];
rz(-1.4297952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098786548) q[0];
sx q[0];
rz(-1.0142066) q[0];
sx q[0];
rz(-0.3996398) q[0];
x q[1];
rz(-1.3892608) q[2];
sx q[2];
rz(-1.8971271) q[2];
sx q[2];
rz(1.2943314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2473879) q[1];
sx q[1];
rz(-2.6726843) q[1];
sx q[1];
rz(-0.50007485) q[1];
x q[2];
rz(-2.4938857) q[3];
sx q[3];
rz(-2.9053024) q[3];
sx q[3];
rz(-1.3767729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8810205) q[2];
sx q[2];
rz(-1.842097) q[2];
sx q[2];
rz(2.1232429) q[2];
rz(-1.5492505) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(-2.3840267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26579478) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(2.1739668) q[0];
rz(2.8014917) q[1];
sx q[1];
rz(-2.3022771) q[1];
sx q[1];
rz(-1.0338773) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63003016) q[0];
sx q[0];
rz(-1.017319) q[0];
sx q[0];
rz(2.9787977) q[0];
rz(-pi) q[1];
rz(1.1393093) q[2];
sx q[2];
rz(-0.98879204) q[2];
sx q[2];
rz(-1.3912569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6714194) q[1];
sx q[1];
rz(-0.88693383) q[1];
sx q[1];
rz(-2.629423) q[1];
rz(-pi) q[2];
rz(2.9176209) q[3];
sx q[3];
rz(-1.7143608) q[3];
sx q[3];
rz(-2.3825213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0391482) q[2];
sx q[2];
rz(-1.6951963) q[2];
sx q[2];
rz(3.1026802) q[2];
rz(-0.14032042) q[3];
sx q[3];
rz(-3.0571627) q[3];
sx q[3];
rz(1.8261568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4843531) q[0];
sx q[0];
rz(-2.1359279) q[0];
sx q[0];
rz(-1.2580309) q[0];
rz(2.925442) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(2.7770538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98436873) q[0];
sx q[0];
rz(-1.8175392) q[0];
sx q[0];
rz(2.0329727) q[0];
x q[1];
rz(2.0272354) q[2];
sx q[2];
rz(-2.3203712) q[2];
sx q[2];
rz(0.26320339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8748624) q[1];
sx q[1];
rz(-1.5577496) q[1];
sx q[1];
rz(3.0860098) q[1];
rz(1.6031426) q[3];
sx q[3];
rz(-1.5129287) q[3];
sx q[3];
rz(-0.75266121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2898499) q[2];
sx q[2];
rz(-1.2052636) q[2];
sx q[2];
rz(-2.5002948) q[2];
rz(-1.4801721) q[3];
sx q[3];
rz(-1.2484173) q[3];
sx q[3];
rz(-2.4680468) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73364532) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(2.1570878) q[1];
sx q[1];
rz(-1.4823109) q[1];
sx q[1];
rz(-1.0479814) q[1];
rz(0.11646902) q[2];
sx q[2];
rz(-1.7670198) q[2];
sx q[2];
rz(2.7311538) q[2];
rz(-2.409561) q[3];
sx q[3];
rz(-0.66301262) q[3];
sx q[3];
rz(0.62779203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
