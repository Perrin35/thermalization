OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72579223) q[0];
sx q[0];
rz(-1.1644726) q[0];
sx q[0];
rz(1.8902984) q[0];
rz(2.8407821) q[1];
sx q[1];
rz(-1.7164813) q[1];
sx q[1];
rz(0.11597522) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.011256) q[0];
sx q[0];
rz(-2.7948871) q[0];
sx q[0];
rz(-1.1954855) q[0];
x q[1];
rz(1.0865583) q[2];
sx q[2];
rz(-1.3645483) q[2];
sx q[2];
rz(0.95824403) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8172215) q[1];
sx q[1];
rz(-1.6332375) q[1];
sx q[1];
rz(2.3425779) q[1];
rz(-1.4316145) q[3];
sx q[3];
rz(-1.0324036) q[3];
sx q[3];
rz(-1.11275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3731602) q[2];
sx q[2];
rz(-3.0765522) q[2];
sx q[2];
rz(-1.6981286) q[2];
rz(-0.69624919) q[3];
sx q[3];
rz(-1.5158451) q[3];
sx q[3];
rz(2.7878917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732805) q[0];
sx q[0];
rz(-0.074957632) q[0];
sx q[0];
rz(0.36175501) q[0];
rz(1.3842281) q[1];
sx q[1];
rz(-2.7499966) q[1];
sx q[1];
rz(2.1143544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3310259) q[0];
sx q[0];
rz(-2.2510656) q[0];
sx q[0];
rz(2.9085338) q[0];
rz(-pi) q[1];
rz(-0.79723704) q[2];
sx q[2];
rz(-1.452871) q[2];
sx q[2];
rz(0.069534289) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9421922) q[1];
sx q[1];
rz(-2.2052599) q[1];
sx q[1];
rz(-0.30997194) q[1];
rz(-0.22543474) q[3];
sx q[3];
rz(-1.6525558) q[3];
sx q[3];
rz(-1.1573302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1343753) q[2];
sx q[2];
rz(-2.5773039) q[2];
sx q[2];
rz(-2.1103653) q[2];
rz(0.68650308) q[3];
sx q[3];
rz(-1.7347696) q[3];
sx q[3];
rz(0.321872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6512076) q[0];
sx q[0];
rz(-1.729916) q[0];
sx q[0];
rz(-1.7899293) q[0];
rz(1.6199813) q[1];
sx q[1];
rz(-0.23623513) q[1];
sx q[1];
rz(-1.9550386) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4410858) q[0];
sx q[0];
rz(-1.1916627) q[0];
sx q[0];
rz(-0.4613045) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.988838) q[2];
sx q[2];
rz(-0.4604333) q[2];
sx q[2];
rz(2.3178315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3061593) q[1];
sx q[1];
rz(-3.1305538) q[1];
sx q[1];
rz(-0.61223642) q[1];
rz(-1.2159186) q[3];
sx q[3];
rz(-2.4198341) q[3];
sx q[3];
rz(2.1510368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1658907) q[2];
sx q[2];
rz(-1.2624792) q[2];
sx q[2];
rz(0.28124896) q[2];
rz(0.37505546) q[3];
sx q[3];
rz(-0.91383362) q[3];
sx q[3];
rz(-2.8121172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220045) q[0];
sx q[0];
rz(-2.3276038) q[0];
sx q[0];
rz(2.8170012) q[0];
rz(-1.548467) q[1];
sx q[1];
rz(-0.32061583) q[1];
sx q[1];
rz(0.89210192) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88282132) q[0];
sx q[0];
rz(-0.020892512) q[0];
sx q[0];
rz(-2.5740795) q[0];
rz(-pi) q[1];
rz(0.11406819) q[2];
sx q[2];
rz(-2.4270396) q[2];
sx q[2];
rz(0.29334208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2617787) q[1];
sx q[1];
rz(-2.752465) q[1];
sx q[1];
rz(-2.3107982) q[1];
x q[2];
rz(-0.033644955) q[3];
sx q[3];
rz(-2.1639898) q[3];
sx q[3];
rz(2.193424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0804245) q[2];
sx q[2];
rz(-1.3865546) q[2];
sx q[2];
rz(2.5648153) q[2];
rz(-0.12442496) q[3];
sx q[3];
rz(-0.7500698) q[3];
sx q[3];
rz(-3.0456544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4633923) q[0];
sx q[0];
rz(-2.7233349) q[0];
sx q[0];
rz(2.8611355) q[0];
rz(0.29817835) q[1];
sx q[1];
rz(-0.069998048) q[1];
sx q[1];
rz(0.15204522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2346245) q[0];
sx q[0];
rz(-1.5432602) q[0];
sx q[0];
rz(-1.6236134) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0481669) q[2];
sx q[2];
rz(-1.657147) q[2];
sx q[2];
rz(-2.5180343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3144532) q[1];
sx q[1];
rz(-1.9681853) q[1];
sx q[1];
rz(2.023815) q[1];
x q[2];
rz(1.5492113) q[3];
sx q[3];
rz(-0.38281554) q[3];
sx q[3];
rz(0.60071731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8449479) q[2];
sx q[2];
rz(-1.351202) q[2];
sx q[2];
rz(-0.039999261) q[2];
rz(-0.1376888) q[3];
sx q[3];
rz(-2.6929893) q[3];
sx q[3];
rz(2.3562446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26507759) q[0];
sx q[0];
rz(-3.0410933) q[0];
sx q[0];
rz(2.5608089) q[0];
rz(2.4526217) q[1];
sx q[1];
rz(-3.0035778) q[1];
sx q[1];
rz(-1.8438011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15928024) q[0];
sx q[0];
rz(-1.1508688) q[0];
sx q[0];
rz(-2.6830656) q[0];
rz(-pi) q[1];
rz(-1.6766729) q[2];
sx q[2];
rz(-2.5953802) q[2];
sx q[2];
rz(0.27413163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93377289) q[1];
sx q[1];
rz(-1.1807654) q[1];
sx q[1];
rz(1.8092501) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59249119) q[3];
sx q[3];
rz(-1.1969229) q[3];
sx q[3];
rz(-0.44523525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7350404) q[2];
sx q[2];
rz(-1.8031969) q[2];
sx q[2];
rz(1.5470541) q[2];
rz(0.43153396) q[3];
sx q[3];
rz(-1.1003234) q[3];
sx q[3];
rz(0.6507473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887017) q[0];
sx q[0];
rz(-3.1259013) q[0];
sx q[0];
rz(0.59590644) q[0];
rz(1.3504922) q[1];
sx q[1];
rz(-2.6563783) q[1];
sx q[1];
rz(2.5686666) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030339) q[0];
sx q[0];
rz(-1.5679533) q[0];
sx q[0];
rz(-0.009441998) q[0];
rz(1.2946044) q[2];
sx q[2];
rz(-1.9825299) q[2];
sx q[2];
rz(0.48897435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.044989) q[1];
sx q[1];
rz(-1.2536238) q[1];
sx q[1];
rz(-2.5326292) q[1];
x q[2];
rz(-1.9455276) q[3];
sx q[3];
rz(-0.84991377) q[3];
sx q[3];
rz(-1.8955613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4493745) q[2];
sx q[2];
rz(-0.91392475) q[2];
sx q[2];
rz(-2.2268028) q[2];
rz(-1.6080914) q[3];
sx q[3];
rz(-1.1398818) q[3];
sx q[3];
rz(2.1515414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6221301) q[0];
sx q[0];
rz(-1.4254445) q[0];
sx q[0];
rz(3.07716) q[0];
rz(2.875115) q[1];
sx q[1];
rz(-0.12424145) q[1];
sx q[1];
rz(-1.2640094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5910931) q[0];
sx q[0];
rz(-1.2471218) q[0];
sx q[0];
rz(-1.2901911) q[0];
rz(-0.70049501) q[2];
sx q[2];
rz(-1.2681172) q[2];
sx q[2];
rz(1.8526015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19169584) q[1];
sx q[1];
rz(-2.765871) q[1];
sx q[1];
rz(-0.57447489) q[1];
x q[2];
rz(0.76412075) q[3];
sx q[3];
rz(-2.4001038) q[3];
sx q[3];
rz(-1.4457857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5233351) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(-2.0972283) q[2];
rz(0.70762819) q[3];
sx q[3];
rz(-0.39053598) q[3];
sx q[3];
rz(-1.0467168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7758961) q[0];
sx q[0];
rz(-1.5713659) q[0];
sx q[0];
rz(-0.43029341) q[0];
rz(-0.73768342) q[1];
sx q[1];
rz(-0.084641181) q[1];
sx q[1];
rz(0.35668361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5045387) q[0];
sx q[0];
rz(-1.2794727) q[0];
sx q[0];
rz(-1.6700406) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8356933) q[2];
sx q[2];
rz(-1.5415493) q[2];
sx q[2];
rz(-2.4550329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.051239911) q[1];
sx q[1];
rz(-1.9899217) q[1];
sx q[1];
rz(2.5263949) q[1];
rz(-pi) q[2];
rz(-1.6453708) q[3];
sx q[3];
rz(-1.4514203) q[3];
sx q[3];
rz(2.125691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3714527) q[2];
sx q[2];
rz(-2.3944103) q[2];
sx q[2];
rz(-0.39539567) q[2];
rz(-1.4622408) q[3];
sx q[3];
rz(-1.7489) q[3];
sx q[3];
rz(3.1126378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28703654) q[0];
sx q[0];
rz(-0.77571464) q[0];
sx q[0];
rz(0.70242822) q[0];
rz(-2.9632945) q[1];
sx q[1];
rz(-2.4874004) q[1];
sx q[1];
rz(1.5231232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.130647) q[0];
sx q[0];
rz(-2.6551882) q[0];
sx q[0];
rz(2.0964699) q[0];
x q[1];
rz(1.4728925) q[2];
sx q[2];
rz(-1.5563097) q[2];
sx q[2];
rz(-1.9928428) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13924456) q[1];
sx q[1];
rz(-1.4301456) q[1];
sx q[1];
rz(0.5690677) q[1];
rz(-0.63796343) q[3];
sx q[3];
rz(-2.5416982) q[3];
sx q[3];
rz(-1.1092345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.496333) q[2];
sx q[2];
rz(-0.96480227) q[2];
sx q[2];
rz(-0.94950914) q[2];
rz(-1.1358787) q[3];
sx q[3];
rz(-0.94616008) q[3];
sx q[3];
rz(-2.5778263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5180494) q[0];
sx q[0];
rz(-1.8907056) q[0];
sx q[0];
rz(-0.84778669) q[0];
rz(-1.4354979) q[1];
sx q[1];
rz(-1.86263) q[1];
sx q[1];
rz(-2.7772171) q[1];
rz(2.2709104) q[2];
sx q[2];
rz(-1.489594) q[2];
sx q[2];
rz(1.7832047) q[2];
rz(1.8783384) q[3];
sx q[3];
rz(-1.6292385) q[3];
sx q[3];
rz(1.0077976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
