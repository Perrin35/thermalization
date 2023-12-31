OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(-0.0012794415) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89844184) q[0];
sx q[0];
rz(-1.4282707) q[0];
sx q[0];
rz(1.5205163) q[0];
x q[1];
rz(-2.7156419) q[2];
sx q[2];
rz(-0.89086878) q[2];
sx q[2];
rz(1.3198864) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4349164) q[1];
sx q[1];
rz(-1.7129363) q[1];
sx q[1];
rz(-0.040277004) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48304708) q[3];
sx q[3];
rz(-2.8222198) q[3];
sx q[3];
rz(1.3787624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.3886064) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4202704) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401706) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(-2.7313822) q[0];
rz(-1.3865115) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(-2.266303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83982044) q[1];
sx q[1];
rz(-0.41026792) q[1];
sx q[1];
rz(-2.6778585) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1249173) q[3];
sx q[3];
rz(-1.8801873) q[3];
sx q[3];
rz(2.4966937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9440207) q[0];
sx q[0];
rz(-0.70957843) q[0];
sx q[0];
rz(-2.5463085) q[0];
rz(2.5580102) q[2];
sx q[2];
rz(-0.43791134) q[2];
sx q[2];
rz(-0.32354087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2521378) q[1];
sx q[1];
rz(-0.368202) q[1];
sx q[1];
rz(-0.77106573) q[1];
rz(-2.5123185) q[3];
sx q[3];
rz(-2.4089703) q[3];
sx q[3];
rz(-0.34793684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-0.11985699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.170341) q[0];
sx q[0];
rz(-0.60702885) q[0];
sx q[0];
rz(-1.6334565) q[0];
x q[1];
rz(1.4255964) q[2];
sx q[2];
rz(-1.2135266) q[2];
sx q[2];
rz(1.6574588) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2377492) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(3.1234427) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5424764) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(0.37483254) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4362674) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(1.5572085) q[0];
x q[1];
rz(0.24057062) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(-1.0742998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1357437) q[1];
sx q[1];
rz(-1.6372576) q[1];
sx q[1];
rz(-0.60484109) q[1];
x q[2];
rz(0.92091839) q[3];
sx q[3];
rz(-0.71483597) q[3];
sx q[3];
rz(-2.9361847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.7317584) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(-0.05274996) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71125644) q[0];
sx q[0];
rz(-0.38861409) q[0];
sx q[0];
rz(-1.3769763) q[0];
rz(-pi) q[1];
rz(1.5724206) q[2];
sx q[2];
rz(-0.44518984) q[2];
sx q[2];
rz(2.0507398) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3546238) q[1];
sx q[1];
rz(-1.0458535) q[1];
sx q[1];
rz(-2.900219) q[1];
rz(-1.5105641) q[3];
sx q[3];
rz(-2.4619953) q[3];
sx q[3];
rz(-0.70590245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.3247103) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727407) q[0];
sx q[0];
rz(-1.4434837) q[0];
sx q[0];
rz(-1.5936113) q[0];
x q[1];
rz(1.0075188) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(-1.8958467) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1454989) q[1];
sx q[1];
rz(-2.5719574) q[1];
sx q[1];
rz(-2.3949404) q[1];
rz(-2.2277101) q[3];
sx q[3];
rz(-1.8808639) q[3];
sx q[3];
rz(-3.0403746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(-0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(-1.8519648) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136918) q[0];
sx q[0];
rz(-1.1192338) q[0];
sx q[0];
rz(-2.5149462) q[0];
rz(-pi) q[1];
rz(-0.81823924) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(-0.73275369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21266567) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(-0.11282632) q[1];
rz(-2.4786948) q[3];
sx q[3];
rz(-1.7015966) q[3];
sx q[3];
rz(-2.7152284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(-1.0827433) q[2];
rz(3.0454214) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-2.8589378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3960421) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(2.1349483) q[0];
x q[1];
rz(-2.5731509) q[2];
sx q[2];
rz(-0.38707765) q[2];
sx q[2];
rz(2.8708411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0375367) q[1];
sx q[1];
rz(-1.9754585) q[1];
sx q[1];
rz(1.3794823) q[1];
x q[2];
rz(-2.2563124) q[3];
sx q[3];
rz(-1.4797398) q[3];
sx q[3];
rz(-1.0202927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-0.58399502) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0081351) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(1.0159147) q[0];
rz(0.85068662) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(-1.016664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0360003) q[1];
sx q[1];
rz(-0.75165527) q[1];
sx q[1];
rz(2.4113301) q[1];
rz(-pi) q[2];
rz(2.3144249) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(2.8005023) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-2.5088359) q[2];
sx q[2];
rz(-0.7706332) q[2];
sx q[2];
rz(2.3494233) q[2];
rz(-1.3689465) q[3];
sx q[3];
rz(-1.9484083) q[3];
sx q[3];
rz(1.2231135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
