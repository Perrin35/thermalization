OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(4.5067956) q[0];
sx q[0];
rz(11.542008) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5380733) q[0];
sx q[0];
rz(-2.3581714) q[0];
sx q[0];
rz(0.51622434) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65806234) q[2];
sx q[2];
rz(-1.0339289) q[2];
sx q[2];
rz(-1.8368349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.484326) q[1];
sx q[1];
rz(-0.53521672) q[1];
sx q[1];
rz(0.80429299) q[1];
rz(-0.3124247) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(-1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99012016) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(2.5734076) q[0];
x q[1];
rz(2.3677164) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(-2.3351923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84320074) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(-0.028224736) q[1];
rz(-pi) q[2];
rz(-0.81289566) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(1.2113435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(-2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-2.5684165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5432376) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(-1.7217365) q[0];
rz(1.2531883) q[2];
sx q[2];
rz(-2.3177958) q[2];
sx q[2];
rz(-2.2044646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55858559) q[1];
sx q[1];
rz(-1.9531986) q[1];
sx q[1];
rz(-1.8797727) q[1];
rz(-pi) q[2];
rz(-2.7335448) q[3];
sx q[3];
rz(-1.1432262) q[3];
sx q[3];
rz(-1.6392631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(-2.0129054) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4578611) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(1.213221) q[0];
x q[1];
rz(-3.1231819) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(-1.6056431) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1945222) q[1];
sx q[1];
rz(-0.94049373) q[1];
sx q[1];
rz(-1.4428201) q[1];
x q[2];
rz(-1.4705212) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5495816) q[0];
sx q[0];
rz(-1.3322543) q[0];
sx q[0];
rz(0.15079389) q[0];
rz(-pi) q[1];
rz(0.95894496) q[2];
sx q[2];
rz(-2.358837) q[2];
sx q[2];
rz(-0.81629717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7534605) q[1];
sx q[1];
rz(-1.7963444) q[1];
sx q[1];
rz(1.7727586) q[1];
rz(-2.0586117) q[3];
sx q[3];
rz(-1.2851464) q[3];
sx q[3];
rz(0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(-2.5573964) q[0];
rz(1.9110514) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(3.0029283) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5748782) q[0];
sx q[0];
rz(-3.1209271) q[0];
sx q[0];
rz(-0.1621577) q[0];
x q[1];
rz(0.3615985) q[2];
sx q[2];
rz(-2.9271759) q[2];
sx q[2];
rz(-2.8117361) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0181959) q[1];
sx q[1];
rz(-2.1235848) q[1];
sx q[1];
rz(2.4559896) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3994282) q[3];
sx q[3];
rz(-1.2483178) q[3];
sx q[3];
rz(-0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038135197) q[0];
sx q[0];
rz(-2.3893444) q[0];
sx q[0];
rz(-2.8318846) q[0];
rz(-pi) q[1];
rz(-2.1937218) q[2];
sx q[2];
rz(-1.9588753) q[2];
sx q[2];
rz(0.58448234) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1075322) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(2.2570616) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1950486) q[3];
sx q[3];
rz(-2.4157899) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(3.0016622) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0807954) q[0];
sx q[0];
rz(-1.1757438) q[0];
sx q[0];
rz(-3.0271544) q[0];
x q[1];
rz(-2.2419937) q[2];
sx q[2];
rz(-1.1102144) q[2];
sx q[2];
rz(-0.34924289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8184549) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(-0.44169359) q[1];
rz(-pi) q[2];
rz(0.017459004) q[3];
sx q[3];
rz(-0.73988065) q[3];
sx q[3];
rz(-0.12714126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.63697469) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(-1.3508266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352167) q[0];
sx q[0];
rz(-1.750964) q[0];
sx q[0];
rz(-3.0336477) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58198858) q[2];
sx q[2];
rz(-2.4798923) q[2];
sx q[2];
rz(2.2616507) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13547922) q[1];
sx q[1];
rz(-1.0891501) q[1];
sx q[1];
rz(2.3856132) q[1];
x q[2];
rz(0.25165598) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(-2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(2.3012327) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47537739) q[2];
sx q[2];
rz(-1.4350841) q[2];
sx q[2];
rz(-3.140608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.039779546) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(-1.5323557) q[1];
x q[2];
rz(-1.4320847) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(2.9885938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13835779) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(2.6976363) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(1.8757204) q[2];
sx q[2];
rz(-3.0921428) q[2];
sx q[2];
rz(-2.5947239) q[2];
rz(0.98106445) q[3];
sx q[3];
rz(-0.55064252) q[3];
sx q[3];
rz(-2.8513089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
