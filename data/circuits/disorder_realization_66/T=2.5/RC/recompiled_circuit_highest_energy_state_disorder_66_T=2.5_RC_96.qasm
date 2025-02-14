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
rz(2.1386327) q[0];
sx q[0];
rz(-0.47337368) q[0];
sx q[0];
rz(-1.0869429) q[0];
rz(2.3776157) q[1];
sx q[1];
rz(-2.6978701) q[1];
sx q[1];
rz(0.8313764) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.000769) q[0];
sx q[0];
rz(-2.0481264) q[0];
sx q[0];
rz(0.60388375) q[0];
rz(-pi) q[1];
rz(1.126312) q[2];
sx q[2];
rz(-2.630027) q[2];
sx q[2];
rz(-1.4064521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8352141) q[1];
sx q[1];
rz(-1.6740546) q[1];
sx q[1];
rz(1.4500965) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7870907) q[3];
sx q[3];
rz(-3.0648764) q[3];
sx q[3];
rz(1.8767954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5272556) q[2];
sx q[2];
rz(-1.8680806) q[2];
sx q[2];
rz(2.7296076) q[2];
rz(2.8968503) q[3];
sx q[3];
rz(-0.63260806) q[3];
sx q[3];
rz(-2.8542724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7016542) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(-0.42020759) q[0];
rz(0.48928753) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(1.9583826) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.228738) q[0];
sx q[0];
rz(-1.1381554) q[0];
sx q[0];
rz(0.88254024) q[0];
rz(-0.31416201) q[2];
sx q[2];
rz(-1.8793686) q[2];
sx q[2];
rz(-1.4663638) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.816531) q[1];
sx q[1];
rz(-1.9516212) q[1];
sx q[1];
rz(2.517609) q[1];
x q[2];
rz(-0.11759956) q[3];
sx q[3];
rz(-1.3983852) q[3];
sx q[3];
rz(-0.93061479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6861787) q[2];
sx q[2];
rz(-1.4718461) q[2];
sx q[2];
rz(-3.0012644) q[2];
rz(-2.8651107) q[3];
sx q[3];
rz(-2.3639207) q[3];
sx q[3];
rz(-2.8592143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.33271933) q[0];
sx q[0];
rz(-2.7820899) q[0];
sx q[0];
rz(-1.7669539) q[0];
rz(-2.2287492) q[1];
sx q[1];
rz(-0.54250598) q[1];
sx q[1];
rz(-2.1106145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08374005) q[0];
sx q[0];
rz(-0.3763323) q[0];
sx q[0];
rz(-2.0876838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39856492) q[2];
sx q[2];
rz(-1.0462648) q[2];
sx q[2];
rz(2.2664859) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1304969) q[1];
sx q[1];
rz(-1.5807932) q[1];
sx q[1];
rz(1.1586055) q[1];
rz(-pi) q[2];
rz(2.8171478) q[3];
sx q[3];
rz(-1.7534755) q[3];
sx q[3];
rz(-0.35445359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.012152) q[2];
sx q[2];
rz(-2.0871711) q[2];
sx q[2];
rz(-2.5035109) q[2];
rz(3.0408527) q[3];
sx q[3];
rz(-1.0120069) q[3];
sx q[3];
rz(-0.49873275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8345555) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(-1.3822973) q[0];
rz(-2.8221829) q[1];
sx q[1];
rz(-1.6326135) q[1];
sx q[1];
rz(2.3841948) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1019615) q[0];
sx q[0];
rz(-1.5622458) q[0];
sx q[0];
rz(0.074683384) q[0];
x q[1];
rz(-1.6998197) q[2];
sx q[2];
rz(-1.7119383) q[2];
sx q[2];
rz(-0.8069969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.08257016) q[1];
sx q[1];
rz(-2.3772514) q[1];
sx q[1];
rz(-0.90103006) q[1];
rz(-2.0946461) q[3];
sx q[3];
rz(-2.110384) q[3];
sx q[3];
rz(-1.3803409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3697529) q[2];
sx q[2];
rz(-0.87551337) q[2];
sx q[2];
rz(0.80779752) q[2];
rz(-1.7317023) q[3];
sx q[3];
rz(-1.2597224) q[3];
sx q[3];
rz(-0.65214777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659502) q[0];
sx q[0];
rz(-0.37794161) q[0];
sx q[0];
rz(-0.98528969) q[0];
rz(1.4957042) q[1];
sx q[1];
rz(-2.0031877) q[1];
sx q[1];
rz(0.92934242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7791361) q[0];
sx q[0];
rz(-0.99186388) q[0];
sx q[0];
rz(-1.2499362) q[0];
x q[1];
rz(1.9787637) q[2];
sx q[2];
rz(-1.0024602) q[2];
sx q[2];
rz(-1.6108244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1869609) q[1];
sx q[1];
rz(-1.7699047) q[1];
sx q[1];
rz(1.0331927) q[1];
rz(-pi) q[2];
rz(-2.1319904) q[3];
sx q[3];
rz(-1.582904) q[3];
sx q[3];
rz(-2.5502513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0890961) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(1.0864786) q[2];
rz(-0.9440445) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(2.0199147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97359598) q[0];
sx q[0];
rz(-2.7041628) q[0];
sx q[0];
rz(0.22015372) q[0];
rz(-1.6379697) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(2.5880609) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65992113) q[0];
sx q[0];
rz(-1.4559457) q[0];
sx q[0];
rz(-0.19821313) q[0];
rz(0.36350162) q[2];
sx q[2];
rz(-1.7365626) q[2];
sx q[2];
rz(1.7004418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.625861) q[1];
sx q[1];
rz(-1.5222094) q[1];
sx q[1];
rz(-0.69339417) q[1];
rz(0.27724482) q[3];
sx q[3];
rz(-2.436815) q[3];
sx q[3];
rz(0.98990209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0452051) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(0.20827797) q[2];
rz(1.1194718) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(-2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49523062) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(2.0679423) q[0];
rz(0.13058361) q[1];
sx q[1];
rz(-1.5675631) q[1];
sx q[1];
rz(1.0098339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4376917) q[0];
sx q[0];
rz(-0.44095518) q[0];
sx q[0];
rz(2.6954014) q[0];
x q[1];
rz(-0.45744894) q[2];
sx q[2];
rz(-2.9764796) q[2];
sx q[2];
rz(-0.091006309) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6679935) q[1];
sx q[1];
rz(-0.58071857) q[1];
sx q[1];
rz(-2.9734427) q[1];
x q[2];
rz(0.44949406) q[3];
sx q[3];
rz(-0.71555863) q[3];
sx q[3];
rz(2.4982832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4449571) q[2];
sx q[2];
rz(-2.770165) q[2];
sx q[2];
rz(2.8182287) q[2];
rz(0.67445406) q[3];
sx q[3];
rz(-0.89265299) q[3];
sx q[3];
rz(-3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12482878) q[0];
sx q[0];
rz(-0.49377307) q[0];
sx q[0];
rz(2.0523409) q[0];
rz(-0.59843868) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(0.79900297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4161637) q[0];
sx q[0];
rz(-1.2667286) q[0];
sx q[0];
rz(1.2771525) q[0];
rz(-3.0472894) q[2];
sx q[2];
rz(-2.2348197) q[2];
sx q[2];
rz(1.9349328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.028513718) q[1];
sx q[1];
rz(-1.7195819) q[1];
sx q[1];
rz(0.30537511) q[1];
rz(-pi) q[2];
rz(0.24183065) q[3];
sx q[3];
rz(-1.5760836) q[3];
sx q[3];
rz(1.9092803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4716855) q[2];
sx q[2];
rz(-2.6354463) q[2];
sx q[2];
rz(1.252582) q[2];
rz(-1.0505098) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85225409) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(0.5156714) q[0];
rz(-1.6853261) q[1];
sx q[1];
rz(-1.8003502) q[1];
sx q[1];
rz(-0.32404831) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0607325) q[0];
sx q[0];
rz(-1.3172842) q[0];
sx q[0];
rz(-0.41187615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.222175) q[2];
sx q[2];
rz(-1.2044665) q[2];
sx q[2];
rz(-1.9969783) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16897136) q[1];
sx q[1];
rz(-1.7105127) q[1];
sx q[1];
rz(2.0129544) q[1];
rz(-pi) q[2];
rz(0.079599722) q[3];
sx q[3];
rz(-1.0145502) q[3];
sx q[3];
rz(1.85301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1590165) q[2];
sx q[2];
rz(-2.0197208) q[2];
sx q[2];
rz(0.70875657) q[2];
rz(0.7835663) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(1.6243352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5358955) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(2.9851483) q[0];
rz(-2.5018196) q[1];
sx q[1];
rz(-0.6548869) q[1];
sx q[1];
rz(-0.99008647) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36562409) q[0];
sx q[0];
rz(-0.98593119) q[0];
sx q[0];
rz(-0.82230277) q[0];
rz(-pi) q[1];
rz(-0.089870139) q[2];
sx q[2];
rz(-1.9192358) q[2];
sx q[2];
rz(-0.2350829) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8955947) q[1];
sx q[1];
rz(-0.57064547) q[1];
sx q[1];
rz(0.24773189) q[1];
x q[2];
rz(-0.20179468) q[3];
sx q[3];
rz(-1.4425689) q[3];
sx q[3];
rz(-1.9596254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5797552) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(-0.61335316) q[2];
rz(0.26398811) q[3];
sx q[3];
rz(-1.8050906) q[3];
sx q[3];
rz(-2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.092125208) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(0.77656219) q[2];
sx q[2];
rz(-2.6040417) q[2];
sx q[2];
rz(-0.14234409) q[2];
rz(0.61870863) q[3];
sx q[3];
rz(-1.5386464) q[3];
sx q[3];
rz(1.9428941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
