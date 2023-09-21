OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(-1.5919332) q[1];
sx q[1];
rz(-3.0631493) q[1];
sx q[1];
rz(0.6426386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519878) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(-2.6023988) q[0];
x q[1];
rz(2.6334555) q[2];
sx q[2];
rz(-1.7388441) q[2];
sx q[2];
rz(0.44067581) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3712284) q[1];
sx q[1];
rz(-2.9712354) q[1];
sx q[1];
rz(-0.78070663) q[1];
rz(-1.60728) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(2.8521188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47444433) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(1.2791963) q[2];
rz(-0.71875087) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-2.1610778) q[0];
rz(2.9837043) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(0.78871361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65608998) q[0];
sx q[0];
rz(-2.9733109) q[0];
sx q[0];
rz(-2.3165354) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12077232) q[2];
sx q[2];
rz(-0.25794068) q[2];
sx q[2];
rz(2.068012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1566369) q[1];
sx q[1];
rz(-2.5002694) q[1];
sx q[1];
rz(-2.6167469) q[1];
x q[2];
rz(0.85487811) q[3];
sx q[3];
rz(-0.93920556) q[3];
sx q[3];
rz(2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(0.186084) q[2];
rz(0.6535334) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9300951) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-2.511456) q[0];
rz(-0.12763003) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(0.72174597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8080224) q[0];
sx q[0];
rz(-0.70735065) q[0];
sx q[0];
rz(0.39444123) q[0];
rz(-pi) q[1];
rz(2.0314991) q[2];
sx q[2];
rz(-2.2113872) q[2];
sx q[2];
rz(-3.1415423) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7194781) q[1];
sx q[1];
rz(-1.1269224) q[1];
sx q[1];
rz(1.5281954) q[1];
rz(-pi) q[2];
rz(2.1980522) q[3];
sx q[3];
rz(-1.846608) q[3];
sx q[3];
rz(2.4081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-2.2325113) q[2];
rz(-0.51820731) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5575314) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-0.50023729) q[1];
sx q[1];
rz(2.3775878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15427904) q[0];
sx q[0];
rz(-0.89218441) q[0];
sx q[0];
rz(0.016088386) q[0];
rz(0.64455428) q[2];
sx q[2];
rz(-1.6719712) q[2];
sx q[2];
rz(-0.25601706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2149787) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(-2.742393) q[1];
rz(2.0352827) q[3];
sx q[3];
rz(-2.1811857) q[3];
sx q[3];
rz(0.27184091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(0.53304535) q[2];
rz(0.26238966) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9168636) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(1.8977144) q[0];
rz(0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(-2.7382543) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3161702) q[0];
sx q[0];
rz(-2.8916725) q[0];
sx q[0];
rz(2.4106246) q[0];
x q[1];
rz(2.7062347) q[2];
sx q[2];
rz(-0.73409664) q[2];
sx q[2];
rz(1.3370607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4407318) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(2.6101019) q[1];
rz(-pi) q[2];
rz(-2.1020528) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(1.8464551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(-2.3590951) q[2];
rz(-1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(-3.0420711) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(0.0064370357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2086439) q[0];
sx q[0];
rz(-2.3017578) q[0];
sx q[0];
rz(-1.2293925) q[0];
x q[1];
rz(1.8338649) q[2];
sx q[2];
rz(-0.38882133) q[2];
sx q[2];
rz(3.1281272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7240119) q[1];
sx q[1];
rz(-1.15861) q[1];
sx q[1];
rz(-2.8860303) q[1];
x q[2];
rz(0.84900093) q[3];
sx q[3];
rz(-0.95522049) q[3];
sx q[3];
rz(-1.701223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36879888) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-2.2951365) q[2];
rz(-0.99772292) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(3.085882) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6877277) q[0];
sx q[0];
rz(-1.9181607) q[0];
sx q[0];
rz(-2.0485282) q[0];
x q[1];
rz(-2.0918526) q[2];
sx q[2];
rz(-1.1546635) q[2];
sx q[2];
rz(-3.1396438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22735587) q[1];
sx q[1];
rz(-2.1319234) q[1];
sx q[1];
rz(-0.44740541) q[1];
x q[2];
rz(-2.3485687) q[3];
sx q[3];
rz(-2.6245955) q[3];
sx q[3];
rz(-1.6263863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0760076) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-2.356142) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(0.41539645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3128368) q[0];
sx q[0];
rz(-3.0506595) q[0];
sx q[0];
rz(1.998741) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6255105) q[0];
sx q[0];
rz(-1.372638) q[0];
sx q[0];
rz(2.0546191) q[0];
x q[1];
rz(1.1025238) q[2];
sx q[2];
rz(-1.6981914) q[2];
sx q[2];
rz(-2.1786736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5059698) q[1];
sx q[1];
rz(-1.1715679) q[1];
sx q[1];
rz(-1.3969621) q[1];
rz(-1.0981406) q[3];
sx q[3];
rz(-1.6695767) q[3];
sx q[3];
rz(-3.0702101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(-1.9454983) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(-0.51913613) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25871823) q[0];
sx q[0];
rz(-1.9472194) q[0];
sx q[0];
rz(-0.14278485) q[0];
rz(-pi) q[1];
rz(1.3349418) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(0.73667919) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0874333) q[1];
sx q[1];
rz(-0.79257733) q[1];
sx q[1];
rz(1.4448318) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5095652) q[3];
sx q[3];
rz(-0.29187376) q[3];
sx q[3];
rz(-2.8458418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-3.1006151) q[2];
rz(0.86769062) q[3];
sx q[3];
rz(-0.65338457) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7460019) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(1.4279667) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(-1.6428927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258543) q[0];
sx q[0];
rz(-3.0634974) q[0];
sx q[0];
rz(-1.8321091) q[0];
rz(-pi) q[1];
rz(-3.0936436) q[2];
sx q[2];
rz(-1.4158632) q[2];
sx q[2];
rz(-2.5028554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5012706) q[1];
sx q[1];
rz(-0.81201279) q[1];
sx q[1];
rz(2.0444319) q[1];
rz(1.6937709) q[3];
sx q[3];
rz(-1.3696559) q[3];
sx q[3];
rz(0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(2.995058) q[2];
rz(-2.3274029) q[3];
sx q[3];
rz(-2.345293) q[3];
sx q[3];
rz(2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2232589) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(1.2659484) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(-1.3314432) q[2];
sx q[2];
rz(-1.6806921) q[2];
sx q[2];
rz(0.42721911) q[2];
rz(0.90445789) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];