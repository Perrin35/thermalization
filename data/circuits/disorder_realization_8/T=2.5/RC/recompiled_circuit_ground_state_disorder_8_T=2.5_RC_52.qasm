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
rz(-0.30081055) q[1];
sx q[1];
rz(-1.4251113) q[1];
sx q[1];
rz(-0.11597522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2271954) q[0];
sx q[0];
rz(-1.6956788) q[0];
sx q[0];
rz(1.2465121) q[0];
rz(-pi) q[1];
rz(-1.9931727) q[2];
sx q[2];
rz(-2.6184888) q[2];
sx q[2];
rz(0.98382271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.955817) q[1];
sx q[1];
rz(-0.80090947) q[1];
sx q[1];
rz(3.0545727) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7099781) q[3];
sx q[3];
rz(-2.109189) q[3];
sx q[3];
rz(2.0288426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7684324) q[2];
sx q[2];
rz(-0.065040437) q[2];
sx q[2];
rz(-1.4434641) q[2];
rz(-0.69624919) q[3];
sx q[3];
rz(-1.5158451) q[3];
sx q[3];
rz(2.7878917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732805) q[0];
sx q[0];
rz(-3.066635) q[0];
sx q[0];
rz(2.7798376) q[0];
rz(-1.7573645) q[1];
sx q[1];
rz(-0.39159602) q[1];
sx q[1];
rz(1.0272383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8105668) q[0];
sx q[0];
rz(-2.2510656) q[0];
sx q[0];
rz(2.9085338) q[0];
rz(-2.9774819) q[2];
sx q[2];
rz(-0.80399738) q[2];
sx q[2];
rz(1.6156593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9421922) q[1];
sx q[1];
rz(-0.93633274) q[1];
sx q[1];
rz(0.30997194) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22543474) q[3];
sx q[3];
rz(-1.6525558) q[3];
sx q[3];
rz(-1.1573302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0072173) q[2];
sx q[2];
rz(-2.5773039) q[2];
sx q[2];
rz(-2.1103653) q[2];
rz(-2.4550896) q[3];
sx q[3];
rz(-1.7347696) q[3];
sx q[3];
rz(0.321872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-1.4903851) q[0];
sx q[0];
rz(-1.4116766) q[0];
sx q[0];
rz(1.3516634) q[0];
rz(-1.5216113) q[1];
sx q[1];
rz(-2.9053575) q[1];
sx q[1];
rz(-1.1865541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4410858) q[0];
sx q[0];
rz(-1.94993) q[0];
sx q[0];
rz(2.6802882) q[0];
x q[1];
rz(2.988838) q[2];
sx q[2];
rz(-2.6811594) q[2];
sx q[2];
rz(2.3178315) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83543336) q[1];
sx q[1];
rz(-3.1305538) q[1];
sx q[1];
rz(0.61223642) q[1];
rz(-0.29680829) q[3];
sx q[3];
rz(-2.2389484) q[3];
sx q[3];
rz(2.6096024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1658907) q[2];
sx q[2];
rz(-1.8791135) q[2];
sx q[2];
rz(-0.28124896) q[2];
rz(-2.7665372) q[3];
sx q[3];
rz(-2.227759) q[3];
sx q[3];
rz(-0.32947549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.3220045) q[0];
sx q[0];
rz(-0.81398886) q[0];
sx q[0];
rz(0.32459146) q[0];
rz(-1.548467) q[1];
sx q[1];
rz(-0.32061583) q[1];
sx q[1];
rz(0.89210192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4504334) q[0];
sx q[0];
rz(-1.5884134) q[0];
sx q[0];
rz(-1.582028) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71132724) q[2];
sx q[2];
rz(-1.4961424) q[2];
sx q[2];
rz(-1.1911281) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1010236) q[1];
sx q[1];
rz(-1.8547579) q[1];
sx q[1];
rz(2.8718487) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1079477) q[3];
sx q[3];
rz(-2.1639898) q[3];
sx q[3];
rz(-2.193424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.061168142) q[2];
sx q[2];
rz(-1.3865546) q[2];
sx q[2];
rz(-2.5648153) q[2];
rz(0.12442496) q[3];
sx q[3];
rz(-0.7500698) q[3];
sx q[3];
rz(3.0456544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4633923) q[0];
sx q[0];
rz(-2.7233349) q[0];
sx q[0];
rz(-2.8611355) q[0];
rz(2.8434143) q[1];
sx q[1];
rz(-3.0715946) q[1];
sx q[1];
rz(0.15204522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66528377) q[0];
sx q[0];
rz(-1.6235933) q[0];
sx q[0];
rz(0.027574551) q[0];
rz(-pi) q[1];
rz(2.0934257) q[2];
sx q[2];
rz(-1.4844456) q[2];
sx q[2];
rz(0.62355838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92805198) q[1];
sx q[1];
rz(-0.593363) q[1];
sx q[1];
rz(-0.80634547) q[1];
rz(-pi) q[2];
rz(-1.5492113) q[3];
sx q[3];
rz(-0.38281554) q[3];
sx q[3];
rz(-0.60071731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8449479) q[2];
sx q[2];
rz(-1.7903906) q[2];
sx q[2];
rz(-3.1015934) q[2];
rz(0.1376888) q[3];
sx q[3];
rz(-0.44860336) q[3];
sx q[3];
rz(-0.78534809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8765151) q[0];
sx q[0];
rz(-0.10049937) q[0];
sx q[0];
rz(2.5608089) q[0];
rz(2.4526217) q[1];
sx q[1];
rz(-0.13801485) q[1];
sx q[1];
rz(-1.2977915) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9286683) q[0];
sx q[0];
rz(-1.1547233) q[0];
sx q[0];
rz(-1.1088158) q[0];
x q[1];
rz(-1.0270732) q[2];
sx q[2];
rz(-1.5158733) q[2];
sx q[2];
rz(-1.3872272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7766469) q[1];
sx q[1];
rz(-2.6876535) q[1];
sx q[1];
rz(-2.620082) q[1];
rz(1.1290324) q[3];
sx q[3];
rz(-1.0241177) q[3];
sx q[3];
rz(-2.2571079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40655228) q[2];
sx q[2];
rz(-1.3383957) q[2];
sx q[2];
rz(1.5470541) q[2];
rz(-2.7100587) q[3];
sx q[3];
rz(-1.1003234) q[3];
sx q[3];
rz(0.6507473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65289098) q[0];
sx q[0];
rz(-0.015691375) q[0];
sx q[0];
rz(0.59590644) q[0];
rz(-1.7911004) q[1];
sx q[1];
rz(-0.48521438) q[1];
sx q[1];
rz(-2.5686666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3086711) q[0];
sx q[0];
rz(-3.1317319) q[0];
sx q[0];
rz(2.8491151) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7155034) q[2];
sx q[2];
rz(-1.3182148) q[2];
sx q[2];
rz(1.1947699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73997073) q[1];
sx q[1];
rz(-2.1453806) q[1];
sx q[1];
rz(-1.9514685) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7468729) q[3];
sx q[3];
rz(-2.3449247) q[3];
sx q[3];
rz(1.358169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4493745) q[2];
sx q[2];
rz(-2.2276679) q[2];
sx q[2];
rz(2.2268028) q[2];
rz(1.6080914) q[3];
sx q[3];
rz(-2.0017109) q[3];
sx q[3];
rz(-0.99005121) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51946259) q[0];
sx q[0];
rz(-1.4254445) q[0];
sx q[0];
rz(0.064432681) q[0];
rz(-0.26647767) q[1];
sx q[1];
rz(-3.0173512) q[1];
sx q[1];
rz(1.2640094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071113599) q[0];
sx q[0];
rz(-1.3051254) q[0];
sx q[0];
rz(2.8056953) q[0];
rz(-pi) q[1];
rz(-0.70049501) q[2];
sx q[2];
rz(-1.2681172) q[2];
sx q[2];
rz(-1.2889912) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79962423) q[1];
sx q[1];
rz(-1.2576629) q[1];
sx q[1];
rz(-1.3596441) q[1];
rz(-pi) q[2];
rz(-2.1355992) q[3];
sx q[3];
rz(-2.0801615) q[3];
sx q[3];
rz(-0.53087528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5233351) q[2];
sx q[2];
rz(-1.9674415) q[2];
sx q[2];
rz(-1.0443643) q[2];
rz(-2.4339645) q[3];
sx q[3];
rz(-2.7510567) q[3];
sx q[3];
rz(1.0467168) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7758961) q[0];
sx q[0];
rz(-1.5702268) q[0];
sx q[0];
rz(2.7112992) q[0];
rz(2.4039092) q[1];
sx q[1];
rz(-3.0569515) q[1];
sx q[1];
rz(2.784909) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094848763) q[0];
sx q[0];
rz(-1.665846) q[0];
sx q[0];
rz(-2.8489091) q[0];
x q[1];
rz(1.4595152) q[2];
sx q[2];
rz(-2.8751237) q[2];
sx q[2];
rz(-2.1499718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7995924) q[1];
sx q[1];
rz(-2.1260602) q[1];
sx q[1];
rz(2.0702122) q[1];
x q[2];
rz(1.6453708) q[3];
sx q[3];
rz(-1.6901724) q[3];
sx q[3];
rz(2.125691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3714527) q[2];
sx q[2];
rz(-0.74718237) q[2];
sx q[2];
rz(-2.746197) q[2];
rz(-1.4622408) q[3];
sx q[3];
rz(-1.7489) q[3];
sx q[3];
rz(3.1126378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8545561) q[0];
sx q[0];
rz(-2.365878) q[0];
sx q[0];
rz(-2.4391644) q[0];
rz(-0.17829819) q[1];
sx q[1];
rz(-2.4874004) q[1];
sx q[1];
rz(-1.5231232) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010945646) q[0];
sx q[0];
rz(-2.6551882) q[0];
sx q[0];
rz(-2.0964699) q[0];
x q[1];
rz(-3.1270364) q[2];
sx q[2];
rz(-1.4729028) q[2];
sx q[2];
rz(2.720969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6206086) q[1];
sx q[1];
rz(-2.1335601) q[1];
sx q[1];
rz(-1.7373134) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5036292) q[3];
sx q[3];
rz(-0.59989446) q[3];
sx q[3];
rz(-1.1092345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.496333) q[2];
sx q[2];
rz(-0.96480227) q[2];
sx q[2];
rz(0.94950914) q[2];
rz(-1.1358787) q[3];
sx q[3];
rz(-0.94616008) q[3];
sx q[3];
rz(0.56376636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5180494) q[0];
sx q[0];
rz(-1.2508871) q[0];
sx q[0];
rz(2.293806) q[0];
rz(1.7060948) q[1];
sx q[1];
rz(-1.86263) q[1];
sx q[1];
rz(-2.7772171) q[1];
rz(1.6964396) q[2];
sx q[2];
rz(-0.70401618) q[2];
sx q[2];
rz(-2.8332016) q[2];
rz(-1.7617211) q[3];
sx q[3];
rz(-0.31287258) q[3];
sx q[3];
rz(2.3967299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
