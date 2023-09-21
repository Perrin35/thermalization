OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773168) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(-3.0741865) q[0];
rz(-pi) q[1];
rz(-1.921466) q[2];
sx q[2];
rz(-1.3436683) q[2];
sx q[2];
rz(2.7820058) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.623917) q[1];
sx q[1];
rz(-2.3530934) q[1];
sx q[1];
rz(2.4967525) q[1];
x q[2];
rz(-2.1360374) q[3];
sx q[3];
rz(-1.2473277) q[3];
sx q[3];
rz(0.934787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.6764486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10509051) q[0];
sx q[0];
rz(-2.4173792) q[0];
sx q[0];
rz(-0.0037395517) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25603489) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(-1.7325967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.022545594) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(-2.4382298) q[1];
rz(-1.915669) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(-1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-2.3410472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313989) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(-1.9171159) q[0];
x q[1];
rz(0.82528798) q[2];
sx q[2];
rz(-2.2824259) q[2];
sx q[2];
rz(-2.3724144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87677466) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(-3.1289711) q[1];
rz(0.3018467) q[3];
sx q[3];
rz(-2.2491124) q[3];
sx q[3];
rz(0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.7819972) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0264498) q[0];
sx q[0];
rz(-2.0949674) q[0];
sx q[0];
rz(1.3036222) q[0];
rz(-0.98767878) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(1.093986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3738721) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(0.14212455) q[1];
rz(1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(-0.15743263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.00502914) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78445804) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(-2.5255192) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9833897) q[2];
sx q[2];
rz(-1.8100097) q[2];
sx q[2];
rz(1.4520793) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2770734) q[1];
sx q[1];
rz(-2.652087) q[1];
sx q[1];
rz(0.31788748) q[1];
x q[2];
rz(-1.6793628) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(-1.997228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(2.9763124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80711354) q[0];
sx q[0];
rz(-1.665859) q[0];
sx q[0];
rz(1.6001742) q[0];
rz(-pi) q[1];
rz(0.68768244) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(2.8962367) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2416934) q[1];
sx q[1];
rz(-2.4447933) q[1];
sx q[1];
rz(0.65710575) q[1];
rz(-pi) q[2];
rz(1.1214439) q[3];
sx q[3];
rz(-1.3580139) q[3];
sx q[3];
rz(-1.3495812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(-0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38024494) q[0];
sx q[0];
rz(-1.5878829) q[0];
sx q[0];
rz(-1.6197617) q[0];
rz(-2.1963504) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(1.2298825) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2939824) q[1];
sx q[1];
rz(-0.81969205) q[1];
sx q[1];
rz(-1.3241029) q[1];
x q[2];
rz(2.001708) q[3];
sx q[3];
rz(-2.0322554) q[3];
sx q[3];
rz(2.0737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(2.2198026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066276) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(2.448003) q[0];
x q[1];
rz(0.83306649) q[2];
sx q[2];
rz(-1.7852011) q[2];
sx q[2];
rz(0.094878541) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7505155) q[1];
sx q[1];
rz(-2.8120496) q[1];
sx q[1];
rz(2.4604172) q[1];
rz(-0.17955762) q[3];
sx q[3];
rz(-2.2710544) q[3];
sx q[3];
rz(-1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(0.67970651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25526991) q[0];
sx q[0];
rz(-2.5810044) q[0];
sx q[0];
rz(2.7532996) q[0];
x q[1];
rz(0.87551261) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(-1.8347486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6913773) q[1];
sx q[1];
rz(-0.9141578) q[1];
sx q[1];
rz(-1.1378098) q[1];
x q[2];
rz(2.4217442) q[3];
sx q[3];
rz(-0.82092972) q[3];
sx q[3];
rz(1.6380701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(-2.4662468) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.1766599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424073) q[0];
sx q[0];
rz(-0.17051324) q[0];
sx q[0];
rz(1.2810345) q[0];
rz(-0.93675905) q[2];
sx q[2];
rz(-1.1144131) q[2];
sx q[2];
rz(0.87181834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6224222) q[1];
sx q[1];
rz(-1.4960438) q[1];
sx q[1];
rz(1.012411) q[1];
rz(-pi) q[2];
rz(-0.15575274) q[3];
sx q[3];
rz(-2.6253346) q[3];
sx q[3];
rz(-1.932365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7174299) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(2.7453616) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-2.4049315) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(2.7908294) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];