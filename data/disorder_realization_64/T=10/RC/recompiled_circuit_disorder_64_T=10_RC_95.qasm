OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0639122) q[0];
sx q[0];
rz(-0.67718107) q[0];
sx q[0];
rz(2.1190686) q[0];
rz(-pi) q[1];
rz(1.7982499) q[2];
sx q[2];
rz(-1.4706352) q[2];
sx q[2];
rz(1.7286466) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.761258) q[1];
sx q[1];
rz(-0.77612703) q[1];
sx q[1];
rz(-2.9556729) q[1];
rz(-pi) q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(2.3232834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(0.60418207) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(-1.0634134) q[0];
rz(2.2564607) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(0.0016454776) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5981818) q[0];
sx q[0];
rz(-0.046292154) q[0];
sx q[0];
rz(-1.8020736) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64237853) q[2];
sx q[2];
rz(-1.3628236) q[2];
sx q[2];
rz(-3.1285398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.21316321) q[1];
sx q[1];
rz(-2.3191116) q[1];
sx q[1];
rz(0.33555062) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9886338) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(-2.8305588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(3.1365085) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(1.089383) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82747805) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(0.88721888) q[0];
rz(-pi) q[1];
rz(1.1459648) q[2];
sx q[2];
rz(-1.0152738) q[2];
sx q[2];
rz(2.5512763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.032420302) q[1];
sx q[1];
rz(-1.65616) q[1];
sx q[1];
rz(-1.9568155) q[1];
rz(-1.009922) q[3];
sx q[3];
rz(-2.2787333) q[3];
sx q[3];
rz(0.55299711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0638782) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(-0.88469488) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4041876) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(0.23097158) q[0];
rz(0.14881046) q[2];
sx q[2];
rz(-2.4850922) q[2];
sx q[2];
rz(1.3615001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6872014) q[1];
sx q[1];
rz(-1.193421) q[1];
sx q[1];
rz(-2.6881933) q[1];
rz(2.8344645) q[3];
sx q[3];
rz(-1.9046475) q[3];
sx q[3];
rz(3.0577554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.7255406) q[0];
rz(0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(1.0353154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028776289) q[0];
sx q[0];
rz(-2.007764) q[0];
sx q[0];
rz(-0.90941888) q[0];
rz(-pi) q[1];
rz(1.5315227) q[2];
sx q[2];
rz(-2.0619259) q[2];
sx q[2];
rz(-2.479535) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(-3.0461237) q[1];
x q[2];
rz(-1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(-2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(-2.9343658) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743054) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(-0.62367237) q[0];
rz(-pi) q[1];
rz(-1.1057165) q[2];
sx q[2];
rz(-2.1207223) q[2];
sx q[2];
rz(1.9351026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44240272) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(1.3139903) q[1];
rz(-pi) q[2];
rz(-1.3979785) q[3];
sx q[3];
rz(-2.380905) q[3];
sx q[3];
rz(0.93483227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-2.4198789) q[2];
rz(-1.2747814) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-0.73202837) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7043982) q[0];
sx q[0];
rz(-1.6219553) q[0];
sx q[0];
rz(0.41914661) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8740011) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(0.50819699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3953494) q[1];
sx q[1];
rz(-0.33848539) q[1];
sx q[1];
rz(0.45792087) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52280207) q[3];
sx q[3];
rz(-1.8409981) q[3];
sx q[3];
rz(3.0772046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(0.83731246) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(2.1829139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230969) q[0];
sx q[0];
rz(-0.76508689) q[0];
sx q[0];
rz(-2.3575248) q[0];
rz(-pi) q[1];
rz(1.1218698) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(-1.6758855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9194591) q[1];
sx q[1];
rz(-0.53680116) q[1];
sx q[1];
rz(-1.3608576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58416768) q[3];
sx q[3];
rz(-1.6730047) q[3];
sx q[3];
rz(-0.16971961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(-0.46554309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532928) q[0];
sx q[0];
rz(-0.71209891) q[0];
sx q[0];
rz(3.0112991) q[0];
rz(-pi) q[1];
rz(1.2836254) q[2];
sx q[2];
rz(-0.42771491) q[2];
sx q[2];
rz(0.62921333) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8852639) q[1];
sx q[1];
rz(-1.1210821) q[1];
sx q[1];
rz(-1.1880258) q[1];
x q[2];
rz(-2.8903928) q[3];
sx q[3];
rz(-1.2217055) q[3];
sx q[3];
rz(0.72124764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(-2.8695316) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(2.5591992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84401417) q[0];
sx q[0];
rz(-1.3676924) q[0];
sx q[0];
rz(-0.88589478) q[0];
rz(-1.0803797) q[2];
sx q[2];
rz(-1.4368125) q[2];
sx q[2];
rz(1.6756563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6873684) q[1];
sx q[1];
rz(-0.2174938) q[1];
sx q[1];
rz(0.99517676) q[1];
rz(-3.0155229) q[3];
sx q[3];
rz(-0.21890103) q[3];
sx q[3];
rz(2.253988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(1.8021884) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(-1.8680686) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(2.6639832) q[3];
sx q[3];
rz(-2.0255247) q[3];
sx q[3];
rz(-0.075102641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
