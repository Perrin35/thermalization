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
rz(-0.79987502) q[0];
sx q[0];
rz(2.3246111) q[0];
sx q[0];
rz(9.8819879) q[0];
rz(0.41181052) q[1];
sx q[1];
rz(-1.5195941) q[1];
sx q[1];
rz(-2.5817459) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074926404) q[0];
sx q[0];
rz(-1.1242965) q[0];
sx q[0];
rz(-1.2139113) q[0];
rz(-pi) q[1];
rz(1.0798321) q[2];
sx q[2];
rz(-2.5837971) q[2];
sx q[2];
rz(1.3809134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76418855) q[1];
sx q[1];
rz(-1.2930096) q[1];
sx q[1];
rz(-2.6242816) q[1];
rz(0.22038898) q[3];
sx q[3];
rz(-1.427009) q[3];
sx q[3];
rz(-1.7124364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7655699) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(-1.8454856) q[2];
rz(0.62093312) q[3];
sx q[3];
rz(-0.66578484) q[3];
sx q[3];
rz(-0.92500979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45944443) q[0];
sx q[0];
rz(-2.5542673) q[0];
sx q[0];
rz(2.4272954) q[0];
rz(2.4192877) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(1.3246271) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52377578) q[0];
sx q[0];
rz(-1.1259698) q[0];
sx q[0];
rz(2.8033048) q[0];
rz(-2.4087853) q[2];
sx q[2];
rz(-1.6116953) q[2];
sx q[2];
rz(-2.3529904) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2683619) q[1];
sx q[1];
rz(-1.2030081) q[1];
sx q[1];
rz(1.4002383) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7198661) q[3];
sx q[3];
rz(-0.6908292) q[3];
sx q[3];
rz(0.2836424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.033919949) q[2];
sx q[2];
rz(-1.718797) q[2];
sx q[2];
rz(2.1495492) q[2];
rz(-2.3658559) q[3];
sx q[3];
rz(-2.3596767) q[3];
sx q[3];
rz(1.0824664) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7139605) q[0];
sx q[0];
rz(-1.7949224) q[0];
sx q[0];
rz(2.2295075) q[0];
rz(0.38463792) q[1];
sx q[1];
rz(-0.96133989) q[1];
sx q[1];
rz(-0.49547637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3289483) q[0];
sx q[0];
rz(-1.5248796) q[0];
sx q[0];
rz(-1.4981235) q[0];
x q[1];
rz(-2.936061) q[2];
sx q[2];
rz(-2.3816817) q[2];
sx q[2];
rz(-0.854137) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66049536) q[1];
sx q[1];
rz(-2.4156791) q[1];
sx q[1];
rz(0.076151163) q[1];
rz(2.1927102) q[3];
sx q[3];
rz(-1.5601741) q[3];
sx q[3];
rz(-2.0257906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89874011) q[2];
sx q[2];
rz(-1.1246559) q[2];
sx q[2];
rz(0.43453547) q[2];
rz(-0.09856002) q[3];
sx q[3];
rz(-1.7193272) q[3];
sx q[3];
rz(0.77384531) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3944405) q[0];
sx q[0];
rz(-2.9065865) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(-2.0647743) q[1];
sx q[1];
rz(-1.7348758) q[1];
sx q[1];
rz(1.9487618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6859602) q[0];
sx q[0];
rz(-1.3979027) q[0];
sx q[0];
rz(2.8894823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13612408) q[2];
sx q[2];
rz(-1.4413712) q[2];
sx q[2];
rz(2.4905175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0620246) q[1];
sx q[1];
rz(-2.8753198) q[1];
sx q[1];
rz(2.5063199) q[1];
rz(-pi) q[2];
x q[2];
rz(1.773473) q[3];
sx q[3];
rz(-0.87966387) q[3];
sx q[3];
rz(-0.99032573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8863135) q[2];
sx q[2];
rz(-0.89526075) q[2];
sx q[2];
rz(-2.8423584) q[2];
rz(0.29087654) q[3];
sx q[3];
rz(-1.2521005) q[3];
sx q[3];
rz(2.4860184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84151477) q[0];
sx q[0];
rz(-0.97989196) q[0];
sx q[0];
rz(3.063524) q[0];
rz(2.1878751) q[1];
sx q[1];
rz(-0.51906145) q[1];
sx q[1];
rz(1.567522) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6758976) q[0];
sx q[0];
rz(-2.2392096) q[0];
sx q[0];
rz(-0.23082478) q[0];
rz(-1.2260796) q[2];
sx q[2];
rz(-2.2376752) q[2];
sx q[2];
rz(-0.67701159) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5213274) q[1];
sx q[1];
rz(-0.6129188) q[1];
sx q[1];
rz(-1.6271724) q[1];
x q[2];
rz(0.6930954) q[3];
sx q[3];
rz(-1.6247107) q[3];
sx q[3];
rz(2.0283606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2786402) q[2];
sx q[2];
rz(-0.43206698) q[2];
sx q[2];
rz(2.8734109) q[2];
rz(-0.48926085) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(1.6930273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0972524) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(-0.46446717) q[0];
rz(0.52344549) q[1];
sx q[1];
rz(-0.76952666) q[1];
sx q[1];
rz(2.4109667) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084186144) q[0];
sx q[0];
rz(-2.3825186) q[0];
sx q[0];
rz(-1.7763863) q[0];
rz(-pi) q[1];
rz(2.4695005) q[2];
sx q[2];
rz(-1.8037829) q[2];
sx q[2];
rz(-1.9547878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82416049) q[1];
sx q[1];
rz(-0.80866122) q[1];
sx q[1];
rz(2.5416994) q[1];
rz(-pi) q[2];
rz(-0.2517638) q[3];
sx q[3];
rz(-1.442896) q[3];
sx q[3];
rz(0.75899103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0475433) q[2];
sx q[2];
rz(-2.0619679) q[2];
sx q[2];
rz(0.13761061) q[2];
rz(1.1329457) q[3];
sx q[3];
rz(-1.9652941) q[3];
sx q[3];
rz(-1.8784116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15261821) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(3.1076987) q[0];
rz(2.7821817) q[1];
sx q[1];
rz(-1.5243328) q[1];
sx q[1];
rz(-0.83438897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9890404) q[0];
sx q[0];
rz(-1.1942399) q[0];
sx q[0];
rz(1.9045619) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5691146) q[2];
sx q[2];
rz(-0.77478638) q[2];
sx q[2];
rz(0.090473378) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8270478) q[1];
sx q[1];
rz(-0.27092182) q[1];
sx q[1];
rz(1.9344058) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6336176) q[3];
sx q[3];
rz(-2.2414226) q[3];
sx q[3];
rz(-2.2986029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72839165) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(1.0580074) q[2];
rz(-2.6672065) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(1.0559731) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742842) q[0];
sx q[0];
rz(-1.625076) q[0];
sx q[0];
rz(1.7838595) q[0];
rz(2.7108497) q[1];
sx q[1];
rz(-1.6669225) q[1];
sx q[1];
rz(2.6745083) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3586836) q[0];
sx q[0];
rz(-0.043178044) q[0];
sx q[0];
rz(-1.6048649) q[0];
x q[1];
rz(-1.378304) q[2];
sx q[2];
rz(-2.2644832) q[2];
sx q[2];
rz(0.81499962) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7469646) q[1];
sx q[1];
rz(-2.10022) q[1];
sx q[1];
rz(-1.549333) q[1];
rz(1.6428493) q[3];
sx q[3];
rz(-2.5024838) q[3];
sx q[3];
rz(-3.1187583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0507386) q[2];
sx q[2];
rz(-0.41768062) q[2];
sx q[2];
rz(-1.0806855) q[2];
rz(-0.68197) q[3];
sx q[3];
rz(-0.8050279) q[3];
sx q[3];
rz(2.4431156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68898106) q[0];
sx q[0];
rz(-2.6230951) q[0];
sx q[0];
rz(0.80775753) q[0];
rz(-1.0001596) q[1];
sx q[1];
rz(-2.2554485) q[1];
sx q[1];
rz(0.79661405) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8971491) q[0];
sx q[0];
rz(-1.5287491) q[0];
sx q[0];
rz(-1.6680662) q[0];
rz(-pi) q[1];
rz(-2.5162637) q[2];
sx q[2];
rz(-1.3742067) q[2];
sx q[2];
rz(-2.1670053) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1353242) q[1];
sx q[1];
rz(-0.43244967) q[1];
sx q[1];
rz(1.9073553) q[1];
rz(-pi) q[2];
rz(1.2563989) q[3];
sx q[3];
rz(-1.8125497) q[3];
sx q[3];
rz(2.7890361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5920068) q[2];
sx q[2];
rz(-1.0909811) q[2];
sx q[2];
rz(-2.8263212) q[2];
rz(1.573805) q[3];
sx q[3];
rz(-1.9939634) q[3];
sx q[3];
rz(2.018759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0386117) q[0];
sx q[0];
rz(-0.6655612) q[0];
sx q[0];
rz(2.3200206) q[0];
rz(-0.483825) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(-2.9169567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6100138) q[0];
sx q[0];
rz(-1.2481127) q[0];
sx q[0];
rz(-0.15845297) q[0];
rz(2.5098652) q[2];
sx q[2];
rz(-1.016482) q[2];
sx q[2];
rz(2.4265223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7901296) q[1];
sx q[1];
rz(-0.7143414) q[1];
sx q[1];
rz(-1.2495561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12489071) q[3];
sx q[3];
rz(-1.7833059) q[3];
sx q[3];
rz(0.24617174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75767526) q[2];
sx q[2];
rz(-3.0249247) q[2];
sx q[2];
rz(-0.095001027) q[2];
rz(1.654024) q[3];
sx q[3];
rz(-0.58949685) q[3];
sx q[3];
rz(-0.76475638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8944396) q[0];
sx q[0];
rz(-0.40774397) q[0];
sx q[0];
rz(2.8751873) q[0];
rz(-15/(8*pi)) q[1];
sx q[1];
rz(-1.4529556) q[1];
sx q[1];
rz(-1.6642889) q[1];
rz(2.0043787) q[2];
sx q[2];
rz(-2.4162393) q[2];
sx q[2];
rz(-2.2888714) q[2];
rz(3.1367875) q[3];
sx q[3];
rz(-2.7660696) q[3];
sx q[3];
rz(0.39733359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
