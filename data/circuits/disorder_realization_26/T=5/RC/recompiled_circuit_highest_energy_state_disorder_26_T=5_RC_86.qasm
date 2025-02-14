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
rz(-2.969279) q[0];
sx q[0];
rz(-0.20792374) q[0];
sx q[0];
rz(-1.2907668) q[0];
rz(0.15200226) q[1];
sx q[1];
rz(-2.4988889) q[1];
sx q[1];
rz(0.42833498) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9865532) q[0];
sx q[0];
rz(-1.6800092) q[0];
sx q[0];
rz(0.10515736) q[0];
rz(-pi) q[1];
rz(-0.67240388) q[2];
sx q[2];
rz(-1.8032296) q[2];
sx q[2];
rz(-1.4990136) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0103156) q[1];
sx q[1];
rz(-2.4020264) q[1];
sx q[1];
rz(-3.0775978) q[1];
x q[2];
rz(2.400757) q[3];
sx q[3];
rz(-2.8489698) q[3];
sx q[3];
rz(3.0081841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9316445) q[2];
sx q[2];
rz(-0.95982176) q[2];
sx q[2];
rz(-1.9317365) q[2];
rz(-3.0620388) q[3];
sx q[3];
rz(-2.7456561) q[3];
sx q[3];
rz(2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88456589) q[0];
sx q[0];
rz(-2.631812) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(0.89029038) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(-0.30337897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.035159) q[0];
sx q[0];
rz(-1.037853) q[0];
sx q[0];
rz(2.919801) q[0];
rz(0.25945865) q[2];
sx q[2];
rz(-2.0739809) q[2];
sx q[2];
rz(0.34215701) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1002101) q[1];
sx q[1];
rz(-0.1100556) q[1];
sx q[1];
rz(0.96918126) q[1];
rz(-pi) q[2];
rz(-0.56191617) q[3];
sx q[3];
rz(-2.4059835) q[3];
sx q[3];
rz(3.0829142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38637912) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(2.4066822) q[2];
rz(2.2928061) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(2.7809704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.3291149) q[0];
sx q[0];
rz(-0.37549967) q[0];
sx q[0];
rz(0.078711674) q[0];
rz(-0.82376897) q[1];
sx q[1];
rz(-1.0853826) q[1];
sx q[1];
rz(-1.2944006) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98662739) q[0];
sx q[0];
rz(-1.5091578) q[0];
sx q[0];
rz(-3.0810647) q[0];
rz(0.95125385) q[2];
sx q[2];
rz(-1.143537) q[2];
sx q[2];
rz(0.45742971) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0261198) q[1];
sx q[1];
rz(-1.8515665) q[1];
sx q[1];
rz(1.6927648) q[1];
x q[2];
rz(-2.2749316) q[3];
sx q[3];
rz(-1.7539644) q[3];
sx q[3];
rz(-0.66732348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89982975) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(0.81643334) q[2];
rz(-1.624931) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(0.70037705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.65275943) q[0];
sx q[0];
rz(-0.82010287) q[0];
sx q[0];
rz(-0.56831992) q[0];
rz(1.9991416) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(1.4521339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84684337) q[0];
sx q[0];
rz(-0.94943014) q[0];
sx q[0];
rz(-1.8445548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4770437) q[2];
sx q[2];
rz(-0.20186587) q[2];
sx q[2];
rz(-2.9143726) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2842015) q[1];
sx q[1];
rz(-2.439154) q[1];
sx q[1];
rz(1.0739011) q[1];
rz(2.9630757) q[3];
sx q[3];
rz(-1.9048637) q[3];
sx q[3];
rz(-2.8802383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6292754) q[2];
sx q[2];
rz(-0.53048152) q[2];
sx q[2];
rz(3.0647035) q[2];
rz(-0.39938375) q[3];
sx q[3];
rz(-2.2279584) q[3];
sx q[3];
rz(-2.5975749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15453108) q[0];
sx q[0];
rz(-2.7847544) q[0];
sx q[0];
rz(-0.29990184) q[0];
rz(1.1497644) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(2.6531175) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91575275) q[0];
sx q[0];
rz(-1.7560648) q[0];
sx q[0];
rz(0.033323296) q[0];
x q[1];
rz(2.0769507) q[2];
sx q[2];
rz(-2.1035668) q[2];
sx q[2];
rz(1.8440994) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27389535) q[1];
sx q[1];
rz(-2.5576943) q[1];
sx q[1];
rz(0.03429596) q[1];
x q[2];
rz(0.3533479) q[3];
sx q[3];
rz(-2.4690557) q[3];
sx q[3];
rz(-1.8272682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38146314) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(-1.4786973) q[2];
rz(1.1100769) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(0.64322513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396624) q[0];
sx q[0];
rz(-1.1818385) q[0];
sx q[0];
rz(2.8231743) q[0];
rz(3.1069801) q[1];
sx q[1];
rz(-0.56762677) q[1];
sx q[1];
rz(-2.6908223) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.634003) q[0];
sx q[0];
rz(-1.5849772) q[0];
sx q[0];
rz(-2.0812278) q[0];
rz(-pi) q[1];
rz(-2.320596) q[2];
sx q[2];
rz(-2.2911173) q[2];
sx q[2];
rz(-2.3483417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5391716) q[1];
sx q[1];
rz(-2.1694467) q[1];
sx q[1];
rz(-0.4954229) q[1];
x q[2];
rz(-1.9616021) q[3];
sx q[3];
rz(-0.79422659) q[3];
sx q[3];
rz(3.1193747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.047711756) q[2];
sx q[2];
rz(-2.9517951) q[2];
sx q[2];
rz(-3.0799358) q[2];
rz(3.0430074) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-0.70257598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.3848569) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(3.1016438) q[0];
rz(-0.7705676) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(-2.9133453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9271999) q[0];
sx q[0];
rz(-1.4935555) q[0];
sx q[0];
rz(1.6595675) q[0];
rz(-pi) q[1];
x q[1];
rz(1.575861) q[2];
sx q[2];
rz(-2.8354079) q[2];
sx q[2];
rz(1.926946) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3769987) q[1];
sx q[1];
rz(-1.6239163) q[1];
sx q[1];
rz(0.082100987) q[1];
x q[2];
rz(2.1114717) q[3];
sx q[3];
rz(-2.432593) q[3];
sx q[3];
rz(1.0741155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.063529) q[2];
sx q[2];
rz(-1.0588131) q[2];
sx q[2];
rz(2.883319) q[2];
rz(-2.9410948) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(2.958278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16266009) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-2.8343416) q[0];
rz(-1.2112674) q[1];
sx q[1];
rz(-0.4937506) q[1];
sx q[1];
rz(0.58120751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8756105) q[0];
sx q[0];
rz(-1.5395682) q[0];
sx q[0];
rz(-1.9755548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70971428) q[2];
sx q[2];
rz(-1.4462556) q[2];
sx q[2];
rz(1.0106076) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8381184) q[1];
sx q[1];
rz(-1.3708769) q[1];
sx q[1];
rz(2.3568627) q[1];
rz(-1.2138661) q[3];
sx q[3];
rz(-1.0634616) q[3];
sx q[3];
rz(-1.1905673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6792949) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(0.031166859) q[2];
rz(0.22552414) q[3];
sx q[3];
rz(-1.7799107) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088260055) q[0];
sx q[0];
rz(-0.044476155) q[0];
sx q[0];
rz(-0.70892507) q[0];
rz(0.21253474) q[1];
sx q[1];
rz(-1.0147076) q[1];
sx q[1];
rz(-2.2841891) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1266076) q[0];
sx q[0];
rz(-1.7005655) q[0];
sx q[0];
rz(0.0051104498) q[0];
rz(-pi) q[1];
rz(-0.51841684) q[2];
sx q[2];
rz(-2.5441493) q[2];
sx q[2];
rz(3.0068676) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.227945) q[1];
sx q[1];
rz(-2.1008607) q[1];
sx q[1];
rz(-2.0258521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69714947) q[3];
sx q[3];
rz(-1.4065885) q[3];
sx q[3];
rz(-2.177161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4400441) q[2];
sx q[2];
rz(-0.35228071) q[2];
sx q[2];
rz(-2.3360543) q[2];
rz(-0.37197477) q[3];
sx q[3];
rz(-1.493908) q[3];
sx q[3];
rz(0.39723799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032967903) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(-0.97880542) q[0];
rz(0.72273123) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(-0.61789787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2513951) q[0];
sx q[0];
rz(-2.3171436) q[0];
sx q[0];
rz(-1.9491117) q[0];
rz(1.0523791) q[2];
sx q[2];
rz(-2.1627277) q[2];
sx q[2];
rz(3.0822138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97014654) q[1];
sx q[1];
rz(-1.3916124) q[1];
sx q[1];
rz(0.1712562) q[1];
rz(-pi) q[2];
rz(0.083656351) q[3];
sx q[3];
rz(-1.4865051) q[3];
sx q[3];
rz(0.82814661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2021947) q[2];
sx q[2];
rz(-1.7859744) q[2];
sx q[2];
rz(-0.00016577684) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-1.0060468) q[3];
sx q[3];
rz(2.5463026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7547739) q[0];
sx q[0];
rz(-1.5712354) q[0];
sx q[0];
rz(-1.5729217) q[0];
rz(-1.8008925) q[1];
sx q[1];
rz(-2.0410213) q[1];
sx q[1];
rz(-1.6165728) q[1];
rz(-2.6877838) q[2];
sx q[2];
rz(-2.1593675) q[2];
sx q[2];
rz(-0.16271954) q[2];
rz(1.9289005) q[3];
sx q[3];
rz(-1.1546338) q[3];
sx q[3];
rz(3.0178937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
