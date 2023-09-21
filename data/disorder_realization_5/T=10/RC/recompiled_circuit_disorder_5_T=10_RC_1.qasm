OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(4.7441626) q[0];
sx q[0];
rz(6.8466853) q[0];
sx q[0];
rz(9.8817649) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5029992) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(0.55708416) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0511191) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(2.5475516) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1162122) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(1.2799353) q[1];
x q[2];
rz(-1.0702707) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(-2.4430032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0322545) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.6363268) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(-2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(2.3388458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4883603) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(1.0351719) q[0];
rz(-pi) q[1];
rz(-3.0375508) q[2];
sx q[2];
rz(-1.0696971) q[2];
sx q[2];
rz(0.83545557) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74233494) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(1.6595112) q[1];
rz(-1.4071776) q[3];
sx q[3];
rz(-0.67577261) q[3];
sx q[3];
rz(0.27392745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(-1.4260028) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(-0.081993016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(-2.8564575) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(-1.8018988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08683603) q[0];
sx q[0];
rz(-1.9924379) q[0];
sx q[0];
rz(2.3479168) q[0];
rz(0.3091829) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(0.14112976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6321286) q[1];
sx q[1];
rz(-0.2466991) q[1];
sx q[1];
rz(2.2798666) q[1];
x q[2];
rz(-2.6684746) q[3];
sx q[3];
rz(-2.3948673) q[3];
sx q[3];
rz(-2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(0.66657153) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4857805) q[0];
sx q[0];
rz(-2.862145) q[0];
sx q[0];
rz(0.73894545) q[0];
rz(-2.8104086) q[2];
sx q[2];
rz(-0.74197717) q[2];
sx q[2];
rz(-0.77510288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(2.3198747) q[1];
rz(-pi) q[2];
rz(-2.9922819) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(0.028545054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2864705) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(-2.4812223) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(-1.1594835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5975208) q[0];
sx q[0];
rz(-0.87771767) q[0];
sx q[0];
rz(-1.9835299) q[0];
x q[1];
rz(-2.9858066) q[2];
sx q[2];
rz(-2.163379) q[2];
sx q[2];
rz(-0.32183811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78427181) q[1];
sx q[1];
rz(-1.1457774) q[1];
sx q[1];
rz(0.97370633) q[1];
x q[2];
rz(-2.5043143) q[3];
sx q[3];
rz(-1.1253469) q[3];
sx q[3];
rz(0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7065457) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(1.2716028) q[2];
rz(-1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-0.10087068) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6504165) q[0];
sx q[0];
rz(-0.6498148) q[0];
sx q[0];
rz(-0.25214904) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2783373) q[2];
sx q[2];
rz(-2.0403701) q[2];
sx q[2];
rz(-1.8243607) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(1.3376065) q[3];
sx q[3];
rz(-1.9662074) q[3];
sx q[3];
rz(2.6449441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-2.2163088) q[2];
rz(-0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(-0.68774736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(1.2290918) q[0];
rz(1.8887927) q[2];
sx q[2];
rz(-1.3557634) q[2];
sx q[2];
rz(3.0766578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93393275) q[1];
sx q[1];
rz(-1.4414409) q[1];
sx q[1];
rz(1.594286) q[1];
rz(3.1256166) q[3];
sx q[3];
rz(-1.0025347) q[3];
sx q[3];
rz(2.361134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-2.1607416) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(-1.942873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6262561) q[0];
sx q[0];
rz(-1.6394098) q[0];
sx q[0];
rz(-0.4347883) q[0];
x q[1];
rz(-3.0555658) q[2];
sx q[2];
rz(-1.2089529) q[2];
sx q[2];
rz(-1.3081683) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3938013) q[1];
sx q[1];
rz(-1.0683021) q[1];
sx q[1];
rz(1.6277114) q[1];
x q[2];
rz(-1.3098605) q[3];
sx q[3];
rz(-1.121215) q[3];
sx q[3];
rz(-0.98228067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(-2.3525227) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25710042) q[0];
sx q[0];
rz(-1.677209) q[0];
sx q[0];
rz(-1.5110925) q[0];
rz(1.1912212) q[2];
sx q[2];
rz(-1.2588725) q[2];
sx q[2];
rz(1.003007) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0578627) q[1];
sx q[1];
rz(-1.607556) q[1];
sx q[1];
rz(-2.4576393) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4890355) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.062722) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(2.771647) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18098022) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(1.6880077) q[0];
rz(-pi) q[1];
rz(0.83689883) q[2];
sx q[2];
rz(-0.95015929) q[2];
sx q[2];
rz(-1.9377973) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93563423) q[1];
sx q[1];
rz(-1.4050254) q[1];
sx q[1];
rz(-2.2531855) q[1];
rz(3.1179908) q[3];
sx q[3];
rz(-0.58348237) q[3];
sx q[3];
rz(-2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(-0.20239057) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(2.7783685) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(1.9831052) q[2];
sx q[2];
rz(-1.9956797) q[2];
sx q[2];
rz(3.0394625) q[2];
rz(0.17584569) q[3];
sx q[3];
rz(-2.1538215) q[3];
sx q[3];
rz(0.62721696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
