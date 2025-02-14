OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(-1.1385318) q[0];
sx q[0];
rz(1.0983374) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(2.0055973) q[1];
sx q[1];
rz(10.027167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9838969) q[0];
sx q[0];
rz(-2.0838782) q[0];
sx q[0];
rz(-1.0668355) q[0];
x q[1];
rz(1.6260398) q[2];
sx q[2];
rz(-2.2836402) q[2];
sx q[2];
rz(-1.6875658) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.019004025) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(0.12855512) q[1];
rz(0.037558719) q[3];
sx q[3];
rz(-1.9989487) q[3];
sx q[3];
rz(1.1809837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7263553) q[2];
sx q[2];
rz(-2.0550315) q[2];
sx q[2];
rz(1.497867) q[2];
rz(3.030153) q[3];
sx q[3];
rz(-1.8569511) q[3];
sx q[3];
rz(2.1845412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9621256) q[0];
sx q[0];
rz(-2.2158556) q[0];
sx q[0];
rz(-0.16394462) q[0];
rz(2.1666849) q[1];
sx q[1];
rz(-0.68760741) q[1];
sx q[1];
rz(1.499739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83368084) q[0];
sx q[0];
rz(-1.8215067) q[0];
sx q[0];
rz(2.3019019) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4663839) q[2];
sx q[2];
rz(-2.2134502) q[2];
sx q[2];
rz(2.2810138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1999535) q[1];
sx q[1];
rz(-2.2806045) q[1];
sx q[1];
rz(-1.468424) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9256066) q[3];
sx q[3];
rz(-1.8895928) q[3];
sx q[3];
rz(0.91451281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6770596) q[2];
sx q[2];
rz(-2.634282) q[2];
sx q[2];
rz(-0.23059174) q[2];
rz(-2.4584037) q[3];
sx q[3];
rz(-2.4468827) q[3];
sx q[3];
rz(2.3901239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7815919) q[0];
sx q[0];
rz(-1.3104023) q[0];
sx q[0];
rz(-0.42136425) q[0];
rz(1.6257809) q[1];
sx q[1];
rz(-0.49585626) q[1];
sx q[1];
rz(-0.96744195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7712647) q[0];
sx q[0];
rz(-1.8974842) q[0];
sx q[0];
rz(2.8730434) q[0];
x q[1];
rz(-1.9052299) q[2];
sx q[2];
rz(-2.6504271) q[2];
sx q[2];
rz(-0.14310357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5109668) q[1];
sx q[1];
rz(-0.31019638) q[1];
sx q[1];
rz(-1.0794182) q[1];
rz(-pi) q[2];
rz(-2.580242) q[3];
sx q[3];
rz(-1.2863942) q[3];
sx q[3];
rz(-1.4903414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4030054) q[2];
sx q[2];
rz(-2.1295197) q[2];
sx q[2];
rz(0.15815132) q[2];
rz(-0.75913298) q[3];
sx q[3];
rz(-1.682155) q[3];
sx q[3];
rz(0.37091836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89449969) q[0];
sx q[0];
rz(-2.5573754) q[0];
sx q[0];
rz(1.1868813) q[0];
rz(-3.0184556) q[1];
sx q[1];
rz(-1.5712527) q[1];
sx q[1];
rz(-1.8298967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2037012) q[0];
sx q[0];
rz(-0.16898705) q[0];
sx q[0];
rz(1.0615055) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6475296) q[2];
sx q[2];
rz(-1.6729681) q[2];
sx q[2];
rz(-0.42021423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9667184) q[1];
sx q[1];
rz(-2.5481564) q[1];
sx q[1];
rz(-2.3834641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5667241) q[3];
sx q[3];
rz(-1.5713672) q[3];
sx q[3];
rz(-0.79889311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19795236) q[2];
sx q[2];
rz(-1.8215048) q[2];
sx q[2];
rz(-2.2894335) q[2];
rz(-0.85396829) q[3];
sx q[3];
rz(-1.9210509) q[3];
sx q[3];
rz(-2.061969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78344408) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(-0.070601687) q[0];
rz(-1.0380896) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(-2.3322754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22915086) q[0];
sx q[0];
rz(-0.99780592) q[0];
sx q[0];
rz(1.4848764) q[0];
rz(0.74639928) q[2];
sx q[2];
rz(-1.047102) q[2];
sx q[2];
rz(-2.6735898) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5018511) q[1];
sx q[1];
rz(-1.0305522) q[1];
sx q[1];
rz(-1.8748724) q[1];
x q[2];
rz(-1.7820939) q[3];
sx q[3];
rz(-0.70970067) q[3];
sx q[3];
rz(-3.078408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.9522379) q[2];
rz(3.0009624) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(0.3895337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73914948) q[0];
sx q[0];
rz(-0.91366714) q[0];
sx q[0];
rz(1.6074578) q[0];
rz(-0.88273478) q[1];
sx q[1];
rz(-0.84461707) q[1];
sx q[1];
rz(0.61558634) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.405571) q[0];
sx q[0];
rz(-1.6628436) q[0];
sx q[0];
rz(0.39393257) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2698055) q[2];
sx q[2];
rz(-2.1004371) q[2];
sx q[2];
rz(2.3123534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2193377) q[1];
sx q[1];
rz(-0.23592792) q[1];
sx q[1];
rz(-0.94249819) q[1];
rz(-pi) q[2];
rz(2.1911931) q[3];
sx q[3];
rz(-1.5801893) q[3];
sx q[3];
rz(0.91419807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9756056) q[2];
sx q[2];
rz(-0.96419445) q[2];
sx q[2];
rz(1.1791505) q[2];
rz(2.9278582) q[3];
sx q[3];
rz(-1.7515747) q[3];
sx q[3];
rz(0.24553044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8971276) q[0];
sx q[0];
rz(-1.1470969) q[0];
sx q[0];
rz(-3.097528) q[0];
rz(-2.7524718) q[1];
sx q[1];
rz(-0.48946425) q[1];
sx q[1];
rz(2.3997831) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866068) q[0];
sx q[0];
rz(-1.2076723) q[0];
sx q[0];
rz(1.2753049) q[0];
rz(-pi) q[1];
rz(2.7905383) q[2];
sx q[2];
rz(-0.85509713) q[2];
sx q[2];
rz(2.9911016) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85838137) q[1];
sx q[1];
rz(-1.8923396) q[1];
sx q[1];
rz(1.47846) q[1];
rz(0.23854511) q[3];
sx q[3];
rz(-2.425634) q[3];
sx q[3];
rz(1.3357287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47702181) q[2];
sx q[2];
rz(-1.7844113) q[2];
sx q[2];
rz(0.50194293) q[2];
rz(-1.4255514) q[3];
sx q[3];
rz(-0.52949667) q[3];
sx q[3];
rz(-2.3486923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75941706) q[0];
sx q[0];
rz(-1.4074396) q[0];
sx q[0];
rz(-0.39696804) q[0];
rz(-1.601903) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(-2.4698965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26123372) q[0];
sx q[0];
rz(-1.6971998) q[0];
sx q[0];
rz(-1.1221627) q[0];
rz(-2.7717088) q[2];
sx q[2];
rz(-1.1544168) q[2];
sx q[2];
rz(0.86341349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7377598) q[1];
sx q[1];
rz(-1.8813526) q[1];
sx q[1];
rz(-2.4871189) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3109452) q[3];
sx q[3];
rz(-1.1201356) q[3];
sx q[3];
rz(-0.78307952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4063065) q[2];
sx q[2];
rz(-0.51484171) q[2];
sx q[2];
rz(1.7740645) q[2];
rz(2.1788518) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(2.4875557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3108567) q[0];
sx q[0];
rz(-1.5850569) q[0];
sx q[0];
rz(2.8209525) q[0];
rz(0.93797183) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(1.4042312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9228678) q[0];
sx q[0];
rz(-1.467112) q[0];
sx q[0];
rz(1.4669048) q[0];
rz(-0.088555468) q[2];
sx q[2];
rz(-2.8635396) q[2];
sx q[2];
rz(-2.0012358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0241037) q[1];
sx q[1];
rz(-1.5908594) q[1];
sx q[1];
rz(1.5320918) q[1];
rz(-pi) q[2];
rz(2.8912084) q[3];
sx q[3];
rz(-0.77810771) q[3];
sx q[3];
rz(-0.66533662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9414901) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(0.74083677) q[2];
rz(2.0790993) q[3];
sx q[3];
rz(-1.1081568) q[3];
sx q[3];
rz(-1.9443996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2049388) q[0];
sx q[0];
rz(-0.42977253) q[0];
sx q[0];
rz(0.75576654) q[0];
rz(1.3142122) q[1];
sx q[1];
rz(-1.6270437) q[1];
sx q[1];
rz(-0.72602138) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7814815) q[0];
sx q[0];
rz(-1.6051634) q[0];
sx q[0];
rz(2.3590259) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2560202) q[2];
sx q[2];
rz(-1.6513593) q[2];
sx q[2];
rz(-1.1159814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0418813) q[1];
sx q[1];
rz(-1.2234294) q[1];
sx q[1];
rz(-0.16176407) q[1];
rz(-pi) q[2];
rz(1.1216702) q[3];
sx q[3];
rz(-1.7054929) q[3];
sx q[3];
rz(1.0659892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3296335) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(2.192396) q[2];
rz(0.1720998) q[3];
sx q[3];
rz(-2.1642919) q[3];
sx q[3];
rz(2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127472) q[0];
sx q[0];
rz(-1.5929359) q[0];
sx q[0];
rz(1.5449217) q[0];
rz(1.159809) q[1];
sx q[1];
rz(-1.8198967) q[1];
sx q[1];
rz(1.5130704) q[1];
rz(-1.415435) q[2];
sx q[2];
rz(-1.0280357) q[2];
sx q[2];
rz(0.24161777) q[2];
rz(0.070318182) q[3];
sx q[3];
rz(-2.0587772) q[3];
sx q[3];
rz(-1.8785431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
