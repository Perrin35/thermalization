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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(0.79467839) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(-0.95212189) q[1];
sx q[1];
rz(-1.9917537) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2993752) q[0];
sx q[0];
rz(-1.8646075) q[0];
sx q[0];
rz(-3.1042645) q[0];
rz(2.8500184) q[2];
sx q[2];
rz(-0.85735029) q[2];
sx q[2];
rz(-2.6387362) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.036630298) q[1];
sx q[1];
rz(-1.399002) q[1];
sx q[1];
rz(2.7378332) q[1];
rz(-1.3142228) q[3];
sx q[3];
rz(-2.024641) q[3];
sx q[3];
rz(0.75405771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6940234) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(1.0688952) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(-1.4890891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25817961) q[0];
sx q[0];
rz(-3.1386107) q[0];
sx q[0];
rz(-0.51134837) q[0];
rz(1.6343575) q[1];
sx q[1];
rz(-1.8664482) q[1];
sx q[1];
rz(1.7587597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6561474) q[0];
sx q[0];
rz(-1.9535747) q[0];
sx q[0];
rz(2.9188927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34592705) q[2];
sx q[2];
rz(-0.9285766) q[2];
sx q[2];
rz(-1.8373035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.4310337) q[1];
sx q[1];
rz(-1.1712499) q[1];
sx q[1];
rz(0.33627681) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17937998) q[3];
sx q[3];
rz(-0.53032833) q[3];
sx q[3];
rz(1.3190003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4649268) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-0.10291544) q[2];
rz(-2.0788705) q[3];
sx q[3];
rz(-1.1954185) q[3];
sx q[3];
rz(-3.0265871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17387834) q[0];
sx q[0];
rz(-2.102484) q[0];
sx q[0];
rz(0.1097196) q[0];
rz(-2.8166855) q[1];
sx q[1];
rz(-1.1512681) q[1];
sx q[1];
rz(0.5330162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.652063) q[0];
sx q[0];
rz(-0.32885131) q[0];
sx q[0];
rz(2.7016599) q[0];
rz(-2.7092878) q[2];
sx q[2];
rz(-0.057261618) q[2];
sx q[2];
rz(2.6309516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17688454) q[1];
sx q[1];
rz(-1.0057276) q[1];
sx q[1];
rz(-0.76660291) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0267065) q[3];
sx q[3];
rz(-2.3264936) q[3];
sx q[3];
rz(0.86642735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32536062) q[2];
sx q[2];
rz(-1.3282789) q[2];
sx q[2];
rz(-1.6858961) q[2];
rz(0.085518941) q[3];
sx q[3];
rz(-1.0695846) q[3];
sx q[3];
rz(2.1948309) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38086712) q[0];
sx q[0];
rz(-2.4873698) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(-1.5707312) q[1];
sx q[1];
rz(-1.9941565) q[1];
sx q[1];
rz(-0.042472366) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077350294) q[0];
sx q[0];
rz(-2.0127065) q[0];
sx q[0];
rz(1.3674221) q[0];
x q[1];
rz(1.3971523) q[2];
sx q[2];
rz(-3.0499027) q[2];
sx q[2];
rz(-1.1806115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3610246) q[1];
sx q[1];
rz(-0.79611049) q[1];
sx q[1];
rz(-1.8382422) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9080584) q[3];
sx q[3];
rz(-2.326205) q[3];
sx q[3];
rz(2.5366207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.7652721) q[2];
sx q[2];
rz(-2.8975272) q[2];
rz(2.3501588) q[3];
sx q[3];
rz(-1.8297198) q[3];
sx q[3];
rz(0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0303845) q[0];
sx q[0];
rz(-0.1143488) q[0];
sx q[0];
rz(-2.9682888) q[0];
rz(2.2068842) q[1];
sx q[1];
rz(-0.84988958) q[1];
sx q[1];
rz(-2.8876143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.105071) q[0];
sx q[0];
rz(-1.6833841) q[0];
sx q[0];
rz(2.3238411) q[0];
rz(-1.9333657) q[2];
sx q[2];
rz(-1.0075648) q[2];
sx q[2];
rz(-2.8993487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9137302) q[1];
sx q[1];
rz(-2.263952) q[1];
sx q[1];
rz(1.597803) q[1];
rz(-3.0766904) q[3];
sx q[3];
rz(-1.9577454) q[3];
sx q[3];
rz(2.153819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0737334) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(-1.7013288) q[2];
rz(2.4654147) q[3];
sx q[3];
rz(-1.22236) q[3];
sx q[3];
rz(3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9192231) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(-2.4134912) q[0];
rz(1.9758196) q[1];
sx q[1];
rz(-0.89465061) q[1];
sx q[1];
rz(2.5808835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3325719) q[0];
sx q[0];
rz(-0.14223465) q[0];
sx q[0];
rz(0.47551544) q[0];
x q[1];
rz(1.4523023) q[2];
sx q[2];
rz(-1.0241707) q[2];
sx q[2];
rz(2.7989863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5821895) q[1];
sx q[1];
rz(-1.3546556) q[1];
sx q[1];
rz(-1.1690421) q[1];
rz(-1.5076163) q[3];
sx q[3];
rz(-1.4734905) q[3];
sx q[3];
rz(-2.7246876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5460983) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(-1.072849) q[2];
rz(0.32235518) q[3];
sx q[3];
rz(-0.41283804) q[3];
sx q[3];
rz(2.7505007) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748465) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(-2.7333976) q[0];
rz(2.8152668) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(1.4761285) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832907) q[0];
sx q[0];
rz(-1.5314449) q[0];
sx q[0];
rz(1.395306) q[0];
x q[1];
rz(0.23952837) q[2];
sx q[2];
rz(-1.0205185) q[2];
sx q[2];
rz(0.36877647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5567814) q[1];
sx q[1];
rz(-1.4540392) q[1];
sx q[1];
rz(0.054971855) q[1];
rz(-1.4011151) q[3];
sx q[3];
rz(-1.3417923) q[3];
sx q[3];
rz(1.9226216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25492111) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(-2.2754748) q[2];
rz(-2.4050889) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(2.2542663) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12164584) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(-1.2787) q[0];
rz(2.1389029) q[1];
sx q[1];
rz(-1.2607144) q[1];
sx q[1];
rz(-2.6967646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072413) q[0];
sx q[0];
rz(-1.0486516) q[0];
sx q[0];
rz(-0.17052167) q[0];
rz(0.72250267) q[2];
sx q[2];
rz(-1.2151657) q[2];
sx q[2];
rz(2.1788545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7463716) q[1];
sx q[1];
rz(-2.3276797) q[1];
sx q[1];
rz(2.9774352) q[1];
rz(-2.5505742) q[3];
sx q[3];
rz(-2.6962387) q[3];
sx q[3];
rz(1.2184008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43850809) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-2.6355696) q[2];
rz(-2.1972726) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-0.73616141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29189062) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(0.010183656) q[0];
rz(0.61095515) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(-3.0593061) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9708389) q[0];
sx q[0];
rz(-2.0505095) q[0];
sx q[0];
rz(-1.1520349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.250483) q[2];
sx q[2];
rz(-1.751305) q[2];
sx q[2];
rz(1.0918822) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.391606) q[1];
sx q[1];
rz(-0.22400236) q[1];
sx q[1];
rz(-0.27952607) q[1];
rz(-pi) q[2];
rz(-0.66466596) q[3];
sx q[3];
rz(-1.7418523) q[3];
sx q[3];
rz(-1.6534352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33132195) q[2];
sx q[2];
rz(-1.6931345) q[2];
sx q[2];
rz(0.70915478) q[2];
rz(1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(-0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3191147) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(-2.5122232) q[0];
rz(1.4156226) q[1];
sx q[1];
rz(-1.6547838) q[1];
sx q[1];
rz(-1.1782882) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.685235) q[0];
sx q[0];
rz(-2.166232) q[0];
sx q[0];
rz(2.6490982) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0731931) q[2];
sx q[2];
rz(-1.2594688) q[2];
sx q[2];
rz(-2.5466998) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7892742) q[1];
sx q[1];
rz(-1.4613918) q[1];
sx q[1];
rz(-2.3477702) q[1];
rz(-pi) q[2];
rz(-2.516516) q[3];
sx q[3];
rz(-0.82603329) q[3];
sx q[3];
rz(2.1476098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87469953) q[2];
sx q[2];
rz(-1.5723672) q[2];
sx q[2];
rz(2.3229522) q[2];
rz(1.6308174) q[3];
sx q[3];
rz(-1.8418334) q[3];
sx q[3];
rz(-2.7394133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4928987) q[0];
sx q[0];
rz(-2.0317827) q[0];
sx q[0];
rz(-1.6775525) q[0];
rz(1.875444) q[1];
sx q[1];
rz(-1.9910973) q[1];
sx q[1];
rz(1.533351) q[1];
rz(-0.89546236) q[2];
sx q[2];
rz(-1.3668899) q[2];
sx q[2];
rz(-1.5820506) q[2];
rz(0.18465445) q[3];
sx q[3];
rz(-1.6600556) q[3];
sx q[3];
rz(-0.32515812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
