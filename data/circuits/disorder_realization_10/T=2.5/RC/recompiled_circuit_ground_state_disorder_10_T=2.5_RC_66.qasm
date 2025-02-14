OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(2.7326512) q[0];
sx q[0];
rz(10.41806) q[0];
rz(-3.2886796) q[1];
sx q[1];
rz(4.1380652) q[1];
sx q[1];
rz(11.148718) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.111747) q[0];
sx q[0];
rz(-0.8518712) q[0];
sx q[0];
rz(1.0875888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70320798) q[2];
sx q[2];
rz(-2.0173912) q[2];
sx q[2];
rz(0.27189068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1163535) q[1];
sx q[1];
rz(-2.1903746) q[1];
sx q[1];
rz(1.2337633) q[1];
rz(-pi) q[2];
rz(-1.8090114) q[3];
sx q[3];
rz(-0.138126) q[3];
sx q[3];
rz(-0.81616831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19109569) q[2];
sx q[2];
rz(-1.82093) q[2];
sx q[2];
rz(2.0868059) q[2];
rz(2.6705006) q[3];
sx q[3];
rz(-1.4548917) q[3];
sx q[3];
rz(2.2733222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489814) q[0];
sx q[0];
rz(-1.6845717) q[0];
sx q[0];
rz(-0.23496041) q[0];
rz(-1.8670392) q[1];
sx q[1];
rz(-1.3095368) q[1];
sx q[1];
rz(2.4539006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0339081) q[0];
sx q[0];
rz(-1.8982197) q[0];
sx q[0];
rz(1.4998687) q[0];
rz(-pi) q[1];
rz(-2.657452) q[2];
sx q[2];
rz(-2.2233367) q[2];
sx q[2];
rz(2.2732041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60581698) q[1];
sx q[1];
rz(-1.5763842) q[1];
sx q[1];
rz(1.5833012) q[1];
rz(-0.092472381) q[3];
sx q[3];
rz(-2.1824129) q[3];
sx q[3];
rz(-3.0870952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.71622744) q[2];
sx q[2];
rz(-1.3024412) q[2];
sx q[2];
rz(-2.9322374) q[2];
rz(0.9730722) q[3];
sx q[3];
rz(-0.89997411) q[3];
sx q[3];
rz(-0.59276855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8415602) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(-2.9969065) q[0];
rz(-1.2380098) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(0.091781052) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9301355) q[0];
sx q[0];
rz(-1.2263023) q[0];
sx q[0];
rz(-0.99665595) q[0];
rz(-1.4745643) q[2];
sx q[2];
rz(-1.6394168) q[2];
sx q[2];
rz(0.13185908) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8615177) q[1];
sx q[1];
rz(-2.8535378) q[1];
sx q[1];
rz(0.83863284) q[1];
rz(-pi) q[2];
rz(2.0030267) q[3];
sx q[3];
rz(-1.5097268) q[3];
sx q[3];
rz(1.40738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38516513) q[2];
sx q[2];
rz(-1.5084718) q[2];
sx q[2];
rz(0.37919322) q[2];
rz(-2.3181629) q[3];
sx q[3];
rz(-0.40612602) q[3];
sx q[3];
rz(0.2894952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8150811) q[0];
sx q[0];
rz(-2.264475) q[0];
sx q[0];
rz(0.99863482) q[0];
rz(-0.48209349) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(1.3179717) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4559193) q[0];
sx q[0];
rz(-0.64252526) q[0];
sx q[0];
rz(0.57149096) q[0];
rz(-0.58546328) q[2];
sx q[2];
rz(-1.778086) q[2];
sx q[2];
rz(2.6475865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1554679) q[1];
sx q[1];
rz(-1.7805702) q[1];
sx q[1];
rz(-0.54297691) q[1];
rz(-pi) q[2];
rz(-0.79341268) q[3];
sx q[3];
rz(-1.1856836) q[3];
sx q[3];
rz(-2.957893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2202997) q[2];
sx q[2];
rz(-2.2396294) q[2];
sx q[2];
rz(-2.656142) q[2];
rz(-0.38393936) q[3];
sx q[3];
rz(-1.5860312) q[3];
sx q[3];
rz(-0.78401047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601198) q[0];
sx q[0];
rz(-0.10426846) q[0];
sx q[0];
rz(-0.01165788) q[0];
rz(0.13721379) q[1];
sx q[1];
rz(-1.5533181) q[1];
sx q[1];
rz(-1.3020017) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6377651) q[0];
sx q[0];
rz(-1.3877467) q[0];
sx q[0];
rz(1.6812912) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2574591) q[2];
sx q[2];
rz(-2.4120286) q[2];
sx q[2];
rz(1.8232249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3350388) q[1];
sx q[1];
rz(-1.8825899) q[1];
sx q[1];
rz(-3.0278518) q[1];
x q[2];
rz(-0.37837677) q[3];
sx q[3];
rz(-1.9073745) q[3];
sx q[3];
rz(-0.61016309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7234708) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(2.8483025) q[2];
rz(0.063118525) q[3];
sx q[3];
rz(-1.2226353) q[3];
sx q[3];
rz(0.89490923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1191331) q[0];
sx q[0];
rz(-2.8526511) q[0];
sx q[0];
rz(-0.038507842) q[0];
rz(1.0844082) q[1];
sx q[1];
rz(-2.7190828) q[1];
sx q[1];
rz(-2.1748621) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9733676) q[0];
sx q[0];
rz(-1.7666139) q[0];
sx q[0];
rz(1.7652579) q[0];
rz(1.837908) q[2];
sx q[2];
rz(-1.5582286) q[2];
sx q[2];
rz(1.723701) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37361426) q[1];
sx q[1];
rz(-1.9100338) q[1];
sx q[1];
rz(1.148073) q[1];
rz(-pi) q[2];
rz(1.2472263) q[3];
sx q[3];
rz(-0.67646356) q[3];
sx q[3];
rz(2.885779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66095573) q[2];
sx q[2];
rz(-0.2350685) q[2];
sx q[2];
rz(2.1601775) q[2];
rz(1.6437982) q[3];
sx q[3];
rz(-1.2621597) q[3];
sx q[3];
rz(1.5952544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82866955) q[0];
sx q[0];
rz(-2.4607615) q[0];
sx q[0];
rz(2.8269826) q[0];
rz(-0.18868748) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(2.6729118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22671902) q[0];
sx q[0];
rz(-1.3283936) q[0];
sx q[0];
rz(-1.7049432) q[0];
rz(1.738564) q[2];
sx q[2];
rz(-1.5689384) q[2];
sx q[2];
rz(2.6986966) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6275121) q[1];
sx q[1];
rz(-2.2805164) q[1];
sx q[1];
rz(-0.40436097) q[1];
rz(-0.58705892) q[3];
sx q[3];
rz(-0.54577845) q[3];
sx q[3];
rz(-3.0445638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20588747) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(-2.1708798) q[2];
rz(-0.5361706) q[3];
sx q[3];
rz(-1.2090809) q[3];
sx q[3];
rz(0.71603388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69990528) q[0];
sx q[0];
rz(-2.5869885) q[0];
sx q[0];
rz(1.4458789) q[0];
rz(1.0384167) q[1];
sx q[1];
rz(-1.1035792) q[1];
sx q[1];
rz(0.081092484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2807407) q[0];
sx q[0];
rz(-0.74374712) q[0];
sx q[0];
rz(-2.70983) q[0];
rz(-pi) q[1];
rz(-2.7620188) q[2];
sx q[2];
rz(-2.4502769) q[2];
sx q[2];
rz(-2.8065681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3861919) q[1];
sx q[1];
rz(-0.31495783) q[1];
sx q[1];
rz(-1.9506728) q[1];
rz(-pi) q[2];
rz(-1.1015973) q[3];
sx q[3];
rz(-1.0476255) q[3];
sx q[3];
rz(-1.8947471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1991835) q[2];
sx q[2];
rz(-0.77745357) q[2];
sx q[2];
rz(-0.24277631) q[2];
rz(1.5593922) q[3];
sx q[3];
rz(-0.42583164) q[3];
sx q[3];
rz(-1.2529713) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9026069) q[0];
sx q[0];
rz(-1.0317529) q[0];
sx q[0];
rz(1.2549988) q[0];
rz(1.7651419) q[1];
sx q[1];
rz(-2.1170728) q[1];
sx q[1];
rz(0.66744101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4723929) q[0];
sx q[0];
rz(-2.071864) q[0];
sx q[0];
rz(-2.7612812) q[0];
rz(-pi) q[1];
rz(1.5282643) q[2];
sx q[2];
rz(-2.0669524) q[2];
sx q[2];
rz(-2.4537078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.094634) q[1];
sx q[1];
rz(-2.486645) q[1];
sx q[1];
rz(1.0818693) q[1];
rz(0.38102229) q[3];
sx q[3];
rz(-1.8748611) q[3];
sx q[3];
rz(-1.1883433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61975512) q[2];
sx q[2];
rz(-0.73885584) q[2];
sx q[2];
rz(0.45836207) q[2];
rz(1.6058263) q[3];
sx q[3];
rz(-1.0777487) q[3];
sx q[3];
rz(2.2569236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89852029) q[0];
sx q[0];
rz(-0.5721108) q[0];
sx q[0];
rz(2.7761053) q[0];
rz(-2.1416523) q[1];
sx q[1];
rz(-1.4391856) q[1];
sx q[1];
rz(-1.4568636) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92577584) q[0];
sx q[0];
rz(-2.121637) q[0];
sx q[0];
rz(-2.9047545) q[0];
x q[1];
rz(3.0295154) q[2];
sx q[2];
rz(-1.3044895) q[2];
sx q[2];
rz(0.2923686) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7813606) q[1];
sx q[1];
rz(-1.7443027) q[1];
sx q[1];
rz(2.3886015) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0356009) q[3];
sx q[3];
rz(-0.92886954) q[3];
sx q[3];
rz(-1.2479051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0008056) q[2];
sx q[2];
rz(-1.3212589) q[2];
sx q[2];
rz(-1.0055536) q[2];
rz(-0.65070659) q[3];
sx q[3];
rz(-1.7020099) q[3];
sx q[3];
rz(0.85056359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093813048) q[0];
sx q[0];
rz(-0.78421264) q[0];
sx q[0];
rz(-1.0017851) q[0];
rz(1.7306937) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(3.0244556) q[2];
sx q[2];
rz(-1.2436554) q[2];
sx q[2];
rz(0.40851442) q[2];
rz(-3.0372314) q[3];
sx q[3];
rz(-1.6634533) q[3];
sx q[3];
rz(-2.8692393) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
