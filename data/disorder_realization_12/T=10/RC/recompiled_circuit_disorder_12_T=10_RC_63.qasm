OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(0.88820052) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22904299) q[0];
sx q[0];
rz(-2.9648844) q[0];
sx q[0];
rz(-2.8852709) q[0];
x q[1];
rz(1.8445831) q[2];
sx q[2];
rz(-2.5878776) q[2];
sx q[2];
rz(0.48534976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2759309) q[1];
sx q[1];
rz(-1.4038424) q[1];
sx q[1];
rz(1.5702412) q[1];
rz(-pi) q[2];
rz(-1.5376484) q[3];
sx q[3];
rz(-2.2691233) q[3];
sx q[3];
rz(2.2693199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(-0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(3.0156946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1032216) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(-1.0612556) q[0];
rz(-pi) q[1];
rz(0.6217896) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(0.34678005) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1344578) q[1];
sx q[1];
rz(-1.4440698) q[1];
sx q[1];
rz(0.030220672) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010766518) q[3];
sx q[3];
rz(-0.84732238) q[3];
sx q[3];
rz(-0.41124287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(-0.38802567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.1598587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0214329) q[0];
sx q[0];
rz(-0.5172356) q[0];
sx q[0];
rz(-1.73818) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94295393) q[2];
sx q[2];
rz(-1.1477594) q[2];
sx q[2];
rz(1.9726276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4574036) q[1];
sx q[1];
rz(-1.1154798) q[1];
sx q[1];
rz(1.1475569) q[1];
rz(-pi) q[2];
x q[2];
rz(0.056315259) q[3];
sx q[3];
rz(-1.0259797) q[3];
sx q[3];
rz(-0.69308263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(-2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9765587) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(1.6595586) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(-2.8569417) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92589256) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(2.0490993) q[0];
rz(-pi) q[1];
rz(2.9361211) q[2];
sx q[2];
rz(-1.8528419) q[2];
sx q[2];
rz(0.13946433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22420158) q[1];
sx q[1];
rz(-0.14897878) q[1];
sx q[1];
rz(-1.0148744) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1553667) q[3];
sx q[3];
rz(-0.75891906) q[3];
sx q[3];
rz(2.330247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.231679) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652997) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(1.9891706) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(-2.2036536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4612761) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(1.3374431) q[0];
rz(-pi) q[1];
rz(-2.3100501) q[2];
sx q[2];
rz(-0.70280308) q[2];
sx q[2];
rz(1.0520832) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93463072) q[1];
sx q[1];
rz(-0.043255581) q[1];
sx q[1];
rz(2.2224777) q[1];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.218759) q[3];
sx q[3];
rz(0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1753297) q[2];
sx q[2];
rz(-2.9523409) q[2];
sx q[2];
rz(-2.1002634) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(-0.56466651) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(0.55647892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51150409) q[0];
sx q[0];
rz(-0.43404365) q[0];
sx q[0];
rz(-0.50891288) q[0];
rz(-0.85530497) q[2];
sx q[2];
rz(-1.4391293) q[2];
sx q[2];
rz(1.635074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6186657) q[1];
sx q[1];
rz(-1.3892801) q[1];
sx q[1];
rz(-1.8261441) q[1];
rz(-pi) q[2];
rz(2.3658386) q[3];
sx q[3];
rz(-1.5343101) q[3];
sx q[3];
rz(3.0690103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-3.0986837) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(0.50061217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0450889) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(-1.7486497) q[0];
rz(-pi) q[1];
rz(2.826564) q[2];
sx q[2];
rz(-2.0909967) q[2];
sx q[2];
rz(-0.5859642) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40655207) q[1];
sx q[1];
rz(-0.41932377) q[1];
sx q[1];
rz(-0.29962824) q[1];
rz(-pi) q[2];
rz(-1.9686437) q[3];
sx q[3];
rz(-0.85507353) q[3];
sx q[3];
rz(2.0778823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5381955) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-2.4659757) q[2];
rz(-2.7006941) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.0764517) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.4656461) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.081799) q[0];
sx q[0];
rz(-1.1104715) q[0];
sx q[0];
rz(2.9914809) q[0];
rz(1.4618116) q[2];
sx q[2];
rz(-1.3853067) q[2];
sx q[2];
rz(-0.89599228) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9626179) q[1];
sx q[1];
rz(-2.0638421) q[1];
sx q[1];
rz(-0.38260539) q[1];
x q[2];
rz(-0.74387868) q[3];
sx q[3];
rz(-0.79259593) q[3];
sx q[3];
rz(0.63046968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(1.6652997) q[2];
rz(-2.2680797) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(0.73721686) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(2.2163056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6773274) q[0];
sx q[0];
rz(-0.87809169) q[0];
sx q[0];
rz(-0.6662743) q[0];
rz(-2.972702) q[2];
sx q[2];
rz(-0.67018159) q[2];
sx q[2];
rz(0.18667135) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8268938) q[1];
sx q[1];
rz(-1.7539382) q[1];
sx q[1];
rz(1.4383184) q[1];
x q[2];
rz(0.85261811) q[3];
sx q[3];
rz(-1.8947621) q[3];
sx q[3];
rz(0.94911239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.2596624) q[0];
rz(-2.8758077) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(0.75751799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.950338) q[0];
sx q[0];
rz(-1.7636824) q[0];
sx q[0];
rz(0.14001503) q[0];
x q[1];
rz(-0.9805571) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(-0.57934258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1265035) q[1];
sx q[1];
rz(-1.9941829) q[1];
sx q[1];
rz(-2.8979315) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71000368) q[3];
sx q[3];
rz(-1.5922976) q[3];
sx q[3];
rz(-1.3083252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1511128) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-0.63888597) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4929852) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(-1.5007301) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-0.63207788) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(-2.0414447) q[3];
sx q[3];
rz(-0.63334076) q[3];
sx q[3];
rz(0.28950194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
