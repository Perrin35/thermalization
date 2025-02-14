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
rz(0.49864545) q[0];
sx q[0];
rz(3.6874229) q[0];
sx q[0];
rz(9.3064718) q[0];
rz(-2.2220597) q[1];
sx q[1];
rz(-1.3032721) q[1];
sx q[1];
rz(1.2143171) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7987503) q[0];
sx q[0];
rz(-1.3590706) q[0];
sx q[0];
rz(2.0131074) q[0];
rz(-pi) q[1];
rz(0.75096424) q[2];
sx q[2];
rz(-0.77901024) q[2];
sx q[2];
rz(0.89391237) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8439629) q[1];
sx q[1];
rz(-1.7073892) q[1];
sx q[1];
rz(0.94874391) q[1];
rz(-pi) q[2];
rz(1.6475347) q[3];
sx q[3];
rz(-2.8676448) q[3];
sx q[3];
rz(-0.63989598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97293568) q[2];
sx q[2];
rz(-2.0106222) q[2];
sx q[2];
rz(-1.0966148) q[2];
rz(0.24122572) q[3];
sx q[3];
rz(-1.7719519) q[3];
sx q[3];
rz(0.84478861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0519401) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(-0.3845149) q[0];
rz(-2.8334726) q[1];
sx q[1];
rz(-2.6607951) q[1];
sx q[1];
rz(-2.5164129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97323167) q[0];
sx q[0];
rz(-0.28795469) q[0];
sx q[0];
rz(-3.0298427) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7169508) q[2];
sx q[2];
rz(-1.6670818) q[2];
sx q[2];
rz(-1.7801628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3846085) q[1];
sx q[1];
rz(-2.217727) q[1];
sx q[1];
rz(-1.4642843) q[1];
rz(2.1891648) q[3];
sx q[3];
rz(-1.1606207) q[3];
sx q[3];
rz(0.71632121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.008931) q[2];
sx q[2];
rz(-0.93175685) q[2];
sx q[2];
rz(0.30161944) q[2];
rz(2.610142) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6279491) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(-1.1847786) q[0];
rz(-2.9966677) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(1.5941934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14159053) q[0];
sx q[0];
rz(-0.98391082) q[0];
sx q[0];
rz(-2.1678655) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9504531) q[2];
sx q[2];
rz(-1.4002424) q[2];
sx q[2];
rz(-1.7315239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6559534) q[1];
sx q[1];
rz(-1.9536634) q[1];
sx q[1];
rz(-0.98321557) q[1];
rz(-pi) q[2];
rz(2.4800042) q[3];
sx q[3];
rz(-0.49822712) q[3];
sx q[3];
rz(3.1116886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3323815) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(-1.3261718) q[2];
rz(-0.91340804) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(-2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42295414) q[0];
sx q[0];
rz(-0.47684968) q[0];
sx q[0];
rz(1.3612716) q[0];
rz(2.4286229) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(-2.4580809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.303377) q[0];
sx q[0];
rz(-1.5808754) q[0];
sx q[0];
rz(0.99752183) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2095378) q[2];
sx q[2];
rz(-1.352013) q[2];
sx q[2];
rz(2.0444586) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5665413) q[1];
sx q[1];
rz(-1.05541) q[1];
sx q[1];
rz(2.383286) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6708507) q[3];
sx q[3];
rz(-1.81395) q[3];
sx q[3];
rz(-0.2167165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5061364) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(-2.6551841) q[2];
rz(-0.5168612) q[3];
sx q[3];
rz(-1.1812527) q[3];
sx q[3];
rz(1.0335056) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085623398) q[0];
sx q[0];
rz(-1.509868) q[0];
sx q[0];
rz(1.7686718) q[0];
rz(-1.4068475) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(-2.6645606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1005362) q[0];
sx q[0];
rz(-1.1541799) q[0];
sx q[0];
rz(-0.21505298) q[0];
rz(-pi) q[1];
rz(0.3933752) q[2];
sx q[2];
rz(-2.057339) q[2];
sx q[2];
rz(1.3855307) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.841694) q[1];
sx q[1];
rz(-1.8704789) q[1];
sx q[1];
rz(1.0503286) q[1];
rz(-pi) q[2];
rz(2.414235) q[3];
sx q[3];
rz(-2.2148956) q[3];
sx q[3];
rz(-1.7880881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0341952) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(-2.6743215) q[2];
rz(-0.16921903) q[3];
sx q[3];
rz(-1.4110112) q[3];
sx q[3];
rz(0.30538487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071103) q[0];
sx q[0];
rz(-2.2619673) q[0];
sx q[0];
rz(0.18950263) q[0];
rz(-2.8684008) q[1];
sx q[1];
rz(-1.0739645) q[1];
sx q[1];
rz(-0.80064076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36291262) q[0];
sx q[0];
rz(-1.3375459) q[0];
sx q[0];
rz(-3.0972079) q[0];
x q[1];
rz(-2.8932299) q[2];
sx q[2];
rz(-1.6585322) q[2];
sx q[2];
rz(-1.1540292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86294293) q[1];
sx q[1];
rz(-0.23506308) q[1];
sx q[1];
rz(2.1549015) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74108776) q[3];
sx q[3];
rz(-0.68276486) q[3];
sx q[3];
rz(0.53214754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2022986) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(0.536971) q[2];
rz(-1.8592853) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.110431) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.8102616) q[0];
rz(1.7000465) q[1];
sx q[1];
rz(-1.6621637) q[1];
sx q[1];
rz(-1.9190681) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13200296) q[0];
sx q[0];
rz(-1.1545657) q[0];
sx q[0];
rz(-0.25968174) q[0];
rz(3.020633) q[2];
sx q[2];
rz(-2.5865062) q[2];
sx q[2];
rz(2.4400389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68798399) q[1];
sx q[1];
rz(-2.2609432) q[1];
sx q[1];
rz(-0.28120561) q[1];
rz(-1.8172713) q[3];
sx q[3];
rz(-0.45650864) q[3];
sx q[3];
rz(0.64303165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20087251) q[2];
sx q[2];
rz(-2.6437289) q[2];
sx q[2];
rz(-2.8087924) q[2];
rz(0.29916549) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(-0.021473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0126295) q[0];
sx q[0];
rz(-2.3733932) q[0];
sx q[0];
rz(2.8666038) q[0];
rz(3.0601652) q[1];
sx q[1];
rz(-2.3402201) q[1];
sx q[1];
rz(-1.8191232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0596432) q[0];
sx q[0];
rz(-2.9026051) q[0];
sx q[0];
rz(1.3168524) q[0];
x q[1];
rz(2.7733388) q[2];
sx q[2];
rz(-2.2241631) q[2];
sx q[2];
rz(0.39619941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9857707) q[1];
sx q[1];
rz(-0.81822526) q[1];
sx q[1];
rz(2.6236412) q[1];
rz(-pi) q[2];
rz(2.2127719) q[3];
sx q[3];
rz(-0.95212338) q[3];
sx q[3];
rz(-2.259038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4203809) q[2];
sx q[2];
rz(-0.12004852) q[2];
sx q[2];
rz(1.4018641) q[2];
rz(1.2841355) q[3];
sx q[3];
rz(-1.9522342) q[3];
sx q[3];
rz(0.0349667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5818802) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(0.48932073) q[0];
rz(3.1210461) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(0.70125854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4160181) q[0];
sx q[0];
rz(-1.3301992) q[0];
sx q[0];
rz(-0.98576905) q[0];
rz(-pi) q[1];
x q[1];
rz(0.047766165) q[2];
sx q[2];
rz(-0.86776185) q[2];
sx q[2];
rz(-0.060279708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32966954) q[1];
sx q[1];
rz(-2.3282302) q[1];
sx q[1];
rz(1.2438891) q[1];
rz(-2.4550426) q[3];
sx q[3];
rz(-1.1995763) q[3];
sx q[3];
rz(-2.8438501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4150841) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(-2.6611967) q[2];
rz(1.4091617) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3163863) q[0];
sx q[0];
rz(-1.2403064) q[0];
sx q[0];
rz(2.7301042) q[0];
rz(-1.8972137) q[1];
sx q[1];
rz(-1.139541) q[1];
sx q[1];
rz(2.8172475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.024687) q[0];
sx q[0];
rz(-1.9431264) q[0];
sx q[0];
rz(0.056554746) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3408808) q[2];
sx q[2];
rz(-2.6691648) q[2];
sx q[2];
rz(-0.97324449) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0592546) q[1];
sx q[1];
rz(-1.6538701) q[1];
sx q[1];
rz(-2.1359813) q[1];
rz(-pi) q[2];
rz(-1.5839229) q[3];
sx q[3];
rz(-1.4495951) q[3];
sx q[3];
rz(-1.2016202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1856498) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-0.11697098) q[2];
rz(2.1864435) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(-2.5999033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87880001) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(-1.1204489) q[1];
sx q[1];
rz(-2.3255377) q[1];
sx q[1];
rz(-1.9768523) q[1];
rz(-2.929959) q[2];
sx q[2];
rz(-1.0734049) q[2];
sx q[2];
rz(2.887101) q[2];
rz(1.4799812) q[3];
sx q[3];
rz(-2.869893) q[3];
sx q[3];
rz(0.13520959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
