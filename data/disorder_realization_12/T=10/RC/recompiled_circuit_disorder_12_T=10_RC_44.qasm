OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(1.4847423) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(-2.2533921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9125497) q[0];
sx q[0];
rz(-2.9648844) q[0];
sx q[0];
rz(-0.25632174) q[0];
rz(-pi) q[1];
rz(2.1076803) q[2];
sx q[2];
rz(-1.4281338) q[2];
sx q[2];
rz(0.85096525) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2759309) q[1];
sx q[1];
rz(-1.4038424) q[1];
sx q[1];
rz(1.5702412) q[1];
x q[2];
rz(2.4429951) q[3];
sx q[3];
rz(-1.5961831) q[3];
sx q[3];
rz(0.7198402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-2.0236012) q[2];
rz(2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730826) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(-1.2163935) q[1];
sx q[1];
rz(-1.1647859) q[1];
sx q[1];
rz(-3.0156946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038371041) q[0];
sx q[0];
rz(-1.9403337) q[0];
sx q[0];
rz(-1.0612556) q[0];
rz(0.89181487) q[2];
sx q[2];
rz(-1.0620772) q[2];
sx q[2];
rz(1.5430792) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1344578) q[1];
sx q[1];
rz(-1.6975228) q[1];
sx q[1];
rz(0.030220672) q[1];
rz(-0.010766518) q[3];
sx q[3];
rz(-2.2942703) q[3];
sx q[3];
rz(-0.41124287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(3.0575867) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(-1.9817339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1201598) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(1.73818) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50767501) q[2];
sx q[2];
rz(-2.1360364) q[2];
sx q[2];
rz(2.450168) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.68418903) q[1];
sx q[1];
rz(-1.1154798) q[1];
sx q[1];
rz(1.9940358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0252762) q[3];
sx q[3];
rz(-1.5226411) q[3];
sx q[3];
rz(0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40538654) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(-0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9765587) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(-0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(-0.28465095) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4581504) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(3.0577858) q[0];
rz(-0.95732032) q[2];
sx q[2];
rz(-0.34733221) q[2];
sx q[2];
rz(0.50328244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33672562) q[1];
sx q[1];
rz(-1.6972099) q[1];
sx q[1];
rz(0.079041914) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.986226) q[3];
sx q[3];
rz(-0.75891906) q[3];
sx q[3];
rz(-0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9099137) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(2.5740734) q[2];
rz(2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-1.8150785) q[0];
rz(1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(0.93793905) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4834578) q[0];
sx q[0];
rz(-1.5489274) q[0];
sx q[0];
rz(1.4786426) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51860923) q[2];
sx q[2];
rz(-1.0728288) q[2];
sx q[2];
rz(-1.2155611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93463072) q[1];
sx q[1];
rz(-0.043255581) q[1];
sx q[1];
rz(-2.2224777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(2.6854533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(-2.5842216) q[0];
rz(-2.5769261) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-2.5851137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508931) q[0];
sx q[0];
rz(-1.3644344) q[0];
sx q[0];
rz(-0.38462374) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7700023) q[2];
sx q[2];
rz(-2.4161985) q[2];
sx q[2];
rz(-2.9273916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6528875) q[1];
sx q[1];
rz(-2.8294551) q[1];
sx q[1];
rz(-2.1991792) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5197093) q[3];
sx q[3];
rz(-2.3458977) q[3];
sx q[3];
rz(1.6076128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(-1.5103643) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0790134) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(0.50061217) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635054) q[0];
sx q[0];
rz(-0.18112077) q[0];
sx q[0];
rz(1.7630793) q[0];
rz(-pi) q[1];
x q[1];
rz(2.826564) q[2];
sx q[2];
rz(-2.0909967) q[2];
sx q[2];
rz(2.5556285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.171247) q[1];
sx q[1];
rz(-1.7016181) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.385419) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(-2.365436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5381955) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(0.44089857) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(-0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(1.4656461) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059793652) q[0];
sx q[0];
rz(-1.1104715) q[0];
sx q[0];
rz(0.1501118) q[0];
x q[1];
rz(2.9550214) q[2];
sx q[2];
rz(-1.6779043) q[2];
sx q[2];
rz(0.69498108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88361909) q[1];
sx q[1];
rz(-2.5273364) q[1];
sx q[1];
rz(2.1780464) q[1];
x q[2];
rz(2.500324) q[3];
sx q[3];
rz(-2.0740168) q[3];
sx q[3];
rz(0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(1.476293) q[2];
rz(-0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(-0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575689) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-2.2163056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46426526) q[0];
sx q[0];
rz(-2.263501) q[0];
sx q[0];
rz(-2.4753184) q[0];
x q[1];
rz(2.4783752) q[2];
sx q[2];
rz(-1.466201) q[2];
sx q[2];
rz(-1.6246206) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9097594) q[1];
sx q[1];
rz(-1.4405466) q[1];
sx q[1];
rz(-2.9568683) q[1];
x q[2];
rz(-2.7221189) q[3];
sx q[3];
rz(-2.2443716) q[3];
sx q[3];
rz(-0.89299612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(2.1378689) q[2];
rz(3.051493) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039624778) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(2.3840747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19125464) q[0];
sx q[0];
rz(-1.7636824) q[0];
sx q[0];
rz(0.14001503) q[0];
rz(-pi) q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(-2.5622501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6706558) q[1];
sx q[1];
rz(-2.6568012) q[1];
sx q[1];
rz(2.0623341) q[1];
rz(-pi) q[2];
rz(-2.431589) q[3];
sx q[3];
rz(-1.5922976) q[3];
sx q[3];
rz(1.8332675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6486075) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(-1.6408625) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(0.63207788) q[2];
sx q[2];
rz(-0.52543228) q[2];
sx q[2];
rz(0.57731522) q[2];
rz(0.32140857) q[3];
sx q[3];
rz(-1.0151498) q[3];
sx q[3];
rz(-2.2890454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];