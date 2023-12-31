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
rz(-1.6568503) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(0.88820052) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0892544) q[0];
sx q[0];
rz(-1.5262145) q[0];
sx q[0];
rz(-2.9705439) q[0];
x q[1];
rz(-2.9759679) q[2];
sx q[2];
rz(-2.1016444) q[2];
sx q[2];
rz(-0.80425516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8623212) q[1];
sx q[1];
rz(-0.16695484) q[1];
sx q[1];
rz(3.1382986) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6039443) q[3];
sx q[3];
rz(-2.2691233) q[3];
sx q[3];
rz(2.2693199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.858294) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(2.0236012) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730826) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(0.65482393) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-0.12589802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1032216) q[0];
sx q[0];
rz(-1.9403337) q[0];
sx q[0];
rz(-1.0612556) q[0];
x q[1];
rz(-0.84463859) q[2];
sx q[2];
rz(-0.8234878) q[2];
sx q[2];
rz(0.57074947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44015917) q[1];
sx q[1];
rz(-1.5408181) q[1];
sx q[1];
rz(-1.6975801) q[1];
rz(-0.010766518) q[3];
sx q[3];
rz(-2.2942703) q[3];
sx q[3];
rz(2.7303498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(-0.38802567) q[2];
rz(1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(0.079332381) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(-1.1598587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1201598) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(1.4034127) q[0];
rz(-2.6339176) q[2];
sx q[2];
rz(-2.1360364) q[2];
sx q[2];
rz(2.450168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11324727) q[1];
sx q[1];
rz(-2.5302437) q[1];
sx q[1];
rz(0.69797413) q[1];
rz(2.1163164) q[3];
sx q[3];
rz(-1.6189515) q[3];
sx q[3];
rz(0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40538654) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(3.0333056) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9765587) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(-2.8569417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2157001) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(2.0490993) q[0];
x q[1];
rz(0.20547159) q[2];
sx q[2];
rz(-1.2887508) q[2];
sx q[2];
rz(-3.0021283) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22420158) q[1];
sx q[1];
rz(-2.9926139) q[1];
sx q[1];
rz(-1.0148744) q[1];
rz(-pi) q[2];
rz(-2.1553667) q[3];
sx q[3];
rz(-2.3826736) q[3];
sx q[3];
rz(-0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(2.5740734) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(1.1119941) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.652997) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(-1.3265142) q[0];
rz(1.9891706) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(0.93793905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468231) q[0];
sx q[0];
rz(-1.4786647) q[0];
sx q[0];
rz(-0.021962086) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51860923) q[2];
sx q[2];
rz(-2.0687639) q[2];
sx q[2];
rz(1.2155611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93463072) q[1];
sx q[1];
rz(-0.043255581) q[1];
sx q[1];
rz(-0.91911493) q[1];
x q[2];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.218759) q[3];
sx q[3];
rz(-2.6854533) q[3];
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
rz(-1.3025618) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(-1.8235824) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0628478) q[0];
sx q[0];
rz(-1.1947462) q[0];
sx q[0];
rz(1.7929121) q[0];
x q[1];
rz(-1.7700023) q[2];
sx q[2];
rz(-2.4161985) q[2];
sx q[2];
rz(-2.9273916) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1408128) q[1];
sx q[1];
rz(-1.3197348) q[1];
sx q[1];
rz(-2.954133) q[1];
rz(-pi) q[2];
rz(-2.3658386) q[3];
sx q[3];
rz(-1.5343101) q[3];
sx q[3];
rz(0.072582399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.640655) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(2.6409805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4780873) q[0];
sx q[0];
rz(-2.9604719) q[0];
sx q[0];
rz(1.3785133) q[0];
x q[1];
rz(-1.074965) q[2];
sx q[2];
rz(-0.60056409) q[2];
sx q[2];
rz(-1.9753089) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2523633) q[1];
sx q[1];
rz(-1.4503308) q[1];
sx q[1];
rz(-2.7389588) q[1];
rz(-pi) q[2];
rz(-1.9686437) q[3];
sx q[3];
rz(-2.2865191) q[3];
sx q[3];
rz(-2.0778823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(-0.44089857) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(-0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-1.0674397) q[0];
rz(2.7087129) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.4656461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.081799) q[0];
sx q[0];
rz(-2.0311211) q[0];
sx q[0];
rz(-0.1501118) q[0];
x q[1];
rz(1.6797811) q[2];
sx q[2];
rz(-1.7562859) q[2];
sx q[2];
rz(2.2456004) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17897478) q[1];
sx q[1];
rz(-2.0638421) q[1];
sx q[1];
rz(0.38260539) q[1];
rz(-0.74387868) q[3];
sx q[3];
rz(-2.3489967) q[3];
sx q[3];
rz(2.511123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33891588) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(1.476293) q[2];
rz(2.2680797) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(-0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2575689) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(-0.73721686) q[0];
rz(0.018741477) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-0.92528701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46426526) q[0];
sx q[0];
rz(-0.87809169) q[0];
sx q[0];
rz(-2.4753184) q[0];
x q[1];
rz(-0.16889062) q[2];
sx q[2];
rz(-2.4714111) q[2];
sx q[2];
rz(0.18667135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9463897) q[1];
sx q[1];
rz(-2.9159947) q[1];
sx q[1];
rz(-0.619508) q[1];
rz(-0.85261811) q[3];
sx q[3];
rz(-1.8947621) q[3];
sx q[3];
rz(2.1924803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(-2.3840747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82523358) q[0];
sx q[0];
rz(-0.23783437) q[0];
sx q[0];
rz(0.95038484) q[0];
rz(-2.9847758) q[2];
sx q[2];
rz(-0.98625253) q[2];
sx q[2];
rz(2.2371694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4841008) q[1];
sx q[1];
rz(-1.7925646) q[1];
sx q[1];
rz(-2.0055254) q[1];
rz(2.431589) q[3];
sx q[3];
rz(-1.5492951) q[3];
sx q[3];
rz(-1.3083252) q[3];
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
rz(-0.63888597) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(-2.1209774) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6486075) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(-1.5007301) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-1.2407672) q[2];
sx q[2];
rz(-1.9874265) q[2];
sx q[2];
rz(-0.12513587) q[2];
rz(2.8201841) q[3];
sx q[3];
rz(-2.1264429) q[3];
sx q[3];
rz(0.85254729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
