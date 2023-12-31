OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(-0.43963471) q[0];
sx q[0];
rz(3.0602732) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62277943) q[0];
sx q[0];
rz(-2.8840027) q[0];
sx q[0];
rz(1.326484) q[0];
rz(-pi) q[1];
rz(-1.2230258) q[2];
sx q[2];
rz(-2.7070621) q[2];
sx q[2];
rz(-0.44473106) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1374955) q[1];
sx q[1];
rz(-1.6988519) q[1];
sx q[1];
rz(2.0069684) q[1];
rz(2.8738535) q[3];
sx q[3];
rz(-2.8149238) q[3];
sx q[3];
rz(0.79145811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2540934) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(0.19533531) q[0];
rz(-2.7665566) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(2.9017752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46562425) q[0];
sx q[0];
rz(-1.3249319) q[0];
sx q[0];
rz(-1.4093536) q[0];
rz(-pi) q[1];
rz(0.2113091) q[2];
sx q[2];
rz(-1.276187) q[2];
sx q[2];
rz(-2.353637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1344192) q[1];
sx q[1];
rz(-1.5181784) q[1];
sx q[1];
rz(-2.3514071) q[1];
x q[2];
rz(-2.2858743) q[3];
sx q[3];
rz(-0.021589605) q[3];
sx q[3];
rz(1.3947226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6372765) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(-1.8117388) q[2];
rz(-1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(-0.66816107) q[0];
rz(1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(-1.0659165) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33879963) q[0];
sx q[0];
rz(-0.25164139) q[0];
sx q[0];
rz(1.7358857) q[0];
rz(-1.4270646) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(1.8260173) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8733858) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
x q[2];
rz(0.12001868) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(-1.1413871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.187591) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-0.032547396) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.4867841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.810881) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(2.7964554) q[0];
x q[1];
rz(-0.60947946) q[2];
sx q[2];
rz(-1.8595427) q[2];
sx q[2];
rz(-0.8537054) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50258892) q[1];
sx q[1];
rz(-1.2856312) q[1];
sx q[1];
rz(1.4147948) q[1];
rz(-pi) q[2];
rz(1.5781457) q[3];
sx q[3];
rz(-0.74439936) q[3];
sx q[3];
rz(-0.053754427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(2.0452943) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0356045) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.460357) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(0.016074093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2747297) q[0];
sx q[0];
rz(-0.99039536) q[0];
sx q[0];
rz(-0.6443364) q[0];
x q[1];
rz(1.5889421) q[2];
sx q[2];
rz(-0.43076736) q[2];
sx q[2];
rz(-0.22242966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88735089) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(0.90708797) q[1];
rz(1.5752951) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(-0.11210657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(0.62189046) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(1.2247359) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82299267) q[0];
sx q[0];
rz(-2.4791414) q[0];
sx q[0];
rz(-1.2160765) q[0];
rz(2.9333026) q[2];
sx q[2];
rz(-1.1366833) q[2];
sx q[2];
rz(-1.8722033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9310301) q[1];
sx q[1];
rz(-1.6258874) q[1];
sx q[1];
rz(2.1767666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.46035) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(-0.16050592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(0.58037037) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(2.4408128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66877767) q[0];
sx q[0];
rz(-0.46572177) q[0];
sx q[0];
rz(-1.3186243) q[0];
rz(-pi) q[1];
rz(-1.6875661) q[2];
sx q[2];
rz(-1.8883369) q[2];
sx q[2];
rz(1.5460528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80220561) q[1];
sx q[1];
rz(-1.4805668) q[1];
sx q[1];
rz(-2.1369364) q[1];
rz(-pi) q[2];
rz(-3.0573781) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.7250852) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(-0.8917121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72043334) q[0];
sx q[0];
rz(-1.7283068) q[0];
sx q[0];
rz(0.79173761) q[0];
rz(-2.0217998) q[2];
sx q[2];
rz(-2.1037256) q[2];
sx q[2];
rz(1.230513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8430431) q[1];
sx q[1];
rz(-1.9798568) q[1];
sx q[1];
rz(-0.72960735) q[1];
x q[2];
rz(1.9705087) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(-1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(-2.8083535) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.3649712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38395912) q[0];
sx q[0];
rz(-2.3968292) q[0];
sx q[0];
rz(-0.061116771) q[0];
rz(1.981609) q[2];
sx q[2];
rz(-1.4940726) q[2];
sx q[2];
rz(3.0551747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2862128) q[1];
sx q[1];
rz(-2.3533822) q[1];
sx q[1];
rz(2.6469995) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19946675) q[3];
sx q[3];
rz(-1.27928) q[3];
sx q[3];
rz(-2.9150972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(0.96819425) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4614842) q[0];
sx q[0];
rz(-0.52163863) q[0];
sx q[0];
rz(2.3011544) q[0];
rz(-pi) q[1];
rz(-2.7231611) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(0.07428169) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26607516) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(-2.7060899) q[1];
rz(-1.8929385) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-0.9957046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(-2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-0.33404074) q[2];
sx q[2];
rz(-0.77790778) q[2];
sx q[2];
rz(1.6890656) q[2];
rz(1.0352186) q[3];
sx q[3];
rz(-2.513701) q[3];
sx q[3];
rz(2.7660478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
