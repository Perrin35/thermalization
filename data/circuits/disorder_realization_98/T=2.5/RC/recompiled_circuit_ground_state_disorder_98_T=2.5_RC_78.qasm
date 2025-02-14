OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(-0.21259354) q[0];
sx q[0];
rz(0.73417443) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(-0.34595481) q[1];
sx q[1];
rz(-0.024756519) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9182565) q[0];
sx q[0];
rz(-1.6331216) q[0];
sx q[0];
rz(-1.614109) q[0];
x q[1];
rz(-1.8046298) q[2];
sx q[2];
rz(-1.4604092) q[2];
sx q[2];
rz(1.8821007) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7721482) q[1];
sx q[1];
rz(-1.2719063) q[1];
sx q[1];
rz(1.5159831) q[1];
rz(1.7832028) q[3];
sx q[3];
rz(-1.5161361) q[3];
sx q[3];
rz(2.2509991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2572702) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(-0.30882588) q[2];
rz(2.3158) q[3];
sx q[3];
rz(-1.0079931) q[3];
sx q[3];
rz(-0.76396137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.17077133) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(-2.1139076) q[0];
rz(0.81614256) q[1];
sx q[1];
rz(-1.1183389) q[1];
sx q[1];
rz(2.3291086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.621429) q[0];
sx q[0];
rz(-1.7296711) q[0];
sx q[0];
rz(-1.2170674) q[0];
rz(0.37336739) q[2];
sx q[2];
rz(-1.5575711) q[2];
sx q[2];
rz(0.91886273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3025488) q[1];
sx q[1];
rz(-1.4329646) q[1];
sx q[1];
rz(2.0853569) q[1];
x q[2];
rz(0.3489209) q[3];
sx q[3];
rz(-1.2349718) q[3];
sx q[3];
rz(0.15062697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5673148) q[2];
sx q[2];
rz(-1.4894166) q[2];
sx q[2];
rz(0.91892773) q[2];
rz(-0.53141665) q[3];
sx q[3];
rz(-2.0960977) q[3];
sx q[3];
rz(0.04118583) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4384005) q[0];
sx q[0];
rz(-1.0696609) q[0];
sx q[0];
rz(-1.137314) q[0];
rz(-1.3905585) q[1];
sx q[1];
rz(-0.5568234) q[1];
sx q[1];
rz(-2.4086319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8370286) q[0];
sx q[0];
rz(-2.6246855) q[0];
sx q[0];
rz(-2.5190973) q[0];
rz(-pi) q[1];
rz(0.82436136) q[2];
sx q[2];
rz(-0.75468735) q[2];
sx q[2];
rz(3.1122308) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3324686) q[1];
sx q[1];
rz(-2.4753503) q[1];
sx q[1];
rz(0.82102832) q[1];
rz(-pi) q[2];
rz(-0.0065782733) q[3];
sx q[3];
rz(-2.0195144) q[3];
sx q[3];
rz(2.6980163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0040032337) q[2];
sx q[2];
rz(-1.4121476) q[2];
sx q[2];
rz(-2.1027193) q[2];
rz(-2.1393356) q[3];
sx q[3];
rz(-1.301731) q[3];
sx q[3];
rz(2.8514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20659474) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(-0.90718734) q[0];
rz(-0.15779237) q[1];
sx q[1];
rz(-2.1267499) q[1];
sx q[1];
rz(1.8603604) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8195551) q[0];
sx q[0];
rz(-2.0896308) q[0];
sx q[0];
rz(0.71625336) q[0];
rz(-pi) q[1];
rz(-0.8550941) q[2];
sx q[2];
rz(-1.8786598) q[2];
sx q[2];
rz(1.5941053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7042027) q[1];
sx q[1];
rz(-1.9949732) q[1];
sx q[1];
rz(2.4553039) q[1];
rz(0.90003711) q[3];
sx q[3];
rz(-1.0484015) q[3];
sx q[3];
rz(-2.1417422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24718757) q[2];
sx q[2];
rz(-2.2053568) q[2];
sx q[2];
rz(1.2897162) q[2];
rz(-1.3442518) q[3];
sx q[3];
rz(-1.684609) q[3];
sx q[3];
rz(-0.70014203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808788) q[0];
sx q[0];
rz(-0.49837708) q[0];
sx q[0];
rz(2.1233249) q[0];
rz(-2.0385888) q[1];
sx q[1];
rz(-2.4296727) q[1];
sx q[1];
rz(1.8011372) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.074891) q[0];
sx q[0];
rz(-1.929787) q[0];
sx q[0];
rz(-0.51407459) q[0];
x q[1];
rz(1.4730155) q[2];
sx q[2];
rz(-1.9488397) q[2];
sx q[2];
rz(-2.5387272) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.75234883) q[1];
sx q[1];
rz(-1.6977786) q[1];
sx q[1];
rz(-1.4402706) q[1];
rz(-pi) q[2];
rz(3.0839447) q[3];
sx q[3];
rz(-1.0312005) q[3];
sx q[3];
rz(0.12069139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(-1.0020024) q[2];
rz(2.4501948) q[3];
sx q[3];
rz(-1.9611497) q[3];
sx q[3];
rz(1.6208167) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6028676) q[0];
sx q[0];
rz(-2.0724917) q[0];
sx q[0];
rz(2.5163203) q[0];
rz(1.9484776) q[1];
sx q[1];
rz(-2.4517877) q[1];
sx q[1];
rz(-0.56732059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7869572) q[0];
sx q[0];
rz(-0.79878317) q[0];
sx q[0];
rz(1.3243933) q[0];
rz(-pi) q[1];
rz(0.82639931) q[2];
sx q[2];
rz(-1.5763088) q[2];
sx q[2];
rz(2.7720087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7116136) q[1];
sx q[1];
rz(-2.1890854) q[1];
sx q[1];
rz(0.031876335) q[1];
rz(-1.7501168) q[3];
sx q[3];
rz(-1.3973425) q[3];
sx q[3];
rz(-1.7141527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83711964) q[2];
sx q[2];
rz(-1.8581055) q[2];
sx q[2];
rz(-1.5597957) q[2];
rz(0.61257735) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(-2.9858203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0953858) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(2.0147391) q[0];
rz(-1.2696179) q[1];
sx q[1];
rz(-1.9536628) q[1];
sx q[1];
rz(-0.52106214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4241735) q[0];
sx q[0];
rz(-1.9195286) q[0];
sx q[0];
rz(-0.86535378) q[0];
rz(1.2418141) q[2];
sx q[2];
rz(-1.5397203) q[2];
sx q[2];
rz(2.6779793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.55406694) q[1];
sx q[1];
rz(-0.8272285) q[1];
sx q[1];
rz(0.83283333) q[1];
x q[2];
rz(0.54102202) q[3];
sx q[3];
rz(-1.8086026) q[3];
sx q[3];
rz(1.0353119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0236987) q[2];
sx q[2];
rz(-1.7774899) q[2];
sx q[2];
rz(-1.2269616) q[2];
rz(-1.6795233) q[3];
sx q[3];
rz(-1.6242124) q[3];
sx q[3];
rz(0.44155651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284477) q[0];
sx q[0];
rz(-1.0297091) q[0];
sx q[0];
rz(-0.5471158) q[0];
rz(0.088134915) q[1];
sx q[1];
rz(-1.3476177) q[1];
sx q[1];
rz(0.42253447) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2937909) q[0];
sx q[0];
rz(-1.592684) q[0];
sx q[0];
rz(1.2897911) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5498075) q[2];
sx q[2];
rz(-0.6066423) q[2];
sx q[2];
rz(0.085651407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48848885) q[1];
sx q[1];
rz(-3.0526027) q[1];
sx q[1];
rz(2.2527534) q[1];
rz(-pi) q[2];
rz(-2.2219031) q[3];
sx q[3];
rz(-0.22907478) q[3];
sx q[3];
rz(0.13815115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2253458) q[2];
sx q[2];
rz(-1.6371181) q[2];
sx q[2];
rz(-0.28727356) q[2];
rz(-2.5987127) q[3];
sx q[3];
rz(-0.77016872) q[3];
sx q[3];
rz(-3.0834901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420237) q[0];
sx q[0];
rz(-0.68411198) q[0];
sx q[0];
rz(2.3642819) q[0];
rz(2.2151392) q[1];
sx q[1];
rz(-1.1421721) q[1];
sx q[1];
rz(1.3949589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3789916) q[0];
sx q[0];
rz(-1.9585573) q[0];
sx q[0];
rz(-2.7975525) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59916536) q[2];
sx q[2];
rz(-1.7788506) q[2];
sx q[2];
rz(-2.4924202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0277176) q[1];
sx q[1];
rz(-2.0248955) q[1];
sx q[1];
rz(-2.7732631) q[1];
x q[2];
rz(-0.93392196) q[3];
sx q[3];
rz(-1.0524629) q[3];
sx q[3];
rz(-2.3766488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7159783) q[2];
sx q[2];
rz(-1.3724816) q[2];
sx q[2];
rz(-1.0905637) q[2];
rz(-2.1770832) q[3];
sx q[3];
rz(-1.0114074) q[3];
sx q[3];
rz(1.2151659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57508093) q[0];
sx q[0];
rz(-0.33841857) q[0];
sx q[0];
rz(1.6315208) q[0];
rz(-0.034320023) q[1];
sx q[1];
rz(-1.3395373) q[1];
sx q[1];
rz(0.43201772) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2902108) q[0];
sx q[0];
rz(-2.599596) q[0];
sx q[0];
rz(1.4567514) q[0];
rz(-pi) q[1];
rz(1.2197184) q[2];
sx q[2];
rz(-1.4939927) q[2];
sx q[2];
rz(-1.5157645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68709438) q[1];
sx q[1];
rz(-1.7155572) q[1];
sx q[1];
rz(-1.471038) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98938359) q[3];
sx q[3];
rz(-1.3477147) q[3];
sx q[3];
rz(-1.3815961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4711275) q[2];
sx q[2];
rz(-1.9894783) q[2];
sx q[2];
rz(-2.067789) q[2];
rz(0.18887575) q[3];
sx q[3];
rz(-2.5329068) q[3];
sx q[3];
rz(-2.7591738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0976681) q[0];
sx q[0];
rz(-2.8207939) q[0];
sx q[0];
rz(2.4378142) q[0];
rz(1.7882998) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(-2.399171) q[2];
sx q[2];
rz(-2.0423642) q[2];
sx q[2];
rz(-2.5004417) q[2];
rz(-2.682229) q[3];
sx q[3];
rz(-1.50926) q[3];
sx q[3];
rz(2.0987233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
