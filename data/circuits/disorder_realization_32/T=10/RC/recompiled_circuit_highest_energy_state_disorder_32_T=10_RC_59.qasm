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
rz(1.2122756) q[0];
sx q[0];
rz(-2.0097998) q[0];
sx q[0];
rz(1.3728859) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(4.5166587) q[1];
sx q[1];
rz(5.4969129) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4117811) q[0];
sx q[0];
rz(-1.3656817) q[0];
sx q[0];
rz(1.1954444) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2462772) q[2];
sx q[2];
rz(-2.3626872) q[2];
sx q[2];
rz(-1.9112183) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3481625) q[1];
sx q[1];
rz(-1.5409267) q[1];
sx q[1];
rz(1.2814327) q[1];
rz(-2.1699675) q[3];
sx q[3];
rz(-1.2279945) q[3];
sx q[3];
rz(-0.179571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11975153) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(1.2500259) q[2];
rz(-0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(-1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.64330548) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(0.071320891) q[0];
rz(-1.9885063) q[1];
sx q[1];
rz(-2.3801897) q[1];
sx q[1];
rz(-2.6128795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708955) q[0];
sx q[0];
rz(-2.5636854) q[0];
sx q[0];
rz(0.96630567) q[0];
rz(2.9527093) q[2];
sx q[2];
rz(-1.3391277) q[2];
sx q[2];
rz(-2.846039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15859998) q[1];
sx q[1];
rz(-2.0865177) q[1];
sx q[1];
rz(0.92293585) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9170179) q[3];
sx q[3];
rz(-1.16515) q[3];
sx q[3];
rz(0.39164513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(-0.028623494) q[2];
rz(-2.7875767) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(-2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(0.89135528) q[0];
rz(-2.640653) q[1];
sx q[1];
rz(-1.5762065) q[1];
sx q[1];
rz(-3.0827789) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0552154) q[0];
sx q[0];
rz(-3.0432467) q[0];
sx q[0];
rz(-0.50758657) q[0];
rz(-pi) q[1];
rz(2.8562653) q[2];
sx q[2];
rz(-1.0614191) q[2];
sx q[2];
rz(0.90778186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9998628) q[1];
sx q[1];
rz(-2.0577621) q[1];
sx q[1];
rz(-0.75090639) q[1];
rz(-0.18043442) q[3];
sx q[3];
rz(-1.5984319) q[3];
sx q[3];
rz(0.35865006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3567051) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(-2.9050262) q[2];
rz(2.2208354) q[3];
sx q[3];
rz(-1.0509138) q[3];
sx q[3];
rz(-1.7304272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22223602) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(-2.2663569) q[0];
rz(1.596176) q[1];
sx q[1];
rz(-0.86219209) q[1];
sx q[1];
rz(-0.045305591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52462477) q[0];
sx q[0];
rz(-3.0201206) q[0];
sx q[0];
rz(-1.0018906) q[0];
rz(0.50778054) q[2];
sx q[2];
rz(-2.1740401) q[2];
sx q[2];
rz(-2.809462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26805763) q[1];
sx q[1];
rz(-2.2378439) q[1];
sx q[1];
rz(3.0309543) q[1];
x q[2];
rz(1.8604467) q[3];
sx q[3];
rz(-2.5389606) q[3];
sx q[3];
rz(-1.8528191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8594325) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(-1.947594) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(-2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0328261) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(2.4591675) q[0];
rz(1.2884864) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.3312181) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2829943) q[0];
sx q[0];
rz(-2.4147075) q[0];
sx q[0];
rz(1.5461393) q[0];
rz(2.8814828) q[2];
sx q[2];
rz(-1.9528198) q[2];
sx q[2];
rz(-1.1191238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9115548) q[1];
sx q[1];
rz(-1.4304407) q[1];
sx q[1];
rz(-0.36169238) q[1];
x q[2];
rz(-0.87792895) q[3];
sx q[3];
rz(-2.4932541) q[3];
sx q[3];
rz(0.058365783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1149301) q[2];
sx q[2];
rz(-0.59938359) q[2];
sx q[2];
rz(-0.66693532) q[2];
rz(0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(-2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126467) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(0.70372787) q[0];
rz(-2.9723889) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(0.15883787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7089091) q[0];
sx q[0];
rz(-1.428621) q[0];
sx q[0];
rz(-3.0199023) q[0];
rz(-0.67262291) q[2];
sx q[2];
rz(-0.19572283) q[2];
sx q[2];
rz(-1.7253523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4602214) q[1];
sx q[1];
rz(-1.7029783) q[1];
sx q[1];
rz(-0.48355196) q[1];
rz(-pi) q[2];
rz(0.30384003) q[3];
sx q[3];
rz(-2.4812077) q[3];
sx q[3];
rz(0.25750289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3682897) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(2.2840195) q[2];
rz(-1.1693906) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17539772) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(-0.71863693) q[0];
rz(2.8596558) q[1];
sx q[1];
rz(-1.3749296) q[1];
sx q[1];
rz(-2.4729572) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91808337) q[0];
sx q[0];
rz(-1.9094206) q[0];
sx q[0];
rz(0.91135599) q[0];
rz(-pi) q[1];
rz(1.093542) q[2];
sx q[2];
rz(-0.68533939) q[2];
sx q[2];
rz(2.6743741) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.201828) q[1];
sx q[1];
rz(-1.9376905) q[1];
sx q[1];
rz(-2.7903627) q[1];
rz(-pi) q[2];
rz(-2.5103522) q[3];
sx q[3];
rz(-1.5774014) q[3];
sx q[3];
rz(-0.94830482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4947027) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(-2.1051712) q[2];
rz(0.63452619) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46102872) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(-2.4155937) q[0];
rz(1.9793824) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(-1.9451709) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1370173) q[0];
sx q[0];
rz(-1.5992303) q[0];
sx q[0];
rz(-1.397055) q[0];
rz(-1.7199207) q[2];
sx q[2];
rz(-0.31055488) q[2];
sx q[2];
rz(0.20937411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58652523) q[1];
sx q[1];
rz(-1.0690926) q[1];
sx q[1];
rz(2.1074159) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4287339) q[3];
sx q[3];
rz(-2.4247243) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90765816) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(-2.2080803) q[2];
rz(2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(-0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.9762978) q[0];
sx q[0];
rz(-0.10460654) q[0];
sx q[0];
rz(-0.053255178) q[0];
rz(2.8540197) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(2.8823749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142562) q[0];
sx q[0];
rz(-1.5467318) q[0];
sx q[0];
rz(-1.5177478) q[0];
rz(-pi) q[1];
rz(1.0561814) q[2];
sx q[2];
rz(-0.26517235) q[2];
sx q[2];
rz(-2.3725703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0269912) q[1];
sx q[1];
rz(-0.82102126) q[1];
sx q[1];
rz(0.084084078) q[1];
x q[2];
rz(-0.76595347) q[3];
sx q[3];
rz(-1.6807091) q[3];
sx q[3];
rz(0.069086941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58664924) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(2.9602236) q[2];
rz(2.7465316) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(-2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3996537) q[0];
sx q[0];
rz(-1.5927277) q[0];
sx q[0];
rz(-0.3821061) q[0];
rz(0.29640472) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(1.2688676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4110376) q[0];
sx q[0];
rz(-2.0417622) q[0];
sx q[0];
rz(2.6593326) q[0];
x q[1];
rz(-0.13258719) q[2];
sx q[2];
rz(-0.90182038) q[2];
sx q[2];
rz(-0.35945177) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8812534) q[1];
sx q[1];
rz(-1.2311544) q[1];
sx q[1];
rz(-0.022625523) q[1];
rz(2.6878854) q[3];
sx q[3];
rz(-1.4611011) q[3];
sx q[3];
rz(-3.0806993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84539139) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(0.38983795) q[2];
rz(1.2975533) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(1.3564431) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(-0.46089725) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(1.9080347) q[3];
sx q[3];
rz(-0.81546406) q[3];
sx q[3];
rz(-2.245695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
