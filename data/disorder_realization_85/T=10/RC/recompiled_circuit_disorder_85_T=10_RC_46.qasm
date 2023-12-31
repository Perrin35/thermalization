OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(-2.2493275) q[0];
x q[1];
rz(1.0797834) q[2];
sx q[2];
rz(-0.91677374) q[2];
sx q[2];
rz(0.71066463) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10378097) q[1];
sx q[1];
rz(-1.722464) q[1];
sx q[1];
rz(-1.2396461) q[1];
rz(-1.6269496) q[3];
sx q[3];
rz(-2.3893642) q[3];
sx q[3];
rz(-2.5889531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(0.98891813) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724021) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-2.4172799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(-2.2749167) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.486859) q[2];
sx q[2];
rz(-1.3720023) q[2];
sx q[2];
rz(-2.4059911) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.69246768) q[1];
sx q[1];
rz(-1.7763224) q[1];
sx q[1];
rz(-1.9301901) q[1];
x q[2];
rz(-0.69085391) q[3];
sx q[3];
rz(-0.30403954) q[3];
sx q[3];
rz(-3.0512878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(0.36188564) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.2190855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65788668) q[0];
sx q[0];
rz(-1.3077967) q[0];
sx q[0];
rz(-1.9348683) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8692603) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(1.9321835) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29875007) q[1];
sx q[1];
rz(-2.0992273) q[1];
sx q[1];
rz(0.28038402) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3027906) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-2.5915495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7200155) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.4455618) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(0.51825994) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(2.3148361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8985734) q[0];
sx q[0];
rz(-1.1613701) q[0];
sx q[0];
rz(-0.55222521) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7051758) q[2];
sx q[2];
rz(-2.3339286) q[2];
sx q[2];
rz(-0.32854167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7552232) q[1];
sx q[1];
rz(-1.9519023) q[1];
sx q[1];
rz(1.8783046) q[1];
rz(-pi) q[2];
rz(1.8825674) q[3];
sx q[3];
rz(-2.2707553) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(1.1258874) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-2.856423) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.973424) q[0];
sx q[0];
rz(-2.0452721) q[0];
sx q[0];
rz(3.0939177) q[0];
rz(1.1826913) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(-0.91458048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3030745) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(1.5694373) q[1];
rz(-pi) q[2];
rz(0.22484803) q[3];
sx q[3];
rz(-2.7307011) q[3];
sx q[3];
rz(-1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(-2.6021393) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(1.3670115) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052664) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(-1.1195539) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96732803) q[2];
sx q[2];
rz(-2.2852995) q[2];
sx q[2];
rz(-0.29010233) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65946603) q[1];
sx q[1];
rz(-2.8409993) q[1];
sx q[1];
rz(-1.2947004) q[1];
rz(-0.35297024) q[3];
sx q[3];
rz(-1.4735231) q[3];
sx q[3];
rz(-0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(-2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0872333) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(-1.5646294) q[0];
rz(-pi) q[1];
rz(-0.41202338) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(0.7574946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2504187) q[1];
sx q[1];
rz(-0.95953566) q[1];
sx q[1];
rz(0.82959081) q[1];
x q[2];
rz(2.6550754) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-0.61075413) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(2.4954605) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883811) q[0];
sx q[0];
rz(-1.4697945) q[0];
sx q[0];
rz(1.7779011) q[0];
rz(-1.7392776) q[2];
sx q[2];
rz(-1.104276) q[2];
sx q[2];
rz(-2.434935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38478002) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(-0.53201075) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(-2.3074647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(2.6867552) q[2];
rz(2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(-2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-0.79750693) q[0];
rz(-0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(-0.10841766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90005504) q[0];
sx q[0];
rz(-2.2336322) q[0];
sx q[0];
rz(-0.073297757) q[0];
x q[1];
rz(1.8234532) q[2];
sx q[2];
rz(-0.18441072) q[2];
sx q[2];
rz(0.84469634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81356955) q[1];
sx q[1];
rz(-1.4038101) q[1];
sx q[1];
rz(-1.7896673) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86650642) q[3];
sx q[3];
rz(-1.4276854) q[3];
sx q[3];
rz(1.0673616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(-0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(3.0864339) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3653152) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(-2.9676653) q[0];
x q[1];
rz(2.8137384) q[2];
sx q[2];
rz(-0.58460669) q[2];
sx q[2];
rz(0.84529982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2436279) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(0.90378739) q[1];
rz(-1.1147898) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(1.959521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-0.49017635) q[2];
rz(-3.0040719) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(-0.31221496) q[2];
sx q[2];
rz(-1.6630465) q[2];
sx q[2];
rz(1.9065471) q[2];
rz(-0.18110885) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
