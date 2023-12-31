OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(0.18145951) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(1.0531309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7979413) q[0];
sx q[0];
rz(-0.5125674) q[0];
sx q[0];
rz(-1.1287862) q[0];
x q[1];
rz(1.2744781) q[2];
sx q[2];
rz(-0.9689435) q[2];
sx q[2];
rz(-1.559343) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6932106) q[1];
sx q[1];
rz(-2.3464077) q[1];
sx q[1];
rz(2.8661212) q[1];
x q[2];
rz(-0.16206046) q[3];
sx q[3];
rz(-1.9063623) q[3];
sx q[3];
rz(-0.1048564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(-1.7738316) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(0.59511551) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(1.6275303) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9428974) q[0];
sx q[0];
rz(-1.9778344) q[0];
sx q[0];
rz(-2.1172949) q[0];
x q[1];
rz(2.3081231) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(-1.2534864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82799998) q[1];
sx q[1];
rz(-0.26121059) q[1];
sx q[1];
rz(-2.2952495) q[1];
rz(-2.0239003) q[3];
sx q[3];
rz(-0.11370224) q[3];
sx q[3];
rz(-0.5773471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(2.3454323) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(-0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650836) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(-2.8308716) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(-1.6832738) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761988) q[0];
sx q[0];
rz(-1.5237234) q[0];
sx q[0];
rz(-1.7182299) q[0];
rz(2.4887423) q[2];
sx q[2];
rz(-2.0712426) q[2];
sx q[2];
rz(0.99393883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4255193) q[1];
sx q[1];
rz(-0.41285535) q[1];
sx q[1];
rz(-0.13456657) q[1];
rz(-0.31483105) q[3];
sx q[3];
rz(-1.5545087) q[3];
sx q[3];
rz(-0.44486886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1594499) q[2];
sx q[2];
rz(-2.2950256) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-2.5722752) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(-1.9212978) q[1];
sx q[1];
rz(-1.9492457) q[1];
sx q[1];
rz(1.6569998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.785448) q[0];
sx q[0];
rz(-1.6011642) q[0];
sx q[0];
rz(-1.5201475) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2414615) q[2];
sx q[2];
rz(-2.6000458) q[2];
sx q[2];
rz(-0.75377084) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.93814072) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(0.52312619) q[1];
rz(-0.72317601) q[3];
sx q[3];
rz(-1.3190862) q[3];
sx q[3];
rz(0.46079208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-1.0162639) q[2];
rz(-1.7381564) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.50773412) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(1.9790861) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(0.17366017) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0807053) q[0];
sx q[0];
rz(-1.6433435) q[0];
sx q[0];
rz(1.5990431) q[0];
rz(-2.6661751) q[2];
sx q[2];
rz(-2.7042275) q[2];
sx q[2];
rz(1.6944483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34827161) q[1];
sx q[1];
rz(-1.4655359) q[1];
sx q[1];
rz(-3.1131016) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8097314) q[3];
sx q[3];
rz(-1.7680941) q[3];
sx q[3];
rz(-2.8062537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(-2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(-2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.1059234) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.1434198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15018806) q[0];
sx q[0];
rz(-0.17230573) q[0];
sx q[0];
rz(2.6323787) q[0];
rz(2.5383699) q[2];
sx q[2];
rz(-1.8458741) q[2];
sx q[2];
rz(-3.0657257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.032635078) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(2.0202548) q[1];
x q[2];
rz(1.0295168) q[3];
sx q[3];
rz(-1.8347077) q[3];
sx q[3];
rz(-1.908386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64289552) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(1.5303622) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(-2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(-1.4402333) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92111174) q[0];
sx q[0];
rz(-2.8163914) q[0];
sx q[0];
rz(3.0082506) q[0];
x q[1];
rz(1.557017) q[2];
sx q[2];
rz(-0.85859495) q[2];
sx q[2];
rz(2.6028002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.95722317) q[1];
sx q[1];
rz(-1.8436699) q[1];
sx q[1];
rz(-0.89783122) q[1];
rz(-0.26182884) q[3];
sx q[3];
rz(-1.4187519) q[3];
sx q[3];
rz(0.1482299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40522727) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-0.48842946) q[2];
rz(1.8296261) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.2088998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8457444) q[0];
sx q[0];
rz(-2.309572) q[0];
sx q[0];
rz(0.31006281) q[0];
rz(-pi) q[1];
rz(-2.9572763) q[2];
sx q[2];
rz(-1.8039928) q[2];
sx q[2];
rz(1.2554982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2541131) q[1];
sx q[1];
rz(-1.5548318) q[1];
sx q[1];
rz(-1.4404526) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62187059) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(-0.57884502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.97380012) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(0.31162509) q[0];
rz(-0.82178003) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(-1.5100381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841543) q[0];
sx q[0];
rz(-0.3785924) q[0];
sx q[0];
rz(0.5666825) q[0];
rz(-pi) q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(2.0575112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4842802) q[1];
sx q[1];
rz(-1.3687951) q[1];
sx q[1];
rz(0.61985086) q[1];
x q[2];
rz(2.908913) q[3];
sx q[3];
rz(-1.3921229) q[3];
sx q[3];
rz(-2.132638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8682378) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-0.81595016) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8907392) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(0.2302641) q[0];
rz(0.62581217) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(-0.65840107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31357665) q[0];
sx q[0];
rz(-2.1470634) q[0];
sx q[0];
rz(2.1616031) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.621643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.876993) q[1];
sx q[1];
rz(-2.0473285) q[1];
sx q[1];
rz(2.9546253) q[1];
rz(-pi) q[2];
rz(2.5913521) q[3];
sx q[3];
rz(-0.60866683) q[3];
sx q[3];
rz(-2.254571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(-2.8038483) q[2];
rz(1.0234458) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.9256928) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(-1.2383923) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(-0.94454371) q[2];
sx q[2];
rz(-1.7164451) q[2];
sx q[2];
rz(-3.0838983) q[2];
rz(-0.099048793) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
