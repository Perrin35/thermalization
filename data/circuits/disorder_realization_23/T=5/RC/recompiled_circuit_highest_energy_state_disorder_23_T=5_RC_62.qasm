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
rz(0.55573207) q[0];
sx q[0];
rz(-1.8620123) q[0];
sx q[0];
rz(-0.32787856) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(2.1305003) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0450789) q[0];
sx q[0];
rz(-2.156684) q[0];
sx q[0];
rz(-1.7631329) q[0];
x q[1];
rz(1.3921803) q[2];
sx q[2];
rz(-1.1682604) q[2];
sx q[2];
rz(-0.91031633) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15926898) q[1];
sx q[1];
rz(-1.2071274) q[1];
sx q[1];
rz(-0.77467847) q[1];
x q[2];
rz(-2.903721) q[3];
sx q[3];
rz(-1.9783881) q[3];
sx q[3];
rz(-1.067124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27935394) q[2];
sx q[2];
rz(-2.3591159) q[2];
sx q[2];
rz(1.5581101) q[2];
rz(0.33485788) q[3];
sx q[3];
rz(-1.0840651) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8849477) q[0];
sx q[0];
rz(-2.5769951) q[0];
sx q[0];
rz(-2.8096492) q[0];
rz(-0.360082) q[1];
sx q[1];
rz(-1.304345) q[1];
sx q[1];
rz(-0.28775451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.677414) q[0];
sx q[0];
rz(-2.8197643) q[0];
sx q[0];
rz(-2.4610956) q[0];
rz(-1.2505202) q[2];
sx q[2];
rz(-0.8647635) q[2];
sx q[2];
rz(1.6635739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8043878) q[1];
sx q[1];
rz(-1.8998977) q[1];
sx q[1];
rz(1.9371402) q[1];
x q[2];
rz(-0.61057456) q[3];
sx q[3];
rz(-0.27974162) q[3];
sx q[3];
rz(0.66204643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.687279) q[2];
sx q[2];
rz(-1.6087029) q[2];
sx q[2];
rz(1.7506556) q[2];
rz(-2.2823997) q[3];
sx q[3];
rz(-1.8682559) q[3];
sx q[3];
rz(0.19101645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39621064) q[0];
sx q[0];
rz(-1.5353545) q[0];
sx q[0];
rz(-2.9111653) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-1.0942065) q[1];
sx q[1];
rz(-0.28883019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.13076) q[0];
sx q[0];
rz(-1.1625746) q[0];
sx q[0];
rz(-3.0768763) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46858139) q[2];
sx q[2];
rz(-1.6372674) q[2];
sx q[2];
rz(2.2867212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9625712) q[1];
sx q[1];
rz(-1.6238535) q[1];
sx q[1];
rz(-1.1607443) q[1];
rz(2.9656319) q[3];
sx q[3];
rz(-1.7909414) q[3];
sx q[3];
rz(-2.1484321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1032224) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(3.1008516) q[2];
rz(0.1013969) q[3];
sx q[3];
rz(-1.8056168) q[3];
sx q[3];
rz(2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.7539702) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(-0.23400865) q[0];
rz(-2.9435844) q[1];
sx q[1];
rz(-1.0375236) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957243) q[0];
sx q[0];
rz(-1.0875174) q[0];
sx q[0];
rz(-3.0935174) q[0];
rz(-pi) q[1];
rz(2.9806541) q[2];
sx q[2];
rz(-0.51921028) q[2];
sx q[2];
rz(-0.96566654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3819921) q[1];
sx q[1];
rz(-0.97642094) q[1];
sx q[1];
rz(2.4305953) q[1];
x q[2];
rz(1.181482) q[3];
sx q[3];
rz(-1.5411756) q[3];
sx q[3];
rz(-2.9148341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6733751) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(-0.77486983) q[2];
rz(-1.8153927) q[3];
sx q[3];
rz(-1.9522791) q[3];
sx q[3];
rz(0.51138043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.73959094) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(0.032489754) q[0];
rz(-1.912311) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(3.0585739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0262526) q[0];
sx q[0];
rz(-1.7949545) q[0];
sx q[0];
rz(0.46333939) q[0];
rz(-pi) q[1];
rz(-0.26814383) q[2];
sx q[2];
rz(-1.4073155) q[2];
sx q[2];
rz(0.016506052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8889035) q[1];
sx q[1];
rz(-1.1071883) q[1];
sx q[1];
rz(2.3867334) q[1];
x q[2];
rz(-1.5918077) q[3];
sx q[3];
rz(-0.63696955) q[3];
sx q[3];
rz(1.6629499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3518389) q[2];
sx q[2];
rz(-0.80078501) q[2];
sx q[2];
rz(0.16072533) q[2];
rz(1.1927346) q[3];
sx q[3];
rz(-2.9294117) q[3];
sx q[3];
rz(-0.76782697) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6650894) q[0];
sx q[0];
rz(-1.1192717) q[0];
sx q[0];
rz(2.8506668) q[0];
rz(-1.7581958) q[1];
sx q[1];
rz(-1.29888) q[1];
sx q[1];
rz(2.6417522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4668149) q[0];
sx q[0];
rz(-1.5272015) q[0];
sx q[0];
rz(-1.6131667) q[0];
x q[1];
rz(-0.84264522) q[2];
sx q[2];
rz(-1.1367186) q[2];
sx q[2];
rz(-1.7378716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19491296) q[1];
sx q[1];
rz(-0.99338594) q[1];
sx q[1];
rz(2.1140079) q[1];
rz(0.53414102) q[3];
sx q[3];
rz(-0.93399601) q[3];
sx q[3];
rz(2.4780688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9342039) q[2];
sx q[2];
rz(-0.61177212) q[2];
sx q[2];
rz(-2.440051) q[2];
rz(0.97405854) q[3];
sx q[3];
rz(-0.8166703) q[3];
sx q[3];
rz(-0.90132236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262254) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(2.4819964) q[0];
rz(-1.9024128) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(-1.0741796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.488193) q[0];
sx q[0];
rz(-1.5852514) q[0];
sx q[0];
rz(-0.01343347) q[0];
rz(-pi) q[1];
rz(-1.1848106) q[2];
sx q[2];
rz(-1.1710656) q[2];
sx q[2];
rz(1.0077493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2905137) q[1];
sx q[1];
rz(-2.004262) q[1];
sx q[1];
rz(2.9032533) q[1];
rz(-pi) q[2];
rz(0.45589186) q[3];
sx q[3];
rz(-0.8178725) q[3];
sx q[3];
rz(-0.65480937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2913975) q[2];
sx q[2];
rz(-0.75866282) q[2];
sx q[2];
rz(-2.5568753) q[2];
rz(2.0753453) q[3];
sx q[3];
rz(-1.9138347) q[3];
sx q[3];
rz(2.7339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(-2.6608652) q[0];
rz(1.9302543) q[1];
sx q[1];
rz(-1.5296661) q[1];
sx q[1];
rz(-0.617625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077972875) q[0];
sx q[0];
rz(-1.5665861) q[0];
sx q[0];
rz(1.8959037) q[0];
x q[1];
rz(-2.5343393) q[2];
sx q[2];
rz(-2.0130081) q[2];
sx q[2];
rz(-2.2970125) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1943975) q[1];
sx q[1];
rz(-0.36809599) q[1];
sx q[1];
rz(0.93666623) q[1];
x q[2];
rz(0.17980735) q[3];
sx q[3];
rz(-2.9463065) q[3];
sx q[3];
rz(0.92474588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9628613) q[2];
sx q[2];
rz(-1.7469254) q[2];
sx q[2];
rz(-2.2507131) q[2];
rz(-2.1122872) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.7903098) q[0];
sx q[0];
rz(-2.2064378) q[0];
sx q[0];
rz(1.9679605) q[0];
rz(1.952518) q[1];
sx q[1];
rz(-0.74780858) q[1];
sx q[1];
rz(2.660451) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089398459) q[0];
sx q[0];
rz(-1.3476831) q[0];
sx q[0];
rz(2.5714178) q[0];
rz(-pi) q[1];
rz(2.5817211) q[2];
sx q[2];
rz(-1.6563043) q[2];
sx q[2];
rz(-2.1712803) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12901017) q[1];
sx q[1];
rz(-2.5222188) q[1];
sx q[1];
rz(1.315809) q[1];
rz(-pi) q[2];
rz(1.0403544) q[3];
sx q[3];
rz(-1.181349) q[3];
sx q[3];
rz(-2.3688909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9159307) q[2];
sx q[2];
rz(-1.91232) q[2];
sx q[2];
rz(-1.6602328) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(-1.2272629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34761053) q[0];
sx q[0];
rz(-2.1078258) q[0];
sx q[0];
rz(-0.069742918) q[0];
rz(-2.8522885) q[1];
sx q[1];
rz(-1.2846839) q[1];
sx q[1];
rz(-2.0893673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69221321) q[0];
sx q[0];
rz(-1.2468296) q[0];
sx q[0];
rz(0.37925668) q[0];
x q[1];
rz(-1.7830816) q[2];
sx q[2];
rz(-0.67321482) q[2];
sx q[2];
rz(-2.5517983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94967647) q[1];
sx q[1];
rz(-1.1947617) q[1];
sx q[1];
rz(-2.7533349) q[1];
x q[2];
rz(0.77731384) q[3];
sx q[3];
rz(-0.56985006) q[3];
sx q[3];
rz(-0.22138813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59794402) q[2];
sx q[2];
rz(-1.9660549) q[2];
sx q[2];
rz(-3.0808466) q[2];
rz(2.8525823) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(-1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.65171) q[0];
sx q[0];
rz(-1.6429506) q[0];
sx q[0];
rz(2.0605675) q[0];
rz(-2.091485) q[1];
sx q[1];
rz(-1.0625912) q[1];
sx q[1];
rz(-1.2407632) q[1];
rz(-1.8619887) q[2];
sx q[2];
rz(-0.48710258) q[2];
sx q[2];
rz(-1.0500548) q[2];
rz(-2.1331187) q[3];
sx q[3];
rz(-2.6431966) q[3];
sx q[3];
rz(0.24503844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
