OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(-1.2712103) q[0];
rz(-2.864569) q[1];
sx q[1];
rz(-2.6695873) q[1];
sx q[1];
rz(-0.0013874887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822617) q[0];
sx q[0];
rz(-1.9388655) q[0];
sx q[0];
rz(-2.9666535) q[0];
rz(-pi) q[1];
rz(0.087287993) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(2.0729614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.672294) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(1.4908355) q[1];
x q[2];
rz(2.1625159) q[3];
sx q[3];
rz(-2.7032529) q[3];
sx q[3];
rz(-1.0867659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-2.7367676) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(0.81545365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97074189) q[0];
sx q[0];
rz(-0.65119699) q[0];
sx q[0];
rz(-0.09911508) q[0];
rz(-2.8380727) q[2];
sx q[2];
rz(-2.3887861) q[2];
sx q[2];
rz(2.7941861) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8784139) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(-3.0221699) q[1];
rz(1.3482413) q[3];
sx q[3];
rz(-0.9369623) q[3];
sx q[3];
rz(2.8702877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-0.99951807) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.2721491) q[0];
sx q[0];
rz(0.15655984) q[0];
rz(-pi) q[1];
rz(0.83061647) q[2];
sx q[2];
rz(-0.81095552) q[2];
sx q[2];
rz(0.22731552) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1332902) q[1];
sx q[1];
rz(-1.3160719) q[1];
sx q[1];
rz(-2.6363274) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7129094) q[3];
sx q[3];
rz(-0.51321533) q[3];
sx q[3];
rz(1.7538479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.1497315) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(-0.91039175) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-0.27483637) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4764458) q[0];
sx q[0];
rz(-0.11867141) q[0];
sx q[0];
rz(-0.84400405) q[0];
x q[1];
rz(-0.0015953548) q[2];
sx q[2];
rz(-2.0312107) q[2];
sx q[2];
rz(-2.7574725) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7819314) q[1];
sx q[1];
rz(-1.1210103) q[1];
sx q[1];
rz(0.91314544) q[1];
x q[2];
rz(-2.2673625) q[3];
sx q[3];
rz(-1.4733553) q[3];
sx q[3];
rz(-1.8271354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(-1.516974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6277916) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(2.6898726) q[0];
x q[1];
rz(0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.5460207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6301873) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(1.1955111) q[1];
rz(-pi) q[2];
rz(-0.77002854) q[3];
sx q[3];
rz(-2.3838245) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(-1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.8776241) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(2.8009159) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5013803) q[0];
sx q[0];
rz(-3.0684154) q[0];
sx q[0];
rz(-2.0206547) q[0];
rz(-pi) q[1];
rz(-0.31008115) q[2];
sx q[2];
rz(-2.3611464) q[2];
sx q[2];
rz(-0.23852894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1250455) q[1];
sx q[1];
rz(-2.6066337) q[1];
sx q[1];
rz(2.6782481) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0601222) q[3];
sx q[3];
rz(-1.193207) q[3];
sx q[3];
rz(-1.9469572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(2.6766434) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5349605) q[0];
sx q[0];
rz(-2.6323942) q[0];
sx q[0];
rz(1.1688933) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49352383) q[2];
sx q[2];
rz(-1.0527305) q[2];
sx q[2];
rz(1.7871737) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3073687) q[1];
sx q[1];
rz(-2.8982179) q[1];
sx q[1];
rz(-2.6878396) q[1];
rz(1.5246478) q[3];
sx q[3];
rz(-1.5449636) q[3];
sx q[3];
rz(1.3290562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(-2.8542744) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(0.28433329) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-0.078358738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6440455) q[0];
sx q[0];
rz(-2.3345778) q[0];
sx q[0];
rz(0.09597309) q[0];
rz(-pi) q[1];
rz(-0.23711726) q[2];
sx q[2];
rz(-1.8097005) q[2];
sx q[2];
rz(2.4799926) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8790508) q[1];
sx q[1];
rz(-1.3971551) q[1];
sx q[1];
rz(-0.91211984) q[1];
rz(-1.7275229) q[3];
sx q[3];
rz(-1.5369475) q[3];
sx q[3];
rz(-1.0364929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(2.890214) q[2];
rz(-2.5583983) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(-2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-2.267568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286205) q[0];
sx q[0];
rz(-0.7681094) q[0];
sx q[0];
rz(0.012499768) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85272312) q[2];
sx q[2];
rz(-2.6211779) q[2];
sx q[2];
rz(0.66327099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2487138) q[1];
sx q[1];
rz(-1.4993748) q[1];
sx q[1];
rz(2.0715269) q[1];
rz(-2.5328818) q[3];
sx q[3];
rz(-1.2864283) q[3];
sx q[3];
rz(-1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(-1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(-1.9627409) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783979) q[0];
sx q[0];
rz(-0.34391719) q[0];
sx q[0];
rz(0.92349903) q[0];
rz(-2.5095021) q[2];
sx q[2];
rz(-2.9846016) q[2];
sx q[2];
rz(1.1319515) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57943343) q[1];
sx q[1];
rz(-0.65083083) q[1];
sx q[1];
rz(1.2594373) q[1];
x q[2];
rz(-1.7970656) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(2.2422092) q[2];
rz(-0.87456885) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(-3.1065431) q[2];
sx q[2];
rz(-1.1309584) q[2];
sx q[2];
rz(-2.1396648) q[2];
rz(1.7189797) q[3];
sx q[3];
rz(-1.8437181) q[3];
sx q[3];
rz(-3.0317882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];