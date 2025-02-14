OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73622048) q[0];
sx q[0];
rz(3.1625746) q[0];
sx q[0];
rz(7.4641135) q[0];
rz(2.6612072) q[1];
sx q[1];
rz(-1.1552224) q[1];
sx q[1];
rz(0.46673271) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8202756) q[0];
sx q[0];
rz(-1.2877687) q[0];
sx q[0];
rz(-0.77518605) q[0];
rz(-pi) q[1];
rz(2.0601022) q[2];
sx q[2];
rz(-2.0976411) q[2];
sx q[2];
rz(-2.26109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4157652) q[1];
sx q[1];
rz(-1.5752324) q[1];
sx q[1];
rz(1.8091101) q[1];
rz(-pi) q[2];
rz(-2.4634697) q[3];
sx q[3];
rz(-2.0270542) q[3];
sx q[3];
rz(-0.18517906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0878318) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(-1.9161179) q[2];
rz(-2.7671704) q[3];
sx q[3];
rz(-2.008308) q[3];
sx q[3];
rz(-2.5882914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.18749172) q[0];
sx q[0];
rz(-2.2611389) q[0];
sx q[0];
rz(-2.6075897) q[0];
rz(-2.7502637) q[1];
sx q[1];
rz(-2.6983039) q[1];
sx q[1];
rz(2.2543529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3852291) q[0];
sx q[0];
rz(-1.8779715) q[0];
sx q[0];
rz(0.36184575) q[0];
rz(-2.3848955) q[2];
sx q[2];
rz(-1.1771508) q[2];
sx q[2];
rz(-1.4243038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5767867) q[1];
sx q[1];
rz(-2.8999355) q[1];
sx q[1];
rz(2.5741626) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0530472) q[3];
sx q[3];
rz(-1.8705006) q[3];
sx q[3];
rz(-0.24417711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54681626) q[2];
sx q[2];
rz(-1.3204601) q[2];
sx q[2];
rz(-2.7776264) q[2];
rz(-1.4173896) q[3];
sx q[3];
rz(-2.851749) q[3];
sx q[3];
rz(-0.83436203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2367547) q[0];
sx q[0];
rz(-0.098243864) q[0];
sx q[0];
rz(0.33016095) q[0];
rz(-2.5579021) q[1];
sx q[1];
rz(-0.81283641) q[1];
sx q[1];
rz(2.6469753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9512105) q[0];
sx q[0];
rz(-1.8415762) q[0];
sx q[0];
rz(0.27373224) q[0];
x q[1];
rz(-2.0614659) q[2];
sx q[2];
rz(-1.3441836) q[2];
sx q[2];
rz(-1.8267711) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6909852) q[1];
sx q[1];
rz(-2.4276884) q[1];
sx q[1];
rz(0.83414487) q[1];
rz(-pi) q[2];
rz(2.1547537) q[3];
sx q[3];
rz(-1.47746) q[3];
sx q[3];
rz(-3.0350041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3459449) q[2];
sx q[2];
rz(-2.4730885) q[2];
sx q[2];
rz(3.0188959) q[2];
rz(3.0986541) q[3];
sx q[3];
rz(-1.0791082) q[3];
sx q[3];
rz(0.19449657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.73404679) q[0];
sx q[0];
rz(-2.6341944) q[0];
sx q[0];
rz(-0.88974446) q[0];
rz(-0.45097688) q[1];
sx q[1];
rz(-0.86800066) q[1];
sx q[1];
rz(2.6873592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4047336) q[0];
sx q[0];
rz(-1.6772207) q[0];
sx q[0];
rz(2.7989796) q[0];
rz(-pi) q[1];
rz(0.003963917) q[2];
sx q[2];
rz(-0.42734738) q[2];
sx q[2];
rz(1.3766118) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.046084) q[1];
sx q[1];
rz(-2.3254776) q[1];
sx q[1];
rz(2.7296394) q[1];
rz(1.4934214) q[3];
sx q[3];
rz(-1.2168125) q[3];
sx q[3];
rz(3.1396239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1526327) q[2];
sx q[2];
rz(-2.5071414) q[2];
sx q[2];
rz(-0.75759849) q[2];
rz(-2.9518413) q[3];
sx q[3];
rz(-1.3187783) q[3];
sx q[3];
rz(0.54722133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106638) q[0];
sx q[0];
rz(-0.78721109) q[0];
sx q[0];
rz(1.2934562) q[0];
rz(0.58213195) q[1];
sx q[1];
rz(-1.4385834) q[1];
sx q[1];
rz(2.9621946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6027733) q[0];
sx q[0];
rz(-1.3673377) q[0];
sx q[0];
rz(-0.18360965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30443983) q[2];
sx q[2];
rz(-0.6938664) q[2];
sx q[2];
rz(-2.2562698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7779667) q[1];
sx q[1];
rz(-2.1898664) q[1];
sx q[1];
rz(-2.9113063) q[1];
x q[2];
rz(0.50197451) q[3];
sx q[3];
rz(-0.93794146) q[3];
sx q[3];
rz(-2.8090734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3680215) q[2];
sx q[2];
rz(-0.37006912) q[2];
sx q[2];
rz(2.5717112) q[2];
rz(0.44024769) q[3];
sx q[3];
rz(-1.2443685) q[3];
sx q[3];
rz(-0.35278916) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84949795) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(2.0000892) q[0];
rz(1.796272) q[1];
sx q[1];
rz(-2.1514386) q[1];
sx q[1];
rz(-0.17123953) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9123403) q[0];
sx q[0];
rz(-1.0919337) q[0];
sx q[0];
rz(-1.6845869) q[0];
rz(-pi) q[1];
rz(-2.6944567) q[2];
sx q[2];
rz(-1.0903745) q[2];
sx q[2];
rz(0.26408476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0192193) q[1];
sx q[1];
rz(-1.7369441) q[1];
sx q[1];
rz(-0.56143399) q[1];
x q[2];
rz(0.01252894) q[3];
sx q[3];
rz(-1.0141148) q[3];
sx q[3];
rz(-1.0906681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3693927) q[2];
sx q[2];
rz(-2.2056396) q[2];
sx q[2];
rz(-0.13747036) q[2];
rz(1.2482268) q[3];
sx q[3];
rz(-2.5745001) q[3];
sx q[3];
rz(1.9198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084994706) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(-1.6531264) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.4005125) q[1];
sx q[1];
rz(1.8280425) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1689455) q[0];
sx q[0];
rz(-0.15014507) q[0];
sx q[0];
rz(-2.1779968) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.684217) q[2];
sx q[2];
rz(-1.3807502) q[2];
sx q[2];
rz(2.6148877) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0848799) q[1];
sx q[1];
rz(-0.093357714) q[1];
sx q[1];
rz(2.8356524) q[1];
rz(-3.0902052) q[3];
sx q[3];
rz(-0.87552658) q[3];
sx q[3];
rz(2.0688847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5588348) q[2];
sx q[2];
rz(-1.02966) q[2];
sx q[2];
rz(0.47490698) q[2];
rz(2.9389935) q[3];
sx q[3];
rz(-0.83297268) q[3];
sx q[3];
rz(1.1420265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.1982034) q[0];
sx q[0];
rz(-0.13309637) q[0];
sx q[0];
rz(2.8357847) q[0];
rz(-0.43931475) q[1];
sx q[1];
rz(-0.82249928) q[1];
sx q[1];
rz(1.3154715) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1409802) q[0];
sx q[0];
rz(-2.4580553) q[0];
sx q[0];
rz(-0.017869259) q[0];
rz(-pi) q[1];
rz(-1.1313296) q[2];
sx q[2];
rz(-0.60053289) q[2];
sx q[2];
rz(3.0366922) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2478359) q[1];
sx q[1];
rz(-0.79036056) q[1];
sx q[1];
rz(1.4001346) q[1];
rz(-pi) q[2];
rz(0.99268861) q[3];
sx q[3];
rz(-0.6738014) q[3];
sx q[3];
rz(2.0836308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4149912) q[2];
sx q[2];
rz(-1.5433658) q[2];
sx q[2];
rz(2.7169363) q[2];
rz(-0.031938227) q[3];
sx q[3];
rz(-1.6274803) q[3];
sx q[3];
rz(-0.64938515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53751078) q[0];
sx q[0];
rz(-0.6518971) q[0];
sx q[0];
rz(-2.5867468) q[0];
rz(2.8765053) q[1];
sx q[1];
rz(-1.4684497) q[1];
sx q[1];
rz(-1.9404985) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4557649) q[0];
sx q[0];
rz(-1.0170833) q[0];
sx q[0];
rz(2.035805) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3987819) q[2];
sx q[2];
rz(-1.2988678) q[2];
sx q[2];
rz(-3.0957785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9399413) q[1];
sx q[1];
rz(-1.2229511) q[1];
sx q[1];
rz(-1.6099206) q[1];
x q[2];
rz(3.115377) q[3];
sx q[3];
rz(-0.75002065) q[3];
sx q[3];
rz(1.2312082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.889664) q[2];
sx q[2];
rz(-2.1519075) q[2];
sx q[2];
rz(-2.2970693) q[2];
rz(1.1318413) q[3];
sx q[3];
rz(-2.2364538) q[3];
sx q[3];
rz(-1.7822781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1691386) q[0];
sx q[0];
rz(-2.7511981) q[0];
sx q[0];
rz(-1.9688695) q[0];
rz(2.4523465) q[1];
sx q[1];
rz(-1.0968364) q[1];
sx q[1];
rz(-2.8776339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64484064) q[0];
sx q[0];
rz(-0.94366108) q[0];
sx q[0];
rz(-2.1261313) q[0];
rz(-2.1446635) q[2];
sx q[2];
rz(-2.3341114) q[2];
sx q[2];
rz(-1.6107744) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3766438) q[1];
sx q[1];
rz(-1.3851435) q[1];
sx q[1];
rz(1.8913881) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92688074) q[3];
sx q[3];
rz(-1.4729233) q[3];
sx q[3];
rz(0.23668469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3802203) q[2];
sx q[2];
rz(-1.994588) q[2];
sx q[2];
rz(2.0210361) q[2];
rz(-1.606696) q[3];
sx q[3];
rz(-0.94529072) q[3];
sx q[3];
rz(-1.631261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(3.0705538) q[0];
sx q[0];
rz(-0.65503913) q[0];
sx q[0];
rz(-0.10792637) q[0];
rz(-2.4212266) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(2.2137523) q[2];
sx q[2];
rz(-0.8561047) q[2];
sx q[2];
rz(2.718953) q[2];
rz(1.4955487) q[3];
sx q[3];
rz(-0.9552707) q[3];
sx q[3];
rz(0.09494119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
