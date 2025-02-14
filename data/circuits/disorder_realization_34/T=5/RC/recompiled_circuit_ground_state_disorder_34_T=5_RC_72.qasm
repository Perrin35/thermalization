OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(-0.11556927) q[1];
sx q[1];
rz(2.5791383) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90949517) q[0];
sx q[0];
rz(-1.2065071) q[0];
sx q[0];
rz(-3.0845736) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8461343) q[2];
sx q[2];
rz(-0.52342969) q[2];
sx q[2];
rz(-1.1688423) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3823279) q[1];
sx q[1];
rz(-0.60324962) q[1];
sx q[1];
rz(-0.77253491) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1645557) q[3];
sx q[3];
rz(-1.2066168) q[3];
sx q[3];
rz(1.4573569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2354551) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(0.79929024) q[2];
rz(-2.6702787) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(2.2057064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039463194) q[0];
sx q[0];
rz(-1.8164182) q[0];
sx q[0];
rz(-1.348173) q[0];
rz(1.7680291) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.9893533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9236167) q[0];
sx q[0];
rz(-2.2647175) q[0];
sx q[0];
rz(2.7569735) q[0];
rz(0.67518465) q[2];
sx q[2];
rz(-2.084888) q[2];
sx q[2];
rz(2.6509283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79346953) q[1];
sx q[1];
rz(-1.7382227) q[1];
sx q[1];
rz(0.69116418) q[1];
rz(-0.51898414) q[3];
sx q[3];
rz(-2.0618084) q[3];
sx q[3];
rz(1.4903013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79933244) q[2];
sx q[2];
rz(-0.41877425) q[2];
sx q[2];
rz(2.2104134) q[2];
rz(-0.10654199) q[3];
sx q[3];
rz(-2.0276766) q[3];
sx q[3];
rz(-2.54134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0837285) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(-0.51783836) q[0];
rz(2.2528516) q[1];
sx q[1];
rz(-0.70777142) q[1];
sx q[1];
rz(-0.4471561) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5260552) q[0];
sx q[0];
rz(-1.3096022) q[0];
sx q[0];
rz(1.5047856) q[0];
rz(1.0165527) q[2];
sx q[2];
rz(-2.2187761) q[2];
sx q[2];
rz(-0.78314834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67990002) q[1];
sx q[1];
rz(-2.5352806) q[1];
sx q[1];
rz(-1.8099422) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6154061) q[3];
sx q[3];
rz(-1.8136214) q[3];
sx q[3];
rz(3.1265885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7670224) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(-0.53263295) q[2];
rz(-2.8940708) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(-1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6240876) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(-0.10661539) q[0];
rz(0.34635776) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(-1.1603629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713523) q[0];
sx q[0];
rz(-1.3246857) q[0];
sx q[0];
rz(-1.3918124) q[0];
rz(-0.82616634) q[2];
sx q[2];
rz(-1.0001422) q[2];
sx q[2];
rz(-2.014239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9306158) q[1];
sx q[1];
rz(-1.2158356) q[1];
sx q[1];
rz(2.8205754) q[1];
rz(-2.0959693) q[3];
sx q[3];
rz(-1.9672638) q[3];
sx q[3];
rz(0.1016271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13566636) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(0.076210991) q[2];
rz(0.58586079) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(-1.3425672) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0214486) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(2.7914877) q[0];
rz(0.95589751) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(1.9116481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1317539) q[0];
sx q[0];
rz(-1.3309877) q[0];
sx q[0];
rz(1.1935704) q[0];
rz(1.7857871) q[2];
sx q[2];
rz(-0.16915288) q[2];
sx q[2];
rz(-2.7403938) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0240682) q[1];
sx q[1];
rz(-1.9000016) q[1];
sx q[1];
rz(2.9419521) q[1];
rz(2.941972) q[3];
sx q[3];
rz(-2.1484882) q[3];
sx q[3];
rz(-0.22838372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.6492856) q[2];
rz(2.7367075) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(0.83612061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(0.58240044) q[0];
rz(-2.4549585) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(-1.2581717) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95520407) q[0];
sx q[0];
rz(-2.0058647) q[0];
sx q[0];
rz(-3.0032936) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.011884) q[2];
sx q[2];
rz(-0.61043533) q[2];
sx q[2];
rz(-1.287078) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88731474) q[1];
sx q[1];
rz(-1.6266903) q[1];
sx q[1];
rz(-0.069737597) q[1];
rz(-pi) q[2];
rz(1.889713) q[3];
sx q[3];
rz(-2.5237759) q[3];
sx q[3];
rz(-1.8998673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.9932237) q[2];
sx q[2];
rz(-1.0180391) q[2];
rz(-0.24614075) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(-1.6233981) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4356284) q[0];
sx q[0];
rz(-0.80393296) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(-2.9980581) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-0.20763436) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5136826) q[0];
sx q[0];
rz(-1.2259036) q[0];
sx q[0];
rz(2.3570127) q[0];
rz(-1.064038) q[2];
sx q[2];
rz(-0.64219507) q[2];
sx q[2];
rz(-0.80201496) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5112111) q[1];
sx q[1];
rz(-2.2464754) q[1];
sx q[1];
rz(2.8065794) q[1];
rz(-2.8945699) q[3];
sx q[3];
rz(-2.7678856) q[3];
sx q[3];
rz(1.5242239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2988854) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(-2.5353954) q[2];
rz(0.68743622) q[3];
sx q[3];
rz(-2.0076553) q[3];
sx q[3];
rz(2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48924482) q[0];
sx q[0];
rz(-3.1135961) q[0];
sx q[0];
rz(-1.0580753) q[0];
rz(3.0319013) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.6995957) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3198411) q[0];
sx q[0];
rz(-1.1761001) q[0];
sx q[0];
rz(1.9072397) q[0];
rz(1.2953561) q[2];
sx q[2];
rz(-2.0552962) q[2];
sx q[2];
rz(2.7000008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.86275872) q[1];
sx q[1];
rz(-2.2087084) q[1];
sx q[1];
rz(2.5274123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85956802) q[3];
sx q[3];
rz(-2.7894434) q[3];
sx q[3];
rz(1.5052049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66120061) q[2];
sx q[2];
rz(-2.9471687) q[2];
sx q[2];
rz(-1.2083758) q[2];
rz(2.4783573) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709568) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(-3.048625) q[0];
rz(-1.8611106) q[1];
sx q[1];
rz(-0.72851506) q[1];
sx q[1];
rz(-3.0063937) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8672392) q[0];
sx q[0];
rz(-1.5177736) q[0];
sx q[0];
rz(-1.6163325) q[0];
rz(-1.7241571) q[2];
sx q[2];
rz(-1.8259154) q[2];
sx q[2];
rz(1.1197787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0307903) q[1];
sx q[1];
rz(-0.25486481) q[1];
sx q[1];
rz(0.084666208) q[1];
rz(-pi) q[2];
rz(-0.074196176) q[3];
sx q[3];
rz(-1.7230986) q[3];
sx q[3];
rz(-1.1049113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4325503) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(-2.3700628) q[2];
rz(2.6500474) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(-2.5468723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58525697) q[0];
sx q[0];
rz(-0.62194967) q[0];
sx q[0];
rz(-1.1248032) q[0];
rz(1.8428165) q[1];
sx q[1];
rz(-2.5262084) q[1];
sx q[1];
rz(-0.76464701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9287531) q[0];
sx q[0];
rz(-3.0055586) q[0];
sx q[0];
rz(-2.4555951) q[0];
x q[1];
rz(-0.9773639) q[2];
sx q[2];
rz(-2.368401) q[2];
sx q[2];
rz(0.14466454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5717126) q[1];
sx q[1];
rz(-0.24001828) q[1];
sx q[1];
rz(-1.0102788) q[1];
rz(-pi) q[2];
rz(-1.9493136) q[3];
sx q[3];
rz(-2.0928185) q[3];
sx q[3];
rz(1.4498364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7654968) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(-2.9343904) q[2];
rz(0.97992212) q[3];
sx q[3];
rz(-2.0082974) q[3];
sx q[3];
rz(2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021066396) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(2.6407241) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(1.6899213) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(-1.408314) q[3];
sx q[3];
rz(-1.8058895) q[3];
sx q[3];
rz(-1.4598082) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
