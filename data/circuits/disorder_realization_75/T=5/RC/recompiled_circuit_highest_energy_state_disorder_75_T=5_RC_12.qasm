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
rz(-1.8259814) q[0];
sx q[0];
rz(-2.313518) q[0];
sx q[0];
rz(0.84870422) q[0];
rz(-1.6500213) q[1];
sx q[1];
rz(-2.4529011) q[1];
sx q[1];
rz(-0.42660776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72683452) q[0];
sx q[0];
rz(-1.4953185) q[0];
sx q[0];
rz(3.0407136) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0399917) q[2];
sx q[2];
rz(-1.4590605) q[2];
sx q[2];
rz(-0.32411175) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.035543899) q[1];
sx q[1];
rz(-1.1263761) q[1];
sx q[1];
rz(-1.6515093) q[1];
rz(-pi) q[2];
rz(0.95152609) q[3];
sx q[3];
rz(-2.2689156) q[3];
sx q[3];
rz(-1.2941456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9244869) q[2];
sx q[2];
rz(-1.7542398) q[2];
sx q[2];
rz(-0.039487751) q[2];
rz(2.110179) q[3];
sx q[3];
rz(-2.1125427) q[3];
sx q[3];
rz(1.2379117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51204387) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(-1.398983) q[0];
rz(1.5139187) q[1];
sx q[1];
rz(-0.84295034) q[1];
sx q[1];
rz(1.4167851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48081452) q[0];
sx q[0];
rz(-2.3153458) q[0];
sx q[0];
rz(-1.2673402) q[0];
rz(-pi) q[1];
rz(-1.8292839) q[2];
sx q[2];
rz(-0.80855364) q[2];
sx q[2];
rz(-2.0214579) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21525225) q[1];
sx q[1];
rz(-1.4344341) q[1];
sx q[1];
rz(0.3991613) q[1];
rz(1.5228769) q[3];
sx q[3];
rz(-1.9780434) q[3];
sx q[3];
rz(1.6919235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9713355) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(-0.7106759) q[2];
rz(-0.014287861) q[3];
sx q[3];
rz(-0.24737869) q[3];
sx q[3];
rz(1.3361196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9209442) q[0];
sx q[0];
rz(-2.1823688) q[0];
sx q[0];
rz(0.4162108) q[0];
rz(-1.8572218) q[1];
sx q[1];
rz(-1.4764079) q[1];
sx q[1];
rz(-0.31825569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6775433) q[0];
sx q[0];
rz(-0.47992063) q[0];
sx q[0];
rz(2.2024758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4893583) q[2];
sx q[2];
rz(-0.95931897) q[2];
sx q[2];
rz(0.91064168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7075338) q[1];
sx q[1];
rz(-2.4726082) q[1];
sx q[1];
rz(2.9218052) q[1];
rz(0.89983799) q[3];
sx q[3];
rz(-2.7551015) q[3];
sx q[3];
rz(-0.46903601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2661813) q[2];
sx q[2];
rz(-0.77989945) q[2];
sx q[2];
rz(1.2904588) q[2];
rz(0.053704638) q[3];
sx q[3];
rz(-1.8219681) q[3];
sx q[3];
rz(1.7558338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5829492) q[0];
sx q[0];
rz(-2.4531893) q[0];
sx q[0];
rz(2.131856) q[0];
rz(-0.91521493) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(-2.7161652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260507) q[0];
sx q[0];
rz(-0.73662355) q[0];
sx q[0];
rz(-1.3213586) q[0];
rz(0.82773005) q[2];
sx q[2];
rz(-1.6249738) q[2];
sx q[2];
rz(-0.82967134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2293416) q[1];
sx q[1];
rz(-2.3315363) q[1];
sx q[1];
rz(0.99397578) q[1];
x q[2];
rz(0.0056267319) q[3];
sx q[3];
rz(-1.4079908) q[3];
sx q[3];
rz(-0.73782792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2414744) q[2];
sx q[2];
rz(-1.5510617) q[2];
sx q[2];
rz(3.0305064) q[2];
rz(0.19242081) q[3];
sx q[3];
rz(-0.40168732) q[3];
sx q[3];
rz(-0.69453159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38839328) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(2.2764192) q[0];
rz(0.49258891) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(1.385484) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139112) q[0];
sx q[0];
rz(-1.4340785) q[0];
sx q[0];
rz(0.78929957) q[0];
rz(-pi) q[1];
rz(1.1329012) q[2];
sx q[2];
rz(-0.81134331) q[2];
sx q[2];
rz(-2.0591625) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4069479) q[1];
sx q[1];
rz(-2.4666767) q[1];
sx q[1];
rz(-0.65924834) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4597857) q[3];
sx q[3];
rz(-2.1958133) q[3];
sx q[3];
rz(-0.30407076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8869141) q[2];
sx q[2];
rz(-0.46923894) q[2];
sx q[2];
rz(-0.706642) q[2];
rz(-2.8142269) q[3];
sx q[3];
rz(-1.8492161) q[3];
sx q[3];
rz(-2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8364328) q[0];
sx q[0];
rz(-1.3494455) q[0];
sx q[0];
rz(-2.7850372) q[0];
rz(-0.56218475) q[1];
sx q[1];
rz(-2.0088582) q[1];
sx q[1];
rz(-2.1551989) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0537181) q[0];
sx q[0];
rz(-1.0759504) q[0];
sx q[0];
rz(2.3943564) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7964727) q[2];
sx q[2];
rz(-1.9202931) q[2];
sx q[2];
rz(-1.0687168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2960259) q[1];
sx q[1];
rz(-0.81261501) q[1];
sx q[1];
rz(1.5350828) q[1];
rz(-pi) q[2];
rz(1.3273193) q[3];
sx q[3];
rz(-1.7469353) q[3];
sx q[3];
rz(1.0276664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76546136) q[2];
sx q[2];
rz(-2.4262846) q[2];
sx q[2];
rz(2.4410655) q[2];
rz(2.1721407) q[3];
sx q[3];
rz(-2.1078096) q[3];
sx q[3];
rz(-0.9303003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9483865) q[0];
sx q[0];
rz(-2.8128615) q[0];
sx q[0];
rz(-0.1513924) q[0];
rz(-1.6311749) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(-1.4600533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23510012) q[0];
sx q[0];
rz(-2.0564046) q[0];
sx q[0];
rz(-0.67295815) q[0];
rz(1.0437772) q[2];
sx q[2];
rz(-1.2536546) q[2];
sx q[2];
rz(-0.70995599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20138473) q[1];
sx q[1];
rz(-1.0321069) q[1];
sx q[1];
rz(-0.45392764) q[1];
rz(-pi) q[2];
rz(-2.9956508) q[3];
sx q[3];
rz(-1.3781007) q[3];
sx q[3];
rz(-1.3789498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90525642) q[2];
sx q[2];
rz(-0.6826123) q[2];
sx q[2];
rz(-2.7899817) q[2];
rz(-0.44175276) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(2.9898804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041895954) q[0];
sx q[0];
rz(-1.8459039) q[0];
sx q[0];
rz(-1.9805441) q[0];
rz(-0.64758045) q[1];
sx q[1];
rz(-1.7627962) q[1];
sx q[1];
rz(0.82239282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8063342) q[0];
sx q[0];
rz(-1.0653235) q[0];
sx q[0];
rz(-0.45787207) q[0];
rz(-pi) q[1];
rz(2.0089125) q[2];
sx q[2];
rz(-2.2575976) q[2];
sx q[2];
rz(-0.39114726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5288189) q[1];
sx q[1];
rz(-0.83136035) q[1];
sx q[1];
rz(-2.7465613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3052182) q[3];
sx q[3];
rz(-0.55787239) q[3];
sx q[3];
rz(-1.316586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9245727) q[2];
sx q[2];
rz(-1.1447039) q[2];
sx q[2];
rz(-1.3059957) q[2];
rz(-2.4235587) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(-3.0927299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7056535) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(0.017024592) q[0];
rz(-1.1964993) q[1];
sx q[1];
rz(-0.39533177) q[1];
sx q[1];
rz(-2.8065525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8117479) q[0];
sx q[0];
rz(-0.56808725) q[0];
sx q[0];
rz(-1.5519867) q[0];
rz(-pi) q[1];
rz(0.64249383) q[2];
sx q[2];
rz(-1.9064925) q[2];
sx q[2];
rz(3.1106126) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2372053) q[1];
sx q[1];
rz(-0.6977302) q[1];
sx q[1];
rz(-1.8664788) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0552956) q[3];
sx q[3];
rz(-2.006975) q[3];
sx q[3];
rz(2.0791813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42387858) q[2];
sx q[2];
rz(-0.7462036) q[2];
sx q[2];
rz(-2.8625873) q[2];
rz(1.772607) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(-2.2122808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245215) q[0];
sx q[0];
rz(-0.14685024) q[0];
sx q[0];
rz(-0.76706925) q[0];
rz(-2.7686139) q[1];
sx q[1];
rz(-1.4992799) q[1];
sx q[1];
rz(-2.9708718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6304008) q[0];
sx q[0];
rz(-0.83626473) q[0];
sx q[0];
rz(-1.6096576) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35358785) q[2];
sx q[2];
rz(-2.7025526) q[2];
sx q[2];
rz(-0.072364213) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77826277) q[1];
sx q[1];
rz(-2.7365541) q[1];
sx q[1];
rz(-2.2447849) q[1];
x q[2];
rz(-1.7307698) q[3];
sx q[3];
rz(-2.0334646) q[3];
sx q[3];
rz(1.8335952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3501222) q[2];
sx q[2];
rz(-2.5779724) q[2];
sx q[2];
rz(-0.0027837022) q[2];
rz(2.2400098) q[3];
sx q[3];
rz(-2.1193347) q[3];
sx q[3];
rz(-2.3367052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0881385) q[0];
sx q[0];
rz(-0.32129856) q[0];
sx q[0];
rz(-1.970485) q[0];
rz(-2.3114655) q[1];
sx q[1];
rz(-1.8142038) q[1];
sx q[1];
rz(1.289191) q[1];
rz(-2.1813761) q[2];
sx q[2];
rz(-1.1398906) q[2];
sx q[2];
rz(0.7791145) q[2];
rz(-2.3385901) q[3];
sx q[3];
rz(-1.0754801) q[3];
sx q[3];
rz(0.51343244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
