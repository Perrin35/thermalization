OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17833248) q[0];
sx q[0];
rz(-1.4890716) q[0];
sx q[0];
rz(2.2464377) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(-1.4562343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.668894) q[0];
sx q[0];
rz(-1.5905252) q[0];
sx q[0];
rz(1.701645) q[0];
rz(-pi) q[1];
rz(0.99631359) q[2];
sx q[2];
rz(-0.62553863) q[2];
sx q[2];
rz(-1.1686981) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97988843) q[1];
sx q[1];
rz(-2.4474505) q[1];
sx q[1];
rz(0.23692268) q[1];
rz(-pi) q[2];
rz(2.8217444) q[3];
sx q[3];
rz(-1.6228798) q[3];
sx q[3];
rz(-0.44714123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(1.1520977) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(2.6706085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7171213) q[0];
sx q[0];
rz(-2.0524128) q[0];
sx q[0];
rz(0.061901285) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23439622) q[2];
sx q[2];
rz(-2.1550551) q[2];
sx q[2];
rz(2.5258979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4091332) q[1];
sx q[1];
rz(-1.0832548) q[1];
sx q[1];
rz(0.061846102) q[1];
rz(-pi) q[2];
rz(0.69555517) q[3];
sx q[3];
rz(-1.5091981) q[3];
sx q[3];
rz(2.0640645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-2.7022865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96734756) q[0];
sx q[0];
rz(-3.0558911) q[0];
sx q[0];
rz(1.9232737) q[0];
rz(1.6224242) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(2.1215631) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5369536) q[1];
sx q[1];
rz(-0.96452689) q[1];
sx q[1];
rz(-0.21142516) q[1];
rz(2.1268232) q[3];
sx q[3];
rz(-1.1960256) q[3];
sx q[3];
rz(2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(-0.20047323) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4865206) q[0];
sx q[0];
rz(-2.1694896) q[0];
sx q[0];
rz(0.18874164) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34747296) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(2.4328872) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17393915) q[1];
sx q[1];
rz(-1.1974918) q[1];
sx q[1];
rz(1.8147545) q[1];
rz(-pi) q[2];
rz(0.61769684) q[3];
sx q[3];
rz(-3.0129637) q[3];
sx q[3];
rz(-2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(0.74329174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76737228) q[0];
sx q[0];
rz(-1.7905856) q[0];
sx q[0];
rz(2.2348997) q[0];
x q[1];
rz(-1.4666124) q[2];
sx q[2];
rz(-1.0717234) q[2];
sx q[2];
rz(-0.68945976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6965461) q[1];
sx q[1];
rz(-1.3475218) q[1];
sx q[1];
rz(2.9162507) q[1];
rz(-pi) q[2];
rz(-0.27637847) q[3];
sx q[3];
rz(-0.82377269) q[3];
sx q[3];
rz(-0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(2.8732079) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(2.6348689) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-2.0862897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366093) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(-1.2654632) q[0];
rz(-pi) q[1];
rz(-2.6199117) q[2];
sx q[2];
rz(-0.98419596) q[2];
sx q[2];
rz(-0.66643836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.036990449) q[1];
sx q[1];
rz(-1.9676419) q[1];
sx q[1];
rz(-0.87581046) q[1];
rz(-pi) q[2];
rz(0.7430195) q[3];
sx q[3];
rz(-1.395874) q[3];
sx q[3];
rz(-0.35374853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(-1.4121217) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360146) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(1.3789603) q[0];
rz(-0.022106604) q[2];
sx q[2];
rz(-1.7709641) q[2];
sx q[2];
rz(0.18891639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7635203) q[1];
sx q[1];
rz(-2.1280648) q[1];
sx q[1];
rz(-2.4850363) q[1];
rz(-pi) q[2];
rz(2.1366828) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(0.79808455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(2.3708564) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.27281) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.9708721) q[0];
rz(-0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.2957113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9840178) q[0];
sx q[0];
rz(-0.57849738) q[0];
sx q[0];
rz(-2.5111141) q[0];
rz(-pi) q[1];
rz(0.92408085) q[2];
sx q[2];
rz(-2.7087822) q[2];
sx q[2];
rz(2.2518287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6756729) q[1];
sx q[1];
rz(-2.188176) q[1];
sx q[1];
rz(-0.27523756) q[1];
rz(3.1136884) q[3];
sx q[3];
rz(-1.1301665) q[3];
sx q[3];
rz(2.2925216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(2.4639159) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(-2.303404) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7064338) q[0];
sx q[0];
rz(-0.43276946) q[0];
sx q[0];
rz(-0.5612527) q[0];
rz(-1.6423177) q[2];
sx q[2];
rz(-2.3841249) q[2];
sx q[2];
rz(1.0345936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2182525) q[1];
sx q[1];
rz(-0.76949161) q[1];
sx q[1];
rz(-1.283265) q[1];
rz(-pi) q[2];
rz(2.8348654) q[3];
sx q[3];
rz(-0.2704276) q[3];
sx q[3];
rz(0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(2.1789815) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(0.094749711) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(2.9097897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0203665) q[0];
sx q[0];
rz(-1.0640642) q[0];
sx q[0];
rz(0.10911848) q[0];
rz(-pi) q[1];
rz(-0.78328697) q[2];
sx q[2];
rz(-0.33595339) q[2];
sx q[2];
rz(0.17620262) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17813645) q[1];
sx q[1];
rz(-1.8354715) q[1];
sx q[1];
rz(-1.4241649) q[1];
x q[2];
rz(-1.8292571) q[3];
sx q[3];
rz(-2.4371394) q[3];
sx q[3];
rz(-2.7750912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647472) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.8112524) q[2];
sx q[2];
rz(-1.8869055) q[2];
sx q[2];
rz(-1.4958924) q[2];
rz(-0.10235056) q[3];
sx q[3];
rz(-0.1889189) q[3];
sx q[3];
rz(1.2253996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
