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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(2.8454236) q[1];
sx q[1];
rz(12.256395) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5910019) q[0];
sx q[0];
rz(-2.7059116) q[0];
sx q[0];
rz(0.26623078) q[0];
rz(-pi) q[1];
x q[1];
rz(1.593179) q[2];
sx q[2];
rz(-1.7760008) q[2];
sx q[2];
rz(-0.6485282) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1811318) q[1];
sx q[1];
rz(-0.38552654) q[1];
sx q[1];
rz(-2.1147644) q[1];
rz(1.3499851) q[3];
sx q[3];
rz(-2.3901) q[3];
sx q[3];
rz(-1.8210851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4787204) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(3.0459246) q[2];
rz(-3.0293448) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(-1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9005168) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(2.526793) q[0];
rz(-0.88848937) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(-2.6420171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083778) q[0];
sx q[0];
rz(-1.3691804) q[0];
sx q[0];
rz(-1.2382384) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26878727) q[2];
sx q[2];
rz(-0.54020451) q[2];
sx q[2];
rz(0.42551431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.080729) q[1];
sx q[1];
rz(-1.7434967) q[1];
sx q[1];
rz(-0.91367803) q[1];
x q[2];
rz(-2.1533215) q[3];
sx q[3];
rz(-1.2856021) q[3];
sx q[3];
rz(0.56166431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(-0.020817967) q[2];
rz(1.2402395) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(1.9564691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6388539) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(0.33367208) q[0];
rz(-0.73792136) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(-0.12942448) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724587) q[0];
sx q[0];
rz(-0.6967623) q[0];
sx q[0];
rz(1.8201226) q[0];
rz(2.1652075) q[2];
sx q[2];
rz(-2.1529641) q[2];
sx q[2];
rz(-0.79603031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5065352) q[1];
sx q[1];
rz(-1.660907) q[1];
sx q[1];
rz(-1.6930626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3593216) q[3];
sx q[3];
rz(-1.9214298) q[3];
sx q[3];
rz(2.3334437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0176598) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(-3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081886) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(0.56513894) q[0];
rz(0.42916974) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(-2.3133004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.013731) q[0];
sx q[0];
rz(-0.7452226) q[0];
sx q[0];
rz(-3.0237314) q[0];
rz(-2.056869) q[2];
sx q[2];
rz(-1.8258494) q[2];
sx q[2];
rz(-1.4866831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.284621) q[1];
sx q[1];
rz(-1.0880252) q[1];
sx q[1];
rz(1.8716145) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8935562) q[3];
sx q[3];
rz(-2.0741077) q[3];
sx q[3];
rz(2.4865884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4565178) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(2.2020014) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(-0.00069869839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(1.8121207) q[0];
rz(-1.9895408) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(-0.00024814127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329929) q[0];
sx q[0];
rz(-1.8155351) q[0];
sx q[0];
rz(-0.18969638) q[0];
rz(2.1049961) q[2];
sx q[2];
rz(-1.4956258) q[2];
sx q[2];
rz(0.49803621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9226886) q[1];
sx q[1];
rz(-2.5210025) q[1];
sx q[1];
rz(0.94938486) q[1];
rz(-2.4100634) q[3];
sx q[3];
rz(-1.5168744) q[3];
sx q[3];
rz(1.9526854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1382711) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(3.1357583) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(-2.0148923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(-1.6538612) q[0];
rz(0.030390175) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(-0.94246513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4543991) q[0];
sx q[0];
rz(-1.8158578) q[0];
sx q[0];
rz(1.7646199) q[0];
rz(0.29877383) q[2];
sx q[2];
rz(-1.3232097) q[2];
sx q[2];
rz(-2.5556759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9760175) q[1];
sx q[1];
rz(-1.2846795) q[1];
sx q[1];
rz(2.0598434) q[1];
rz(-pi) q[2];
rz(-1.0932572) q[3];
sx q[3];
rz(-2.6160588) q[3];
sx q[3];
rz(-0.20697396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59948644) q[2];
sx q[2];
rz(-2.3607871) q[2];
sx q[2];
rz(1.3055118) q[2];
rz(2.5944338) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-0.80593306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9611573) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(-3.0410774) q[0];
rz(0.48249498) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(2.2844792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470951) q[0];
sx q[0];
rz(-2.408354) q[0];
sx q[0];
rz(0.49219699) q[0];
rz(-pi) q[1];
rz(-1.6202507) q[2];
sx q[2];
rz(-0.92453814) q[2];
sx q[2];
rz(-2.6924999) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7776612) q[1];
sx q[1];
rz(-0.45501935) q[1];
sx q[1];
rz(0.48952405) q[1];
x q[2];
rz(-1.1825652) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(1.5807815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7414005) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(0.3024438) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7804467) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(-3.1106023) q[0];
rz(2.567645) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(1.2773638) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63420682) q[0];
sx q[0];
rz(-1.9045826) q[0];
sx q[0];
rz(-2.7884363) q[0];
rz(2.0691643) q[2];
sx q[2];
rz(-1.3565836) q[2];
sx q[2];
rz(-1.6139849) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.053169202) q[1];
sx q[1];
rz(-0.49234566) q[1];
sx q[1];
rz(1.7264992) q[1];
rz(-1.7788309) q[3];
sx q[3];
rz(-2.8505508) q[3];
sx q[3];
rz(0.42359776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.10451) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(0.48745298) q[2];
rz(2.4449352) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(-0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0697407) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(-2.790614) q[0];
rz(-0.17414302) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(1.7399656) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6232672) q[0];
sx q[0];
rz(-1.4484157) q[0];
sx q[0];
rz(2.8920679) q[0];
x q[1];
rz(2.6675111) q[2];
sx q[2];
rz(-0.80898058) q[2];
sx q[2];
rz(-2.7681729) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3199215) q[1];
sx q[1];
rz(-2.955759) q[1];
sx q[1];
rz(-1.7286848) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97934874) q[3];
sx q[3];
rz(-2.2205847) q[3];
sx q[3];
rz(-1.5486919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(2.5049211) q[2];
rz(-1.1104256) q[3];
sx q[3];
rz(-1.5444376) q[3];
sx q[3];
rz(0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(3.1029347) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-2.0775332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87133133) q[0];
sx q[0];
rz(-3.1360021) q[0];
sx q[0];
rz(-2.5020775) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27711192) q[2];
sx q[2];
rz(-1.849035) q[2];
sx q[2];
rz(-2.3269175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0088996) q[1];
sx q[1];
rz(-1.6311967) q[1];
sx q[1];
rz(1.6113043) q[1];
rz(-0.45952103) q[3];
sx q[3];
rz(-1.5395313) q[3];
sx q[3];
rz(0.36324901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(-0.074020298) q[2];
rz(-2.7600539) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(2.7874302) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030180177) q[0];
sx q[0];
rz(-1.8288061) q[0];
sx q[0];
rz(-0.70963138) q[0];
rz(2.5297655) q[1];
sx q[1];
rz(-0.71129967) q[1];
sx q[1];
rz(1.5536972) q[1];
rz(-1.189497) q[2];
sx q[2];
rz(-2.0443889) q[2];
sx q[2];
rz(-2.2257795) q[2];
rz(-2.6388219) q[3];
sx q[3];
rz(-0.81450653) q[3];
sx q[3];
rz(-0.52342879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
