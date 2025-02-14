OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.81808972) q[0];
sx q[0];
rz(-0.40677318) q[0];
sx q[0];
rz(0.50954252) q[0];
rz(2.8264363) q[1];
sx q[1];
rz(-0.1935614) q[1];
sx q[1];
rz(-1.0055746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5339342) q[0];
sx q[0];
rz(-0.21398057) q[0];
sx q[0];
rz(-2.7817552) q[0];
rz(2.2174913) q[2];
sx q[2];
rz(-1.8515203) q[2];
sx q[2];
rz(1.4626056) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.352658) q[1];
sx q[1];
rz(-2.1772258) q[1];
sx q[1];
rz(-2.2878617) q[1];
rz(-pi) q[2];
rz(-0.96989743) q[3];
sx q[3];
rz(-0.88774046) q[3];
sx q[3];
rz(-1.1163464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2931508) q[2];
sx q[2];
rz(-1.8895431) q[2];
sx q[2];
rz(-2.0330009) q[2];
rz(1.1340002) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(3.0602684) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9950681) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(-0.31196892) q[0];
rz(-0.23090714) q[1];
sx q[1];
rz(-2.0377908) q[1];
sx q[1];
rz(-1.1280967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474079) q[0];
sx q[0];
rz(-2.0387406) q[0];
sx q[0];
rz(-2.9293961) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0064193) q[2];
sx q[2];
rz(-0.95016236) q[2];
sx q[2];
rz(-1.0301318) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.038874168) q[1];
sx q[1];
rz(-0.46157122) q[1];
sx q[1];
rz(0.25074236) q[1];
rz(2.6873782) q[3];
sx q[3];
rz(-0.81791211) q[3];
sx q[3];
rz(1.0354158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6850623) q[2];
sx q[2];
rz(-1.637746) q[2];
sx q[2];
rz(-3.077363) q[2];
rz(1.8188933) q[3];
sx q[3];
rz(-0.6367681) q[3];
sx q[3];
rz(0.88671154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927758) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(1.2491666) q[0];
rz(2.8104172) q[1];
sx q[1];
rz(-2.0024029) q[1];
sx q[1];
rz(-1.9452728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0977444) q[0];
sx q[0];
rz(-0.00034887408) q[0];
sx q[0];
rz(-2.1972138) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.021042391) q[2];
sx q[2];
rz(-2.8940563) q[2];
sx q[2];
rz(-0.68138441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0374245) q[1];
sx q[1];
rz(-1.0598039) q[1];
sx q[1];
rz(-2.4591706) q[1];
rz(0.75342859) q[3];
sx q[3];
rz(-1.2810858) q[3];
sx q[3];
rz(-0.10623194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5195878) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(2.2288442) q[2];
rz(-1.7031472) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(-1.8057757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66252935) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(2.4712439) q[0];
rz(-0.66728512) q[1];
sx q[1];
rz(-1.0083464) q[1];
sx q[1];
rz(2.537312) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42603618) q[0];
sx q[0];
rz(-0.73764387) q[0];
sx q[0];
rz(3.0771653) q[0];
x q[1];
rz(1.4154345) q[2];
sx q[2];
rz(-1.1929907) q[2];
sx q[2];
rz(0.31766674) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6091023) q[1];
sx q[1];
rz(-2.8081642) q[1];
sx q[1];
rz(-0.87581265) q[1];
rz(-pi) q[2];
rz(1.8205582) q[3];
sx q[3];
rz(-0.74640025) q[3];
sx q[3];
rz(1.6370893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1112572) q[2];
sx q[2];
rz(-1.5735441) q[2];
sx q[2];
rz(0.37720171) q[2];
rz(0.12353573) q[3];
sx q[3];
rz(-1.433452) q[3];
sx q[3];
rz(-1.7326573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443026) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(-2.8422728) q[0];
rz(2.0206644) q[1];
sx q[1];
rz(-1.2052373) q[1];
sx q[1];
rz(2.1790806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4967921) q[0];
sx q[0];
rz(-1.776665) q[0];
sx q[0];
rz(-1.0394761) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6766432) q[2];
sx q[2];
rz(-2.1735149) q[2];
sx q[2];
rz(0.51747396) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0244151) q[1];
sx q[1];
rz(-1.0512253) q[1];
sx q[1];
rz(-2.5409798) q[1];
rz(-pi) q[2];
rz(0.59904091) q[3];
sx q[3];
rz(-2.0192462) q[3];
sx q[3];
rz(-0.34235172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0472497) q[2];
sx q[2];
rz(-2.3361358) q[2];
sx q[2];
rz(2.1979525) q[2];
rz(-1.0654248) q[3];
sx q[3];
rz(-1.3506972) q[3];
sx q[3];
rz(1.0984727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0672673) q[0];
sx q[0];
rz(-1.4370947) q[0];
sx q[0];
rz(-0.0083010439) q[0];
rz(2.6240194) q[1];
sx q[1];
rz(-2.4868496) q[1];
sx q[1];
rz(1.5483206) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8682713) q[0];
sx q[0];
rz(-1.5467321) q[0];
sx q[0];
rz(1.9217743) q[0];
rz(-0.48262934) q[2];
sx q[2];
rz(-2.5833774) q[2];
sx q[2];
rz(-1.7354154) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.949985) q[1];
sx q[1];
rz(-1.9298501) q[1];
sx q[1];
rz(2.784347) q[1];
x q[2];
rz(0.41554873) q[3];
sx q[3];
rz(-0.99667785) q[3];
sx q[3];
rz(-0.94438636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8325309) q[2];
sx q[2];
rz(-1.6522202) q[2];
sx q[2];
rz(-1.8050516) q[2];
rz(-3.0858223) q[3];
sx q[3];
rz(-2.1499108) q[3];
sx q[3];
rz(-2.706004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.5830773) q[0];
sx q[0];
rz(-2.6641088) q[0];
sx q[0];
rz(0.68786311) q[0];
rz(-1.4631924) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(-0.24737839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8198422) q[0];
sx q[0];
rz(-0.30748707) q[0];
sx q[0];
rz(2.5006258) q[0];
rz(-pi) q[1];
rz(0.74784634) q[2];
sx q[2];
rz(-0.81163844) q[2];
sx q[2];
rz(2.5876665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.069889594) q[1];
sx q[1];
rz(-2.1462198) q[1];
sx q[1];
rz(-2.6399122) q[1];
x q[2];
rz(2.0603397) q[3];
sx q[3];
rz(-2.8791028) q[3];
sx q[3];
rz(2.9406527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4443724) q[2];
sx q[2];
rz(-1.3345382) q[2];
sx q[2];
rz(-2.7174301) q[2];
rz(1.0423202) q[3];
sx q[3];
rz(-0.76516953) q[3];
sx q[3];
rz(2.585129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0031849) q[0];
sx q[0];
rz(-0.98588949) q[0];
sx q[0];
rz(2.5307122) q[0];
rz(0.25018397) q[1];
sx q[1];
rz(-1.7411722) q[1];
sx q[1];
rz(0.47666034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.902001) q[0];
sx q[0];
rz(-1.6098813) q[0];
sx q[0];
rz(-2.0530967) q[0];
x q[1];
rz(-1.119461) q[2];
sx q[2];
rz(-2.9381466) q[2];
sx q[2];
rz(2.1687393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6338773) q[1];
sx q[1];
rz(-2.4832279) q[1];
sx q[1];
rz(-2.4509199) q[1];
rz(-pi) q[2];
rz(-2.1564756) q[3];
sx q[3];
rz(-1.9352311) q[3];
sx q[3];
rz(1.5353312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1499947) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(-0.67982137) q[2];
rz(-0.46197915) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(-1.5010887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944549) q[0];
sx q[0];
rz(-1.3752022) q[0];
sx q[0];
rz(0.15400259) q[0];
rz(0.18140659) q[1];
sx q[1];
rz(-1.0702952) q[1];
sx q[1];
rz(3.0214686) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2699566) q[0];
sx q[0];
rz(-1.9737195) q[0];
sx q[0];
rz(-0.064489207) q[0];
rz(1.3072877) q[2];
sx q[2];
rz(-1.3804091) q[2];
sx q[2];
rz(-2.9817977) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1644205) q[1];
sx q[1];
rz(-1.394632) q[1];
sx q[1];
rz(3.0422496) q[1];
x q[2];
rz(2.613853) q[3];
sx q[3];
rz(-1.6956009) q[3];
sx q[3];
rz(0.11858701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7141562) q[2];
sx q[2];
rz(-1.3849881) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(0.51042026) q[3];
sx q[3];
rz(-0.84158689) q[3];
sx q[3];
rz(-2.4418805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089461483) q[0];
sx q[0];
rz(-2.1928146) q[0];
sx q[0];
rz(-0.38598886) q[0];
rz(-1.5486859) q[1];
sx q[1];
rz(-2.6451151) q[1];
sx q[1];
rz(-1.5492424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-1.5423466) q[0];
sx q[0];
rz(1.7398119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1524197) q[2];
sx q[2];
rz(-2.1206822) q[2];
sx q[2];
rz(1.3370507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1116699) q[1];
sx q[1];
rz(-0.1055461) q[1];
sx q[1];
rz(1.8284069) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.066927197) q[3];
sx q[3];
rz(-1.4176344) q[3];
sx q[3];
rz(0.39439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3647032) q[2];
sx q[2];
rz(-2.2025755) q[2];
sx q[2];
rz(1.6949863) q[2];
rz(1.0363091) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(-3.115263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.885289) q[0];
sx q[0];
rz(-2.1623609) q[0];
sx q[0];
rz(1.4872861) q[0];
rz(-2.9215095) q[1];
sx q[1];
rz(-2.0368123) q[1];
sx q[1];
rz(-1.4727551) q[1];
rz(1.2575116) q[2];
sx q[2];
rz(-2.2844615) q[2];
sx q[2];
rz(1.5649038) q[2];
rz(-2.5840633) q[3];
sx q[3];
rz(-0.57303187) q[3];
sx q[3];
rz(-1.8106885) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
