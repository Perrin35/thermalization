OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3235029) q[0];
sx q[0];
rz(-2.7348195) q[0];
sx q[0];
rz(-0.50954252) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(-2.9480313) q[1];
sx q[1];
rz(-2.136018) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5339342) q[0];
sx q[0];
rz(-0.21398057) q[0];
sx q[0];
rz(-2.7817552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.124473) q[2];
sx q[2];
rz(-0.69689489) q[2];
sx q[2];
rz(0.45988894) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.352658) q[1];
sx q[1];
rz(-0.96436687) q[1];
sx q[1];
rz(-2.2878617) q[1];
rz(-pi) q[2];
rz(0.77862424) q[3];
sx q[3];
rz(-2.0247685) q[3];
sx q[3];
rz(-3.0954645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2931508) q[2];
sx q[2];
rz(-1.2520496) q[2];
sx q[2];
rz(-2.0330009) q[2];
rz(-2.0075924) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(3.0602684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9950681) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(0.31196892) q[0];
rz(0.23090714) q[1];
sx q[1];
rz(-1.1038019) q[1];
sx q[1];
rz(2.013496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9474079) q[0];
sx q[0];
rz(-2.0387406) q[0];
sx q[0];
rz(2.9293961) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6083469) q[2];
sx q[2];
rz(-2.4002078) q[2];
sx q[2];
rz(0.35517755) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7573537) q[1];
sx q[1];
rz(-1.4600672) q[1];
sx q[1];
rz(2.6925681) q[1];
x q[2];
rz(2.3771068) q[3];
sx q[3];
rz(-1.2448881) q[3];
sx q[3];
rz(-0.21316646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6850623) q[2];
sx q[2];
rz(-1.637746) q[2];
sx q[2];
rz(0.064229639) q[2];
rz(-1.3226994) q[3];
sx q[3];
rz(-0.6367681) q[3];
sx q[3];
rz(-2.2548811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927758) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(-1.2491666) q[0];
rz(2.8104172) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(-1.1963199) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0977444) q[0];
sx q[0];
rz(-3.1412438) q[0];
sx q[0];
rz(2.1972138) q[0];
x q[1];
rz(-1.565479) q[2];
sx q[2];
rz(-1.8182767) q[2];
sx q[2];
rz(2.4385045) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0374245) q[1];
sx q[1];
rz(-1.0598039) q[1];
sx q[1];
rz(-0.68242208) q[1];
rz(-pi) q[2];
rz(-0.41091856) q[3];
sx q[3];
rz(-2.344729) q[3];
sx q[3];
rz(-1.381402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5195878) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(-0.91274846) q[2];
rz(-1.7031472) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(1.335817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790633) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(-2.4712439) q[0];
rz(0.66728512) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(-0.60428062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9491228) q[0];
sx q[0];
rz(-1.6141103) q[0];
sx q[0];
rz(-0.73661026) q[0];
x q[1];
rz(2.759614) q[2];
sx q[2];
rz(-1.4264709) q[2];
sx q[2];
rz(-1.9461746) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5126219) q[1];
sx q[1];
rz(-1.3596467) q[1];
sx q[1];
rz(-1.3107915) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9168455) q[3];
sx q[3];
rz(-2.2888657) q[3];
sx q[3];
rz(1.8389194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1112572) q[2];
sx q[2];
rz(-1.5680485) q[2];
sx q[2];
rz(0.37720171) q[2];
rz(3.0180569) q[3];
sx q[3];
rz(-1.433452) q[3];
sx q[3];
rz(-1.4089353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1971624) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(0.29931983) q[0];
rz(-1.1209283) q[1];
sx q[1];
rz(-1.9363554) q[1];
sx q[1];
rz(-2.1790806) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4967921) q[0];
sx q[0];
rz(-1.776665) q[0];
sx q[0];
rz(1.0394761) q[0];
rz(-pi) q[1];
x q[1];
rz(2.148284) q[2];
sx q[2];
rz(-0.74314144) q[2];
sx q[2];
rz(-2.9346643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.923063) q[1];
sx q[1];
rz(-2.0835365) q[1];
sx q[1];
rz(2.1770337) q[1];
rz(-pi) q[2];
rz(0.706418) q[3];
sx q[3];
rz(-0.73143857) q[3];
sx q[3];
rz(2.4789435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0943429) q[2];
sx q[2];
rz(-2.3361358) q[2];
sx q[2];
rz(2.1979525) q[2];
rz(1.0654248) q[3];
sx q[3];
rz(-1.3506972) q[3];
sx q[3];
rz(-1.0984727) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074325398) q[0];
sx q[0];
rz(-1.704498) q[0];
sx q[0];
rz(-3.1332916) q[0];
rz(0.51757327) q[1];
sx q[1];
rz(-0.65474302) q[1];
sx q[1];
rz(-1.5932721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36309886) q[0];
sx q[0];
rz(-2.7898247) q[0];
sx q[0];
rz(1.6406873) q[0];
x q[1];
rz(2.6589633) q[2];
sx q[2];
rz(-2.5833774) q[2];
sx q[2];
rz(1.4061773) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.949985) q[1];
sx q[1];
rz(-1.9298501) q[1];
sx q[1];
rz(-2.784347) q[1];
rz(2.1287789) q[3];
sx q[3];
rz(-2.4468385) q[3];
sx q[3];
rz(-1.6266803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8325309) q[2];
sx q[2];
rz(-1.4893724) q[2];
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
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585153) q[0];
sx q[0];
rz(-2.6641088) q[0];
sx q[0];
rz(-0.68786311) q[0];
rz(1.4631924) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(0.24737839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32175049) q[0];
sx q[0];
rz(-0.30748707) q[0];
sx q[0];
rz(-0.64096682) q[0];
rz(-pi) q[1];
rz(2.3937463) q[2];
sx q[2];
rz(-0.81163844) q[2];
sx q[2];
rz(0.55392619) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0717031) q[1];
sx q[1];
rz(-0.99537288) q[1];
sx q[1];
rz(-0.50168049) q[1];
rz(2.0603397) q[3];
sx q[3];
rz(-0.26248989) q[3];
sx q[3];
rz(0.20093991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6972202) q[2];
sx q[2];
rz(-1.3345382) q[2];
sx q[2];
rz(-0.42416254) q[2];
rz(-1.0423202) q[3];
sx q[3];
rz(-0.76516953) q[3];
sx q[3];
rz(-2.585129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.0031849) q[0];
sx q[0];
rz(-2.1557032) q[0];
sx q[0];
rz(0.61088046) q[0];
rz(2.8914087) q[1];
sx q[1];
rz(-1.7411722) q[1];
sx q[1];
rz(-0.47666034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395916) q[0];
sx q[0];
rz(-1.6098813) q[0];
sx q[0];
rz(-2.0530967) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.089739775) q[2];
sx q[2];
rz(-1.3879657) q[2];
sx q[2];
rz(0.51328608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82655079) q[1];
sx q[1];
rz(-2.0619062) q[1];
sx q[1];
rz(-2.0286312) q[1];
rz(-pi) q[2];
rz(2.1564756) q[3];
sx q[3];
rz(-1.9352311) q[3];
sx q[3];
rz(1.6062615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1499947) q[2];
sx q[2];
rz(-2.1074882) q[2];
sx q[2];
rz(0.67982137) q[2];
rz(0.46197915) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(1.5010887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81944549) q[0];
sx q[0];
rz(-1.7663904) q[0];
sx q[0];
rz(0.15400259) q[0];
rz(0.18140659) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(-3.0214686) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2699566) q[0];
sx q[0];
rz(-1.9737195) q[0];
sx q[0];
rz(0.064489207) q[0];
x q[1];
rz(0.93379069) q[2];
sx q[2];
rz(-0.32378886) q[2];
sx q[2];
rz(1.1188913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97717211) q[1];
sx q[1];
rz(-1.394632) q[1];
sx q[1];
rz(3.0422496) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52773962) q[3];
sx q[3];
rz(-1.6956009) q[3];
sx q[3];
rz(0.11858701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7141562) q[2];
sx q[2];
rz(-1.3849881) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(-2.6311724) q[3];
sx q[3];
rz(-2.3000058) q[3];
sx q[3];
rz(-0.69971219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.089461483) q[0];
sx q[0];
rz(-2.1928146) q[0];
sx q[0];
rz(-0.38598886) q[0];
rz(1.5929068) q[1];
sx q[1];
rz(-0.49647757) q[1];
sx q[1];
rz(-1.5923502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693394) q[0];
sx q[0];
rz(-1.5992461) q[0];
sx q[0];
rz(-1.4017808) q[0];
x q[1];
rz(-0.73076325) q[2];
sx q[2];
rz(-0.77789069) q[2];
sx q[2];
rz(0.90532263) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4262271) q[1];
sx q[1];
rz(-1.5976397) q[1];
sx q[1];
rz(-1.6728841) q[1];
x q[2];
rz(1.9796124) q[3];
sx q[3];
rz(-0.16704112) q[3];
sx q[3];
rz(0.80834889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7768895) q[2];
sx q[2];
rz(-0.93901712) q[2];
sx q[2];
rz(1.6949863) q[2];
rz(1.0363091) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(0.026329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885289) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(-0.22008315) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(-0.73843602) q[2];
sx q[2];
rz(-1.3356575) q[2];
sx q[2];
rz(2.9267333) q[2];
rz(2.5840633) q[3];
sx q[3];
rz(-2.5685608) q[3];
sx q[3];
rz(1.3309042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
