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
rz(0.95712823) q[0];
sx q[0];
rz(4.794802) q[0];
sx q[0];
rz(10.573536) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(1.6869071) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5243149) q[0];
sx q[0];
rz(-1.9348659) q[0];
sx q[0];
rz(-0.84366701) q[0];
x q[1];
rz(1.4613011) q[2];
sx q[2];
rz(-0.81297648) q[2];
sx q[2];
rz(2.2535498) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7976746) q[1];
sx q[1];
rz(-2.3189088) q[1];
sx q[1];
rz(1.3901346) q[1];
x q[2];
rz(0.1249073) q[3];
sx q[3];
rz(-2.4797399) q[3];
sx q[3];
rz(-0.13368363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.181432) q[2];
sx q[2];
rz(-2.1655607) q[2];
sx q[2];
rz(1.9350447) q[2];
rz(-1.9801961) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(-1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.96381617) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(-0.97208446) q[0];
rz(-2.8731335) q[1];
sx q[1];
rz(-1.4413036) q[1];
sx q[1];
rz(-2.6867902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3362124) q[0];
sx q[0];
rz(-2.8180362) q[0];
sx q[0];
rz(-1.1624883) q[0];
rz(1.9591397) q[2];
sx q[2];
rz(-1.4253311) q[2];
sx q[2];
rz(2.5487473) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.763096) q[1];
sx q[1];
rz(-0.46751577) q[1];
sx q[1];
rz(3.0276831) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56549325) q[3];
sx q[3];
rz(-1.0988937) q[3];
sx q[3];
rz(-1.3513235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0316169) q[2];
sx q[2];
rz(-2.4133108) q[2];
sx q[2];
rz(0.99962437) q[2];
rz(-3.055618) q[3];
sx q[3];
rz(-1.5558259) q[3];
sx q[3];
rz(3.1302248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.40083945) q[0];
sx q[0];
rz(-1.6735621) q[0];
sx q[0];
rz(2.9174347) q[0];
rz(1.5466746) q[1];
sx q[1];
rz(-1.5983351) q[1];
sx q[1];
rz(1.3067783) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9973361) q[0];
sx q[0];
rz(-1.6929417) q[0];
sx q[0];
rz(0.55295487) q[0];
x q[1];
rz(2.4073657) q[2];
sx q[2];
rz(-0.9817183) q[2];
sx q[2];
rz(-0.40225077) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8301218) q[1];
sx q[1];
rz(-1.9852299) q[1];
sx q[1];
rz(-0.45763514) q[1];
x q[2];
rz(-0.66931386) q[3];
sx q[3];
rz(-1.770854) q[3];
sx q[3];
rz(-1.1049096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.811502) q[2];
sx q[2];
rz(-1.6197438) q[2];
sx q[2];
rz(-1.0203934) q[2];
rz(1.8118106) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(1.6167195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.735585) q[0];
sx q[0];
rz(-1.0574874) q[0];
sx q[0];
rz(-1.1757346) q[0];
rz(-0.27578393) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(0.63794678) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36786554) q[0];
sx q[0];
rz(-1.5380895) q[0];
sx q[0];
rz(-1.5601394) q[0];
x q[1];
rz(-2.5204896) q[2];
sx q[2];
rz(-1.3912462) q[2];
sx q[2];
rz(-0.68815069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0093858) q[1];
sx q[1];
rz(-1.8957545) q[1];
sx q[1];
rz(-0.64909776) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8218587) q[3];
sx q[3];
rz(-1.481062) q[3];
sx q[3];
rz(-3.0872726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5059169) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(-1.4471311) q[2];
rz(0.21008374) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(-0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90784812) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(-2.0727169) q[0];
rz(1.2999889) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(1.2058421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.962709) q[0];
sx q[0];
rz(-0.41085748) q[0];
sx q[0];
rz(-1.7715447) q[0];
rz(-pi) q[1];
rz(1.5832434) q[2];
sx q[2];
rz(-2.4887391) q[2];
sx q[2];
rz(-1.659153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8925807) q[1];
sx q[1];
rz(-2.4155185) q[1];
sx q[1];
rz(1.4145265) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30773325) q[3];
sx q[3];
rz(-2.640297) q[3];
sx q[3];
rz(-3.1395819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24870366) q[2];
sx q[2];
rz(-2.7619669) q[2];
sx q[2];
rz(-3.0730263) q[2];
rz(0.93962234) q[3];
sx q[3];
rz(-1.7328123) q[3];
sx q[3];
rz(1.9127964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762961) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(-0.5214386) q[0];
rz(-1.5154845) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(-0.62561402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2170048) q[0];
sx q[0];
rz(-2.6007605) q[0];
sx q[0];
rz(-1.5103673) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4418774) q[2];
sx q[2];
rz(-1.9471418) q[2];
sx q[2];
rz(-2.3835045) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82078104) q[1];
sx q[1];
rz(-0.50034517) q[1];
sx q[1];
rz(1.149575) q[1];
rz(0.13497495) q[3];
sx q[3];
rz(-2.7162281) q[3];
sx q[3];
rz(-2.7736026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0079415) q[2];
sx q[2];
rz(-2.5456754) q[2];
sx q[2];
rz(3.0360743) q[2];
rz(2.1775235) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(2.156637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725726) q[0];
sx q[0];
rz(-2.9264086) q[0];
sx q[0];
rz(1.4129114) q[0];
rz(0.37881306) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(2.5884195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7460497) q[0];
sx q[0];
rz(-0.93727124) q[0];
sx q[0];
rz(-2.6468011) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4804391) q[2];
sx q[2];
rz(-1.4298265) q[2];
sx q[2];
rz(0.80328926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70771101) q[1];
sx q[1];
rz(-1.6110782) q[1];
sx q[1];
rz(-0.28357419) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91316789) q[3];
sx q[3];
rz(-2.8165157) q[3];
sx q[3];
rz(0.43077786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(-0.41213948) q[2];
rz(-0.66017094) q[3];
sx q[3];
rz(-0.5414525) q[3];
sx q[3];
rz(-0.49303833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93765813) q[0];
sx q[0];
rz(-2.0964607) q[0];
sx q[0];
rz(-0.20189051) q[0];
rz(-2.0837325) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(0.26022628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3500811) q[0];
sx q[0];
rz(-0.66257325) q[0];
sx q[0];
rz(0.64278472) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7041928) q[2];
sx q[2];
rz(-2.045131) q[2];
sx q[2];
rz(2.2651644) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7937745) q[1];
sx q[1];
rz(-1.8500237) q[1];
sx q[1];
rz(-0.3600959) q[1];
rz(-pi) q[2];
rz(1.7395571) q[3];
sx q[3];
rz(-2.2243554) q[3];
sx q[3];
rz(-0.26168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95320931) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(2.5592213) q[2];
rz(0.44238704) q[3];
sx q[3];
rz(-1.9059537) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22170947) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(-1.4269933) q[0];
rz(1.7859979) q[1];
sx q[1];
rz(-1.6537063) q[1];
sx q[1];
rz(1.4367163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6528931) q[0];
sx q[0];
rz(-2.821229) q[0];
sx q[0];
rz(-1.6344749) q[0];
rz(-pi) q[1];
x q[1];
rz(2.605805) q[2];
sx q[2];
rz(-0.54507191) q[2];
sx q[2];
rz(-1.6824012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9475262) q[1];
sx q[1];
rz(-1.8119651) q[1];
sx q[1];
rz(3.126866) q[1];
rz(0.15738486) q[3];
sx q[3];
rz(-2.0920355) q[3];
sx q[3];
rz(0.97053274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23472486) q[2];
sx q[2];
rz(-1.2805254) q[2];
sx q[2];
rz(-1.5391763) q[2];
rz(0.60798821) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(-2.4517945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36668396) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(-0.82345024) q[0];
rz(-2.2460294) q[1];
sx q[1];
rz(-1.3500373) q[1];
sx q[1];
rz(0.38111883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279042) q[0];
sx q[0];
rz(-1.1263337) q[0];
sx q[0];
rz(0.60143394) q[0];
rz(2.7627583) q[2];
sx q[2];
rz(-0.4288579) q[2];
sx q[2];
rz(-2.2414686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9420112) q[1];
sx q[1];
rz(-1.2855654) q[1];
sx q[1];
rz(0.99883325) q[1];
rz(2.1249843) q[3];
sx q[3];
rz(-2.1234406) q[3];
sx q[3];
rz(-1.9345528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8054008) q[2];
sx q[2];
rz(-2.7298584) q[2];
sx q[2];
rz(-2.1545048) q[2];
rz(1.5385212) q[3];
sx q[3];
rz(-2.5542732) q[3];
sx q[3];
rz(0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22454746) q[0];
sx q[0];
rz(-1.2397091) q[0];
sx q[0];
rz(1.5201257) q[0];
rz(0.70980258) q[1];
sx q[1];
rz(-2.5820422) q[1];
sx q[1];
rz(-1.6834264) q[1];
rz(-2.7771525) q[2];
sx q[2];
rz(-1.6779998) q[2];
sx q[2];
rz(-2.6371434) q[2];
rz(0.98593203) q[3];
sx q[3];
rz(-0.91888792) q[3];
sx q[3];
rz(2.1666913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
