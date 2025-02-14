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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(1.992835) q[0];
rz(-2.5454638) q[1];
sx q[1];
rz(-2.7355173) q[1];
sx q[1];
rz(-1.6869071) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8073976) q[0];
sx q[0];
rz(-2.3435623) q[0];
sx q[0];
rz(-1.0502771) q[0];
rz(-pi) q[1];
rz(2.3807736) q[2];
sx q[2];
rz(-1.4913412) q[2];
sx q[2];
rz(2.3834199) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.64950424) q[1];
sx q[1];
rz(-1.4387114) q[1];
sx q[1];
rz(-2.3852939) q[1];
x q[2];
rz(-0.65806453) q[3];
sx q[3];
rz(-1.6474373) q[3];
sx q[3];
rz(1.6057526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9601606) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(1.206548) q[2];
rz(1.9801961) q[3];
sx q[3];
rz(-0.56777081) q[3];
sx q[3];
rz(1.5709706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96381617) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(-2.1695082) q[0];
rz(0.26845911) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(-0.45480248) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3362124) q[0];
sx q[0];
rz(-2.8180362) q[0];
sx q[0];
rz(-1.1624883) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9846086) q[2];
sx q[2];
rz(-1.1867701) q[2];
sx q[2];
rz(2.2228754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5059591) q[1];
sx q[1];
rz(-2.035043) q[1];
sx q[1];
rz(-1.5134769) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0271026) q[3];
sx q[3];
rz(-1.0732526) q[3];
sx q[3];
rz(-0.061391679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1099757) q[2];
sx q[2];
rz(-2.4133108) q[2];
sx q[2];
rz(2.1419683) q[2];
rz(-3.055618) q[3];
sx q[3];
rz(-1.5857668) q[3];
sx q[3];
rz(0.011367817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407532) q[0];
sx q[0];
rz(-1.4680306) q[0];
sx q[0];
rz(0.22415796) q[0];
rz(-1.5466746) q[1];
sx q[1];
rz(-1.5983351) q[1];
sx q[1];
rz(-1.3067783) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9973361) q[0];
sx q[0];
rz(-1.4486509) q[0];
sx q[0];
rz(2.5886378) q[0];
x q[1];
rz(-0.8379015) q[2];
sx q[2];
rz(-0.97992491) q[2];
sx q[2];
rz(-1.5082347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8301218) q[1];
sx q[1];
rz(-1.1563627) q[1];
sx q[1];
rz(-0.45763514) q[1];
rz(2.4722788) q[3];
sx q[3];
rz(-1.3707386) q[3];
sx q[3];
rz(-2.0366831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.811502) q[2];
sx q[2];
rz(-1.5218488) q[2];
sx q[2];
rz(-2.1211993) q[2];
rz(-1.8118106) q[3];
sx q[3];
rz(-1.1251757) q[3];
sx q[3];
rz(-1.5248732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.735585) q[0];
sx q[0];
rz(-1.0574874) q[0];
sx q[0];
rz(-1.965858) q[0];
rz(0.27578393) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(-0.63794678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36786554) q[0];
sx q[0];
rz(-1.5380895) q[0];
sx q[0];
rz(-1.5601394) q[0];
rz(-pi) q[1];
rz(-2.5204896) q[2];
sx q[2];
rz(-1.3912462) q[2];
sx q[2];
rz(2.453442) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6762661) q[1];
sx q[1];
rz(-0.96083896) q[1];
sx q[1];
rz(-1.9708939) q[1];
rz(-1.4762991) q[3];
sx q[3];
rz(-1.8891984) q[3];
sx q[3];
rz(-1.6547784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5059169) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(-1.4471311) q[2];
rz(2.9315089) q[3];
sx q[3];
rz(-0.45322067) q[3];
sx q[3];
rz(2.213018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2337445) q[0];
sx q[0];
rz(-0.15234983) q[0];
sx q[0];
rz(2.0727169) q[0];
rz(1.8416038) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(1.9357505) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.962709) q[0];
sx q[0];
rz(-0.41085748) q[0];
sx q[0];
rz(-1.7715447) q[0];
x q[1];
rz(1.5832434) q[2];
sx q[2];
rz(-2.4887391) q[2];
sx q[2];
rz(1.4824397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5806535) q[1];
sx q[1];
rz(-1.4672797) q[1];
sx q[1];
rz(2.2907881) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4063115) q[3];
sx q[3];
rz(-1.0950673) q[3];
sx q[3];
rz(0.34970899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24870366) q[2];
sx q[2];
rz(-2.7619669) q[2];
sx q[2];
rz(0.068566337) q[2];
rz(-2.2019703) q[3];
sx q[3];
rz(-1.7328123) q[3];
sx q[3];
rz(-1.2287963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.762961) q[0];
sx q[0];
rz(-1.3968422) q[0];
sx q[0];
rz(-0.5214386) q[0];
rz(-1.6261082) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(-2.5159786) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2170048) q[0];
sx q[0];
rz(-2.6007605) q[0];
sx q[0];
rz(1.5103673) q[0];
rz(-pi) q[1];
rz(-2.8270673) q[2];
sx q[2];
rz(-2.7447766) q[2];
sx q[2];
rz(2.0443969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.292899) q[1];
sx q[1];
rz(-2.0239415) q[1];
sx q[1];
rz(-2.921656) q[1];
x q[2];
rz(-0.13497495) q[3];
sx q[3];
rz(-2.7162281) q[3];
sx q[3];
rz(2.7736026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1336512) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(-0.1055183) q[2];
rz(0.96406913) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(0.9849557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725726) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(1.7286812) q[0];
rz(2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(0.55317318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39554292) q[0];
sx q[0];
rz(-0.93727124) q[0];
sx q[0];
rz(-2.6468011) q[0];
x q[1];
rz(3.000053) q[2];
sx q[2];
rz(-1.6602548) q[2];
sx q[2];
rz(-2.3613561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4158226) q[1];
sx q[1];
rz(-2.8552481) q[1];
sx q[1];
rz(2.9985266) q[1];
x q[2];
rz(-0.91316789) q[3];
sx q[3];
rz(-2.8165157) q[3];
sx q[3];
rz(0.43077786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46099123) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(0.41213948) q[2];
rz(2.4814217) q[3];
sx q[3];
rz(-2.6001402) q[3];
sx q[3];
rz(-2.6485543) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93765813) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(-2.9397021) q[0];
rz(-2.0837325) q[1];
sx q[1];
rz(-0.36483929) q[1];
sx q[1];
rz(2.8813664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3500811) q[0];
sx q[0];
rz(-2.4790194) q[0];
sx q[0];
rz(-0.64278472) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.260602) q[2];
sx q[2];
rz(-2.5079003) q[2];
sx q[2];
rz(1.6735759) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7937745) q[1];
sx q[1];
rz(-1.8500237) q[1];
sx q[1];
rz(2.7814968) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21589888) q[3];
sx q[3];
rz(-0.67189081) q[3];
sx q[3];
rz(-0.011530576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95320931) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(-2.5592213) q[2];
rz(-0.44238704) q[3];
sx q[3];
rz(-1.235639) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9198832) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(1.7145994) q[0];
rz(1.3555948) q[1];
sx q[1];
rz(-1.6537063) q[1];
sx q[1];
rz(-1.4367163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021649927) q[0];
sx q[0];
rz(-1.5908373) q[0];
sx q[0];
rz(-1.2510386) q[0];
x q[1];
rz(0.48060178) q[2];
sx q[2];
rz(-1.8386823) q[2];
sx q[2];
rz(0.3581697) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2556527) q[1];
sx q[1];
rz(-2.8999834) q[1];
sx q[1];
rz(1.5109946) q[1];
rz(0.15738486) q[3];
sx q[3];
rz(-1.0495571) q[3];
sx q[3];
rz(-0.97053274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23472486) q[2];
sx q[2];
rz(-1.8610672) q[2];
sx q[2];
rz(-1.5391763) q[2];
rz(-0.60798821) q[3];
sx q[3];
rz(-1.7323078) q[3];
sx q[3];
rz(-2.4517945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749087) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(2.3181424) q[0];
rz(0.89556328) q[1];
sx q[1];
rz(-1.7915553) q[1];
sx q[1];
rz(-0.38111883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6247647) q[0];
sx q[0];
rz(-0.73113686) q[0];
sx q[0];
rz(-0.69964377) q[0];
x q[1];
rz(1.4032768) q[2];
sx q[2];
rz(-1.9674656) q[2];
sx q[2];
rz(1.3127017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5911853) q[1];
sx q[1];
rz(-2.1169615) q[1];
sx q[1];
rz(-0.3355432) q[1];
rz(0.70638871) q[3];
sx q[3];
rz(-2.3800142) q[3];
sx q[3];
rz(1.0672081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8054008) q[2];
sx q[2];
rz(-0.41173428) q[2];
sx q[2];
rz(-2.1545048) q[2];
rz(-1.5385212) q[3];
sx q[3];
rz(-0.58731949) q[3];
sx q[3];
rz(-2.9015598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9170452) q[0];
sx q[0];
rz(-1.9018835) q[0];
sx q[0];
rz(-1.6214669) q[0];
rz(-0.70980258) q[1];
sx q[1];
rz(-0.55955049) q[1];
sx q[1];
rz(1.4581663) q[1];
rz(-2.8483656) q[2];
sx q[2];
rz(-2.7623889) q[2];
sx q[2];
rz(1.8017906) q[2];
rz(-2.1556606) q[3];
sx q[3];
rz(-0.91888792) q[3];
sx q[3];
rz(2.1666913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
