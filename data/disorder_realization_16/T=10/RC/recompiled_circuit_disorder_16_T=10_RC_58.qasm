OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(3.9711877) q[0];
sx q[0];
rz(9.2708099) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(-2.8032803) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4727288) q[0];
sx q[0];
rz(-0.40574408) q[0];
sx q[0];
rz(-2.2292482) q[0];
x q[1];
rz(0.89262427) q[2];
sx q[2];
rz(-1.8632338) q[2];
sx q[2];
rz(0.48722789) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6602064) q[1];
sx q[1];
rz(-1.2848789) q[1];
sx q[1];
rz(1.6260765) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5987414) q[3];
sx q[3];
rz(-2.6290647) q[3];
sx q[3];
rz(0.41748369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14264318) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(1.1738698) q[2];
rz(0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(0.064963438) q[0];
rz(0.57463542) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.2423135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14318289) q[0];
sx q[0];
rz(-1.3230723) q[0];
sx q[0];
rz(-1.7032911) q[0];
rz(-0.93810268) q[2];
sx q[2];
rz(-1.4422851) q[2];
sx q[2];
rz(-3.1046257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6537135) q[1];
sx q[1];
rz(-2.5839845) q[1];
sx q[1];
rz(-2.8370268) q[1];
x q[2];
rz(-1.0817238) q[3];
sx q[3];
rz(-2.2197669) q[3];
sx q[3];
rz(-1.058941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(2.5578965) q[2];
rz(-0.57404533) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(-0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(1.0167936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8903538) q[0];
sx q[0];
rz(-1.4892502) q[0];
sx q[0];
rz(-1.606705) q[0];
rz(-2.617308) q[2];
sx q[2];
rz(-1.1384083) q[2];
sx q[2];
rz(-1.0334894) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.180417) q[1];
sx q[1];
rz(-0.75913402) q[1];
sx q[1];
rz(-1.7708066) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.049116491) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(-1.0055055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-2.0920848) q[2];
rz(-0.63878757) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(-1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(-1.6695492) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-0.24681117) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033337489) q[0];
sx q[0];
rz(-2.4318683) q[0];
sx q[0];
rz(1.3413341) q[0];
x q[1];
rz(-0.39101379) q[2];
sx q[2];
rz(-2.6303929) q[2];
sx q[2];
rz(-0.15399394) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3226763) q[1];
sx q[1];
rz(-0.37441844) q[1];
sx q[1];
rz(-0.14426343) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9645343) q[3];
sx q[3];
rz(-0.65440946) q[3];
sx q[3];
rz(0.27458336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(2.4345496) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(2.2533806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8108869) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(-1.8141618) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9527234) q[0];
sx q[0];
rz(-1.5367322) q[0];
sx q[0];
rz(-2.6012095) q[0];
rz(-0.18647285) q[2];
sx q[2];
rz(-1.7261793) q[2];
sx q[2];
rz(1.965167) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76368139) q[1];
sx q[1];
rz(-1.7421725) q[1];
sx q[1];
rz(0.59314368) q[1];
rz(2.4928717) q[3];
sx q[3];
rz(-1.7626764) q[3];
sx q[3];
rz(-2.3313525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43626943) q[2];
sx q[2];
rz(-2.1537809) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-0.60384408) q[3];
sx q[3];
rz(0.50271547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1335063) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(2.0603518) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(3.016901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69783083) q[0];
sx q[0];
rz(-0.19653453) q[0];
sx q[0];
rz(-0.4561119) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6090392) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(-2.2720624) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5387419) q[1];
sx q[1];
rz(-2.3067143) q[1];
sx q[1];
rz(-1.6470626) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23645225) q[3];
sx q[3];
rz(-1.7835622) q[3];
sx q[3];
rz(0.80602431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51320118) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(-0.46328059) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(2.0708864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28744222) q[0];
sx q[0];
rz(-1.5702015) q[0];
sx q[0];
rz(0.48405148) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35371874) q[2];
sx q[2];
rz(-0.36834799) q[2];
sx q[2];
rz(3.0596717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0106525) q[1];
sx q[1];
rz(-0.7796692) q[1];
sx q[1];
rz(2.9995549) q[1];
rz(1.5504863) q[3];
sx q[3];
rz(-1.9029402) q[3];
sx q[3];
rz(-0.16196812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-0.32361844) q[2];
rz(0.98179022) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(-1.9288829) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-2.1941197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3377209) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(-1.7659811) q[0];
rz(-1.7296373) q[2];
sx q[2];
rz(-1.7454595) q[2];
sx q[2];
rz(-0.58089248) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2944813) q[1];
sx q[1];
rz(-2.1335019) q[1];
sx q[1];
rz(-0.28114762) q[1];
rz(2.794907) q[3];
sx q[3];
rz(-2.3035435) q[3];
sx q[3];
rz(-0.29688641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-2.6718111) q[2];
rz(-1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(2.8619213) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(-1.3289733) q[0];
rz(-0.7912311) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(-0.20283094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348939) q[0];
sx q[0];
rz(-1.5968423) q[0];
sx q[0];
rz(-3.1025725) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22325309) q[2];
sx q[2];
rz(-1.3666144) q[2];
sx q[2];
rz(-2.5399361) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.054246) q[1];
sx q[1];
rz(-3.0294703) q[1];
sx q[1];
rz(1.1152769) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4167452) q[3];
sx q[3];
rz(-2.027958) q[3];
sx q[3];
rz(2.9726213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(-0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(0.18173519) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(-0.28082401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1998394) q[0];
sx q[0];
rz(-1.117525) q[0];
sx q[0];
rz(-2.0487294) q[0];
rz(-pi) q[1];
rz(-0.21364613) q[2];
sx q[2];
rz(-1.3823969) q[2];
sx q[2];
rz(1.633916) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9210789) q[1];
sx q[1];
rz(-1.6069357) q[1];
sx q[1];
rz(-1.6284579) q[1];
rz(-pi) q[2];
x q[2];
rz(2.093408) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(0.89695938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4939209) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(0.87396809) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-0.95332425) q[2];
sx q[2];
rz(-0.8669903) q[2];
sx q[2];
rz(0.16624761) q[2];
rz(-1.1602041) q[3];
sx q[3];
rz(-1.2413597) q[3];
sx q[3];
rz(1.6551457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
