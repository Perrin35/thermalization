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
rz(0.48756227) q[0];
sx q[0];
rz(-0.3806448) q[0];
sx q[0];
rz(-3.0523377) q[0];
rz(2.6788977) q[1];
sx q[1];
rz(-0.2806288) q[1];
sx q[1];
rz(-1.7791003) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69169626) q[0];
sx q[0];
rz(-1.5917771) q[0];
sx q[0];
rz(0.37312656) q[0];
x q[1];
rz(2.9017762) q[2];
sx q[2];
rz(-2.9887298) q[2];
sx q[2];
rz(-0.81365642) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.465724) q[1];
sx q[1];
rz(-1.7935408) q[1];
sx q[1];
rz(2.5506819) q[1];
x q[2];
rz(2.6854679) q[3];
sx q[3];
rz(-1.6509735) q[3];
sx q[3];
rz(-2.2661569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6224299) q[2];
sx q[2];
rz(-0.76202718) q[2];
sx q[2];
rz(-2.2117174) q[2];
rz(-2.8376288) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(0.16873321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794432) q[0];
sx q[0];
rz(-1.5644263) q[0];
sx q[0];
rz(1.2516578) q[0];
rz(-0.21580639) q[1];
sx q[1];
rz(-1.0390751) q[1];
sx q[1];
rz(-0.13914093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5128327) q[0];
sx q[0];
rz(-0.91673512) q[0];
sx q[0];
rz(0.20166986) q[0];
rz(0.47454796) q[2];
sx q[2];
rz(-2.5926771) q[2];
sx q[2];
rz(2.7868556) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3635837) q[1];
sx q[1];
rz(-0.87555779) q[1];
sx q[1];
rz(1.259205) q[1];
rz(-pi) q[2];
rz(-0.82461951) q[3];
sx q[3];
rz(-0.99957217) q[3];
sx q[3];
rz(-2.415641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6824049) q[2];
sx q[2];
rz(-2.6725957) q[2];
sx q[2];
rz(-1.0914618) q[2];
rz(1.546321) q[3];
sx q[3];
rz(-2.0681486) q[3];
sx q[3];
rz(2.6923164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062155398) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(0.011818258) q[0];
rz(0.91105175) q[1];
sx q[1];
rz(-1.3900737) q[1];
sx q[1];
rz(-1.8131088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9244316) q[0];
sx q[0];
rz(-1.7329441) q[0];
sx q[0];
rz(-2.6171706) q[0];
rz(-pi) q[1];
rz(-1.685789) q[2];
sx q[2];
rz(-1.3811688) q[2];
sx q[2];
rz(-0.4624838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6360259) q[1];
sx q[1];
rz(-2.2591002) q[1];
sx q[1];
rz(-1.0632374) q[1];
rz(-0.99081466) q[3];
sx q[3];
rz(-0.13306757) q[3];
sx q[3];
rz(1.7855438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38982424) q[2];
sx q[2];
rz(-1.6168892) q[2];
sx q[2];
rz(0.41979182) q[2];
rz(-1.3587562) q[3];
sx q[3];
rz(-0.32521453) q[3];
sx q[3];
rz(1.1903919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6841986) q[0];
sx q[0];
rz(-2.0036819) q[0];
sx q[0];
rz(2.6275291) q[0];
rz(0.38814107) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(-1.911389) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080175074) q[0];
sx q[0];
rz(-0.40311) q[0];
sx q[0];
rz(-0.70837428) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.94812) q[2];
sx q[2];
rz(-2.2590851) q[2];
sx q[2];
rz(1.6433059) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99707264) q[1];
sx q[1];
rz(-1.2074266) q[1];
sx q[1];
rz(-2.4097869) q[1];
rz(-1.7994653) q[3];
sx q[3];
rz(-2.2284751) q[3];
sx q[3];
rz(-2.2715037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6396883) q[2];
sx q[2];
rz(-1.014726) q[2];
sx q[2];
rz(0.91628966) q[2];
rz(-0.4942975) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(-0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2275527) q[0];
sx q[0];
rz(-0.95714772) q[0];
sx q[0];
rz(0.48602948) q[0];
rz(2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(1.1154307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971216) q[0];
sx q[0];
rz(-2.174518) q[0];
sx q[0];
rz(-1.0869727) q[0];
rz(-pi) q[1];
rz(-0.79430842) q[2];
sx q[2];
rz(-1.8657078) q[2];
sx q[2];
rz(0.25943929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72157114) q[1];
sx q[1];
rz(-0.90769115) q[1];
sx q[1];
rz(-1.0350735) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26392012) q[3];
sx q[3];
rz(-2.0368529) q[3];
sx q[3];
rz(-1.4570683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1439765) q[2];
sx q[2];
rz(-0.3207427) q[2];
sx q[2];
rz(-2.0799267) q[2];
rz(-1.6210506) q[3];
sx q[3];
rz(-1.1989667) q[3];
sx q[3];
rz(2.2475713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0735556) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(-1.2233446) q[0];
rz(2.6262737) q[1];
sx q[1];
rz(-2.4724908) q[1];
sx q[1];
rz(-1.3656778) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0037352) q[0];
sx q[0];
rz(-2.9107339) q[0];
sx q[0];
rz(2.0876711) q[0];
x q[1];
rz(1.3884945) q[2];
sx q[2];
rz(-1.6627835) q[2];
sx q[2];
rz(-1.1223886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9357696) q[1];
sx q[1];
rz(-2.6403815) q[1];
sx q[1];
rz(-2.8233714) q[1];
x q[2];
rz(1.6069769) q[3];
sx q[3];
rz(-1.2734969) q[3];
sx q[3];
rz(1.0042013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.078865163) q[2];
sx q[2];
rz(-1.8026423) q[2];
sx q[2];
rz(-1.8143181) q[2];
rz(-1.4158538) q[3];
sx q[3];
rz(-0.073181987) q[3];
sx q[3];
rz(0.88254005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13407229) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(-0.19675955) q[0];
rz(0.21601954) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(-2.3541727) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2741873) q[0];
sx q[0];
rz(-1.7599154) q[0];
sx q[0];
rz(-1.322079) q[0];
rz(-0.51047275) q[2];
sx q[2];
rz(-1.3101446) q[2];
sx q[2];
rz(-0.54768054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.710121) q[1];
sx q[1];
rz(-2.1619051) q[1];
sx q[1];
rz(1.5159056) q[1];
rz(-pi) q[2];
rz(-1.2489572) q[3];
sx q[3];
rz(-1.7506545) q[3];
sx q[3];
rz(-1.091979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54062033) q[2];
sx q[2];
rz(-1.629402) q[2];
sx q[2];
rz(-2.7579894) q[2];
rz(-2.2331734) q[3];
sx q[3];
rz(-0.81276613) q[3];
sx q[3];
rz(2.9878591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28021321) q[0];
sx q[0];
rz(-0.52556831) q[0];
sx q[0];
rz(2.5049765) q[0];
rz(-2.5909297) q[1];
sx q[1];
rz(-1.5148342) q[1];
sx q[1];
rz(2.6689463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9381959) q[0];
sx q[0];
rz(-2.2455611) q[0];
sx q[0];
rz(1.4569253) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093229712) q[2];
sx q[2];
rz(-1.5179253) q[2];
sx q[2];
rz(2.6950762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0090985) q[1];
sx q[1];
rz(-0.14202296) q[1];
sx q[1];
rz(-2.7743687) q[1];
rz(-pi) q[2];
rz(2.74284) q[3];
sx q[3];
rz(-2.4706783) q[3];
sx q[3];
rz(-0.89445597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4281281) q[2];
sx q[2];
rz(-1.2423923) q[2];
sx q[2];
rz(0.098527519) q[2];
rz(-3.1108917) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(-0.15392412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310265) q[0];
sx q[0];
rz(-0.95777804) q[0];
sx q[0];
rz(-0.73262334) q[0];
rz(-0.0064119617) q[1];
sx q[1];
rz(-2.1156204) q[1];
sx q[1];
rz(1.4220672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61056449) q[0];
sx q[0];
rz(-2.8262666) q[0];
sx q[0];
rz(2.226693) q[0];
x q[1];
rz(-1.2071916) q[2];
sx q[2];
rz(-1.0997314) q[2];
sx q[2];
rz(-0.94719145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26192415) q[1];
sx q[1];
rz(-1.4170483) q[1];
sx q[1];
rz(-1.6052482) q[1];
rz(1.4957756) q[3];
sx q[3];
rz(-2.0891694) q[3];
sx q[3];
rz(-0.013766001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0942568) q[2];
sx q[2];
rz(-1.3553268) q[2];
sx q[2];
rz(2.6195841) q[2];
rz(-0.020307288) q[3];
sx q[3];
rz(-2.4419407) q[3];
sx q[3];
rz(2.844753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8727528) q[0];
sx q[0];
rz(-1.6167384) q[0];
sx q[0];
rz(-1.1312477) q[0];
rz(-2.3616683) q[1];
sx q[1];
rz(-1.2029388) q[1];
sx q[1];
rz(2.6663229) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6492622) q[0];
sx q[0];
rz(-1.9340252) q[0];
sx q[0];
rz(0.80357768) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8653706) q[2];
sx q[2];
rz(-0.68918741) q[2];
sx q[2];
rz(-2.1222494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1058029) q[1];
sx q[1];
rz(-0.33638182) q[1];
sx q[1];
rz(-2.2835284) q[1];
rz(2.4451431) q[3];
sx q[3];
rz(-1.4587194) q[3];
sx q[3];
rz(-0.96986412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.37107006) q[2];
sx q[2];
rz(-1.5658242) q[2];
sx q[2];
rz(-1.9353665) q[2];
rz(1.3953588) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(0.079308184) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64410011) q[0];
sx q[0];
rz(-2.946749) q[0];
sx q[0];
rz(-0.83644833) q[0];
rz(-2.1438228) q[1];
sx q[1];
rz(-1.0918959) q[1];
sx q[1];
rz(0.087654884) q[1];
rz(-3.072425) q[2];
sx q[2];
rz(-1.7887049) q[2];
sx q[2];
rz(-0.6520581) q[2];
rz(0.67881696) q[3];
sx q[3];
rz(-1.4842508) q[3];
sx q[3];
rz(2.1723464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
