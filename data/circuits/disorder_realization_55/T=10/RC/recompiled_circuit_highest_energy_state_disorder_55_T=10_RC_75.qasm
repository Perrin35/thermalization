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
rz(3.6291549) q[0];
sx q[0];
rz(3.5222375) q[0];
sx q[0];
rz(9.335523) q[0];
rz(-0.46269497) q[1];
sx q[1];
rz(-2.8609639) q[1];
sx q[1];
rz(1.7791003) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69169626) q[0];
sx q[0];
rz(-1.5917771) q[0];
sx q[0];
rz(0.37312656) q[0];
rz(-pi) q[1];
rz(2.9930395) q[2];
sx q[2];
rz(-1.6069716) q[2];
sx q[2];
rz(-2.621577) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57682788) q[1];
sx q[1];
rz(-0.62678465) q[1];
sx q[1];
rz(0.38615055) q[1];
rz(-0.18043255) q[3];
sx q[3];
rz(-0.46262925) q[3];
sx q[3];
rz(0.53354665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6224299) q[2];
sx q[2];
rz(-2.3795655) q[2];
sx q[2];
rz(2.2117174) q[2];
rz(2.8376288) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(-0.16873321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794432) q[0];
sx q[0];
rz(-1.5644263) q[0];
sx q[0];
rz(-1.8899348) q[0];
rz(-2.9257863) q[1];
sx q[1];
rz(-1.0390751) q[1];
sx q[1];
rz(-3.0024517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5128327) q[0];
sx q[0];
rz(-0.91673512) q[0];
sx q[0];
rz(-0.20166986) q[0];
x q[1];
rz(-1.2982803) q[2];
sx q[2];
rz(-1.0881698) q[2];
sx q[2];
rz(-2.9543205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77800899) q[1];
sx q[1];
rz(-0.87555779) q[1];
sx q[1];
rz(-1.8823877) q[1];
x q[2];
rz(-0.81273791) q[3];
sx q[3];
rz(-2.2366282) q[3];
sx q[3];
rz(-2.8259371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6824049) q[2];
sx q[2];
rz(-2.6725957) q[2];
sx q[2];
rz(-1.0914618) q[2];
rz(1.5952716) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(2.6923164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794373) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(-0.011818258) q[0];
rz(-0.91105175) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(-1.8131088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2171611) q[0];
sx q[0];
rz(-1.4086485) q[0];
sx q[0];
rz(-2.6171706) q[0];
rz(0.19085796) q[2];
sx q[2];
rz(-1.4578739) q[2];
sx q[2];
rz(1.1300805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78643936) q[1];
sx q[1];
rz(-0.82988534) q[1];
sx q[1];
rz(0.53373465) q[1];
rz(3.0683682) q[3];
sx q[3];
rz(-1.6820046) q[3];
sx q[3];
rz(-1.2014845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7517684) q[2];
sx q[2];
rz(-1.5247034) q[2];
sx q[2];
rz(-2.7218008) q[2];
rz(1.3587562) q[3];
sx q[3];
rz(-2.8163781) q[3];
sx q[3];
rz(-1.9512008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6841986) q[0];
sx q[0];
rz(-2.0036819) q[0];
sx q[0];
rz(-0.5140636) q[0];
rz(0.38814107) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(-1.911389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0614176) q[0];
sx q[0];
rz(-0.40311) q[0];
sx q[0];
rz(-0.70837428) q[0];
rz(-pi) q[1];
rz(1.1934727) q[2];
sx q[2];
rz(-0.88250752) q[2];
sx q[2];
rz(1.4982868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88273222) q[1];
sx q[1];
rz(-2.2453935) q[1];
sx q[1];
rz(2.0432977) q[1];
rz(0.28546412) q[3];
sx q[3];
rz(-0.69068161) q[3];
sx q[3];
rz(2.635298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5019044) q[2];
sx q[2];
rz(-2.1268667) q[2];
sx q[2];
rz(2.225303) q[2];
rz(-2.6472951) q[3];
sx q[3];
rz(-1.1980134) q[3];
sx q[3];
rz(2.3673207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.91404) q[0];
sx q[0];
rz(-2.1844449) q[0];
sx q[0];
rz(0.48602948) q[0];
rz(2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(-2.026162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19768342) q[0];
sx q[0];
rz(-0.75431529) q[0];
sx q[0];
rz(0.59342845) q[0];
rz(-1.979823) q[2];
sx q[2];
rz(-0.81950775) q[2];
sx q[2];
rz(1.5425615) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.045687606) q[1];
sx q[1];
rz(-0.8262016) q[1];
sx q[1];
rz(-0.5788486) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0492184) q[3];
sx q[3];
rz(-0.53077519) q[3];
sx q[3];
rz(-1.1431183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9976161) q[2];
sx q[2];
rz(-2.82085) q[2];
sx q[2];
rz(-2.0799267) q[2];
rz(1.6210506) q[3];
sx q[3];
rz(-1.1989667) q[3];
sx q[3];
rz(-2.2475713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068037085) q[0];
sx q[0];
rz(-0.03980045) q[0];
sx q[0];
rz(1.2233446) q[0];
rz(2.6262737) q[1];
sx q[1];
rz(-2.4724908) q[1];
sx q[1];
rz(1.7759148) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66636234) q[0];
sx q[0];
rz(-1.3705374) q[0];
sx q[0];
rz(3.0259575) q[0];
rz(-3.0480644) q[2];
sx q[2];
rz(-1.7523189) q[2];
sx q[2];
rz(0.43147555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9881043) q[1];
sx q[1];
rz(-1.0968913) q[1];
sx q[1];
rz(1.4010282) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8441098) q[3];
sx q[3];
rz(-1.5362036) q[3];
sx q[3];
rz(-2.5643947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0627275) q[2];
sx q[2];
rz(-1.3389503) q[2];
sx q[2];
rz(1.3272746) q[2];
rz(1.4158538) q[3];
sx q[3];
rz(-0.073181987) q[3];
sx q[3];
rz(2.2590526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13407229) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(0.19675955) q[0];
rz(0.21601954) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(-2.3541727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93381017) q[0];
sx q[0];
rz(-2.8303307) q[0];
sx q[0];
rz(-2.2317076) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51047275) q[2];
sx q[2];
rz(-1.831448) q[2];
sx q[2];
rz(2.5939121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.43147165) q[1];
sx q[1];
rz(-2.1619051) q[1];
sx q[1];
rz(-1.5159056) q[1];
x q[2];
rz(-1.8926354) q[3];
sx q[3];
rz(-1.7506545) q[3];
sx q[3];
rz(1.091979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54062033) q[2];
sx q[2];
rz(-1.5121907) q[2];
sx q[2];
rz(2.7579894) q[2];
rz(-0.90841928) q[3];
sx q[3];
rz(-0.81276613) q[3];
sx q[3];
rz(-2.9878591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28021321) q[0];
sx q[0];
rz(-2.6160243) q[0];
sx q[0];
rz(0.63661611) q[0];
rz(-2.5909297) q[1];
sx q[1];
rz(-1.5148342) q[1];
sx q[1];
rz(-0.47264636) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.11926) q[0];
sx q[0];
rz(-0.68281931) q[0];
sx q[0];
rz(0.14108087) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6246895) q[2];
sx q[2];
rz(-0.10714018) q[2];
sx q[2];
rz(-2.531749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.19785205) q[1];
sx q[1];
rz(-1.6216369) q[1];
sx q[1];
rz(-0.13266854) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2718464) q[3];
sx q[3];
rz(-0.96072324) q[3];
sx q[3];
rz(1.7536193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7134646) q[2];
sx q[2];
rz(-1.8992004) q[2];
sx q[2];
rz(-3.0430651) q[2];
rz(0.030700961) q[3];
sx q[3];
rz(-1.5714329) q[3];
sx q[3];
rz(-2.9876685) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
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
rz(-1.0259722) q[1];
sx q[1];
rz(-1.4220672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32860562) q[0];
sx q[0];
rz(-1.7610794) q[0];
sx q[0];
rz(1.3178131) q[0];
x q[1];
rz(-0.60958013) q[2];
sx q[2];
rz(-0.58660075) q[2];
sx q[2];
rz(-1.644949) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8274424) q[1];
sx q[1];
rz(-1.536751) q[1];
sx q[1];
rz(2.9877547) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4957756) q[3];
sx q[3];
rz(-1.0524233) q[3];
sx q[3];
rz(0.013766001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0473359) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(0.52200851) q[2];
rz(0.020307288) q[3];
sx q[3];
rz(-0.69965196) q[3];
sx q[3];
rz(-0.29683963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2688399) q[0];
sx q[0];
rz(-1.6167384) q[0];
sx q[0];
rz(-2.010345) q[0];
rz(-2.3616683) q[1];
sx q[1];
rz(-1.2029388) q[1];
sx q[1];
rz(-0.47526971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330547) q[0];
sx q[0];
rz(-0.86472874) q[0];
sx q[0];
rz(0.48581328) q[0];
rz(-2.9067801) q[2];
sx q[2];
rz(-2.2249892) q[2];
sx q[2];
rz(-1.393911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0357897) q[1];
sx q[1];
rz(-2.8052108) q[1];
sx q[1];
rz(-0.85806429) q[1];
x q[2];
rz(0.17367878) q[3];
sx q[3];
rz(-2.4376737) q[3];
sx q[3];
rz(-2.6736163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37107006) q[2];
sx q[2];
rz(-1.5658242) q[2];
sx q[2];
rz(1.2062262) q[2];
rz(-1.7462339) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(0.079308184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(0.64410011) q[0];
sx q[0];
rz(-0.19484367) q[0];
sx q[0];
rz(2.3051443) q[0];
rz(0.99776987) q[1];
sx q[1];
rz(-1.0918959) q[1];
sx q[1];
rz(0.087654884) q[1];
rz(-0.069167698) q[2];
sx q[2];
rz(-1.3528878) q[2];
sx q[2];
rz(2.4895346) q[2];
rz(-0.67881696) q[3];
sx q[3];
rz(-1.6573418) q[3];
sx q[3];
rz(-0.96924622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
