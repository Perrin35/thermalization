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
rz(0.08925499) q[0];
rz(2.6788977) q[1];
sx q[1];
rz(-0.2806288) q[1];
sx q[1];
rz(-1.7791003) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69169626) q[0];
sx q[0];
rz(-1.5498156) q[0];
sx q[0];
rz(-0.37312656) q[0];
rz(-pi) q[1];
rz(0.2398165) q[2];
sx q[2];
rz(-2.9887298) q[2];
sx q[2];
rz(-2.3279362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57682788) q[1];
sx q[1];
rz(-2.514808) q[1];
sx q[1];
rz(-0.38615055) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4815349) q[3];
sx q[3];
rz(-1.1162471) q[3];
sx q[3];
rz(0.7346357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6224299) q[2];
sx q[2];
rz(-2.3795655) q[2];
sx q[2];
rz(-0.92987522) q[2];
rz(-0.30396384) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(-0.16873321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794432) q[0];
sx q[0];
rz(-1.5771663) q[0];
sx q[0];
rz(1.8899348) q[0];
rz(-2.9257863) q[1];
sx q[1];
rz(-2.1025175) q[1];
sx q[1];
rz(3.0024517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2073832) q[0];
sx q[0];
rz(-1.4111526) q[0];
sx q[0];
rz(-0.90682323) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49824841) q[2];
sx q[2];
rz(-1.8115269) q[2];
sx q[2];
rz(1.5125076) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5522668) q[1];
sx q[1];
rz(-1.8084452) q[1];
sx q[1];
rz(2.4219805) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3169731) q[3];
sx q[3];
rz(-0.99957217) q[3];
sx q[3];
rz(-2.415641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6824049) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(1.0914618) q[2];
rz(1.546321) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(0.44927621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062155398) q[0];
sx q[0];
rz(-1.9154444) q[0];
sx q[0];
rz(3.1297744) q[0];
rz(0.91105175) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(-1.3284838) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5157943) q[0];
sx q[0];
rz(-0.54667241) q[0];
sx q[0];
rz(2.8258219) q[0];
x q[1];
rz(2.9507347) q[2];
sx q[2];
rz(-1.4578739) q[2];
sx q[2];
rz(2.0115122) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5055667) q[1];
sx q[1];
rz(-0.88249245) q[1];
sx q[1];
rz(2.0783552) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99081466) q[3];
sx q[3];
rz(-3.0085251) q[3];
sx q[3];
rz(1.7855438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38982424) q[2];
sx q[2];
rz(-1.6168892) q[2];
sx q[2];
rz(-2.7218008) q[2];
rz(-1.3587562) q[3];
sx q[3];
rz(-2.8163781) q[3];
sx q[3];
rz(-1.1903919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.457394) q[0];
sx q[0];
rz(-2.0036819) q[0];
sx q[0];
rz(2.6275291) q[0];
rz(2.7534516) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(-1.2302037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.311545) q[0];
sx q[0];
rz(-1.8732949) q[0];
sx q[0];
rz(1.8414458) q[0];
rz(-pi) q[1];
rz(0.72429652) q[2];
sx q[2];
rz(-1.2822552) q[2];
sx q[2];
rz(-2.9674825) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88273222) q[1];
sx q[1];
rz(-2.2453935) q[1];
sx q[1];
rz(1.098295) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67048002) q[3];
sx q[3];
rz(-1.7511715) q[3];
sx q[3];
rz(-0.84202858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5019044) q[2];
sx q[2];
rz(-2.1268667) q[2];
sx q[2];
rz(2.225303) q[2];
rz(2.6472951) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(-0.77427197) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91404) q[0];
sx q[0];
rz(-0.95714772) q[0];
sx q[0];
rz(-0.48602948) q[0];
rz(-2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(-1.1154307) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19768342) q[0];
sx q[0];
rz(-0.75431529) q[0];
sx q[0];
rz(-2.5481642) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40256315) q[2];
sx q[2];
rz(-2.3057115) q[2];
sx q[2];
rz(1.0332359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9420227) q[1];
sx q[1];
rz(-1.1567819) q[1];
sx q[1];
rz(0.7374108) q[1];
rz(-0.26392012) q[3];
sx q[3];
rz(-2.0368529) q[3];
sx q[3];
rz(-1.4570683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9976161) q[2];
sx q[2];
rz(-0.3207427) q[2];
sx q[2];
rz(1.0616659) q[2];
rz(-1.6210506) q[3];
sx q[3];
rz(-1.942626) q[3];
sx q[3];
rz(0.89402136) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068037085) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(1.918248) q[0];
rz(-2.6262737) q[1];
sx q[1];
rz(-2.4724908) q[1];
sx q[1];
rz(-1.7759148) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66636234) q[0];
sx q[0];
rz(-1.7710553) q[0];
sx q[0];
rz(-3.0259575) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3884945) q[2];
sx q[2];
rz(-1.6627835) q[2];
sx q[2];
rz(-1.1223886) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9881043) q[1];
sx q[1];
rz(-2.0447013) q[1];
sx q[1];
rz(-1.4010282) q[1];
rz(2.8441098) q[3];
sx q[3];
rz(-1.5362036) q[3];
sx q[3];
rz(-0.57719798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0627275) q[2];
sx q[2];
rz(-1.8026423) q[2];
sx q[2];
rz(-1.3272746) q[2];
rz(1.7257388) q[3];
sx q[3];
rz(-3.0684107) q[3];
sx q[3];
rz(-0.88254005) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13407229) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(-0.19675955) q[0];
rz(2.9255731) q[1];
sx q[1];
rz(-2.1327503) q[1];
sx q[1];
rz(0.78741995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2741873) q[0];
sx q[0];
rz(-1.7599154) q[0];
sx q[0];
rz(1.8195137) q[0];
x q[1];
rz(-0.51047275) q[2];
sx q[2];
rz(-1.831448) q[2];
sx q[2];
rz(0.54768054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0328788) q[1];
sx q[1];
rz(-1.5252264) q[1];
sx q[1];
rz(2.5497863) q[1];
x q[2];
rz(2.9522252) q[3];
sx q[3];
rz(-1.8872617) q[3];
sx q[3];
rz(-2.6031983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6009723) q[2];
sx q[2];
rz(-1.629402) q[2];
sx q[2];
rz(-2.7579894) q[2];
rz(0.90841928) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(0.15373357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28021321) q[0];
sx q[0];
rz(-0.52556831) q[0];
sx q[0];
rz(-0.63661611) q[0];
rz(2.5909297) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(2.6689463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4387242) q[0];
sx q[0];
rz(-1.6596377) q[0];
sx q[0];
rz(-2.4636561) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51690312) q[2];
sx q[2];
rz(-3.0344525) q[2];
sx q[2];
rz(0.60984367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0090985) q[1];
sx q[1];
rz(-0.14202296) q[1];
sx q[1];
rz(0.36722398) q[1];
rz(1.8697463) q[3];
sx q[3];
rz(-0.96072324) q[3];
sx q[3];
rz(1.3879734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7134646) q[2];
sx q[2];
rz(-1.2423923) q[2];
sx q[2];
rz(3.0430651) q[2];
rz(3.1108917) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(0.15392412) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(0.0064119617) q[1];
sx q[1];
rz(-1.0259722) q[1];
sx q[1];
rz(1.4220672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310282) q[0];
sx q[0];
rz(-0.3153261) q[0];
sx q[0];
rz(-2.226693) q[0];
x q[1];
rz(-0.49894519) q[2];
sx q[2];
rz(-1.2483259) q[2];
sx q[2];
rz(-2.6889963) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8274424) q[1];
sx q[1];
rz(-1.536751) q[1];
sx q[1];
rz(0.15383792) q[1];
rz(-pi) q[2];
rz(-3.0109423) q[3];
sx q[3];
rz(-0.52328309) q[3];
sx q[3];
rz(0.16431683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0942568) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(2.6195841) q[2];
rz(-3.1212854) q[3];
sx q[3];
rz(-0.69965196) q[3];
sx q[3];
rz(-0.29683963) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8727528) q[0];
sx q[0];
rz(-1.5248542) q[0];
sx q[0];
rz(-2.010345) q[0];
rz(0.77992431) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(0.47526971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72544725) q[0];
sx q[0];
rz(-0.83270458) q[0];
sx q[0];
rz(-1.0698143) q[0];
rz(-pi) q[1];
rz(0.90311994) q[2];
sx q[2];
rz(-1.3851056) q[2];
sx q[2];
rz(2.8201495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1058029) q[1];
sx q[1];
rz(-0.33638182) q[1];
sx q[1];
rz(0.85806429) q[1];
rz(-pi) q[2];
rz(-2.4451431) q[3];
sx q[3];
rz(-1.6828732) q[3];
sx q[3];
rz(-0.96986412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37107006) q[2];
sx q[2];
rz(-1.5658242) q[2];
sx q[2];
rz(1.2062262) q[2];
rz(-1.7462339) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(-3.0622845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410011) q[0];
sx q[0];
rz(-0.19484367) q[0];
sx q[0];
rz(2.3051443) q[0];
rz(2.1438228) q[1];
sx q[1];
rz(-2.0496968) q[1];
sx q[1];
rz(-3.0539378) q[1];
rz(-1.2682512) q[2];
sx q[2];
rz(-0.22845636) q[2];
sx q[2];
rz(-0.34194389) q[2];
rz(-1.4597803) q[3];
sx q[3];
rz(-2.2465977) q[3];
sx q[3];
rz(-2.6096596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
