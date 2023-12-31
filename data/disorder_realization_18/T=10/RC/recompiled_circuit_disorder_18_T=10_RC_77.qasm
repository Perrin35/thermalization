OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(-0.95681325) q[0];
sx q[0];
rz(1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0639609) q[0];
sx q[0];
rz(-1.504717) q[0];
sx q[0];
rz(1.328701) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8561864) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(1.0813431) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6822328) q[1];
sx q[1];
rz(-2.0233005) q[1];
sx q[1];
rz(-3.1205936) q[1];
x q[2];
rz(-0.54361312) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.1528667) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(2.2944962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69764187) q[0];
sx q[0];
rz(-1.1753433) q[0];
sx q[0];
rz(-2.6614499) q[0];
rz(-0.77913021) q[2];
sx q[2];
rz(-1.7841633) q[2];
sx q[2];
rz(1.2029635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6716799) q[1];
sx q[1];
rz(-0.54930733) q[1];
sx q[1];
rz(3.0104907) q[1];
rz(-0.3470207) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(1.8300717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.117924) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741094) q[0];
sx q[0];
rz(-1.2615146) q[0];
sx q[0];
rz(1.563619) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2450561) q[2];
sx q[2];
rz(-0.92391652) q[2];
sx q[2];
rz(1.5145472) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69818245) q[1];
sx q[1];
rz(-1.3122307) q[1];
sx q[1];
rz(2.3727388) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7819488) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(-0.8262659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59506455) q[0];
sx q[0];
rz(-2.1273158) q[0];
sx q[0];
rz(-2.3331649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(2.5259279) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77376765) q[1];
sx q[1];
rz(-0.34562472) q[1];
sx q[1];
rz(1.7961545) q[1];
x q[2];
rz(0.63185933) q[3];
sx q[3];
rz(-1.8903036) q[3];
sx q[3];
rz(0.39998049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78695801) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(-1.6164854) q[0];
rz(-pi) q[1];
rz(-0.56474833) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(-2.2861779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51296556) q[1];
sx q[1];
rz(-1.291853) q[1];
sx q[1];
rz(-2.1668424) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9660452) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(-0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(0.0897952) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(-0.18009137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371532) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(0.8888437) q[0];
x q[1];
rz(-2.9526268) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(-1.8408066) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1485968) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-2.6246043) q[1];
rz(-1.5363541) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(0.29355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95773762) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(0.040239008) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979472) q[0];
sx q[0];
rz(-0.25254927) q[0];
sx q[0];
rz(1.9027684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1646541) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(-0.050780642) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8602627) q[1];
sx q[1];
rz(-1.4308235) q[1];
sx q[1];
rz(-0.94957385) q[1];
x q[2];
rz(1.2651029) q[3];
sx q[3];
rz(-1.8743268) q[3];
sx q[3];
rz(-2.2895253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.5054024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818427) q[0];
sx q[0];
rz(-1.8847701) q[0];
sx q[0];
rz(2.7008675) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3178188) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.6549695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(-0.69617747) q[1];
rz(-pi) q[2];
rz(-2.1721341) q[3];
sx q[3];
rz(-1.8441895) q[3];
sx q[3];
rz(-0.45575842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(0.15052477) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(2.495893) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.2649149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3919968) q[0];
sx q[0];
rz(-2.2889334) q[0];
sx q[0];
rz(1.0438265) q[0];
rz(-pi) q[1];
rz(-2.4931156) q[2];
sx q[2];
rz(-0.65995526) q[2];
sx q[2];
rz(0.49396587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4134819) q[1];
sx q[1];
rz(-2.9915447) q[1];
sx q[1];
rz(2.8996182) q[1];
x q[2];
rz(0.42616578) q[3];
sx q[3];
rz(-0.45939547) q[3];
sx q[3];
rz(-0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.4642749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7510371) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(-1.4575973) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9173031) q[2];
sx q[2];
rz(-2.1456492) q[2];
sx q[2];
rz(1.526051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7732685) q[1];
sx q[1];
rz(-1.3001633) q[1];
sx q[1];
rz(-0.86398217) q[1];
rz(-pi) q[2];
rz(-2.0189832) q[3];
sx q[3];
rz(-0.19482329) q[3];
sx q[3];
rz(2.4217055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-2.0586661) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(-0.57337702) q[3];
sx q[3];
rz(-1.2497414) q[3];
sx q[3];
rz(-2.9336815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
