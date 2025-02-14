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
rz(1.8472449) q[0];
sx q[0];
rz(-1.6142774) q[0];
sx q[0];
rz(2.9543028) q[0];
rz(-1.3297431) q[1];
sx q[1];
rz(3.6880479) q[1];
sx q[1];
rz(10.772279) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9634099) q[0];
sx q[0];
rz(-1.4926433) q[0];
sx q[0];
rz(-2.1341538) q[0];
x q[1];
rz(2.0894089) q[2];
sx q[2];
rz(-1.901327) q[2];
sx q[2];
rz(-1.7180819) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.627003) q[1];
sx q[1];
rz(-1.5545067) q[1];
sx q[1];
rz(-0.88912995) q[1];
x q[2];
rz(2.4745117) q[3];
sx q[3];
rz(-0.32235185) q[3];
sx q[3];
rz(-1.9867112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61907855) q[2];
sx q[2];
rz(-1.1828902) q[2];
sx q[2];
rz(0.28156933) q[2];
rz(1.36739) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(-0.27845964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6047769) q[0];
sx q[0];
rz(-2.8474658) q[0];
sx q[0];
rz(-1.3500704) q[0];
rz(3.0916832) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(1.8972338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45634633) q[0];
sx q[0];
rz(-0.35402997) q[0];
sx q[0];
rz(-2.9732605) q[0];
rz(0.18644615) q[2];
sx q[2];
rz(-1.0294339) q[2];
sx q[2];
rz(-0.10554927) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10319732) q[1];
sx q[1];
rz(-1.8382676) q[1];
sx q[1];
rz(-2.9385376) q[1];
rz(-pi) q[2];
rz(2.5998551) q[3];
sx q[3];
rz(-1.1835464) q[3];
sx q[3];
rz(0.70808402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0656978) q[2];
sx q[2];
rz(-1.9588797) q[2];
sx q[2];
rz(0.069570216) q[2];
rz(2.4075497) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(0.61843425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.0760536) q[0];
sx q[0];
rz(-2.9422308) q[0];
sx q[0];
rz(-2.8850436) q[0];
rz(1.8007295) q[1];
sx q[1];
rz(-0.95756617) q[1];
sx q[1];
rz(1.0440913) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279914) q[0];
sx q[0];
rz(-1.2577857) q[0];
sx q[0];
rz(2.9646216) q[0];
x q[1];
rz(1.3568722) q[2];
sx q[2];
rz(-1.7905777) q[2];
sx q[2];
rz(2.8666988) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8736834) q[1];
sx q[1];
rz(-1.2699155) q[1];
sx q[1];
rz(-1.2734702) q[1];
x q[2];
rz(0.99756156) q[3];
sx q[3];
rz(-1.6896393) q[3];
sx q[3];
rz(-3.1021743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23036817) q[2];
sx q[2];
rz(-2.0679097) q[2];
sx q[2];
rz(-2.4252841) q[2];
rz(-2.6835942) q[3];
sx q[3];
rz(-1.4665946) q[3];
sx q[3];
rz(2.9108293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79391795) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(2.4804261) q[0];
rz(1.2541153) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(1.5987781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4371915) q[0];
sx q[0];
rz(-1.3953623) q[0];
sx q[0];
rz(-1.8262922) q[0];
x q[1];
rz(0.91868742) q[2];
sx q[2];
rz(-0.67011681) q[2];
sx q[2];
rz(2.5958259) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.143121) q[1];
sx q[1];
rz(-0.70618343) q[1];
sx q[1];
rz(-0.097340214) q[1];
x q[2];
rz(0.38971735) q[3];
sx q[3];
rz(-1.0709312) q[3];
sx q[3];
rz(-0.12382774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8864112) q[2];
sx q[2];
rz(-2.7913783) q[2];
sx q[2];
rz(1.1345471) q[2];
rz(-0.55401951) q[3];
sx q[3];
rz(-1.7787245) q[3];
sx q[3];
rz(2.5567283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8447113) q[0];
sx q[0];
rz(-3.1394594) q[0];
sx q[0];
rz(-1.1312436) q[0];
rz(-0.058111195) q[1];
sx q[1];
rz(-1.2429712) q[1];
sx q[1];
rz(-1.6910472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7052225) q[0];
sx q[0];
rz(-2.3793325) q[0];
sx q[0];
rz(-0.99135474) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31016953) q[2];
sx q[2];
rz(-2.355994) q[2];
sx q[2];
rz(-2.4082662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1102396) q[1];
sx q[1];
rz(-1.807537) q[1];
sx q[1];
rz(-2.8048022) q[1];
rz(-pi) q[2];
rz(2.2048619) q[3];
sx q[3];
rz(-1.9104244) q[3];
sx q[3];
rz(0.32138164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6821373) q[2];
sx q[2];
rz(-2.5294999) q[2];
sx q[2];
rz(1.3539782) q[2];
rz(2.5891384) q[3];
sx q[3];
rz(-0.8907291) q[3];
sx q[3];
rz(-0.49828211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8950997) q[0];
sx q[0];
rz(-0.95588481) q[0];
sx q[0];
rz(-2.0434875) q[0];
rz(-1.2753298) q[1];
sx q[1];
rz(-1.2397436) q[1];
sx q[1];
rz(-2.1001508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5811036) q[0];
sx q[0];
rz(-2.1745918) q[0];
sx q[0];
rz(0.93224705) q[0];
x q[1];
rz(-1.6026233) q[2];
sx q[2];
rz(-1.4518445) q[2];
sx q[2];
rz(-2.5980189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90064657) q[1];
sx q[1];
rz(-0.52370549) q[1];
sx q[1];
rz(-1.3670735) q[1];
rz(0.84534332) q[3];
sx q[3];
rz(-1.3535557) q[3];
sx q[3];
rz(-0.30423924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8424524) q[2];
sx q[2];
rz(-1.1536529) q[2];
sx q[2];
rz(3.0843206) q[2];
rz(-0.20036571) q[3];
sx q[3];
rz(-2.5860131) q[3];
sx q[3];
rz(-2.949775) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361236) q[0];
sx q[0];
rz(-2.496688) q[0];
sx q[0];
rz(2.7591925) q[0];
rz(-0.96757013) q[1];
sx q[1];
rz(-2.7303374) q[1];
sx q[1];
rz(1.8866906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9318274) q[0];
sx q[0];
rz(-1.5146274) q[0];
sx q[0];
rz(0.61135354) q[0];
rz(-pi) q[1];
rz(-0.81315984) q[2];
sx q[2];
rz(-0.81612464) q[2];
sx q[2];
rz(0.94674158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4285415) q[1];
sx q[1];
rz(-1.9425434) q[1];
sx q[1];
rz(1.1357689) q[1];
rz(-pi) q[2];
rz(-1.2953877) q[3];
sx q[3];
rz(-2.3304238) q[3];
sx q[3];
rz(-1.0397183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3249698) q[2];
sx q[2];
rz(-2.7866456) q[2];
sx q[2];
rz(2.2897282) q[2];
rz(2.8999117) q[3];
sx q[3];
rz(-1.9342187) q[3];
sx q[3];
rz(2.2309301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333176) q[0];
sx q[0];
rz(-0.18263826) q[0];
sx q[0];
rz(1.0688548) q[0];
rz(1.764027) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(-2.4527803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4877164) q[0];
sx q[0];
rz(-2.6453995) q[0];
sx q[0];
rz(2.3997593) q[0];
rz(-pi) q[1];
rz(-2.942932) q[2];
sx q[2];
rz(-0.61214329) q[2];
sx q[2];
rz(-0.3224626) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2548267) q[1];
sx q[1];
rz(-0.51023645) q[1];
sx q[1];
rz(2.4429823) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6121871) q[3];
sx q[3];
rz(-1.7200618) q[3];
sx q[3];
rz(-0.68219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0869202) q[2];
sx q[2];
rz(-0.72267795) q[2];
sx q[2];
rz(1.5642536) q[2];
rz(-0.65586048) q[3];
sx q[3];
rz(-1.7003931) q[3];
sx q[3];
rz(-2.3724469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5250788) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(-2.9922564) q[0];
rz(2.0556045) q[1];
sx q[1];
rz(-2.1859152) q[1];
sx q[1];
rz(-2.65926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59583658) q[0];
sx q[0];
rz(-2.0424423) q[0];
sx q[0];
rz(-0.095998569) q[0];
x q[1];
rz(-1.8033837) q[2];
sx q[2];
rz(-1.5778189) q[2];
sx q[2];
rz(2.1178031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73585549) q[1];
sx q[1];
rz(-2.2453826) q[1];
sx q[1];
rz(1.9231741) q[1];
rz(-pi) q[2];
rz(-1.6685358) q[3];
sx q[3];
rz(-1.1527138) q[3];
sx q[3];
rz(-0.12471499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0444191) q[2];
sx q[2];
rz(-2.3046875) q[2];
sx q[2];
rz(0.4591628) q[2];
rz(-0.098371355) q[3];
sx q[3];
rz(-1.2854853) q[3];
sx q[3];
rz(-1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761053) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(-3.0011445) q[0];
rz(-1.2175995) q[1];
sx q[1];
rz(-2.0001037) q[1];
sx q[1];
rz(-1.5906895) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80181304) q[0];
sx q[0];
rz(-2.411541) q[0];
sx q[0];
rz(0.12224893) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9687499) q[2];
sx q[2];
rz(-1.3982492) q[2];
sx q[2];
rz(2.3177528) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1223476) q[1];
sx q[1];
rz(-0.68210318) q[1];
sx q[1];
rz(-0.42241272) q[1];
rz(-0.012926558) q[3];
sx q[3];
rz(-0.89996808) q[3];
sx q[3];
rz(-2.304043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6259674) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(1.1930126) q[2];
rz(2.6595645) q[3];
sx q[3];
rz(-2.7218282) q[3];
sx q[3];
rz(2.1217864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1358418) q[0];
sx q[0];
rz(-2.0477722) q[0];
sx q[0];
rz(2.3070106) q[0];
rz(0.76345481) q[1];
sx q[1];
rz(-1.3044985) q[1];
sx q[1];
rz(0.34674092) q[1];
rz(-1.687526) q[2];
sx q[2];
rz(-1.5925831) q[2];
sx q[2];
rz(1.5769373) q[2];
rz(-2.5498123) q[3];
sx q[3];
rz(-2.4587678) q[3];
sx q[3];
rz(-3.126161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
