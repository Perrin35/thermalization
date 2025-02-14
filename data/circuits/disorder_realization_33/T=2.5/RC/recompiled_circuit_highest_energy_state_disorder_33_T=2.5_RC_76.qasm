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
rz(-2.9774732) q[0];
sx q[0];
rz(-1.0241221) q[0];
sx q[0];
rz(2.6574988) q[0];
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(0.054952316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.84961) q[0];
sx q[0];
rz(-1.5395482) q[0];
sx q[0];
rz(0.90051767) q[0];
rz(-pi) q[1];
rz(-0.058637549) q[2];
sx q[2];
rz(-1.3212886) q[2];
sx q[2];
rz(-2.3326258) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.37147) q[1];
sx q[1];
rz(-2.9767163) q[1];
sx q[1];
rz(-2.5046964) q[1];
rz(-pi) q[2];
rz(2.0255768) q[3];
sx q[3];
rz(-0.6729799) q[3];
sx q[3];
rz(2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.310828) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(-2.8669299) q[2];
rz(2.2667387) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(-2.8680475) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4817568) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(2.7428395) q[0];
rz(-0.2440456) q[1];
sx q[1];
rz(-1.2363385) q[1];
sx q[1];
rz(2.0358548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6080792) q[0];
sx q[0];
rz(-1.7987804) q[0];
sx q[0];
rz(-0.45324676) q[0];
rz(-pi) q[1];
rz(-0.91715889) q[2];
sx q[2];
rz(-2.0431402) q[2];
sx q[2];
rz(-1.2231959) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0827786) q[1];
sx q[1];
rz(-1.5681439) q[1];
sx q[1];
rz(3.1383731) q[1];
rz(1.3040172) q[3];
sx q[3];
rz(-1.6872129) q[3];
sx q[3];
rz(-2.4698225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(-2.0665118) q[2];
rz(-0.69617802) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(-0.16304326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657846) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(1.3801344) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3724646) q[0];
sx q[0];
rz(-1.6004246) q[0];
sx q[0];
rz(1.6086786) q[0];
rz(-pi) q[1];
rz(1.9335494) q[2];
sx q[2];
rz(-0.13158509) q[2];
sx q[2];
rz(-3.0818444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8113831) q[1];
sx q[1];
rz(-2.7683677) q[1];
sx q[1];
rz(-0.95287816) q[1];
x q[2];
rz(-0.61010078) q[3];
sx q[3];
rz(-1.0786295) q[3];
sx q[3];
rz(0.51612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8850024) q[2];
sx q[2];
rz(-1.8424748) q[2];
sx q[2];
rz(-3.0510862) q[2];
rz(2.6835486) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(-2.0335782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-0.93165398) q[0];
rz(2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(2.1021252) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9477978) q[0];
sx q[0];
rz(-0.76776869) q[0];
sx q[0];
rz(-2.7916273) q[0];
rz(-pi) q[1];
rz(-1.8245071) q[2];
sx q[2];
rz(-1.7895964) q[2];
sx q[2];
rz(-1.4267061) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89371496) q[1];
sx q[1];
rz(-1.3637241) q[1];
sx q[1];
rz(-0.061092579) q[1];
x q[2];
rz(2.5040197) q[3];
sx q[3];
rz(-0.3488144) q[3];
sx q[3];
rz(1.3138527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6197551) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(-2.6452765) q[2];
rz(-2.3450092) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-0.97187483) q[0];
rz(-0.12174363) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(-1.8121388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9864101) q[0];
sx q[0];
rz(-1.2647795) q[0];
sx q[0];
rz(3.0772692) q[0];
rz(-pi) q[1];
rz(2.2586063) q[2];
sx q[2];
rz(-0.60740439) q[2];
sx q[2];
rz(-0.6907874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5962474) q[1];
sx q[1];
rz(-1.1584499) q[1];
sx q[1];
rz(2.418787) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.132344) q[3];
sx q[3];
rz(-1.4917177) q[3];
sx q[3];
rz(1.6759863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(-2.9902048) q[2];
rz(-0.76465145) q[3];
sx q[3];
rz(-2.305856) q[3];
sx q[3];
rz(1.5449272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0695892) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(1.0714916) q[0];
rz(-0.53120652) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(2.5154617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6036442) q[0];
sx q[0];
rz(-1.5660362) q[0];
sx q[0];
rz(-1.5877923) q[0];
rz(2.7886969) q[2];
sx q[2];
rz(-2.5145686) q[2];
sx q[2];
rz(-2.3924946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8529359) q[1];
sx q[1];
rz(-2.1737338) q[1];
sx q[1];
rz(1.1974242) q[1];
x q[2];
rz(-2.4247384) q[3];
sx q[3];
rz(-2.5156913) q[3];
sx q[3];
rz(0.38960534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35817394) q[2];
sx q[2];
rz(-2.5219315) q[2];
sx q[2];
rz(-1.0142856) q[2];
rz(-1.2554393) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(-2.2210806) q[0];
rz(0.11407425) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(-1.1309518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3944954) q[0];
sx q[0];
rz(-0.97309006) q[0];
sx q[0];
rz(1.5639485) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5860687) q[2];
sx q[2];
rz(-1.0706524) q[2];
sx q[2];
rz(-1.9490567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56173813) q[1];
sx q[1];
rz(-0.67719141) q[1];
sx q[1];
rz(-1.3586292) q[1];
x q[2];
rz(-1.4907964) q[3];
sx q[3];
rz(-1.0558075) q[3];
sx q[3];
rz(1.966371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7941234) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(-0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(1.1387775) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(0.29079944) q[0];
rz(0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(-1.1999757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3481846) q[0];
sx q[0];
rz(-0.68773341) q[0];
sx q[0];
rz(2.4853112) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2924215) q[2];
sx q[2];
rz(-2.2361922) q[2];
sx q[2];
rz(2.9862491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0713745) q[1];
sx q[1];
rz(-1.5244251) q[1];
sx q[1];
rz(1.8646311) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7539976) q[3];
sx q[3];
rz(-1.6699413) q[3];
sx q[3];
rz(2.9354036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(1.7928436) q[2];
rz(1.6382943) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1143484) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(-1.0741023) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(-0.46357402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.539042) q[0];
sx q[0];
rz(-0.21115223) q[0];
sx q[0];
rz(-2.0402914) q[0];
x q[1];
rz(2.4899269) q[2];
sx q[2];
rz(-0.32512384) q[2];
sx q[2];
rz(-0.59949694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9755755) q[1];
sx q[1];
rz(-2.3834627) q[1];
sx q[1];
rz(-2.1339971) q[1];
rz(-pi) q[2];
rz(2.998407) q[3];
sx q[3];
rz(-2.5783263) q[3];
sx q[3];
rz(1.3384502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8672436) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(-0.80844936) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(0.43420473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7782068) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(-0.038381902) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-0.29676944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1365741) q[0];
sx q[0];
rz(-1.8474192) q[0];
sx q[0];
rz(1.0199976) q[0];
rz(-pi) q[1];
rz(0.77905853) q[2];
sx q[2];
rz(-0.59528661) q[2];
sx q[2];
rz(-2.1434458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1047639) q[1];
sx q[1];
rz(-1.4955758) q[1];
sx q[1];
rz(2.9192218) q[1];
x q[2];
rz(-2.2567883) q[3];
sx q[3];
rz(-1.8821041) q[3];
sx q[3];
rz(1.6303568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9639637) q[2];
sx q[2];
rz(-2.0202426) q[2];
sx q[2];
rz(-2.5753944) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(1.1894777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(-0.26168564) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(-2.8520201) q[2];
sx q[2];
rz(-0.085554335) q[2];
sx q[2];
rz(2.6680996) q[2];
rz(2.49558) q[3];
sx q[3];
rz(-0.96886841) q[3];
sx q[3];
rz(-2.3652707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
