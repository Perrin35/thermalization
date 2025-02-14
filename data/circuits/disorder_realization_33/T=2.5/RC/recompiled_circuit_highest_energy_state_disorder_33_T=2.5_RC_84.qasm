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
rz(-0.48409387) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(4.7197309) q[1];
sx q[1];
rz(9.3698256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2919827) q[0];
sx q[0];
rz(-1.6020444) q[0];
sx q[0];
rz(-2.241075) q[0];
rz(-3.0829551) q[2];
sx q[2];
rz(-1.8203041) q[2];
sx q[2];
rz(-2.3326258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72803264) q[1];
sx q[1];
rz(-1.4384585) q[1];
sx q[1];
rz(1.4721666) q[1];
x q[2];
rz(2.1922429) q[3];
sx q[3];
rz(-1.2934522) q[3];
sx q[3];
rz(-2.8327033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.310828) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(-0.87485391) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65983588) q[0];
sx q[0];
rz(-2.4517224) q[0];
sx q[0];
rz(0.39875317) q[0];
rz(2.8975471) q[1];
sx q[1];
rz(-1.2363385) q[1];
sx q[1];
rz(2.0358548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.669848) q[0];
sx q[0];
rz(-2.6378184) q[0];
sx q[0];
rz(-2.6543447) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2695838) q[2];
sx q[2];
rz(-0.78561831) q[2];
sx q[2];
rz(-0.18839041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51197375) q[1];
sx q[1];
rz(-1.5675768) q[1];
sx q[1];
rz(-1.5734488) q[1];
rz(0.12064528) q[3];
sx q[3];
rz(-1.3058666) q[3];
sx q[3];
rz(0.93075965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4411321) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(-2.0665118) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(0.16304326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27580801) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.3801344) q[1];
sx q[1];
rz(-0.41788995) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19720896) q[0];
sx q[0];
rz(-1.5329307) q[0];
sx q[0];
rz(0.029649563) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0946629) q[2];
sx q[2];
rz(-1.4478193) q[2];
sx q[2];
rz(0.3058946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9668321) q[1];
sx q[1];
rz(-1.357954) q[1];
sx q[1];
rz(1.2618466) q[1];
x q[2];
rz(-2.1501174) q[3];
sx q[3];
rz(-1.0414755) q[3];
sx q[3];
rz(0.73562276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(3.0510862) q[2];
rz(0.45804405) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8103545) q[0];
sx q[0];
rz(-0.88307035) q[0];
sx q[0];
rz(2.2099387) q[0];
rz(-0.16904198) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(1.0394675) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0213892) q[0];
sx q[0];
rz(-1.3303555) q[0];
sx q[0];
rz(2.405015) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8245071) q[2];
sx q[2];
rz(-1.7895964) q[2];
sx q[2];
rz(-1.7148866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.66450602) q[1];
sx q[1];
rz(-1.6305822) q[1];
sx q[1];
rz(1.7782446) q[1];
rz(-2.5040197) q[3];
sx q[3];
rz(-0.3488144) q[3];
sx q[3];
rz(-1.3138527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52183759) q[2];
sx q[2];
rz(-2.0401185) q[2];
sx q[2];
rz(-2.6452765) q[2];
rz(-0.79658341) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(-1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8747044) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(-0.97187483) q[0];
rz(-0.12174363) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(1.8121388) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15518256) q[0];
sx q[0];
rz(-1.8768132) q[0];
sx q[0];
rz(3.0772692) q[0];
rz(2.0636286) q[2];
sx q[2];
rz(-1.9415641) q[2];
sx q[2];
rz(2.855122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59898224) q[1];
sx q[1];
rz(-0.81331454) q[1];
sx q[1];
rz(-0.58425649) q[1];
rz(1.6869808) q[3];
sx q[3];
rz(-3.0619762) q[3];
sx q[3];
rz(1.5594359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(2.9902048) q[2];
rz(2.3769412) q[3];
sx q[3];
rz(-2.305856) q[3];
sx q[3];
rz(1.5449272) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0695892) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(2.0701011) q[0];
rz(0.53120652) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(-2.5154617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3817978) q[0];
sx q[0];
rz(-0.017649895) q[0];
sx q[0];
rz(1.2977029) q[0];
x q[1];
rz(0.35289571) q[2];
sx q[2];
rz(-2.5145686) q[2];
sx q[2];
rz(-0.74909808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8257013) q[1];
sx q[1];
rz(-2.4448312) q[1];
sx q[1];
rz(2.6543762) q[1];
x q[2];
rz(0.71685426) q[3];
sx q[3];
rz(-0.62590137) q[3];
sx q[3];
rz(2.7519873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7834187) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-2.1273071) q[2];
rz(1.2554393) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(-1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(2.2210806) q[0];
rz(0.11407425) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(1.1309518) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3944954) q[0];
sx q[0];
rz(-2.1685026) q[0];
sx q[0];
rz(1.5776442) q[0];
x q[1];
rz(-0.80318309) q[2];
sx q[2];
rz(-0.72942643) q[2];
sx q[2];
rz(-0.27952295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29207001) q[1];
sx q[1];
rz(-2.2300868) q[1];
sx q[1];
rz(0.16772049) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5163631) q[3];
sx q[3];
rz(-1.5011906) q[3];
sx q[3];
rz(-2.7854837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(2.7327909) q[2];
rz(-0.18320228) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(-1.1387775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(-2.8507932) q[0];
rz(0.94789061) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(1.1999757) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24058293) q[0];
sx q[0];
rz(-1.173061) q[0];
sx q[0];
rz(0.57698864) q[0];
x q[1];
rz(-2.3340204) q[2];
sx q[2];
rz(-1.0243197) q[2];
sx q[2];
rz(-1.228491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34741286) q[1];
sx q[1];
rz(-0.29736667) q[1];
sx q[1];
rz(-1.7296687) q[1];
rz(-pi) q[2];
rz(-3.0407716) q[3];
sx q[3];
rz(-1.3885048) q[3];
sx q[3];
rz(-1.758648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(1.348749) q[2];
rz(1.6382943) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0272442) q[0];
sx q[0];
rz(-0.37186563) q[0];
sx q[0];
rz(2.0674904) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(2.6780186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6025507) q[0];
sx q[0];
rz(-2.9304404) q[0];
sx q[0];
rz(-1.1013012) q[0];
rz(-0.26185449) q[2];
sx q[2];
rz(-1.7657649) q[2];
sx q[2];
rz(-1.5972114) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8817655) q[1];
sx q[1];
rz(-0.9503839) q[1];
sx q[1];
rz(2.6735191) q[1];
rz(1.660668) q[3];
sx q[3];
rz(-1.0139795) q[3];
sx q[3];
rz(-1.1695605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27434906) q[2];
sx q[2];
rz(-2.5090019) q[2];
sx q[2];
rz(0.80844936) q[2];
rz(-0.48458734) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(-3.1032108) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(0.29676944) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39958056) q[0];
sx q[0];
rz(-2.098408) q[0];
sx q[0];
rz(-0.3216089) q[0];
rz(-pi) q[1];
rz(2.3625341) q[2];
sx q[2];
rz(-2.546306) q[2];
sx q[2];
rz(0.99814683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0368288) q[1];
sx q[1];
rz(-1.4955758) q[1];
sx q[1];
rz(0.22237088) q[1];
rz(-1.1007916) q[3];
sx q[3];
rz(-2.3988225) q[3];
sx q[3];
rz(2.7239885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(2.5753944) q[2];
rz(2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(0.26168564) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(1.546312) q[2];
sx q[2];
rz(-1.4888121) q[2];
sx q[2];
rz(-0.76406995) q[2];
rz(-2.2812609) q[3];
sx q[3];
rz(-2.0900149) q[3];
sx q[3];
rz(2.7505977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
