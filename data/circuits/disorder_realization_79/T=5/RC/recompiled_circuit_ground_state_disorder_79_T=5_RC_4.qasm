OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.76685846) q[0];
sx q[0];
rz(-0.91349608) q[0];
sx q[0];
rz(-1.3878393) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(-1.8228276) q[1];
sx q[1];
rz(0.28611046) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090534276) q[0];
sx q[0];
rz(-1.5947475) q[0];
sx q[0];
rz(-2.6811203) q[0];
rz(-pi) q[1];
rz(-2.181796) q[2];
sx q[2];
rz(-1.6382484) q[2];
sx q[2];
rz(0.59997257) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5589511) q[1];
sx q[1];
rz(-2.1182334) q[1];
sx q[1];
rz(-1.2657341) q[1];
rz(2.0116352) q[3];
sx q[3];
rz(-1.5209271) q[3];
sx q[3];
rz(2.4032556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.046595786) q[2];
sx q[2];
rz(-1.6511714) q[2];
sx q[2];
rz(2.4563834) q[2];
rz(1.6738711) q[3];
sx q[3];
rz(-1.0823715) q[3];
sx q[3];
rz(-2.8732324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.0005223) q[0];
sx q[0];
rz(-1.594161) q[0];
sx q[0];
rz(0.94456124) q[0];
rz(3.0682849) q[1];
sx q[1];
rz(-2.3088375) q[1];
sx q[1];
rz(1.8836969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2908173) q[0];
sx q[0];
rz(-0.75277872) q[0];
sx q[0];
rz(-0.87429177) q[0];
x q[1];
rz(2.6075105) q[2];
sx q[2];
rz(-1.3943496) q[2];
sx q[2];
rz(-0.29613972) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2230665) q[1];
sx q[1];
rz(-0.41689532) q[1];
sx q[1];
rz(-0.36347632) q[1];
rz(0.03713921) q[3];
sx q[3];
rz(-1.5243153) q[3];
sx q[3];
rz(-2.6543736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8153317) q[2];
sx q[2];
rz(-1.7504642) q[2];
sx q[2];
rz(2.1168671) q[2];
rz(-0.68638408) q[3];
sx q[3];
rz(-0.73203433) q[3];
sx q[3];
rz(-0.91275233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82655418) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(-2.3231373) q[0];
rz(-1.2089027) q[1];
sx q[1];
rz(-1.6345638) q[1];
sx q[1];
rz(2.2817629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5268974) q[0];
sx q[0];
rz(-2.4458168) q[0];
sx q[0];
rz(1.0403364) q[0];
x q[1];
rz(2.6099714) q[2];
sx q[2];
rz(-0.79382703) q[2];
sx q[2];
rz(0.17228157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2261994) q[1];
sx q[1];
rz(-2.2186154) q[1];
sx q[1];
rz(-2.1658705) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1016261) q[3];
sx q[3];
rz(-0.2519603) q[3];
sx q[3];
rz(-0.89706883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7439421) q[2];
sx q[2];
rz(-0.33471477) q[2];
sx q[2];
rz(-0.115455) q[2];
rz(2.7502821) q[3];
sx q[3];
rz(-2.1152928) q[3];
sx q[3];
rz(1.1439884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35599577) q[0];
sx q[0];
rz(-1.1859256) q[0];
sx q[0];
rz(2.9010229) q[0];
rz(1.3661512) q[1];
sx q[1];
rz(-1.4931449) q[1];
sx q[1];
rz(1.8919401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8842953) q[0];
sx q[0];
rz(-0.69778555) q[0];
sx q[0];
rz(-1.9724413) q[0];
rz(-pi) q[1];
rz(2.1867238) q[2];
sx q[2];
rz(-1.2754585) q[2];
sx q[2];
rz(-1.5952641) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27423672) q[1];
sx q[1];
rz(-2.4411267) q[1];
sx q[1];
rz(0.3697311) q[1];
x q[2];
rz(2.0765188) q[3];
sx q[3];
rz(-2.1902962) q[3];
sx q[3];
rz(0.92403417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8852641) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(-1.3679999) q[2];
rz(-0.3041501) q[3];
sx q[3];
rz(-1.5757898) q[3];
sx q[3];
rz(-0.69653851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91365415) q[0];
sx q[0];
rz(-0.21622394) q[0];
sx q[0];
rz(-2.6044593) q[0];
rz(-2.3917603) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(0.73712635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1340807) q[0];
sx q[0];
rz(-1.5069557) q[0];
sx q[0];
rz(0.095973936) q[0];
rz(-pi) q[1];
rz(0.086949172) q[2];
sx q[2];
rz(-1.3965522) q[2];
sx q[2];
rz(0.50424313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29038844) q[1];
sx q[1];
rz(-1.3121288) q[1];
sx q[1];
rz(-1.438739) q[1];
rz(0.26925663) q[3];
sx q[3];
rz(-0.33447166) q[3];
sx q[3];
rz(-2.0967029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6970814) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(2.9949761) q[2];
rz(1.4491436) q[3];
sx q[3];
rz(-1.2796947) q[3];
sx q[3];
rz(0.52391887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61831063) q[0];
sx q[0];
rz(-1.016541) q[0];
sx q[0];
rz(-1.4215533) q[0];
rz(-1.8824185) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(0.4037942) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84995364) q[0];
sx q[0];
rz(-1.8516415) q[0];
sx q[0];
rz(2.3333059) q[0];
x q[1];
rz(-1.7558891) q[2];
sx q[2];
rz(-1.3178692) q[2];
sx q[2];
rz(1.3651939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48133367) q[1];
sx q[1];
rz(-2.1161882) q[1];
sx q[1];
rz(-0.86472269) q[1];
x q[2];
rz(1.3293224) q[3];
sx q[3];
rz(-0.89697402) q[3];
sx q[3];
rz(-0.50844565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8518565) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(0.16656052) q[2];
rz(1.5038495) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(-3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903776) q[0];
sx q[0];
rz(-2.7523968) q[0];
sx q[0];
rz(-2.26407) q[0];
rz(-2.1794043) q[1];
sx q[1];
rz(-1.1233556) q[1];
sx q[1];
rz(2.7872564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503949) q[0];
sx q[0];
rz(-1.3331067) q[0];
sx q[0];
rz(1.5890122) q[0];
rz(-1.4260068) q[2];
sx q[2];
rz(-1.5476523) q[2];
sx q[2];
rz(1.0095846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-7*pi/15) q[1];
sx q[1];
rz(-2.1051268) q[1];
sx q[1];
rz(-0.85600812) q[1];
rz(-pi) q[2];
rz(1.4952502) q[3];
sx q[3];
rz(-1.0620688) q[3];
sx q[3];
rz(-1.7405603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80498901) q[2];
sx q[2];
rz(-2.4038834) q[2];
sx q[2];
rz(2.0965516) q[2];
rz(-2.7392144) q[3];
sx q[3];
rz(-1.9741524) q[3];
sx q[3];
rz(-1.6654525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104915) q[0];
sx q[0];
rz(-3.0248108) q[0];
sx q[0];
rz(0.85103273) q[0];
rz(-0.35375133) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(-0.85465777) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2609278) q[0];
sx q[0];
rz(-1.3890084) q[0];
sx q[0];
rz(0.28566912) q[0];
x q[1];
rz(1.7206011) q[2];
sx q[2];
rz(-1.8179107) q[2];
sx q[2];
rz(-2.832156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6254297) q[1];
sx q[1];
rz(-0.99402797) q[1];
sx q[1];
rz(-1.7712084) q[1];
rz(-pi) q[2];
rz(2.952997) q[3];
sx q[3];
rz(-0.40939399) q[3];
sx q[3];
rz(2.2379217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19055584) q[2];
sx q[2];
rz(-0.64728105) q[2];
sx q[2];
rz(2.5212042) q[2];
rz(2.9151211) q[3];
sx q[3];
rz(-1.3969996) q[3];
sx q[3];
rz(-0.58102077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1802375) q[0];
sx q[0];
rz(-1.2018452) q[0];
sx q[0];
rz(3.1264547) q[0];
rz(-1.0221647) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(1.9415564) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7040492) q[0];
sx q[0];
rz(-1.747073) q[0];
sx q[0];
rz(-2.2793819) q[0];
rz(-0.087751919) q[2];
sx q[2];
rz(-2.5466036) q[2];
sx q[2];
rz(2.648271) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.697534) q[1];
sx q[1];
rz(-1.5408393) q[1];
sx q[1];
rz(-2.0824964) q[1];
rz(2.6545908) q[3];
sx q[3];
rz(-1.111711) q[3];
sx q[3];
rz(-2.8279357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6156561) q[2];
sx q[2];
rz(-2.2286712) q[2];
sx q[2];
rz(0.23942648) q[2];
rz(3.0827403) q[3];
sx q[3];
rz(-2.3299496) q[3];
sx q[3];
rz(0.27590251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1926159) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(-2.1666727) q[0];
rz(-0.67543593) q[1];
sx q[1];
rz(-0.91921872) q[1];
sx q[1];
rz(0.98178896) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4083207) q[0];
sx q[0];
rz(-1.3106723) q[0];
sx q[0];
rz(0.97113804) q[0];
rz(-1.6390332) q[2];
sx q[2];
rz(-1.1600798) q[2];
sx q[2];
rz(1.0194743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60335873) q[1];
sx q[1];
rz(-0.77419821) q[1];
sx q[1];
rz(-1.9456359) q[1];
x q[2];
rz(1.124106) q[3];
sx q[3];
rz(-1.5367931) q[3];
sx q[3];
rz(1.7336577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56741095) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(-1.9166454) q[2];
rz(0.36568668) q[3];
sx q[3];
rz(-2.1920125) q[3];
sx q[3];
rz(2.001781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097261978) q[0];
sx q[0];
rz(-1.2428357) q[0];
sx q[0];
rz(-1.6999929) q[0];
rz(0.080009566) q[1];
sx q[1];
rz(-1.3792104) q[1];
sx q[1];
rz(3.1389799) q[1];
rz(-0.44543191) q[2];
sx q[2];
rz(-1.4449228) q[2];
sx q[2];
rz(-1.1119643) q[2];
rz(-0.45251485) q[3];
sx q[3];
rz(-1.5913251) q[3];
sx q[3];
rz(-1.6510788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
