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
rz(1.7537533) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(4.4603577) q[1];
sx q[1];
rz(9.7108884) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7095735) q[0];
sx q[0];
rz(-0.46105024) q[0];
sx q[0];
rz(3.0877355) q[0];
rz(-pi) q[1];
rz(0.082290351) q[2];
sx q[2];
rz(-0.9613885) q[2];
sx q[2];
rz(-2.1235958) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58264153) q[1];
sx q[1];
rz(-2.1182334) q[1];
sx q[1];
rz(1.2657341) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.055130868) q[3];
sx q[3];
rz(-1.1305439) q[3];
sx q[3];
rz(-2.3326479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.046595786) q[2];
sx q[2];
rz(-1.6511714) q[2];
sx q[2];
rz(-0.68520927) q[2];
rz(-1.4677216) q[3];
sx q[3];
rz(-1.0823715) q[3];
sx q[3];
rz(-2.8732324) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14107038) q[0];
sx q[0];
rz(-1.5474316) q[0];
sx q[0];
rz(-0.94456124) q[0];
rz(3.0682849) q[1];
sx q[1];
rz(-2.3088375) q[1];
sx q[1];
rz(1.8836969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8507754) q[0];
sx q[0];
rz(-0.75277872) q[0];
sx q[0];
rz(0.87429177) q[0];
x q[1];
rz(0.33690764) q[2];
sx q[2];
rz(-0.55977976) q[2];
sx q[2];
rz(-0.98613662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1546109) q[1];
sx q[1];
rz(-1.715259) q[1];
sx q[1];
rz(-0.39244907) q[1];
x q[2];
rz(3.1044534) q[3];
sx q[3];
rz(-1.6172774) q[3];
sx q[3];
rz(0.4872191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32626095) q[2];
sx q[2];
rz(-1.3911284) q[2];
sx q[2];
rz(-2.1168671) q[2];
rz(2.4552086) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(0.91275233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82655418) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(0.81845534) q[0];
rz(1.2089027) q[1];
sx q[1];
rz(-1.6345638) q[1];
sx q[1];
rz(-2.2817629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53308177) q[0];
sx q[0];
rz(-1.9010592) q[0];
sx q[0];
rz(2.1950366) q[0];
x q[1];
rz(-0.71975885) q[2];
sx q[2];
rz(-1.9406332) q[2];
sx q[2];
rz(1.0074266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8739955) q[1];
sx q[1];
rz(-1.1073879) q[1];
sx q[1];
rz(2.401176) q[1];
x q[2];
rz(-2.0399666) q[3];
sx q[3];
rz(-2.8896324) q[3];
sx q[3];
rz(2.2445238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7439421) q[2];
sx q[2];
rz(-2.8068779) q[2];
sx q[2];
rz(3.0261377) q[2];
rz(-0.3913106) q[3];
sx q[3];
rz(-1.0262998) q[3];
sx q[3];
rz(-1.1439884) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35599577) q[0];
sx q[0];
rz(-1.9556671) q[0];
sx q[0];
rz(0.24056973) q[0];
rz(-1.3661512) q[1];
sx q[1];
rz(-1.6484478) q[1];
sx q[1];
rz(1.8919401) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8842953) q[0];
sx q[0];
rz(-0.69778555) q[0];
sx q[0];
rz(-1.1691514) q[0];
rz(0.35677783) q[2];
sx q[2];
rz(-2.1564004) q[2];
sx q[2];
rz(2.9139589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0084997) q[1];
sx q[1];
rz(-1.805882) q[1];
sx q[1];
rz(0.66621589) q[1];
rz(-pi) q[2];
rz(1.0650738) q[3];
sx q[3];
rz(-0.95129644) q[3];
sx q[3];
rz(0.92403417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8852641) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(-1.7735927) q[2];
rz(2.8374425) q[3];
sx q[3];
rz(-1.5658028) q[3];
sx q[3];
rz(-2.4450541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2279385) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(0.53713334) q[0];
rz(2.3917603) q[1];
sx q[1];
rz(-1.0822108) q[1];
sx q[1];
rz(0.73712635) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1340807) q[0];
sx q[0];
rz(-1.5069557) q[0];
sx q[0];
rz(0.095973936) q[0];
rz(-2.0290899) q[2];
sx q[2];
rz(-2.9470561) q[2];
sx q[2];
rz(3.1032423) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18857107) q[1];
sx q[1];
rz(-0.28975315) q[1];
sx q[1];
rz(-2.6798331) q[1];
x q[2];
rz(-0.26925663) q[3];
sx q[3];
rz(-0.33447166) q[3];
sx q[3];
rz(-1.0448898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4445112) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(0.14661655) q[2];
rz(1.692449) q[3];
sx q[3];
rz(-1.2796947) q[3];
sx q[3];
rz(-0.52391887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61831063) q[0];
sx q[0];
rz(-2.1250516) q[0];
sx q[0];
rz(-1.4215533) q[0];
rz(-1.2591741) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(-0.4037942) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.161927) q[0];
sx q[0];
rz(-0.84507361) q[0];
sx q[0];
rz(2.7620074) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5228259) q[2];
sx q[2];
rz(-2.8293489) q[2];
sx q[2];
rz(2.00756) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5055949) q[1];
sx q[1];
rz(-0.86252102) q[1];
sx q[1];
rz(2.3227033) q[1];
x q[2];
rz(-0.6882235) q[3];
sx q[3];
rz(-1.7587708) q[3];
sx q[3];
rz(-1.926762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8518565) q[2];
sx q[2];
rz(-1.4137665) q[2];
sx q[2];
rz(2.9750321) q[2];
rz(-1.6377431) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(-3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.903776) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(-2.26407) q[0];
rz(-0.96218836) q[1];
sx q[1];
rz(-2.018237) q[1];
sx q[1];
rz(-0.35433623) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0577045) q[0];
sx q[0];
rz(-1.5885) q[0];
sx q[0];
rz(-2.903865) q[0];
x q[1];
rz(-1.4117175) q[2];
sx q[2];
rz(-2.9949778) q[2];
sx q[2];
rz(-2.7377812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-7*pi/15) q[1];
sx q[1];
rz(-1.0364658) q[1];
sx q[1];
rz(-2.2855845) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50994344) q[3];
sx q[3];
rz(-1.5048319) q[3];
sx q[3];
rz(3.0086781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3366036) q[2];
sx q[2];
rz(-2.4038834) q[2];
sx q[2];
rz(2.0965516) q[2];
rz(-2.7392144) q[3];
sx q[3];
rz(-1.9741524) q[3];
sx q[3];
rz(1.4761402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104915) q[0];
sx q[0];
rz(-0.11678188) q[0];
sx q[0];
rz(-2.2905599) q[0];
rz(0.35375133) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(-2.2869349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8806649) q[0];
sx q[0];
rz(-1.3890084) q[0];
sx q[0];
rz(-0.28566912) q[0];
rz(-0.53422569) q[2];
sx q[2];
rz(-2.8534128) q[2];
sx q[2];
rz(-2.2792918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9820199) q[1];
sx q[1];
rz(-2.5347485) q[1];
sx q[1];
rz(2.8446376) q[1];
rz(2.952997) q[3];
sx q[3];
rz(-2.7321987) q[3];
sx q[3];
rz(0.90367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9510368) q[2];
sx q[2];
rz(-0.64728105) q[2];
sx q[2];
rz(0.62038842) q[2];
rz(0.22647151) q[3];
sx q[3];
rz(-1.7445931) q[3];
sx q[3];
rz(2.5605719) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9613551) q[0];
sx q[0];
rz(-1.2018452) q[0];
sx q[0];
rz(0.015137976) q[0];
rz(1.0221647) q[1];
sx q[1];
rz(-0.41003761) q[1];
sx q[1];
rz(1.9415564) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7040492) q[0];
sx q[0];
rz(-1.3945197) q[0];
sx q[0];
rz(0.86221077) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59320088) q[2];
sx q[2];
rz(-1.5216547) q[2];
sx q[2];
rz(1.0047439) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.697534) q[1];
sx q[1];
rz(-1.6007533) q[1];
sx q[1];
rz(-2.0824964) q[1];
rz(0.8127549) q[3];
sx q[3];
rz(-2.4852607) q[3];
sx q[3];
rz(-1.9539273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6156561) q[2];
sx q[2];
rz(-0.91292149) q[2];
sx q[2];
rz(-0.23942648) q[2];
rz(3.0827403) q[3];
sx q[3];
rz(-0.81164304) q[3];
sx q[3];
rz(2.8656901) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9489768) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(0.97491997) q[0];
rz(2.4661567) q[1];
sx q[1];
rz(-0.91921872) q[1];
sx q[1];
rz(0.98178896) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6639121) q[0];
sx q[0];
rz(-0.64723158) q[0];
sx q[0];
rz(-2.0114698) q[0];
x q[1];
rz(-1.6390332) q[2];
sx q[2];
rz(-1.1600798) q[2];
sx q[2];
rz(-2.1221184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0412933) q[1];
sx q[1];
rz(-0.86242341) q[1];
sx q[1];
rz(-0.34379509) q[1];
x q[2];
rz(-2.0174867) q[3];
sx q[3];
rz(-1.5367931) q[3];
sx q[3];
rz(-1.407935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56741095) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(1.2249472) q[2];
rz(0.36568668) q[3];
sx q[3];
rz(-0.94958011) q[3];
sx q[3];
rz(1.1398116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097261978) q[0];
sx q[0];
rz(-1.8987569) q[0];
sx q[0];
rz(1.4415997) q[0];
rz(3.0615831) q[1];
sx q[1];
rz(-1.7623822) q[1];
sx q[1];
rz(-0.0026127876) q[1];
rz(1.4314797) q[2];
sx q[2];
rz(-2.0124544) q[2];
sx q[2];
rz(-2.6228946) q[2];
rz(-0.046924165) q[3];
sx q[3];
rz(-2.6886446) q[3];
sx q[3];
rz(-0.038084134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
