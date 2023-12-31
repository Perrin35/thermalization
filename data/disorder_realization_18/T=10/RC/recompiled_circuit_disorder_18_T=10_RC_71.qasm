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
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0639609) q[0];
sx q[0];
rz(-1.504717) q[0];
sx q[0];
rz(-1.8128916) q[0];
rz(-pi) q[1];
rz(-2.8561864) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(-2.0602496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0393385) q[1];
sx q[1];
rz(-1.551911) q[1];
sx q[1];
rz(-2.0233872) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6712816) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(2.2944962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9053099) q[0];
sx q[0];
rz(-2.5295527) q[0];
sx q[0];
rz(2.406714) q[0];
rz(-pi) q[1];
rz(-1.8663835) q[2];
sx q[2];
rz(-2.327773) q[2];
sx q[2];
rz(-2.5676167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21287316) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(2.5961155) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3470207) q[3];
sx q[3];
rz(-1.9850529) q[3];
sx q[3];
rz(-1.8300717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(-0.71933293) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2704724) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8910599) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(3.1191349) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40076077) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(1.00373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1151162) q[1];
sx q[1];
rz(-0.83362245) q[1];
sx q[1];
rz(-1.923418) q[1];
rz(-pi) q[2];
rz(1.2445883) q[3];
sx q[3];
rz(-2.2766114) q[3];
sx q[3];
rz(-1.0969485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(2.1380651) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-2.3153268) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4429312) q[0];
sx q[0];
rz(-2.1974265) q[0];
sx q[0];
rz(0.71039623) q[0];
x q[1];
rz(-1.516953) q[2];
sx q[2];
rz(-0.87140897) q[2];
sx q[2];
rz(2.2211423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77376765) q[1];
sx q[1];
rz(-0.34562472) q[1];
sx q[1];
rz(-1.3454382) q[1];
rz(-pi) q[2];
rz(-1.9598947) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(-0.94483313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(-2.5240135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78695801) q[0];
sx q[0];
rz(-1.5411396) q[0];
sx q[0];
rz(-1.6164854) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0134301) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(0.15304676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6286271) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(-2.1668424) q[1];
rz(-pi) q[2];
rz(3.0693552) q[3];
sx q[3];
rz(-2.9655955) q[3];
sx q[3];
rz(2.3776059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3072309) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(0.93528549) q[0];
rz(-pi) q[1];
rz(1.8413576) q[2];
sx q[2];
rz(-2.5205043) q[2];
sx q[2];
rz(1.5121216) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.60702) q[1];
sx q[1];
rz(-1.8587451) q[1];
sx q[1];
rz(1.401591) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6052386) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(-0.29355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.9656666) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.4181597) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-3.1013536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38567625) q[0];
sx q[0];
rz(-1.3323116) q[0];
sx q[0];
rz(-0.083906108) q[0];
rz(-pi) q[1];
rz(-1.6778498) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(1.6183311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.28133) q[1];
sx q[1];
rz(-1.4308235) q[1];
sx q[1];
rz(-2.1920188) q[1];
rz(-pi) q[2];
rz(2.8242565) q[3];
sx q[3];
rz(-1.279497) q[3];
sx q[3];
rz(-0.81277646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(2.2824536) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7751986) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(-1.226107) q[0];
x q[1];
rz(-2.3178188) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(-1.6549695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1464403) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(0.050683024) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0309615) q[3];
sx q[3];
rz(-2.4880829) q[3];
sx q[3];
rz(-1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(1.2649149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7495959) q[0];
sx q[0];
rz(-2.2889334) q[0];
sx q[0];
rz(-1.0438265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5876797) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(-1.6162789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7281108) q[1];
sx q[1];
rz(-0.15004798) q[1];
sx q[1];
rz(0.24197443) q[1];
rz(2.7183652) q[3];
sx q[3];
rz(-1.755135) q[3];
sx q[3];
rz(-1.9870027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(2.5850885) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(-2.8961704) q[0];
rz(0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.6773178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3905555) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(-1.4575973) q[0];
x q[1];
rz(-2.1572838) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(3.0629326) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6357395) q[1];
sx q[1];
rz(-0.74843279) q[1];
sx q[1];
rz(1.9745419) q[1];
rz(-pi) q[2];
rz(1.3947992) q[3];
sx q[3];
rz(-1.6547852) q[3];
sx q[3];
rz(1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8650919) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(2.4297498) q[2];
sx q[2];
rz(-1.1887475) q[2];
sx q[2];
rz(2.9954994) q[2];
rz(-2.5682156) q[3];
sx q[3];
rz(-1.8918512) q[3];
sx q[3];
rz(0.20791114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
