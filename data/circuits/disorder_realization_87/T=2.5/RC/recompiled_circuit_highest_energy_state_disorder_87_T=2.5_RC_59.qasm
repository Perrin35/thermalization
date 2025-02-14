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
rz(-0.63646746) q[0];
sx q[0];
rz(-0.31390733) q[0];
sx q[0];
rz(1.3482345) q[0];
rz(-0.95797602) q[1];
sx q[1];
rz(-1.4161243) q[1];
sx q[1];
rz(2.9834566) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414171) q[0];
sx q[0];
rz(-0.99208562) q[0];
sx q[0];
rz(0.35037033) q[0];
x q[1];
rz(1.8991437) q[2];
sx q[2];
rz(-1.2109887) q[2];
sx q[2];
rz(1.7361602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0709405) q[1];
sx q[1];
rz(-1.5717984) q[1];
sx q[1];
rz(1.5654501) q[1];
x q[2];
rz(-2.4918658) q[3];
sx q[3];
rz(-0.5818682) q[3];
sx q[3];
rz(0.054903809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8378975) q[2];
sx q[2];
rz(-2.1968696) q[2];
sx q[2];
rz(2.5486805) q[2];
rz(-2.8494075) q[3];
sx q[3];
rz(-3.1226776) q[3];
sx q[3];
rz(1.7570447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050046571) q[0];
sx q[0];
rz(-0.36588359) q[0];
sx q[0];
rz(-0.094060913) q[0];
rz(1.3919818) q[1];
sx q[1];
rz(-1.6110907) q[1];
sx q[1];
rz(1.403341) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4499324) q[0];
sx q[0];
rz(-0.34087718) q[0];
sx q[0];
rz(0.65861146) q[0];
rz(-0.22461598) q[2];
sx q[2];
rz(-2.086528) q[2];
sx q[2];
rz(1.7040229) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5557068) q[1];
sx q[1];
rz(-0.60190869) q[1];
sx q[1];
rz(-1.5762171) q[1];
rz(2.9280065) q[3];
sx q[3];
rz(-0.41700577) q[3];
sx q[3];
rz(-0.97975376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2469108) q[2];
sx q[2];
rz(-0.44591388) q[2];
sx q[2];
rz(-1.1723588) q[2];
rz(0.42919484) q[3];
sx q[3];
rz(-0.49002886) q[3];
sx q[3];
rz(1.4303077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2369075) q[0];
sx q[0];
rz(-2.560736) q[0];
sx q[0];
rz(0.3592321) q[0];
rz(1.5776186) q[1];
sx q[1];
rz(-0.79505316) q[1];
sx q[1];
rz(-0.94594812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1157784) q[0];
sx q[0];
rz(-0.17318586) q[0];
sx q[0];
rz(-1.0758297) q[0];
x q[1];
rz(-1.8976413) q[2];
sx q[2];
rz(-3.0209401) q[2];
sx q[2];
rz(3.0337703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7316389) q[1];
sx q[1];
rz(-1.8112858) q[1];
sx q[1];
rz(0.13786836) q[1];
x q[2];
rz(2.4427857) q[3];
sx q[3];
rz(-0.10891373) q[3];
sx q[3];
rz(-1.004899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49950162) q[2];
sx q[2];
rz(-1.5311798) q[2];
sx q[2];
rz(-2.0384608) q[2];
rz(-1.057386) q[3];
sx q[3];
rz(-1.5921581) q[3];
sx q[3];
rz(-2.9252083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0353521) q[0];
sx q[0];
rz(-3.048936) q[0];
sx q[0];
rz(-2.4308391) q[0];
rz(0.6657998) q[1];
sx q[1];
rz(-3.1328821) q[1];
sx q[1];
rz(-2.8361368) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9269954) q[0];
sx q[0];
rz(-1.4844795) q[0];
sx q[0];
rz(-1.5284365) q[0];
rz(0.50263672) q[2];
sx q[2];
rz(-1.2900616) q[2];
sx q[2];
rz(-1.0681149) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7240746) q[1];
sx q[1];
rz(-1.4913847) q[1];
sx q[1];
rz(-2.0168138) q[1];
rz(-pi) q[2];
rz(-0.908522) q[3];
sx q[3];
rz(-1.2161331) q[3];
sx q[3];
rz(-0.42820546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1315883) q[2];
sx q[2];
rz(-1.5350716) q[2];
sx q[2];
rz(0.32001495) q[2];
rz(-0.53712505) q[3];
sx q[3];
rz(-0.38083005) q[3];
sx q[3];
rz(2.2297458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.491275) q[0];
sx q[0];
rz(-2.9170051) q[0];
sx q[0];
rz(0.085163072) q[0];
rz(0.59146178) q[1];
sx q[1];
rz(-0.0031091212) q[1];
sx q[1];
rz(-1.9689781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.216476) q[0];
sx q[0];
rz(-1.5763349) q[0];
sx q[0];
rz(-0.001222697) q[0];
rz(-pi) q[1];
rz(1.8109365) q[2];
sx q[2];
rz(-1.2954172) q[2];
sx q[2];
rz(0.20698838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56192815) q[1];
sx q[1];
rz(-0.85612002) q[1];
sx q[1];
rz(0.90335937) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8211604) q[3];
sx q[3];
rz(-0.51051192) q[3];
sx q[3];
rz(1.7701469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11599596) q[2];
sx q[2];
rz(-1.7946578) q[2];
sx q[2];
rz(-1.4541413) q[2];
rz(-1.8986374) q[3];
sx q[3];
rz(-0.60296139) q[3];
sx q[3];
rz(-0.81237826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1989307) q[0];
sx q[0];
rz(-0.16348612) q[0];
sx q[0];
rz(1.7401975) q[0];
rz(2.5247848) q[1];
sx q[1];
rz(-3.1251833) q[1];
sx q[1];
rz(-1.1161463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6907731) q[0];
sx q[0];
rz(-1.4754533) q[0];
sx q[0];
rz(-1.3006163) q[0];
rz(-pi) q[1];
rz(1.1405476) q[2];
sx q[2];
rz(-0.7459085) q[2];
sx q[2];
rz(-0.416428) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68758724) q[1];
sx q[1];
rz(-1.0024888) q[1];
sx q[1];
rz(2.1623728) q[1];
rz(-pi) q[2];
rz(-1.7839637) q[3];
sx q[3];
rz(-1.1305332) q[3];
sx q[3];
rz(-0.21808392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8650032) q[2];
sx q[2];
rz(-1.208655) q[2];
sx q[2];
rz(1.2651944) q[2];
rz(0.64722925) q[3];
sx q[3];
rz(-2.6914458) q[3];
sx q[3];
rz(1.2559206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7728421) q[0];
sx q[0];
rz(-2.5970646) q[0];
sx q[0];
rz(-2.7640589) q[0];
rz(0.15097161) q[1];
sx q[1];
rz(-3.1307463) q[1];
sx q[1];
rz(2.6735701) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0220874) q[0];
sx q[0];
rz(-2.6900253) q[0];
sx q[0];
rz(-3.0877583) q[0];
x q[1];
rz(1.4281316) q[2];
sx q[2];
rz(-2.4141867) q[2];
sx q[2];
rz(-2.3886566) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49761729) q[1];
sx q[1];
rz(-0.68174284) q[1];
sx q[1];
rz(1.4585182) q[1];
rz(1.9028038) q[3];
sx q[3];
rz(-2.598437) q[3];
sx q[3];
rz(-2.4623722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2630792) q[2];
sx q[2];
rz(-0.24363467) q[2];
sx q[2];
rz(1.1070975) q[2];
rz(0.88362068) q[3];
sx q[3];
rz(-0.7928018) q[3];
sx q[3];
rz(-0.71225524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20473075) q[0];
sx q[0];
rz(-2.2624113) q[0];
sx q[0];
rz(0.19464807) q[0];
rz(-1.6814992) q[1];
sx q[1];
rz(-3.1250521) q[1];
sx q[1];
rz(0.3009235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60754362) q[0];
sx q[0];
rz(-0.87042394) q[0];
sx q[0];
rz(2.0714633) q[0];
rz(-pi) q[1];
rz(-2.373502) q[2];
sx q[2];
rz(-0.67736191) q[2];
sx q[2];
rz(-2.0887665) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49773945) q[1];
sx q[1];
rz(-2.5723638) q[1];
sx q[1];
rz(2.7124972) q[1];
rz(-pi) q[2];
rz(-0.31105403) q[3];
sx q[3];
rz(-2.9533953) q[3];
sx q[3];
rz(1.8194906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68227565) q[2];
sx q[2];
rz(-1.514491) q[2];
sx q[2];
rz(-0.21919361) q[2];
rz(-1.3466262) q[3];
sx q[3];
rz(-3.0248088) q[3];
sx q[3];
rz(2.35675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.858736) q[0];
sx q[0];
rz(-2.5114926) q[0];
sx q[0];
rz(-3.1334738) q[0];
rz(-0.66932622) q[1];
sx q[1];
rz(-3.1287584) q[1];
sx q[1];
rz(2.377811) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8618906) q[0];
sx q[0];
rz(-2.6548626) q[0];
sx q[0];
rz(-0.52409117) q[0];
rz(-pi) q[1];
rz(1.6449459) q[2];
sx q[2];
rz(-2.1762245) q[2];
sx q[2];
rz(-2.2260967) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8593401) q[1];
sx q[1];
rz(-1.0251004) q[1];
sx q[1];
rz(0.061435862) q[1];
rz(-pi) q[2];
rz(0.64371806) q[3];
sx q[3];
rz(-0.9399842) q[3];
sx q[3];
rz(-0.35177059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1179489) q[2];
sx q[2];
rz(-3.1274319) q[2];
sx q[2];
rz(-2.1421049) q[2];
rz(-0.1855447) q[3];
sx q[3];
rz(-2.1047635) q[3];
sx q[3];
rz(-0.015983494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9284116) q[0];
sx q[0];
rz(-0.11225926) q[0];
sx q[0];
rz(2.4768594) q[0];
rz(-2.1809273) q[1];
sx q[1];
rz(-3.101109) q[1];
sx q[1];
rz(1.5047081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0069160865) q[0];
sx q[0];
rz(-1.8728372) q[0];
sx q[0];
rz(2.9535314) q[0];
rz(3.0111752) q[2];
sx q[2];
rz(-0.36389029) q[2];
sx q[2];
rz(2.7478349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3469228) q[1];
sx q[1];
rz(-1.1058013) q[1];
sx q[1];
rz(2.5645606) q[1];
rz(-pi) q[2];
rz(0.89791132) q[3];
sx q[3];
rz(-0.6403044) q[3];
sx q[3];
rz(1.698157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8399743) q[2];
sx q[2];
rz(-3.1349389) q[2];
sx q[2];
rz(1.9807695) q[2];
rz(-3.0617359) q[3];
sx q[3];
rz(-3.1294332) q[3];
sx q[3];
rz(1.9399835) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53497159) q[0];
sx q[0];
rz(-1.4749682) q[0];
sx q[0];
rz(-1.5780021) q[0];
rz(3.1040991) q[1];
sx q[1];
rz(-2.9480724) q[1];
sx q[1];
rz(-2.9185157) q[1];
rz(1.846215) q[2];
sx q[2];
rz(-2.7559235) q[2];
sx q[2];
rz(-1.4758103) q[2];
rz(-0.8998733) q[3];
sx q[3];
rz(-1.5597349) q[3];
sx q[3];
rz(-2.8928383) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
