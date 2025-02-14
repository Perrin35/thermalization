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
rz(1.7776547) q[0];
sx q[0];
rz(-0.41843709) q[0];
sx q[0];
rz(-1.3016181) q[0];
rz(0.16297451) q[1];
sx q[1];
rz(-1.4459223) q[1];
sx q[1];
rz(0.016782848) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.421464) q[0];
sx q[0];
rz(-1.4961494) q[0];
sx q[0];
rz(1.2089085) q[0];
rz(-pi) q[1];
rz(1.4799825) q[2];
sx q[2];
rz(-2.9613284) q[2];
sx q[2];
rz(1.6423051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7960733) q[1];
sx q[1];
rz(-1.5658448) q[1];
sx q[1];
rz(-1.5657182) q[1];
x q[2];
rz(2.7990325) q[3];
sx q[3];
rz(-0.16157074) q[3];
sx q[3];
rz(-1.8752961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(1.5622697) q[2];
rz(-0.2114547) q[3];
sx q[3];
rz(-3.1410757) q[3];
sx q[3];
rz(0.16099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.529539) q[0];
sx q[0];
rz(-0.27612975) q[0];
sx q[0];
rz(-1.8147234) q[0];
rz(-0.57948411) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(-0.63900596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20637437) q[0];
sx q[0];
rz(-2.2274096) q[0];
sx q[0];
rz(-0.73120631) q[0];
x q[1];
rz(1.5540358) q[2];
sx q[2];
rz(-1.4498561) q[2];
sx q[2];
rz(3.1230833) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42334786) q[1];
sx q[1];
rz(-3.1210174) q[1];
sx q[1];
rz(-0.59845509) q[1];
rz(-pi) q[2];
rz(2.3611446) q[3];
sx q[3];
rz(-1.6326346) q[3];
sx q[3];
rz(-2.0831747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7063893) q[2];
sx q[2];
rz(-0.13627626) q[2];
sx q[2];
rz(-1.5380247) q[2];
rz(1.5806574) q[3];
sx q[3];
rz(-0.014336421) q[3];
sx q[3];
rz(3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.2494217) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(0.38145915) q[0];
rz(2.4341266) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(1.1245419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1046421) q[0];
sx q[0];
rz(-2.7890887) q[0];
sx q[0];
rz(-0.80926766) q[0];
rz(-pi) q[1];
rz(-3.1156024) q[2];
sx q[2];
rz(-1.455869) q[2];
sx q[2];
rz(0.1564125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4545817) q[1];
sx q[1];
rz(-3.0774766) q[1];
sx q[1];
rz(-0.18973501) q[1];
rz(-1.0550523) q[3];
sx q[3];
rz(-2.0655144) q[3];
sx q[3];
rz(-2.5124541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4418929) q[2];
sx q[2];
rz(-0.012233891) q[2];
sx q[2];
rz(-3.0921248) q[2];
rz(-2.5345645) q[3];
sx q[3];
rz(-0.0012461239) q[3];
sx q[3];
rz(-1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1599051) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(-0.017177563) q[0];
rz(-2.848564) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(-1.5471829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38774224) q[0];
sx q[0];
rz(-2.1374636) q[0];
sx q[0];
rz(2.500534) q[0];
x q[1];
rz(1.3457005) q[2];
sx q[2];
rz(-1.0402816) q[2];
sx q[2];
rz(-0.38166416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80611011) q[1];
sx q[1];
rz(-1.4467738) q[1];
sx q[1];
rz(-1.6308484) q[1];
x q[2];
rz(2.4956216) q[3];
sx q[3];
rz(-1.5785909) q[3];
sx q[3];
rz(2.2199059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25345099) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(0.68224254) q[2];
rz(0.055179723) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(1.8365708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76160112) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(0.4250266) q[0];
rz(1.5399326) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(-2.3262598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3638301) q[0];
sx q[0];
rz(-1.5095995) q[0];
sx q[0];
rz(1.5817002) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4311976) q[2];
sx q[2];
rz(-0.023471467) q[2];
sx q[2];
rz(-1.0271629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7338431) q[1];
sx q[1];
rz(-1.6899278) q[1];
sx q[1];
rz(-1.6581737) q[1];
rz(1.4493222) q[3];
sx q[3];
rz(-1.9125347) q[3];
sx q[3];
rz(-2.4499144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1347947) q[2];
sx q[2];
rz(-0.012601348) q[2];
sx q[2];
rz(-1.4726144) q[2];
rz(0.93004477) q[3];
sx q[3];
rz(-0.01447066) q[3];
sx q[3];
rz(0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8188266) q[0];
sx q[0];
rz(-0.023094026) q[0];
sx q[0];
rz(-1.4267138) q[0];
rz(-0.73000437) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(2.0402562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0520605) q[0];
sx q[0];
rz(-2.2666078) q[0];
sx q[0];
rz(0.3541644) q[0];
rz(-pi) q[1];
rz(2.1997994) q[2];
sx q[2];
rz(-2.859349) q[2];
sx q[2];
rz(-1.4636702) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8423398) q[1];
sx q[1];
rz(-1.4758759) q[1];
sx q[1];
rz(1.4466982) q[1];
rz(-pi) q[2];
rz(-0.032208431) q[3];
sx q[3];
rz(-1.695249) q[3];
sx q[3];
rz(-2.2197753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4250028) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(-1.8343743) q[2];
rz(-2.763125) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(-0.65346658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3017479) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(-0.92754716) q[0];
rz(-1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(-1.6102788) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46857014) q[0];
sx q[0];
rz(-3.0971905) q[0];
sx q[0];
rz(-1.0442249) q[0];
x q[1];
rz(-0.39674098) q[2];
sx q[2];
rz(-1.7974241) q[2];
sx q[2];
rz(2.6643348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.130333) q[1];
sx q[1];
rz(-1.5681055) q[1];
sx q[1];
rz(1.4417159) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0965936) q[3];
sx q[3];
rz(-0.57588314) q[3];
sx q[3];
rz(-0.27920846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7808468) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(-1.2537664) q[2];
rz(0.69416657) q[3];
sx q[3];
rz(-2.4024051) q[3];
sx q[3];
rz(2.9085801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6041782) q[0];
sx q[0];
rz(-1.0056714) q[0];
sx q[0];
rz(-2.1042714) q[0];
rz(-1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(1.6745837) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59975043) q[0];
sx q[0];
rz(-1.6381725) q[0];
sx q[0];
rz(-1.3753433) q[0];
rz(-1.59378) q[2];
sx q[2];
rz(-1.8488374) q[2];
sx q[2];
rz(2.8756623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8407708) q[1];
sx q[1];
rz(-1.5719218) q[1];
sx q[1];
rz(-3.1411533) q[1];
rz(0.45536228) q[3];
sx q[3];
rz(-1.4858559) q[3];
sx q[3];
rz(-0.37173879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1303225) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(-0.043896349) q[2];
rz(-2.6210426) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(1.4588149) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7875824) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(-1.822923) q[0];
rz(-1.4175381) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(1.5444548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70001665) q[0];
sx q[0];
rz(-1.0927943) q[0];
sx q[0];
rz(-1.4411323) q[0];
rz(-0.33716069) q[2];
sx q[2];
rz(-2.4558407) q[2];
sx q[2];
rz(0.27039385) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3785951) q[1];
sx q[1];
rz(-1.7471658) q[1];
sx q[1];
rz(-1.8862944) q[1];
rz(-pi) q[2];
rz(1.0472337) q[3];
sx q[3];
rz(-0.95253836) q[3];
sx q[3];
rz(-2.0766192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80457193) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(-2.9447832) q[2];
rz(-1.9491516) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(0.20096745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8404959) q[0];
sx q[0];
rz(-1.3571285) q[0];
sx q[0];
rz(1.1916196) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-2.4947417) q[1];
sx q[1];
rz(1.5764538) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7708262) q[0];
sx q[0];
rz(-1.2773371) q[0];
sx q[0];
rz(1.7459041) q[0];
rz(2.0247518) q[2];
sx q[2];
rz(-1.7773787) q[2];
sx q[2];
rz(-1.1517186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2597771) q[1];
sx q[1];
rz(-1.5712954) q[1];
sx q[1];
rz(3.1406443) q[1];
x q[2];
rz(-1.4286391) q[3];
sx q[3];
rz(-1.6831493) q[3];
sx q[3];
rz(1.4893099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9578751) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(1.4584165) q[2];
rz(0.030473907) q[3];
sx q[3];
rz(-0.009549791) q[3];
sx q[3];
rz(2.9392346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771582) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(-1.5674113) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(-0.0049736918) q[2];
sx q[2];
rz(-1.4952352) q[2];
sx q[2];
rz(0.16005439) q[2];
rz(-0.81672698) q[3];
sx q[3];
rz(-1.9533659) q[3];
sx q[3];
rz(1.6537651) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
