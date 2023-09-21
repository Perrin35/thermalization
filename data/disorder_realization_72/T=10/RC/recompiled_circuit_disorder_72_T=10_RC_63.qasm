OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7802785) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(1.6367903) q[0];
rz(-2.9688641) q[2];
sx q[2];
rz(-0.85731259) q[2];
sx q[2];
rz(1.2781065) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0608622) q[1];
sx q[1];
rz(-0.84134358) q[1];
sx q[1];
rz(1.7727477) q[1];
x q[2];
rz(-1.3189089) q[3];
sx q[3];
rz(-3.0348274) q[3];
sx q[3];
rz(-1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(-1.0027592) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28264499) q[0];
sx q[0];
rz(-2.1574321) q[0];
sx q[0];
rz(-1.0413917) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3147896) q[2];
sx q[2];
rz(-2.8892592) q[2];
sx q[2];
rz(-1.8581055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7748389) q[1];
sx q[1];
rz(-0.60791053) q[1];
sx q[1];
rz(-1.2971836) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5428883) q[3];
sx q[3];
rz(-1.8667392) q[3];
sx q[3];
rz(0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(2.2580106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2471836) q[0];
sx q[0];
rz(-0.43617019) q[0];
sx q[0];
rz(-0.39788525) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8475464) q[2];
sx q[2];
rz(-1.1883192) q[2];
sx q[2];
rz(1.8665714) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59144943) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(1.8557465) q[1];
rz(-pi) q[2];
x q[2];
rz(2.482445) q[3];
sx q[3];
rz(-1.4775606) q[3];
sx q[3];
rz(2.6574709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(-0.61169949) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54116762) q[0];
sx q[0];
rz(-2.4006002) q[0];
sx q[0];
rz(-0.11359544) q[0];
x q[1];
rz(1.6646531) q[2];
sx q[2];
rz(-1.8641194) q[2];
sx q[2];
rz(-0.23767995) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6039227) q[1];
sx q[1];
rz(-0.99891716) q[1];
sx q[1];
rz(-3.0043383) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9228966) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(2.4179539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(0.27553976) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(0.87669796) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6633776) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(2.1992654) q[0];
x q[1];
rz(-1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(0.63477883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.7261506) q[1];
sx q[1];
rz(-0.032226493) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74113412) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-0.65628624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148465) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2561589) q[0];
sx q[0];
rz(-1.116917) q[0];
sx q[0];
rz(-2.8594349) q[0];
x q[1];
rz(2.6655212) q[2];
sx q[2];
rz(-1.8487612) q[2];
sx q[2];
rz(-1.4897886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.035316) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(1.7229401) q[1];
x q[2];
rz(-0.35950426) q[3];
sx q[3];
rz(-2.5264969) q[3];
sx q[3];
rz(0.86715992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-0.82715183) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9046017) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(0.29445946) q[0];
rz(-1.0625213) q[2];
sx q[2];
rz(-1.695343) q[2];
sx q[2];
rz(-1.6342271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78325242) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(-pi) q[2];
rz(1.8282312) q[3];
sx q[3];
rz(-2.7830887) q[3];
sx q[3];
rz(-0.076971019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-2.9879925) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(-2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(3.1088366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2643124) q[0];
sx q[0];
rz(-2.7240629) q[0];
sx q[0];
rz(-2.2420922) q[0];
rz(2.1342437) q[2];
sx q[2];
rz(-1.1903561) q[2];
sx q[2];
rz(3.0055339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0667463) q[1];
sx q[1];
rz(-0.42145887) q[1];
sx q[1];
rz(-1.7154109) q[1];
rz(-pi) q[2];
rz(-1.1912187) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(2.6759202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1774566) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-0.39045236) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(-0.73293066) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-2.656235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982614) q[0];
sx q[0];
rz(-1.4364463) q[0];
sx q[0];
rz(-1.6249379) q[0];
rz(-pi) q[1];
rz(-1.4114686) q[2];
sx q[2];
rz(-0.97852856) q[2];
sx q[2];
rz(1.9221905) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(1.1361213) q[1];
rz(0.58166196) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(-1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52206802) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-2.5433345) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82829581) q[0];
sx q[0];
rz(-1.7772355) q[0];
sx q[0];
rz(-2.9199827) q[0];
rz(1.7190476) q[2];
sx q[2];
rz(-1.3607927) q[2];
sx q[2];
rz(-1.6809747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2559402) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(0.10140681) q[1];
rz(-pi) q[2];
rz(1.4062905) q[3];
sx q[3];
rz(-1.7269772) q[3];
sx q[3];
rz(2.3567049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-0.65199488) q[2];
rz(-2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(-2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3020637) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(1.6937704) q[2];
sx q[2];
rz(-0.95477827) q[2];
sx q[2];
rz(1.9670602) q[2];
rz(0.89647722) q[3];
sx q[3];
rz(-1.8835246) q[3];
sx q[3];
rz(-2.6261331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];