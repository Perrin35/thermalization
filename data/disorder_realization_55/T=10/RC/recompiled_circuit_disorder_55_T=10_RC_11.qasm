OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(-1.62513) q[0];
sx q[0];
rz(-0.2642785) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(-1.2101313) q[1];
sx q[1];
rz(0.73524737) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7289294) q[0];
sx q[0];
rz(-2.465415) q[0];
sx q[0];
rz(-2.9039608) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60230435) q[2];
sx q[2];
rz(-2.3659083) q[2];
sx q[2];
rz(-1.8204945) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4525675) q[1];
sx q[1];
rz(-1.8568294) q[1];
sx q[1];
rz(-0.16098117) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5427038) q[3];
sx q[3];
rz(-2.1510604) q[3];
sx q[3];
rz(-2.7147646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66951093) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(-2.0377339) q[2];
rz(-1.2708698) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1141777) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(-2.7052178) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(0.26611051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9929745) q[0];
sx q[0];
rz(-0.87478144) q[0];
sx q[0];
rz(-0.50177411) q[0];
rz(1.3219464) q[2];
sx q[2];
rz(-0.58983931) q[2];
sx q[2];
rz(2.313254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3937711) q[1];
sx q[1];
rz(-0.95612477) q[1];
sx q[1];
rz(-2.462903) q[1];
x q[2];
rz(1.3105884) q[3];
sx q[3];
rz(-0.62285138) q[3];
sx q[3];
rz(2.7130896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(-2.3320847) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(-2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.7944638) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1613306) q[0];
sx q[0];
rz(-1.8322924) q[0];
sx q[0];
rz(-2.6677092) q[0];
x q[1];
rz(2.0688829) q[2];
sx q[2];
rz(-0.32215298) q[2];
sx q[2];
rz(1.0108394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.128708) q[1];
sx q[1];
rz(-1.8599659) q[1];
sx q[1];
rz(-2.1916926) q[1];
rz(-pi) q[2];
rz(-0.8637572) q[3];
sx q[3];
rz(-1.45544) q[3];
sx q[3];
rz(-0.74187169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1469664) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.5220801) q[2];
rz(0.26432031) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79214823) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(-2.175892) q[0];
rz(0.72215885) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(0.55975634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4013195) q[0];
sx q[0];
rz(-1.4896605) q[0];
sx q[0];
rz(-2.7149537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83300029) q[2];
sx q[2];
rz(-2.2350395) q[2];
sx q[2];
rz(-2.8766362) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2985845) q[1];
sx q[1];
rz(-2.1857939) q[1];
sx q[1];
rz(3.0857012) q[1];
rz(-pi) q[2];
rz(-2.9066554) q[3];
sx q[3];
rz(-0.80312356) q[3];
sx q[3];
rz(-2.0149751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(-1.4245865) q[2];
rz(0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.150862) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(-2.7752303) q[0];
rz(1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(-0.34367925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3643091) q[0];
sx q[0];
rz(-1.0316348) q[0];
sx q[0];
rz(0.76215141) q[0];
rz(-pi) q[1];
rz(-2.2720488) q[2];
sx q[2];
rz(-2.0336667) q[2];
sx q[2];
rz(-2.9733544) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.9586473) q[1];
sx q[1];
rz(-1.4970386) q[1];
sx q[1];
rz(-1.0720836) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92693365) q[3];
sx q[3];
rz(-1.5537795) q[3];
sx q[3];
rz(3.0683558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3917824) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-0.053744944) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(-3.1164363) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(-0.25973928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1007337) q[0];
sx q[0];
rz(-1.20964) q[0];
sx q[0];
rz(0.33613899) q[0];
rz(-pi) q[1];
rz(-1.514228) q[2];
sx q[2];
rz(-2.8048189) q[2];
sx q[2];
rz(0.77384225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0745084) q[1];
sx q[1];
rz(-1.899252) q[1];
sx q[1];
rz(-0.23898464) q[1];
rz(-pi) q[2];
rz(-2.0201683) q[3];
sx q[3];
rz(-1.5546038) q[3];
sx q[3];
rz(-1.8024973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48866895) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(0.10061131) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1473734) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(2.8314262) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(2.535634) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92354846) q[0];
sx q[0];
rz(-2.179562) q[0];
sx q[0];
rz(-1.8463085) q[0];
rz(-pi) q[1];
rz(0.29166834) q[2];
sx q[2];
rz(-1.3783611) q[2];
sx q[2];
rz(0.33484909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.088674) q[1];
sx q[1];
rz(-2.5175736) q[1];
sx q[1];
rz(-0.96159972) q[1];
rz(-pi) q[2];
rz(-0.53020729) q[3];
sx q[3];
rz(-1.0876417) q[3];
sx q[3];
rz(-0.14164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9014152) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(-2.288738) q[2];
rz(1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119499) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(-1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(-3.0775552) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1916699) q[0];
sx q[0];
rz(-2.4987698) q[0];
sx q[0];
rz(1.886801) q[0];
rz(-pi) q[1];
rz(-1.0649101) q[2];
sx q[2];
rz(-1.3008899) q[2];
sx q[2];
rz(0.050886521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9248326) q[1];
sx q[1];
rz(-1.2158582) q[1];
sx q[1];
rz(1.6096398) q[1];
x q[2];
rz(-1.2588345) q[3];
sx q[3];
rz(-1.6083816) q[3];
sx q[3];
rz(3.104044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.091207592) q[2];
sx q[2];
rz(-2.51077) q[2];
sx q[2];
rz(1.5861661) q[2];
rz(0.88820109) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(0.92938882) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14426194) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(1.7096747) q[0];
rz(0.56888467) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(2.0137537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9040684) q[0];
sx q[0];
rz(-1.5868109) q[0];
sx q[0];
rz(2.9928656) q[0];
x q[1];
rz(0.94294195) q[2];
sx q[2];
rz(-2.0161511) q[2];
sx q[2];
rz(-3.1415591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20967083) q[1];
sx q[1];
rz(-1.5823963) q[1];
sx q[1];
rz(1.329122) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7752152) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(0.60929326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.613712) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(-2.9679427) q[2];
rz(-0.33637834) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(-1.4655112) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-2.5691659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3819645) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(-0.99960021) q[0];
rz(-pi) q[1];
rz(-2.8675251) q[2];
sx q[2];
rz(-2.0147418) q[2];
sx q[2];
rz(-0.62192164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6483557) q[1];
sx q[1];
rz(-1.3322468) q[1];
sx q[1];
rz(-1.7931213) q[1];
rz(-pi) q[2];
rz(2.8045373) q[3];
sx q[3];
rz(-2.4042077) q[3];
sx q[3];
rz(0.22102236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(2.6043716) q[2];
rz(2.0843263) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-2.4479772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-2.6782425) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(3.0436174) q[2];
sx q[2];
rz(-0.91839472) q[2];
sx q[2];
rz(0.67994946) q[2];
rz(-3.0795931) q[3];
sx q[3];
rz(-2.0580895) q[3];
sx q[3];
rz(2.5930391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];