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
rz(1.7293575) q[0];
sx q[0];
rz(-2.5492302) q[0];
sx q[0];
rz(1.2047729) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(-1.187721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76554322) q[0];
sx q[0];
rz(-0.57224579) q[0];
sx q[0];
rz(2.8577198) q[0];
rz(2.4815791) q[2];
sx q[2];
rz(-0.68049351) q[2];
sx q[2];
rz(1.6689491) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.097259911) q[1];
sx q[1];
rz(-2.9625434) q[1];
sx q[1];
rz(1.8392842) q[1];
rz(-1.9378565) q[3];
sx q[3];
rz(-1.4407183) q[3];
sx q[3];
rz(-0.50103044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5997233) q[2];
sx q[2];
rz(-2.5399127) q[2];
sx q[2];
rz(2.0987233) q[2];
rz(0.045182191) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(-1.5359623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0829818) q[0];
sx q[0];
rz(-0.11317145) q[0];
sx q[0];
rz(0.97852069) q[0];
rz(-1.9784031) q[1];
sx q[1];
rz(-2.1470943) q[1];
sx q[1];
rz(1.3226343) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90962553) q[0];
sx q[0];
rz(-2.2920906) q[0];
sx q[0];
rz(-0.90497525) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2873994) q[2];
sx q[2];
rz(-1.3958418) q[2];
sx q[2];
rz(2.8307479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6860477) q[1];
sx q[1];
rz(-1.7864704) q[1];
sx q[1];
rz(0.91051813) q[1];
rz(-pi) q[2];
x q[2];
rz(1.421881) q[3];
sx q[3];
rz(-2.2595539) q[3];
sx q[3];
rz(1.0721579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3375552) q[2];
sx q[2];
rz(-2.940371) q[2];
sx q[2];
rz(-3.1031754) q[2];
rz(1.5086959) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(-1.647515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9134193) q[0];
sx q[0];
rz(-0.67263022) q[0];
sx q[0];
rz(2.9056554) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.6940247) q[1];
sx q[1];
rz(2.6441914) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1318404) q[0];
sx q[0];
rz(-1.3924283) q[0];
sx q[0];
rz(-1.7531036) q[0];
x q[1];
rz(-0.3942986) q[2];
sx q[2];
rz(-1.9579534) q[2];
sx q[2];
rz(0.55768572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5617852) q[1];
sx q[1];
rz(-2.7473463) q[1];
sx q[1];
rz(0.057478776) q[1];
rz(-pi) q[2];
rz(2.892832) q[3];
sx q[3];
rz(-1.8210871) q[3];
sx q[3];
rz(2.8119025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95911372) q[2];
sx q[2];
rz(-1.4664058) q[2];
sx q[2];
rz(-1.9100995) q[2];
rz(1.69151) q[3];
sx q[3];
rz(-1.3030038) q[3];
sx q[3];
rz(0.85795295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0813783) q[0];
sx q[0];
rz(-0.84351081) q[0];
sx q[0];
rz(-0.1928992) q[0];
rz(-1.7861722) q[1];
sx q[1];
rz(-1.5602292) q[1];
sx q[1];
rz(0.17098175) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6714013) q[0];
sx q[0];
rz(-2.6220052) q[0];
sx q[0];
rz(1.4524824) q[0];
rz(-pi) q[1];
rz(0.58329432) q[2];
sx q[2];
rz(-2.3153044) q[2];
sx q[2];
rz(-0.14329958) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4857085) q[1];
sx q[1];
rz(-8*pi/15) q[1];
sx q[1];
rz(2.0018565) q[1];
x q[2];
rz(-0.90733068) q[3];
sx q[3];
rz(-1.7620834) q[3];
sx q[3];
rz(0.80463791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96044797) q[2];
sx q[2];
rz(-1.170271) q[2];
sx q[2];
rz(-2.81874) q[2];
rz(1.7668096) q[3];
sx q[3];
rz(-1.0389453) q[3];
sx q[3];
rz(-0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5906931) q[0];
sx q[0];
rz(-1.0223848) q[0];
sx q[0];
rz(-2.2160231) q[0];
rz(3.0135221) q[1];
sx q[1];
rz(-1.6431199) q[1];
sx q[1];
rz(0.099460348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8444237) q[0];
sx q[0];
rz(-2.4250829) q[0];
sx q[0];
rz(-1.5708718) q[0];
x q[1];
rz(2.3771466) q[2];
sx q[2];
rz(-2.6607155) q[2];
sx q[2];
rz(0.99221992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22597081) q[1];
sx q[1];
rz(-0.56949896) q[1];
sx q[1];
rz(1.9379338) q[1];
rz(1.8811536) q[3];
sx q[3];
rz(-0.25729968) q[3];
sx q[3];
rz(2.8043487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9596935) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(2.7679475) q[2];
rz(-2.2419808) q[3];
sx q[3];
rz(-2.3989232) q[3];
sx q[3];
rz(-2.0067748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23168343) q[0];
sx q[0];
rz(-0.21655701) q[0];
sx q[0];
rz(2.5496971) q[0];
rz(-0.20467155) q[1];
sx q[1];
rz(-0.99212956) q[1];
sx q[1];
rz(3.0904904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0979733) q[0];
sx q[0];
rz(-1.7228925) q[0];
sx q[0];
rz(0.19466227) q[0];
rz(-pi) q[1];
rz(-2.2289071) q[2];
sx q[2];
rz(-1.3649696) q[2];
sx q[2];
rz(-1.4806946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7390427) q[1];
sx q[1];
rz(-1.2616757) q[1];
sx q[1];
rz(1.2039295) q[1];
rz(-0.53364086) q[3];
sx q[3];
rz(-2.009005) q[3];
sx q[3];
rz(2.4166188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0245612) q[2];
sx q[2];
rz(-2.0270963) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(-2.7545605) q[3];
sx q[3];
rz(-2.6670167) q[3];
sx q[3];
rz(0.3064557) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759906) q[0];
sx q[0];
rz(-0.87823534) q[0];
sx q[0];
rz(0.24755092) q[0];
rz(1.4546825) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(-0.036458485) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582583) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(-2.3733632) q[0];
rz(-pi) q[1];
rz(-0.22561947) q[2];
sx q[2];
rz(-0.7157514) q[2];
sx q[2];
rz(0.72474397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4562021) q[1];
sx q[1];
rz(-1.0674607) q[1];
sx q[1];
rz(2.6052942) q[1];
rz(2.2685675) q[3];
sx q[3];
rz(-3.0254078) q[3];
sx q[3];
rz(0.98770638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.054691943) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(-1.7894233) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(-1.1917535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2317113) q[0];
sx q[0];
rz(-0.24869643) q[0];
sx q[0];
rz(-1.6491718) q[0];
rz(0.17852783) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(-2.5915204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9880542) q[0];
sx q[0];
rz(-0.037153553) q[0];
sx q[0];
rz(-2.1376993) q[0];
rz(1.8166601) q[2];
sx q[2];
rz(-2.8865007) q[2];
sx q[2];
rz(-2.6691797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8892563) q[1];
sx q[1];
rz(-1.0358216) q[1];
sx q[1];
rz(1.757213) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1060195) q[3];
sx q[3];
rz(-1.2039935) q[3];
sx q[3];
rz(0.22983009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7940346) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(-0.75322914) q[2];
rz(2.8473162) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(3.0729955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8681965) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(-1.6896601) q[0];
rz(-0.1952576) q[1];
sx q[1];
rz(-1.0804907) q[1];
sx q[1];
rz(1.6498227) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30376745) q[0];
sx q[0];
rz(-1.2653102) q[0];
sx q[0];
rz(-0.49834337) q[0];
rz(0.765018) q[2];
sx q[2];
rz(-1.49053) q[2];
sx q[2];
rz(-1.2470055) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10354708) q[1];
sx q[1];
rz(-2.9730151) q[1];
sx q[1];
rz(1.176531) q[1];
rz(-pi) q[2];
rz(1.5641065) q[3];
sx q[3];
rz(-1.5433528) q[3];
sx q[3];
rz(-0.60718482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1008272) q[2];
sx q[2];
rz(-2.6746174) q[2];
sx q[2];
rz(2.2838498) q[2];
rz(-1.6929251) q[3];
sx q[3];
rz(-1.6489776) q[3];
sx q[3];
rz(-2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(2.1345074) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(2.3585228) q[0];
rz(1.1231517) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(-0.060308594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6864844) q[0];
sx q[0];
rz(-1.2164854) q[0];
sx q[0];
rz(-0.099781009) q[0];
rz(-pi) q[1];
rz(-3.0357828) q[2];
sx q[2];
rz(-0.98710892) q[2];
sx q[2];
rz(-2.2367144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26877182) q[1];
sx q[1];
rz(-0.67292456) q[1];
sx q[1];
rz(-0.10867837) q[1];
x q[2];
rz(-0.29197146) q[3];
sx q[3];
rz(-2.7640016) q[3];
sx q[3];
rz(-0.99536125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15410885) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(1.3101428) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(-1.8346627) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46351984) q[0];
sx q[0];
rz(-2.0250043) q[0];
sx q[0];
rz(2.9153839) q[0];
rz(-0.75640596) q[1];
sx q[1];
rz(-0.49497985) q[1];
sx q[1];
rz(2.3351647) q[1];
rz(-3.1210774) q[2];
sx q[2];
rz(-2.7111369) q[2];
sx q[2];
rz(-2.6028056) q[2];
rz(-1.9202833) q[3];
sx q[3];
rz(-1.7428453) q[3];
sx q[3];
rz(1.3475628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
