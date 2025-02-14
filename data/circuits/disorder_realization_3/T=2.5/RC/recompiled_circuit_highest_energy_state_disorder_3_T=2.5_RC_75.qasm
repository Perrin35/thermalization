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
rz(0.74547493) q[0];
sx q[0];
rz(-0.59362721) q[0];
sx q[0];
rz(-2.797085) q[0];
rz(-1.966882) q[1];
sx q[1];
rz(-0.36829683) q[1];
sx q[1];
rz(-0.60454291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4473576) q[0];
sx q[0];
rz(-1.5120872) q[0];
sx q[0];
rz(0.82253463) q[0];
rz(-pi) q[1];
rz(1.3922607) q[2];
sx q[2];
rz(-1.5024619) q[2];
sx q[2];
rz(-2.9477037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82935059) q[1];
sx q[1];
rz(-1.5728467) q[1];
sx q[1];
rz(-0.41872989) q[1];
rz(-0.97957356) q[3];
sx q[3];
rz(-1.0821078) q[3];
sx q[3];
rz(0.6062909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1623666) q[2];
sx q[2];
rz(-1.484885) q[2];
sx q[2];
rz(1.7171198) q[2];
rz(3.034397) q[3];
sx q[3];
rz(-0.19459477) q[3];
sx q[3];
rz(-0.9596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87218881) q[0];
sx q[0];
rz(-0.93282455) q[0];
sx q[0];
rz(-1.8316356) q[0];
rz(2.682389) q[1];
sx q[1];
rz(-1.867086) q[1];
sx q[1];
rz(-2.8244663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2533644) q[0];
sx q[0];
rz(-1.7087052) q[0];
sx q[0];
rz(-1.6186369) q[0];
rz(-pi) q[1];
rz(-3.0917815) q[2];
sx q[2];
rz(-1.5731817) q[2];
sx q[2];
rz(2.002169) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6739247) q[1];
sx q[1];
rz(-1.6367803) q[1];
sx q[1];
rz(-0.54775796) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7789654) q[3];
sx q[3];
rz(-0.82566264) q[3];
sx q[3];
rz(2.7221808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7121048) q[2];
sx q[2];
rz(-1.8450582) q[2];
sx q[2];
rz(2.3340161) q[2];
rz(-0.56337041) q[3];
sx q[3];
rz(-2.6475776) q[3];
sx q[3];
rz(1.0529244) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6323557) q[0];
sx q[0];
rz(-2.95166) q[0];
sx q[0];
rz(2.307039) q[0];
rz(-2.6368311) q[1];
sx q[1];
rz(-0.36634645) q[1];
sx q[1];
rz(1.4588446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2528238) q[0];
sx q[0];
rz(-1.7703759) q[0];
sx q[0];
rz(-2.5779524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4065532) q[2];
sx q[2];
rz(-2.1872518) q[2];
sx q[2];
rz(0.40555813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1881989) q[1];
sx q[1];
rz(-0.85046221) q[1];
sx q[1];
rz(-1.6745643) q[1];
rz(-0.80648544) q[3];
sx q[3];
rz(-1.1317265) q[3];
sx q[3];
rz(-1.0478549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0316281) q[2];
sx q[2];
rz(-1.4723023) q[2];
sx q[2];
rz(-0.36299452) q[2];
rz(1.1686769) q[3];
sx q[3];
rz(-2.3790338) q[3];
sx q[3];
rz(2.3762083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7211683) q[0];
sx q[0];
rz(-2.6154501) q[0];
sx q[0];
rz(-2.5508733) q[0];
rz(0.72714725) q[1];
sx q[1];
rz(-2.7666481) q[1];
sx q[1];
rz(-0.66853833) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6369978) q[0];
sx q[0];
rz(-1.4894006) q[0];
sx q[0];
rz(0.041725833) q[0];
x q[1];
rz(2.0425778) q[2];
sx q[2];
rz(-2.2792247) q[2];
sx q[2];
rz(2.0076795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3154853) q[1];
sx q[1];
rz(-0.8902227) q[1];
sx q[1];
rz(-2.1904551) q[1];
x q[2];
rz(0.56748029) q[3];
sx q[3];
rz(-2.3154308) q[3];
sx q[3];
rz(1.7989649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10355243) q[2];
sx q[2];
rz(-1.870564) q[2];
sx q[2];
rz(-1.5070149) q[2];
rz(-2.7196837) q[3];
sx q[3];
rz(-2.3814337) q[3];
sx q[3];
rz(2.4222477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.6307395) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(-1.0605633) q[0];
rz(-3.0781436) q[1];
sx q[1];
rz(-2.5109406) q[1];
sx q[1];
rz(-2.5877156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4872525) q[0];
sx q[0];
rz(-0.76566389) q[0];
sx q[0];
rz(0.91057093) q[0];
rz(2.3339726) q[2];
sx q[2];
rz(-0.80349001) q[2];
sx q[2];
rz(2.9832911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0530776) q[1];
sx q[1];
rz(-1.6106399) q[1];
sx q[1];
rz(-0.30327176) q[1];
rz(-3.0054566) q[3];
sx q[3];
rz(-1.8936186) q[3];
sx q[3];
rz(-0.28188595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.267103) q[2];
sx q[2];
rz(-0.59659448) q[2];
sx q[2];
rz(-0.60274094) q[2];
rz(-0.069325773) q[3];
sx q[3];
rz(-1.5963138) q[3];
sx q[3];
rz(2.9749405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.199274) q[0];
sx q[0];
rz(-2.5397904) q[0];
sx q[0];
rz(-2.6875575) q[0];
rz(1.8858689) q[1];
sx q[1];
rz(-0.95971003) q[1];
sx q[1];
rz(-1.4383291) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3011264) q[0];
sx q[0];
rz(-2.2181803) q[0];
sx q[0];
rz(0.50917888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82346114) q[2];
sx q[2];
rz(-2.744906) q[2];
sx q[2];
rz(1.6993864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2809873) q[1];
sx q[1];
rz(-2.1327204) q[1];
sx q[1];
rz(-0.58536305) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0694396) q[3];
sx q[3];
rz(-0.92235288) q[3];
sx q[3];
rz(1.4251777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2890275) q[2];
sx q[2];
rz(-1.270741) q[2];
sx q[2];
rz(-1.3555869) q[2];
rz(1.2191314) q[3];
sx q[3];
rz(-1.6617323) q[3];
sx q[3];
rz(-1.5752972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8236302) q[0];
sx q[0];
rz(-0.83331236) q[0];
sx q[0];
rz(-1.3939567) q[0];
rz(-2.6849003) q[1];
sx q[1];
rz(-2.2629181) q[1];
sx q[1];
rz(1.7535694) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0972892) q[0];
sx q[0];
rz(-2.7345022) q[0];
sx q[0];
rz(1.9096229) q[0];
rz(-pi) q[1];
rz(2.6235103) q[2];
sx q[2];
rz(-2.9112795) q[2];
sx q[2];
rz(1.3467333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4108666) q[1];
sx q[1];
rz(-1.7015839) q[1];
sx q[1];
rz(-2.8128002) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61856602) q[3];
sx q[3];
rz(-1.1453218) q[3];
sx q[3];
rz(-2.8989603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8689416) q[2];
sx q[2];
rz(-0.77241263) q[2];
sx q[2];
rz(2.9504377) q[2];
rz(1.8027421) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(-0.52201456) q[3];
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
rz(-pi) q[0];
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
rz(-0.987748) q[0];
sx q[0];
rz(-2.4697883) q[0];
sx q[0];
rz(-1.9195358) q[0];
rz(-0.98474312) q[1];
sx q[1];
rz(-2.1538815) q[1];
sx q[1];
rz(-0.15880671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.384085) q[0];
sx q[0];
rz(-2.473978) q[0];
sx q[0];
rz(3.0098841) q[0];
rz(-pi) q[1];
rz(0.10649781) q[2];
sx q[2];
rz(-0.42909189) q[2];
sx q[2];
rz(-0.96218357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.040222283) q[1];
sx q[1];
rz(-1.1565546) q[1];
sx q[1];
rz(2.657402) q[1];
x q[2];
rz(1.430159) q[3];
sx q[3];
rz(-1.780073) q[3];
sx q[3];
rz(-1.9491553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.020393546) q[2];
sx q[2];
rz(-2.8195916) q[2];
sx q[2];
rz(-0.89312345) q[2];
rz(0.3802866) q[3];
sx q[3];
rz(-1.2023353) q[3];
sx q[3];
rz(-1.7072385) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0685843) q[0];
sx q[0];
rz(-0.46982729) q[0];
sx q[0];
rz(-1.481886) q[0];
rz(2.1549554) q[1];
sx q[1];
rz(-1.4886798) q[1];
sx q[1];
rz(2.5843487) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.754234) q[0];
sx q[0];
rz(-3.0172303) q[0];
sx q[0];
rz(2.6919305) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2115914) q[2];
sx q[2];
rz(-1.2877727) q[2];
sx q[2];
rz(-2.46794) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4367974) q[1];
sx q[1];
rz(-2.3782502) q[1];
sx q[1];
rz(0.82398681) q[1];
rz(2.557665) q[3];
sx q[3];
rz(-3.1213396) q[3];
sx q[3];
rz(0.7113061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3096932) q[2];
sx q[2];
rz(-0.57277402) q[2];
sx q[2];
rz(2.0834303) q[2];
rz(1.38331) q[3];
sx q[3];
rz(-2.1027095) q[3];
sx q[3];
rz(-0.96341187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6365373) q[0];
sx q[0];
rz(-1.1543244) q[0];
sx q[0];
rz(0.43168798) q[0];
rz(-2.441326) q[1];
sx q[1];
rz(-1.9077178) q[1];
sx q[1];
rz(0.69469992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0330677) q[0];
sx q[0];
rz(-2.0742848) q[0];
sx q[0];
rz(-0.42078542) q[0];
rz(-2.0236597) q[2];
sx q[2];
rz(-1.6647571) q[2];
sx q[2];
rz(-2.8270023) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3251288) q[1];
sx q[1];
rz(-2.0184787) q[1];
sx q[1];
rz(2.1486077) q[1];
rz(-pi) q[2];
rz(1.6142817) q[3];
sx q[3];
rz(-1.3924034) q[3];
sx q[3];
rz(0.94387335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.485864) q[2];
sx q[2];
rz(-2.8317917) q[2];
sx q[2];
rz(-2.7244205) q[2];
rz(-2.396865) q[3];
sx q[3];
rz(-1.797902) q[3];
sx q[3];
rz(-2.1630796) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7930631) q[0];
sx q[0];
rz(-1.2611669) q[0];
sx q[0];
rz(1.4766759) q[0];
rz(1.9229802) q[1];
sx q[1];
rz(-2.0017793) q[1];
sx q[1];
rz(0.58560169) q[1];
rz(-2.2485366) q[2];
sx q[2];
rz(-1.1184659) q[2];
sx q[2];
rz(2.804166) q[2];
rz(0.32998362) q[3];
sx q[3];
rz(-0.97616227) q[3];
sx q[3];
rz(-0.57310692) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
