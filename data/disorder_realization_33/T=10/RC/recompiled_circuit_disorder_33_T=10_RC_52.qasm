OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(-0.27591053) q[0];
sx q[0];
rz(1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64496541) q[0];
sx q[0];
rz(-0.49155203) q[0];
sx q[0];
rz(-2.3636723) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4351032) q[2];
sx q[2];
rz(-0.90105614) q[2];
sx q[2];
rz(-1.1342088) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95333245) q[1];
sx q[1];
rz(-0.42020513) q[1];
sx q[1];
rz(-2.5145867) q[1];
x q[2];
rz(2.0516112) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-0.68457505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-0.21683189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502055) q[0];
sx q[0];
rz(-0.7190401) q[0];
sx q[0];
rz(-1.1262116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47302834) q[2];
sx q[2];
rz(-2.0748667) q[2];
sx q[2];
rz(0.31993983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9456957) q[1];
sx q[1];
rz(-2.2356114) q[1];
sx q[1];
rz(-0.2701387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2651477) q[3];
sx q[3];
rz(-1.0507686) q[3];
sx q[3];
rz(-2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(2.537354) q[0];
rz(1.8151981) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(2.2089829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843825) q[0];
sx q[0];
rz(-2.819448) q[0];
sx q[0];
rz(1.4003217) q[0];
rz(-0.18205299) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(-1.6857266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31308094) q[1];
sx q[1];
rz(-2.5701437) q[1];
sx q[1];
rz(0.22027318) q[1];
rz(-pi) q[2];
rz(0.44585769) q[3];
sx q[3];
rz(-1.0599531) q[3];
sx q[3];
rz(-1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(-2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(-0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.4979699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1456137) q[0];
sx q[0];
rz(-1.0756452) q[0];
sx q[0];
rz(-0.33546319) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.056604071) q[2];
sx q[2];
rz(-1.7610234) q[2];
sx q[2];
rz(2.9375926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30724635) q[1];
sx q[1];
rz(-1.8028959) q[1];
sx q[1];
rz(2.8744065) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46521503) q[3];
sx q[3];
rz(-1.1606996) q[3];
sx q[3];
rz(-2.0801534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(-2.4397819) q[2];
rz(-2.3102405) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
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
rz(2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(1.6197846) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(1.048208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883023) q[0];
sx q[0];
rz(-0.36670812) q[0];
sx q[0];
rz(2.328863) q[0];
x q[1];
rz(0.016096073) q[2];
sx q[2];
rz(-0.65158366) q[2];
sx q[2];
rz(0.96166699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8909) q[1];
sx q[1];
rz(-0.39784583) q[1];
sx q[1];
rz(2.489151) q[1];
rz(-pi) q[2];
rz(2.1861595) q[3];
sx q[3];
rz(-1.2293929) q[3];
sx q[3];
rz(-1.1036901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(0.12983233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(-0.58563389) q[0];
rz(2.1461357) q[2];
sx q[2];
rz(-1.2577004) q[2];
sx q[2];
rz(2.7196333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2643913) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(2.0045723) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54883212) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(0.43743922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(1.6947421) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-2.4553305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670358) q[0];
sx q[0];
rz(-2.0445163) q[0];
sx q[0];
rz(-2.0635598) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5574066) q[2];
sx q[2];
rz(-0.91911941) q[2];
sx q[2];
rz(2.3551031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2346748) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(-0.075637416) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63728441) q[3];
sx q[3];
rz(-2.6316959) q[3];
sx q[3];
rz(-2.2989458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.6954533) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.6400281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088889) q[0];
sx q[0];
rz(-0.49312691) q[0];
sx q[0];
rz(-1.5296442) q[0];
rz(-0.026560606) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(2.314687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0280684) q[1];
sx q[1];
rz(-1.6599732) q[1];
sx q[1];
rz(2.813617) q[1];
x q[2];
rz(0.73236671) q[3];
sx q[3];
rz(-2.3673956) q[3];
sx q[3];
rz(-0.68294169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33655745) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68181224) q[0];
sx q[0];
rz(-2.5373055) q[0];
sx q[0];
rz(-2.2629645) q[0];
x q[1];
rz(-2.4497689) q[2];
sx q[2];
rz(-1.3335388) q[2];
sx q[2];
rz(-2.0123864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6955399) q[1];
sx q[1];
rz(-0.76247588) q[1];
sx q[1];
rz(-1.4941077) q[1];
rz(-pi) q[2];
rz(0.10961253) q[3];
sx q[3];
rz(-1.622756) q[3];
sx q[3];
rz(-2.6250641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-1.055868) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4984109) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-2.9472651) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(2.1077572) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4612504) q[0];
sx q[0];
rz(-2.470812) q[0];
sx q[0];
rz(-1.781342) q[0];
rz(-pi) q[1];
rz(-0.50844426) q[2];
sx q[2];
rz(-2.2061081) q[2];
sx q[2];
rz(2.9490162) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6453213) q[1];
sx q[1];
rz(-1.4693345) q[1];
sx q[1];
rz(0.99781499) q[1];
rz(-pi) q[2];
rz(1.3528353) q[3];
sx q[3];
rz(-0.86943227) q[3];
sx q[3];
rz(2.4234114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(0.6357843) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-1.5244665) q[2];
sx q[2];
rz(-2.5382858) q[2];
sx q[2];
rz(-0.45679191) q[2];
rz(1.2407606) q[3];
sx q[3];
rz(-1.6374554) q[3];
sx q[3];
rz(-1.0367254) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];