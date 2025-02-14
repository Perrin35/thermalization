OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.59392053) q[0];
sx q[0];
rz(-2.5424818) q[0];
sx q[0];
rz(0.17679086) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(3.6511753) q[1];
sx q[1];
rz(9.3378172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6330948) q[0];
sx q[0];
rz(-2.1380391) q[0];
sx q[0];
rz(-1.1192516) q[0];
rz(-pi) q[1];
rz(-0.29524191) q[2];
sx q[2];
rz(-2.549941) q[2];
sx q[2];
rz(1.0585143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0217758) q[1];
sx q[1];
rz(-1.4939185) q[1];
sx q[1];
rz(-1.4426115) q[1];
x q[2];
rz(-2.5218202) q[3];
sx q[3];
rz(-1.3300782) q[3];
sx q[3];
rz(0.71686577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1842492) q[2];
sx q[2];
rz(-0.49746305) q[2];
sx q[2];
rz(0.59172612) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-0.44132909) q[3];
sx q[3];
rz(-2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1205207) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(2.8968107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0140959) q[0];
sx q[0];
rz(-1.2359706) q[0];
sx q[0];
rz(-2.0398519) q[0];
x q[1];
rz(-0.4366283) q[2];
sx q[2];
rz(-2.2712913) q[2];
sx q[2];
rz(-1.3134522) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87585014) q[1];
sx q[1];
rz(-0.85803723) q[1];
sx q[1];
rz(2.4044988) q[1];
rz(-2.6908633) q[3];
sx q[3];
rz(-2.8141032) q[3];
sx q[3];
rz(-2.9779676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.327534) q[2];
sx q[2];
rz(-0.56698292) q[2];
sx q[2];
rz(-2.2954693) q[2];
rz(-0.36492473) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(2.164771) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1035128) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(2.9981151) q[0];
rz(-1.4682651) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-3.0751244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9212355) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(0.14089091) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8861553) q[2];
sx q[2];
rz(-2.0991) q[2];
sx q[2];
rz(-0.87919368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.089515226) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(-3.0245824) q[1];
rz(1.7371561) q[3];
sx q[3];
rz(-0.8373973) q[3];
sx q[3];
rz(-2.211253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0299783) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(0.31271333) q[2];
rz(0.2615658) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(-0.79791445) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0161491) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(-0.087015986) q[0];
rz(1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(0.048390128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66358405) q[0];
sx q[0];
rz(-0.80246325) q[0];
sx q[0];
rz(1.7491231) q[0];
x q[1];
rz(-2.9913285) q[2];
sx q[2];
rz(-1.74904) q[2];
sx q[2];
rz(-1.5706289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9632513) q[1];
sx q[1];
rz(-1.5164485) q[1];
sx q[1];
rz(2.4456294) q[1];
x q[2];
rz(-2.257423) q[3];
sx q[3];
rz(-2.3034952) q[3];
sx q[3];
rz(0.17894408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4309569) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(-0.44130138) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.3683616) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5271673) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(1.6983039) q[1];
sx q[1];
rz(-1.8836421) q[1];
sx q[1];
rz(1.7194933) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2724394) q[0];
sx q[0];
rz(-2.5792988) q[0];
sx q[0];
rz(0.2642239) q[0];
rz(-1.7781813) q[2];
sx q[2];
rz(-1.3355012) q[2];
sx q[2];
rz(-2.4547142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4053012) q[1];
sx q[1];
rz(-1.6561688) q[1];
sx q[1];
rz(3.0194204) q[1];
rz(-pi) q[2];
rz(2.5636219) q[3];
sx q[3];
rz(-0.9909361) q[3];
sx q[3];
rz(2.487929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(-1.5606073) q[2];
rz(-1.3382781) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(-2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018589858) q[0];
sx q[0];
rz(-0.81273166) q[0];
sx q[0];
rz(0.71989584) q[0];
rz(0.24049354) q[1];
sx q[1];
rz(-1.980282) q[1];
sx q[1];
rz(2.8954411) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3703281) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(1.5465897) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42527898) q[2];
sx q[2];
rz(-0.80406791) q[2];
sx q[2];
rz(-0.67415392) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64061058) q[1];
sx q[1];
rz(-2.3965008) q[1];
sx q[1];
rz(0.19265811) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6443374) q[3];
sx q[3];
rz(-2.3905633) q[3];
sx q[3];
rz(2.8145144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64421946) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(-2.6175446) q[3];
sx q[3];
rz(-2.7553813) q[3];
sx q[3];
rz(3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2972357) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(-2.2204087) q[0];
rz(1.4211897) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-0.99501077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667082) q[0];
sx q[0];
rz(-2.6972983) q[0];
sx q[0];
rz(1.0309451) q[0];
x q[1];
rz(-0.065804577) q[2];
sx q[2];
rz(-1.1514246) q[2];
sx q[2];
rz(-2.9247627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2176358) q[1];
sx q[1];
rz(-1.1502277) q[1];
sx q[1];
rz(-2.5687508) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0070442) q[3];
sx q[3];
rz(-1.1863669) q[3];
sx q[3];
rz(0.72641295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2034188) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(0.43756884) q[2];
rz(-2.3214052) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1831128) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(-1.0108277) q[0];
rz(-3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(0.54862499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47893229) q[0];
sx q[0];
rz(-2.9629483) q[0];
sx q[0];
rz(-1.7500072) q[0];
rz(-pi) q[1];
rz(-1.8993511) q[2];
sx q[2];
rz(-1.1588084) q[2];
sx q[2];
rz(3.0871473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1774166) q[1];
sx q[1];
rz(-0.89277041) q[1];
sx q[1];
rz(-1.3864338) q[1];
rz(-pi) q[2];
rz(0.88415481) q[3];
sx q[3];
rz(-0.13152371) q[3];
sx q[3];
rz(0.92708528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34714547) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(0.53317201) q[2];
rz(-2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(2.6326411) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0549523) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(-2.3593498) q[0];
rz(-3.0714463) q[1];
sx q[1];
rz(-0.47824305) q[1];
sx q[1];
rz(2.9152962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5180888) q[0];
sx q[0];
rz(-1.6395578) q[0];
sx q[0];
rz(3.070773) q[0];
rz(-pi) q[1];
rz(-0.89266291) q[2];
sx q[2];
rz(-1.423866) q[2];
sx q[2];
rz(0.28444296) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8574684) q[1];
sx q[1];
rz(-0.54912607) q[1];
sx q[1];
rz(-2.930307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6582011) q[3];
sx q[3];
rz(-2.8897396) q[3];
sx q[3];
rz(0.15628584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(2.8790224) q[3];
sx q[3];
rz(-1.5841443) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(0.18375272) q[0];
rz(2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(-2.5180838) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79270169) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(-1.2280653) q[0];
rz(-0.6208159) q[2];
sx q[2];
rz(-2.0784272) q[2];
sx q[2];
rz(-2.1073282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5851046) q[1];
sx q[1];
rz(-0.10837238) q[1];
sx q[1];
rz(-0.21173112) q[1];
x q[2];
rz(1.9454089) q[3];
sx q[3];
rz(-1.5996029) q[3];
sx q[3];
rz(-2.0458986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3925675) q[2];
sx q[2];
rz(-0.6124658) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(-0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(-2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2035718) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(0.45166311) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(2.6653566) q[2];
sx q[2];
rz(-0.58544896) q[2];
sx q[2];
rz(0.61665012) q[2];
rz(-0.59521159) q[3];
sx q[3];
rz(-2.8122936) q[3];
sx q[3];
rz(-0.82174792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
