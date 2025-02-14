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
rz(0.053057916) q[0];
sx q[0];
rz(0.36202708) q[0];
sx q[0];
rz(10.590635) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(1.4282164) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27487684) q[0];
sx q[0];
rz(-2.6259415) q[0];
sx q[0];
rz(-3.0192119) q[0];
x q[1];
rz(0.25144724) q[2];
sx q[2];
rz(-1.5000952) q[2];
sx q[2];
rz(-1.9858907) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62318351) q[1];
sx q[1];
rz(-1.3819547) q[1];
sx q[1];
rz(-3.0075226) q[1];
x q[2];
rz(0.26784874) q[3];
sx q[3];
rz(-2.043927) q[3];
sx q[3];
rz(2.3426263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6894138) q[2];
sx q[2];
rz(-0.016760085) q[2];
sx q[2];
rz(-2.9407799) q[2];
rz(-2.9928442) q[3];
sx q[3];
rz(-0.0047618682) q[3];
sx q[3];
rz(0.31140056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1397322) q[0];
sx q[0];
rz(-0.59356028) q[0];
sx q[0];
rz(1.0396022) q[0];
rz(-0.014558583) q[1];
sx q[1];
rz(-1.9083551) q[1];
sx q[1];
rz(1.5537517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7015052) q[0];
sx q[0];
rz(-1.0453512) q[0];
sx q[0];
rz(-0.82869014) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.014694728) q[2];
sx q[2];
rz(-1.4966244) q[2];
sx q[2];
rz(-1.6770404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55744104) q[1];
sx q[1];
rz(-1.5366035) q[1];
sx q[1];
rz(-1.8299915) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.513345) q[3];
sx q[3];
rz(-0.37783315) q[3];
sx q[3];
rz(2.6031818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62850922) q[2];
sx q[2];
rz(-1.5336978) q[2];
sx q[2];
rz(1.3879363) q[2];
rz(-1.3758818) q[3];
sx q[3];
rz(-1.0466156) q[3];
sx q[3];
rz(0.29471135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917711) q[0];
sx q[0];
rz(-0.23673683) q[0];
sx q[0];
rz(-2.5340875) q[0];
rz(1.5982184) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(2.1764596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34244363) q[0];
sx q[0];
rz(-1.5502872) q[0];
sx q[0];
rz(0.0040702013) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1943614) q[2];
sx q[2];
rz(-0.5493702) q[2];
sx q[2];
rz(-1.4538764) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9671792) q[1];
sx q[1];
rz(-1.426218) q[1];
sx q[1];
rz(0.036199613) q[1];
rz(-2.2854961) q[3];
sx q[3];
rz(-3.0149547) q[3];
sx q[3];
rz(0.52982012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8049916) q[2];
sx q[2];
rz(-0.6707297) q[2];
sx q[2];
rz(-2.2750308) q[2];
rz(2.0349515) q[3];
sx q[3];
rz(-1.589078) q[3];
sx q[3];
rz(-1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.57257819) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(1.6160075) q[0];
rz(-0.010604803) q[1];
sx q[1];
rz(-0.0037825982) q[1];
sx q[1];
rz(2.3912281) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.01975) q[0];
sx q[0];
rz(-1.4759617) q[0];
sx q[0];
rz(0.66480831) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60549824) q[2];
sx q[2];
rz(-1.7586305) q[2];
sx q[2];
rz(-0.8643291) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8349223) q[1];
sx q[1];
rz(-2.6205898) q[1];
sx q[1];
rz(2.9087594) q[1];
x q[2];
rz(0.37639002) q[3];
sx q[3];
rz(-0.3335267) q[3];
sx q[3];
rz(-3.0390449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41287199) q[2];
sx q[2];
rz(-2.0527288) q[2];
sx q[2];
rz(1.2415761) q[2];
rz(3.1359172) q[3];
sx q[3];
rz(-0.81040502) q[3];
sx q[3];
rz(-2.2976105) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64549696) q[0];
sx q[0];
rz(-3.0912919) q[0];
sx q[0];
rz(-0.92329931) q[0];
rz(-0.80054379) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(2.961535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6430017) q[0];
sx q[0];
rz(-2.8992436) q[0];
sx q[0];
rz(1.7447628) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4550986) q[2];
sx q[2];
rz(-0.1265993) q[2];
sx q[2];
rz(1.1243658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1504768) q[1];
sx q[1];
rz(-1.9228982) q[1];
sx q[1];
rz(-1.8262509) q[1];
rz(-pi) q[2];
rz(1.8534142) q[3];
sx q[3];
rz(-2.7621531) q[3];
sx q[3];
rz(-1.5925549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30183145) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(-1.7054455) q[2];
rz(1.885421) q[3];
sx q[3];
rz(-1.6585645) q[3];
sx q[3];
rz(-3.0459611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489478) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(-2.6982464) q[0];
rz(1.9709142) q[1];
sx q[1];
rz(-3.1405293) q[1];
sx q[1];
rz(-0.43177691) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3928089) q[0];
sx q[0];
rz(-1.8816299) q[0];
sx q[0];
rz(-2.3883467) q[0];
rz(-pi) q[1];
rz(-2.0197996) q[2];
sx q[2];
rz(-0.19059715) q[2];
sx q[2];
rz(-2.7191741) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1237632) q[1];
sx q[1];
rz(-2.4430954) q[1];
sx q[1];
rz(-0.61459728) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8006936) q[3];
sx q[3];
rz(-1.8935003) q[3];
sx q[3];
rz(1.6663807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7410437) q[2];
sx q[2];
rz(-2.2811175) q[2];
sx q[2];
rz(1.7575556) q[2];
rz(-2.4436229) q[3];
sx q[3];
rz(-0.75708404) q[3];
sx q[3];
rz(-0.32148263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8292002) q[0];
sx q[0];
rz(-1.3425403) q[0];
sx q[0];
rz(0.37305748) q[0];
rz(0.27698764) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(-0.75609797) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6269659) q[0];
sx q[0];
rz(-1.1921765) q[0];
sx q[0];
rz(0.65831229) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41474386) q[2];
sx q[2];
rz(-1.2325792) q[2];
sx q[2];
rz(1.2646706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7623118) q[1];
sx q[1];
rz(-1.9729776) q[1];
sx q[1];
rz(2.1625175) q[1];
x q[2];
rz(-0.30310615) q[3];
sx q[3];
rz(-1.167093) q[3];
sx q[3];
rz(1.84246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9296391) q[2];
sx q[2];
rz(-0.55277199) q[2];
sx q[2];
rz(-2.2201404) q[2];
rz(3.0742505) q[3];
sx q[3];
rz(-1.9332956) q[3];
sx q[3];
rz(1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15143722) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(0.13023278) q[0];
rz(0.79365927) q[1];
sx q[1];
rz(-3.1399813) q[1];
sx q[1];
rz(-2.8614955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.139629) q[0];
sx q[0];
rz(-3.0462777) q[0];
sx q[0];
rz(-2.167666) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21799223) q[2];
sx q[2];
rz(-1.483886) q[2];
sx q[2];
rz(0.72959585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4520766) q[1];
sx q[1];
rz(-1.0930287) q[1];
sx q[1];
rz(-2.6803451) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0956826) q[3];
sx q[3];
rz(-2.2244448) q[3];
sx q[3];
rz(0.16330367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0555931) q[2];
sx q[2];
rz(-1.6722101) q[2];
sx q[2];
rz(-0.80297339) q[2];
rz(1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(2.3101961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11237535) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(0.10920864) q[0];
rz(-0.373492) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(2.5901897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5729061) q[0];
sx q[0];
rz(-2.150642) q[0];
sx q[0];
rz(0.93007472) q[0];
x q[1];
rz(2.5481651) q[2];
sx q[2];
rz(-0.22106904) q[2];
sx q[2];
rz(-1.9367557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0903328) q[1];
sx q[1];
rz(-0.96104927) q[1];
sx q[1];
rz(-2.4734797) q[1];
x q[2];
rz(-2.6104796) q[3];
sx q[3];
rz(-2.6817842) q[3];
sx q[3];
rz(2.8961033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.52770829) q[2];
sx q[2];
rz(-1.3091427) q[2];
sx q[2];
rz(-1.8180397) q[2];
rz(-1.2644794) q[3];
sx q[3];
rz(-1.8529961) q[3];
sx q[3];
rz(-0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5453813) q[0];
sx q[0];
rz(-0.63544202) q[0];
sx q[0];
rz(2.4052461) q[0];
rz(-2.9296854) q[1];
sx q[1];
rz(-2.1129463) q[1];
sx q[1];
rz(-1.597499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6074267) q[0];
sx q[0];
rz(-1.5887999) q[0];
sx q[0];
rz(-2.8624363) q[0];
x q[1];
rz(-1.5681463) q[2];
sx q[2];
rz(-1.6011597) q[2];
sx q[2];
rz(-1.7372436) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5527344) q[1];
sx q[1];
rz(-0.97179669) q[1];
sx q[1];
rz(-1.0042648) q[1];
x q[2];
rz(1.8964564) q[3];
sx q[3];
rz(-1.146637) q[3];
sx q[3];
rz(-0.87982925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.836901) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(-1.8677853) q[2];
rz(-1.6972208) q[3];
sx q[3];
rz(-3.0675409) q[3];
sx q[3];
rz(1.4950289) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.973751) q[0];
sx q[0];
rz(-1.5580307) q[0];
sx q[0];
rz(1.8488484) q[0];
rz(-1.6043067) q[1];
sx q[1];
rz(-0.91265408) q[1];
sx q[1];
rz(0.18462054) q[1];
rz(-1.5617076) q[2];
sx q[2];
rz(-1.5307003) q[2];
sx q[2];
rz(-1.8673473) q[2];
rz(0.12243263) q[3];
sx q[3];
rz(-2.1748016) q[3];
sx q[3];
rz(0.55007269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
