OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(2.2428089) q[0];
sx q[0];
rz(8.1221683) q[0];
rz(0.12401914) q[1];
sx q[1];
rz(-1.4065341) q[1];
sx q[1];
rz(-1.6279434) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6684912) q[0];
sx q[0];
rz(-1.2364179) q[0];
sx q[0];
rz(-0.3121434) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0997203) q[2];
sx q[2];
rz(-1.0827107) q[2];
sx q[2];
rz(-1.5630592) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4931655) q[1];
sx q[1];
rz(-1.9110288) q[1];
sx q[1];
rz(-3.0837584) q[1];
x q[2];
rz(-0.94515015) q[3];
sx q[3];
rz(-1.2602451) q[3];
sx q[3];
rz(-2.9356082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0836432) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(-2.7570214) q[2];
rz(0.40258506) q[3];
sx q[3];
rz(-1.7076924) q[3];
sx q[3];
rz(-2.3334077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3835417) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(-2.5480399) q[0];
rz(-1.3213762) q[1];
sx q[1];
rz(-2.0621736) q[1];
sx q[1];
rz(-3.1022601) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44839076) q[0];
sx q[0];
rz(-1.1993919) q[0];
sx q[0];
rz(2.0601963) q[0];
rz(-pi) q[1];
rz(-2.5051475) q[2];
sx q[2];
rz(-1.6628254) q[2];
sx q[2];
rz(-2.4534289) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63168282) q[1];
sx q[1];
rz(-1.0541774) q[1];
sx q[1];
rz(1.6277908) q[1];
rz(-pi) q[2];
rz(1.6013299) q[3];
sx q[3];
rz(-2.1200006) q[3];
sx q[3];
rz(-2.1516678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84052691) q[2];
sx q[2];
rz(-1.4378005) q[2];
sx q[2];
rz(-0.62704101) q[2];
rz(-0.4380694) q[3];
sx q[3];
rz(-1.4984683) q[3];
sx q[3];
rz(-1.1367249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42191926) q[0];
sx q[0];
rz(-1.0872343) q[0];
sx q[0];
rz(1.7701953) q[0];
rz(2.5321391) q[1];
sx q[1];
rz(-2.2513159) q[1];
sx q[1];
rz(-1.2249464) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44642513) q[0];
sx q[0];
rz(-1.37943) q[0];
sx q[0];
rz(1.6686977) q[0];
rz(-2.0535499) q[2];
sx q[2];
rz(-2.4039589) q[2];
sx q[2];
rz(-2.9575153) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65388262) q[1];
sx q[1];
rz(-1.2131897) q[1];
sx q[1];
rz(-2.7951711) q[1];
x q[2];
rz(0.87525778) q[3];
sx q[3];
rz(-1.4039543) q[3];
sx q[3];
rz(2.5341906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3043392) q[2];
sx q[2];
rz(-0.86696583) q[2];
sx q[2];
rz(2.2875817) q[2];
rz(-0.52418661) q[3];
sx q[3];
rz(-1.6322522) q[3];
sx q[3];
rz(-0.9700276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22685856) q[0];
sx q[0];
rz(-2.9840042) q[0];
sx q[0];
rz(-0.7937113) q[0];
rz(-0.0018250068) q[1];
sx q[1];
rz(-2.5853214) q[1];
sx q[1];
rz(0.8108286) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6452519) q[0];
sx q[0];
rz(-0.84028572) q[0];
sx q[0];
rz(-2.8964983) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7063058) q[2];
sx q[2];
rz(-2.1974124) q[2];
sx q[2];
rz(-0.97490189) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5280247) q[1];
sx q[1];
rz(-2.151229) q[1];
sx q[1];
rz(-2.7356262) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.879519) q[3];
sx q[3];
rz(-2.6938872) q[3];
sx q[3];
rz(1.6310504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4578555) q[2];
sx q[2];
rz(-0.33317864) q[2];
sx q[2];
rz(1.6443171) q[2];
rz(-1.3582683) q[3];
sx q[3];
rz(-1.0969011) q[3];
sx q[3];
rz(1.3076521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7233647) q[0];
sx q[0];
rz(-1.6809373) q[0];
sx q[0];
rz(2.8826868) q[0];
rz(-3.018191) q[1];
sx q[1];
rz(-2.5105208) q[1];
sx q[1];
rz(-2.6554328) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82993342) q[0];
sx q[0];
rz(-1.0898542) q[0];
sx q[0];
rz(-2.2416058) q[0];
rz(-pi) q[1];
rz(0.45275432) q[2];
sx q[2];
rz(-1.448302) q[2];
sx q[2];
rz(-1.233135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9100719) q[1];
sx q[1];
rz(-1.2086444) q[1];
sx q[1];
rz(-2.9400184) q[1];
x q[2];
rz(0.87032373) q[3];
sx q[3];
rz(-2.016267) q[3];
sx q[3];
rz(2.0014167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2316124) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(1.7463589) q[2];
rz(2.1868271) q[3];
sx q[3];
rz(-1.3937817) q[3];
sx q[3];
rz(-2.098293) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2024277) q[0];
sx q[0];
rz(-1.2347777) q[0];
sx q[0];
rz(-2.3811316) q[0];
rz(2.3072534) q[1];
sx q[1];
rz(-2.0400679) q[1];
sx q[1];
rz(2.0371425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0898665) q[0];
sx q[0];
rz(-1.7891004) q[0];
sx q[0];
rz(1.0132683) q[0];
rz(-pi) q[1];
rz(2.274674) q[2];
sx q[2];
rz(-2.6353177) q[2];
sx q[2];
rz(0.75870721) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23984662) q[1];
sx q[1];
rz(-0.86565986) q[1];
sx q[1];
rz(2.2589931) q[1];
x q[2];
rz(-2.1451925) q[3];
sx q[3];
rz(-1.8339012) q[3];
sx q[3];
rz(0.97364473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1793648) q[2];
sx q[2];
rz(-2.4374375) q[2];
sx q[2];
rz(-3.1000225) q[2];
rz(-0.19503322) q[3];
sx q[3];
rz(-2.8988367) q[3];
sx q[3];
rz(2.354505) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4389909) q[0];
sx q[0];
rz(-0.85346237) q[0];
sx q[0];
rz(-2.3882197) q[0];
rz(-2.0808749) q[1];
sx q[1];
rz(-0.80442387) q[1];
sx q[1];
rz(-0.20739584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78847892) q[0];
sx q[0];
rz(-0.99035701) q[0];
sx q[0];
rz(0.81932318) q[0];
rz(1.2815525) q[2];
sx q[2];
rz(-2.0194032) q[2];
sx q[2];
rz(1.9236652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3606229) q[1];
sx q[1];
rz(-1.0252254) q[1];
sx q[1];
rz(1.2956133) q[1];
x q[2];
rz(-3.1025299) q[3];
sx q[3];
rz(-1.5785346) q[3];
sx q[3];
rz(1.4259065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9994026) q[2];
sx q[2];
rz(-1.9839857) q[2];
sx q[2];
rz(0.055214971) q[2];
rz(0.84290543) q[3];
sx q[3];
rz(-0.75777811) q[3];
sx q[3];
rz(2.260476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76999369) q[0];
sx q[0];
rz(-0.70796767) q[0];
sx q[0];
rz(-1.9644894) q[0];
rz(-0.87989315) q[1];
sx q[1];
rz(-2.8113007) q[1];
sx q[1];
rz(3.0807307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494722) q[0];
sx q[0];
rz(-0.15639399) q[0];
sx q[0];
rz(-0.91856399) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6598236) q[2];
sx q[2];
rz(-1.2854115) q[2];
sx q[2];
rz(-2.9586359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3538554) q[1];
sx q[1];
rz(-2.1365878) q[1];
sx q[1];
rz(-1.7090399) q[1];
x q[2];
rz(1.3065606) q[3];
sx q[3];
rz(-1.6420393) q[3];
sx q[3];
rz(-0.83282214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0941144) q[2];
sx q[2];
rz(-0.55314174) q[2];
sx q[2];
rz(-1.6806357) q[2];
rz(-0.85584062) q[3];
sx q[3];
rz(-2.1688921) q[3];
sx q[3];
rz(0.87876764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1020553) q[0];
sx q[0];
rz(-2.2886031) q[0];
sx q[0];
rz(3.1374186) q[0];
rz(-1.8904842) q[1];
sx q[1];
rz(-3.0079542) q[1];
sx q[1];
rz(0.68152308) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96723667) q[0];
sx q[0];
rz(-2.3612494) q[0];
sx q[0];
rz(1.1514949) q[0];
rz(-pi) q[1];
rz(2.2987492) q[2];
sx q[2];
rz(-2.3327391) q[2];
sx q[2];
rz(-3.0191772) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3423667) q[1];
sx q[1];
rz(-1.384642) q[1];
sx q[1];
rz(-2.8543945) q[1];
rz(-pi) q[2];
rz(0.32855804) q[3];
sx q[3];
rz(-2.3577981) q[3];
sx q[3];
rz(0.72600049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.556813) q[2];
sx q[2];
rz(-2.6656373) q[2];
sx q[2];
rz(1.3947831) q[2];
rz(-2.7252588) q[3];
sx q[3];
rz(-1.2952015) q[3];
sx q[3];
rz(0.39286119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75547472) q[0];
sx q[0];
rz(-0.32805726) q[0];
sx q[0];
rz(2.1395444) q[0];
rz(2.877032) q[1];
sx q[1];
rz(-1.8160276) q[1];
sx q[1];
rz(-1.4904259) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0984983) q[0];
sx q[0];
rz(-1.4401471) q[0];
sx q[0];
rz(2.5976973) q[0];
x q[1];
rz(0.15204805) q[2];
sx q[2];
rz(-0.81405241) q[2];
sx q[2];
rz(-1.1035827) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7539053) q[1];
sx q[1];
rz(-1.4567944) q[1];
sx q[1];
rz(-0.08367678) q[1];
rz(-pi) q[2];
rz(3.0496554) q[3];
sx q[3];
rz(-1.6306393) q[3];
sx q[3];
rz(-2.6106604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.010290535) q[2];
sx q[2];
rz(-1.4069858) q[2];
sx q[2];
rz(-1.4170125) q[2];
rz(0.33306444) q[3];
sx q[3];
rz(-2.371623) q[3];
sx q[3];
rz(-1.3781579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604816) q[0];
sx q[0];
rz(-1.2662553) q[0];
sx q[0];
rz(1.8929831) q[0];
rz(-0.026451182) q[1];
sx q[1];
rz(-1.5298264) q[1];
sx q[1];
rz(-1.5419921) q[1];
rz(-1.1556861) q[2];
sx q[2];
rz(-0.90183707) q[2];
sx q[2];
rz(3.1169008) q[2];
rz(-1.2534441) q[3];
sx q[3];
rz(-0.98012925) q[3];
sx q[3];
rz(-2.1914225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
