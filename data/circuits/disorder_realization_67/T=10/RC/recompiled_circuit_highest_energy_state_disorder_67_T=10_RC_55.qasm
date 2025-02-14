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
rz(1.500904) q[0];
sx q[0];
rz(-0.7839497) q[0];
sx q[0];
rz(6.1824829) q[0];
rz(-4.5667629) q[1];
sx q[1];
rz(6.9520091) q[1];
sx q[1];
rz(7.8520757) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6467428) q[0];
sx q[0];
rz(-1.6892489) q[0];
sx q[0];
rz(0.54991566) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8694081) q[2];
sx q[2];
rz(-1.6707175) q[2];
sx q[2];
rz(-1.3306432) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71489635) q[1];
sx q[1];
rz(-2.9116803) q[1];
sx q[1];
rz(-0.80674361) q[1];
rz(-pi) q[2];
rz(2.6740736) q[3];
sx q[3];
rz(-2.010502) q[3];
sx q[3];
rz(1.7252418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12681857) q[2];
sx q[2];
rz(-0.36710456) q[2];
sx q[2];
rz(-1.611562) q[2];
rz(-1.642646) q[3];
sx q[3];
rz(-0.25224125) q[3];
sx q[3];
rz(-2.4422395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4636369) q[0];
sx q[0];
rz(-1.7706484) q[0];
sx q[0];
rz(-2.875705) q[0];
rz(1.6193341) q[1];
sx q[1];
rz(-1.8306754) q[1];
sx q[1];
rz(1.6552077) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2280272) q[0];
sx q[0];
rz(-1.5633213) q[0];
sx q[0];
rz(-0.77534239) q[0];
rz(-0.3703385) q[2];
sx q[2];
rz(-2.5458434) q[2];
sx q[2];
rz(2.0496313) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34936541) q[1];
sx q[1];
rz(-1.1046032) q[1];
sx q[1];
rz(1.8833877) q[1];
x q[2];
rz(-2.9488435) q[3];
sx q[3];
rz(-2.0348573) q[3];
sx q[3];
rz(-1.0100067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3050401) q[2];
sx q[2];
rz(-2.063648) q[2];
sx q[2];
rz(3.0744699) q[2];
rz(0.51131311) q[3];
sx q[3];
rz(-1.3498243) q[3];
sx q[3];
rz(-2.2093723) q[3];
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
rz(-pi) q[3];
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
rz(-0.2387282) q[0];
sx q[0];
rz(-1.4816062) q[0];
sx q[0];
rz(2.9929152) q[0];
rz(0.45064926) q[1];
sx q[1];
rz(-1.5573749) q[1];
sx q[1];
rz(-0.91948909) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085657619) q[0];
sx q[0];
rz(-1.3417804) q[0];
sx q[0];
rz(-2.8668501) q[0];
rz(-0.19520031) q[2];
sx q[2];
rz(-2.2487309) q[2];
sx q[2];
rz(-0.90724573) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1295848) q[1];
sx q[1];
rz(-1.9985481) q[1];
sx q[1];
rz(-2.5431551) q[1];
rz(2.6387003) q[3];
sx q[3];
rz(-2.7612895) q[3];
sx q[3];
rz(-1.8675493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86842361) q[2];
sx q[2];
rz(-1.0852852) q[2];
sx q[2];
rz(0.85624179) q[2];
rz(-2.26561) q[3];
sx q[3];
rz(-0.3951422) q[3];
sx q[3];
rz(1.9877079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7885389) q[0];
sx q[0];
rz(-1.5667229) q[0];
sx q[0];
rz(0.62623155) q[0];
rz(-1.4405454) q[1];
sx q[1];
rz(-0.82653058) q[1];
sx q[1];
rz(0.70762077) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6976172) q[0];
sx q[0];
rz(-1.4856021) q[0];
sx q[0];
rz(-2.2966688) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3452937) q[2];
sx q[2];
rz(-0.81419386) q[2];
sx q[2];
rz(-1.0619927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7654082) q[1];
sx q[1];
rz(-1.4896684) q[1];
sx q[1];
rz(2.8721832) q[1];
rz(-pi) q[2];
rz(1.0858363) q[3];
sx q[3];
rz(-1.5991181) q[3];
sx q[3];
rz(1.9076104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6669199) q[2];
sx q[2];
rz(-1.2771353) q[2];
sx q[2];
rz(0.39371583) q[2];
rz(2.3422824) q[3];
sx q[3];
rz(-2.7768713) q[3];
sx q[3];
rz(1.5289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2213152) q[0];
sx q[0];
rz(-2.3096313) q[0];
sx q[0];
rz(-0.078201683) q[0];
rz(-2.4529264) q[1];
sx q[1];
rz(-0.93910256) q[1];
sx q[1];
rz(-2.2130373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16958848) q[0];
sx q[0];
rz(-1.0079591) q[0];
sx q[0];
rz(-1.9408976) q[0];
rz(-pi) q[1];
rz(0.30588116) q[2];
sx q[2];
rz(-1.9123931) q[2];
sx q[2];
rz(2.677013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87009939) q[1];
sx q[1];
rz(-2.6989991) q[1];
sx q[1];
rz(-2.6111994) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5220248) q[3];
sx q[3];
rz(-1.9019526) q[3];
sx q[3];
rz(1.1201657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3277305) q[2];
sx q[2];
rz(-0.50945115) q[2];
sx q[2];
rz(0.65208411) q[2];
rz(0.70513519) q[3];
sx q[3];
rz(-1.4950246) q[3];
sx q[3];
rz(-0.42858538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9014277) q[0];
sx q[0];
rz(-0.21900284) q[0];
sx q[0];
rz(-1.7967615) q[0];
rz(-0.52974686) q[1];
sx q[1];
rz(-1.8579282) q[1];
sx q[1];
rz(-1.8240428) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.844125) q[0];
sx q[0];
rz(-2.560345) q[0];
sx q[0];
rz(1.7183003) q[0];
x q[1];
rz(-2.4286626) q[2];
sx q[2];
rz(-0.33944079) q[2];
sx q[2];
rz(-1.0975099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3403794) q[1];
sx q[1];
rz(-1.4487584) q[1];
sx q[1];
rz(1.9122047) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0798912) q[3];
sx q[3];
rz(-1.7359043) q[3];
sx q[3];
rz(0.77835876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9791457) q[2];
sx q[2];
rz(-1.8628758) q[2];
sx q[2];
rz(2.0549959) q[2];
rz(0.094001683) q[3];
sx q[3];
rz(-1.1007997) q[3];
sx q[3];
rz(-0.47959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78733665) q[0];
sx q[0];
rz(-0.84699637) q[0];
sx q[0];
rz(3.060044) q[0];
rz(-1.5230007) q[1];
sx q[1];
rz(-1.0146419) q[1];
sx q[1];
rz(-2.5340705) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23287671) q[0];
sx q[0];
rz(-1.7542782) q[0];
sx q[0];
rz(0.81277892) q[0];
x q[1];
rz(-2.4685988) q[2];
sx q[2];
rz(-2.5685326) q[2];
sx q[2];
rz(2.429643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1178149) q[1];
sx q[1];
rz(-1.5573606) q[1];
sx q[1];
rz(1.4680844) q[1];
x q[2];
rz(-0.37346249) q[3];
sx q[3];
rz(-1.6420396) q[3];
sx q[3];
rz(0.072865818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.088323204) q[2];
sx q[2];
rz(-2.93556) q[2];
sx q[2];
rz(2.6982809) q[2];
rz(3.0961224) q[3];
sx q[3];
rz(-1.1683522) q[3];
sx q[3];
rz(-0.62062353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670253) q[0];
sx q[0];
rz(-0.98282951) q[0];
sx q[0];
rz(-2.5734651) q[0];
rz(-1.8727632) q[1];
sx q[1];
rz(-1.1548235) q[1];
sx q[1];
rz(-2.5297129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2308919) q[0];
sx q[0];
rz(-0.38815) q[0];
sx q[0];
rz(0.2351454) q[0];
rz(-pi) q[1];
rz(-2.1558495) q[2];
sx q[2];
rz(-0.32801706) q[2];
sx q[2];
rz(2.8190921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7033836) q[1];
sx q[1];
rz(-1.8132231) q[1];
sx q[1];
rz(2.0722778) q[1];
x q[2];
rz(-0.85297008) q[3];
sx q[3];
rz(-2.880734) q[3];
sx q[3];
rz(1.0639695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2511217) q[2];
sx q[2];
rz(-1.8939563) q[2];
sx q[2];
rz(2.9337511) q[2];
rz(3.1244997) q[3];
sx q[3];
rz(-1.0216252) q[3];
sx q[3];
rz(2.0489073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.51181483) q[0];
sx q[0];
rz(-2.9371174) q[0];
sx q[0];
rz(-1.7275607) q[0];
rz(0.066970197) q[1];
sx q[1];
rz(-1.3221062) q[1];
sx q[1];
rz(-0.13967839) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801474) q[0];
sx q[0];
rz(-1.8885888) q[0];
sx q[0];
rz(-1.2387973) q[0];
rz(-pi) q[1];
rz(2.5508444) q[2];
sx q[2];
rz(-0.70066888) q[2];
sx q[2];
rz(0.76393647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.598267) q[1];
sx q[1];
rz(-1.3451335) q[1];
sx q[1];
rz(-1.3545827) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8627251) q[3];
sx q[3];
rz(-0.17660429) q[3];
sx q[3];
rz(0.3738598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4329873) q[2];
sx q[2];
rz(-1.2689509) q[2];
sx q[2];
rz(-2.404786) q[2];
rz(-1.6144729) q[3];
sx q[3];
rz(-1.0090642) q[3];
sx q[3];
rz(1.2705151) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0977741) q[0];
sx q[0];
rz(-2.173824) q[0];
sx q[0];
rz(-1.9737825) q[0];
rz(0.66520005) q[1];
sx q[1];
rz(-1.8635112) q[1];
sx q[1];
rz(0.98242378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3176757) q[0];
sx q[0];
rz(-2.060776) q[0];
sx q[0];
rz(0.22725003) q[0];
rz(-pi) q[1];
rz(-2.5024611) q[2];
sx q[2];
rz(-1.7878727) q[2];
sx q[2];
rz(-0.21737305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.82653) q[1];
sx q[1];
rz(-0.86931397) q[1];
sx q[1];
rz(-1.4884654) q[1];
rz(-pi) q[2];
rz(-1.7299557) q[3];
sx q[3];
rz(-2.7388262) q[3];
sx q[3];
rz(-1.7006766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7925702) q[2];
sx q[2];
rz(-2.4945365) q[2];
sx q[2];
rz(1.5369852) q[2];
rz(2.1215306) q[3];
sx q[3];
rz(-1.6198817) q[3];
sx q[3];
rz(-0.1608688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0237324) q[0];
sx q[0];
rz(-1.5570138) q[0];
sx q[0];
rz(1.8099063) q[0];
rz(2.3125519) q[1];
sx q[1];
rz(-1.6571028) q[1];
sx q[1];
rz(2.170457) q[1];
rz(-1.2185417) q[2];
sx q[2];
rz(-0.92377077) q[2];
sx q[2];
rz(0.9818264) q[2];
rz(-0.42468023) q[3];
sx q[3];
rz(-2.0174572) q[3];
sx q[3];
rz(2.6344217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
