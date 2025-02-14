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
rz(-1.6406887) q[0];
sx q[0];
rz(-2.357643) q[0];
sx q[0];
rz(0.10070237) q[0];
rz(-4.5667629) q[1];
sx q[1];
rz(6.9520091) q[1];
sx q[1];
rz(7.8520757) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9933321) q[0];
sx q[0];
rz(-1.0251704) q[0];
sx q[0];
rz(-1.7094897) q[0];
rz(-pi) q[1];
rz(-2.8694081) q[2];
sx q[2];
rz(-1.6707175) q[2];
sx q[2];
rz(1.8109494) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6493167) q[1];
sx q[1];
rz(-1.4054978) q[1];
sx q[1];
rz(-0.16053546) q[1];
rz(-pi) q[2];
rz(-0.80683462) q[3];
sx q[3];
rz(-0.63043046) q[3];
sx q[3];
rz(-0.85496074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0147741) q[2];
sx q[2];
rz(-2.7744881) q[2];
sx q[2];
rz(1.5300306) q[2];
rz(-1.4989467) q[3];
sx q[3];
rz(-0.25224125) q[3];
sx q[3];
rz(2.4422395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67795578) q[0];
sx q[0];
rz(-1.3709443) q[0];
sx q[0];
rz(-2.875705) q[0];
rz(1.6193341) q[1];
sx q[1];
rz(-1.8306754) q[1];
sx q[1];
rz(1.6552077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9135654) q[0];
sx q[0];
rz(-1.5633213) q[0];
sx q[0];
rz(2.3662503) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7712542) q[2];
sx q[2];
rz(-2.5458434) q[2];
sx q[2];
rz(-1.0919613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.34936541) q[1];
sx q[1];
rz(-2.0369895) q[1];
sx q[1];
rz(-1.258205) q[1];
rz(-pi) q[2];
rz(0.19274917) q[3];
sx q[3];
rz(-2.0348573) q[3];
sx q[3];
rz(-1.0100067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8365525) q[2];
sx q[2];
rz(-2.063648) q[2];
sx q[2];
rz(-3.0744699) q[2];
rz(0.51131311) q[3];
sx q[3];
rz(-1.3498243) q[3];
sx q[3];
rz(-2.2093723) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2387282) q[0];
sx q[0];
rz(-1.6599864) q[0];
sx q[0];
rz(-0.1486775) q[0];
rz(-0.45064926) q[1];
sx q[1];
rz(-1.5573749) q[1];
sx q[1];
rz(0.91948909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7203569) q[0];
sx q[0];
rz(-1.8381869) q[0];
sx q[0];
rz(-1.8084099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2581159) q[2];
sx q[2];
rz(-1.4191437) q[2];
sx q[2];
rz(0.54017457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1053523) q[1];
sx q[1];
rz(-0.72004623) q[1];
sx q[1];
rz(-0.68036637) q[1];
x q[2];
rz(0.50289233) q[3];
sx q[3];
rz(-0.38030312) q[3];
sx q[3];
rz(1.2740434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86842361) q[2];
sx q[2];
rz(-2.0563075) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3530537) q[0];
sx q[0];
rz(-1.5667229) q[0];
sx q[0];
rz(2.5153611) q[0];
rz(-1.7010472) q[1];
sx q[1];
rz(-2.3150621) q[1];
sx q[1];
rz(-2.4339719) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6976172) q[0];
sx q[0];
rz(-1.6559905) q[0];
sx q[0];
rz(0.8449239) q[0];
rz(2.3578016) q[2];
sx q[2];
rz(-1.8194799) q[2];
sx q[2];
rz(2.3907106) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(13/(11*pi)) q[1];
sx q[1];
rz(-1.4896684) q[1];
sx q[1];
rz(2.8721832) q[1];
rz(0.032010664) q[3];
sx q[3];
rz(-1.0860476) q[3];
sx q[3];
rz(2.8196991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6669199) q[2];
sx q[2];
rz(-1.8644574) q[2];
sx q[2];
rz(2.7478768) q[2];
rz(0.79931021) q[3];
sx q[3];
rz(-2.7768713) q[3];
sx q[3];
rz(1.6126311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2213152) q[0];
sx q[0];
rz(-2.3096313) q[0];
sx q[0];
rz(-3.063391) q[0];
rz(2.4529264) q[1];
sx q[1];
rz(-0.93910256) q[1];
sx q[1];
rz(2.2130373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16958848) q[0];
sx q[0];
rz(-1.0079591) q[0];
sx q[0];
rz(1.2006951) q[0];
rz(-1.9276728) q[2];
sx q[2];
rz(-1.8584826) q[2];
sx q[2];
rz(-1.9299802) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2714933) q[1];
sx q[1];
rz(-2.6989991) q[1];
sx q[1];
rz(2.6111994) q[1];
rz(-pi) q[2];
rz(1.171204) q[3];
sx q[3];
rz(-2.152123) q[3];
sx q[3];
rz(2.9188398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81386217) q[2];
sx q[2];
rz(-0.50945115) q[2];
sx q[2];
rz(2.4895085) q[2];
rz(-2.4364575) q[3];
sx q[3];
rz(-1.6465681) q[3];
sx q[3];
rz(-2.7130073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.240165) q[0];
sx q[0];
rz(-0.21900284) q[0];
sx q[0];
rz(-1.3448311) q[0];
rz(-2.6118458) q[1];
sx q[1];
rz(-1.2836645) q[1];
sx q[1];
rz(-1.8240428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12152762) q[0];
sx q[0];
rz(-0.99666599) q[0];
sx q[0];
rz(-0.0962538) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.797768) q[2];
sx q[2];
rz(-1.3161873) q[2];
sx q[2];
rz(1.3019778) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3403794) q[1];
sx q[1];
rz(-1.4487584) q[1];
sx q[1];
rz(1.9122047) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2413792) q[3];
sx q[3];
rz(-0.53295153) q[3];
sx q[3];
rz(0.50607133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.162447) q[2];
sx q[2];
rz(-1.8628758) q[2];
sx q[2];
rz(-2.0549959) q[2];
rz(-0.094001683) q[3];
sx q[3];
rz(-1.1007997) q[3];
sx q[3];
rz(0.47959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78733665) q[0];
sx q[0];
rz(-0.84699637) q[0];
sx q[0];
rz(-0.081548668) q[0];
rz(1.6185919) q[1];
sx q[1];
rz(-1.0146419) q[1];
sx q[1];
rz(-2.5340705) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475246) q[0];
sx q[0];
rz(-2.3660064) q[0];
sx q[0];
rz(-1.8344384) q[0];
rz(-2.6742769) q[2];
sx q[2];
rz(-1.226034) q[2];
sx q[2];
rz(-1.4490102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54840344) q[1];
sx q[1];
rz(-1.4680937) q[1];
sx q[1];
rz(0.013506933) q[1];
x q[2];
rz(-1.4942985) q[3];
sx q[3];
rz(-1.943265) q[3];
sx q[3];
rz(1.6157762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.088323204) q[2];
sx q[2];
rz(-0.20603267) q[2];
sx q[2];
rz(2.6982809) q[2];
rz(-0.045470227) q[3];
sx q[3];
rz(-1.1683522) q[3];
sx q[3];
rz(-0.62062353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670253) q[0];
sx q[0];
rz(-0.98282951) q[0];
sx q[0];
rz(2.5734651) q[0];
rz(1.8727632) q[1];
sx q[1];
rz(-1.1548235) q[1];
sx q[1];
rz(2.5297129) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1639742) q[0];
sx q[0];
rz(-1.9477193) q[0];
sx q[0];
rz(-1.4758171) q[0];
rz(-pi) q[1];
rz(-2.1558495) q[2];
sx q[2];
rz(-0.32801706) q[2];
sx q[2];
rz(2.8190921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5962741) q[1];
sx q[1];
rz(-0.55247149) q[1];
sx q[1];
rz(-1.0956863) q[1];
rz(-2.967784) q[3];
sx q[3];
rz(-1.7663071) q[3];
sx q[3];
rz(1.3427092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.890471) q[2];
sx q[2];
rz(-1.2476363) q[2];
sx q[2];
rz(0.20784155) q[2];
rz(-3.1244997) q[3];
sx q[3];
rz(-1.0216252) q[3];
sx q[3];
rz(-2.0489073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6297778) q[0];
sx q[0];
rz(-2.9371174) q[0];
sx q[0];
rz(-1.7275607) q[0];
rz(-3.0746225) q[1];
sx q[1];
rz(-1.8194865) q[1];
sx q[1];
rz(0.13967839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0020376677) q[0];
sx q[0];
rz(-1.25601) q[0];
sx q[0];
rz(2.806753) q[0];
rz(-2.5508444) q[2];
sx q[2];
rz(-0.70066888) q[2];
sx q[2];
rz(-0.76393647) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5433257) q[1];
sx q[1];
rz(-1.3451335) q[1];
sx q[1];
rz(1.78701) q[1];
rz(1.8627251) q[3];
sx q[3];
rz(-0.17660429) q[3];
sx q[3];
rz(-2.7677329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4329873) q[2];
sx q[2];
rz(-1.8726417) q[2];
sx q[2];
rz(2.404786) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0438185) q[0];
sx q[0];
rz(-2.173824) q[0];
sx q[0];
rz(1.9737825) q[0];
rz(-0.66520005) q[1];
sx q[1];
rz(-1.2780814) q[1];
sx q[1];
rz(-2.1591689) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3176757) q[0];
sx q[0];
rz(-2.060776) q[0];
sx q[0];
rz(-0.22725003) q[0];
rz(-pi) q[1];
rz(2.5024611) q[2];
sx q[2];
rz(-1.7878727) q[2];
sx q[2];
rz(-2.9242196) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9537045) q[1];
sx q[1];
rz(-2.4361103) q[1];
sx q[1];
rz(-0.097037434) q[1];
rz(-3.0741698) q[3];
sx q[3];
rz(-1.1734087) q[3];
sx q[3];
rz(1.8734166) q[3];
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
rz(-1.6046074) q[2];
rz(1.0200621) q[3];
sx q[3];
rz(-1.5217109) q[3];
sx q[3];
rz(-0.1608688) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0237324) q[0];
sx q[0];
rz(-1.5845789) q[0];
sx q[0];
rz(-1.3316863) q[0];
rz(2.3125519) q[1];
sx q[1];
rz(-1.6571028) q[1];
sx q[1];
rz(2.170457) q[1];
rz(-0.6777505) q[2];
sx q[2];
rz(-1.2919147) q[2];
sx q[2];
rz(-0.80703296) q[2];
rz(-0.8603596) q[3];
sx q[3];
rz(-2.5353236) q[3];
sx q[3];
rz(-1.3154588) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
