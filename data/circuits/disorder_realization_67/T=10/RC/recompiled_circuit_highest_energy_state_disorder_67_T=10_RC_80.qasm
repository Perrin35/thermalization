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
rz(-3.0408903) q[0];
rz(-4.5667629) q[1];
sx q[1];
rz(6.9520091) q[1];
sx q[1];
rz(7.8520757) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14826058) q[0];
sx q[0];
rz(-1.0251704) q[0];
sx q[0];
rz(-1.7094897) q[0];
x q[1];
rz(-2.8694081) q[2];
sx q[2];
rz(-1.6707175) q[2];
sx q[2];
rz(-1.3306432) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4922759) q[1];
sx q[1];
rz(-1.4054978) q[1];
sx q[1];
rz(-2.9810572) q[1];
x q[2];
rz(2.334758) q[3];
sx q[3];
rz(-2.5111622) q[3];
sx q[3];
rz(0.85496074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0147741) q[2];
sx q[2];
rz(-2.7744881) q[2];
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
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67795578) q[0];
sx q[0];
rz(-1.3709443) q[0];
sx q[0];
rz(-2.875705) q[0];
rz(1.5222585) q[1];
sx q[1];
rz(-1.8306754) q[1];
sx q[1];
rz(-1.6552077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33544316) q[0];
sx q[0];
rz(-2.3461113) q[0];
sx q[0];
rz(1.5812627) q[0];
x q[1];
rz(1.8114016) q[2];
sx q[2];
rz(-2.1212656) q[2];
sx q[2];
rz(-1.5305504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7759018) q[1];
sx q[1];
rz(-1.8490825) q[1];
sx q[1];
rz(0.48643087) q[1];
rz(-pi) q[2];
rz(-1.2052746) q[3];
sx q[3];
rz(-2.6418002) q[3];
sx q[3];
rz(-2.5427853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3050401) q[2];
sx q[2];
rz(-1.0779447) q[2];
sx q[2];
rz(-3.0744699) q[2];
rz(-0.51131311) q[3];
sx q[3];
rz(-1.7917683) q[3];
sx q[3];
rz(-2.2093723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9028645) q[0];
sx q[0];
rz(-1.6599864) q[0];
sx q[0];
rz(0.1486775) q[0];
rz(2.6909434) q[1];
sx q[1];
rz(-1.5842178) q[1];
sx q[1];
rz(2.2221036) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4212358) q[0];
sx q[0];
rz(-1.8381869) q[0];
sx q[0];
rz(1.3331828) q[0];
rz(-pi) q[1];
rz(1.3344263) q[2];
sx q[2];
rz(-0.70118517) q[2];
sx q[2];
rz(-1.2126306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0362404) q[1];
sx q[1];
rz(-2.4215464) q[1];
sx q[1];
rz(-2.4612263) q[1];
x q[2];
rz(-2.804677) q[3];
sx q[3];
rz(-1.7506699) q[3];
sx q[3];
rz(0.76897393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86842361) q[2];
sx q[2];
rz(-1.0852852) q[2];
sx q[2];
rz(0.85624179) q[2];
rz(-0.87598261) q[3];
sx q[3];
rz(-0.3951422) q[3];
sx q[3];
rz(-1.9877079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-2.3150621) q[1];
sx q[1];
rz(2.4339719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9191878) q[0];
sx q[0];
rz(-2.4116431) q[0];
sx q[0];
rz(-1.6987503) q[0];
x q[1];
rz(-0.78379102) q[2];
sx q[2];
rz(-1.3221127) q[2];
sx q[2];
rz(0.75088203) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.172239) q[1];
sx q[1];
rz(-1.8392977) q[1];
sx q[1];
rz(-1.4866465) q[1];
rz(1.510101) q[3];
sx q[3];
rz(-2.6558721) q[3];
sx q[3];
rz(-2.751089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6669199) q[2];
sx q[2];
rz(-1.8644574) q[2];
sx q[2];
rz(-0.39371583) q[2];
rz(2.3422824) q[3];
sx q[3];
rz(-0.36472133) q[3];
sx q[3];
rz(1.6126311) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2213152) q[0];
sx q[0];
rz(-0.83196139) q[0];
sx q[0];
rz(0.078201683) q[0];
rz(2.4529264) q[1];
sx q[1];
rz(-0.93910256) q[1];
sx q[1];
rz(-0.92855531) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9445223) q[0];
sx q[0];
rz(-1.2598979) q[0];
sx q[0];
rz(-2.5466555) q[0];
rz(-pi) q[1];
rz(-1.2139198) q[2];
sx q[2];
rz(-1.2831101) q[2];
sx q[2];
rz(1.2116125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87009939) q[1];
sx q[1];
rz(-0.44259354) q[1];
sx q[1];
rz(0.53039329) q[1];
rz(2.5220248) q[3];
sx q[3];
rz(-1.2396401) q[3];
sx q[3];
rz(1.1201657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3277305) q[2];
sx q[2];
rz(-0.50945115) q[2];
sx q[2];
rz(2.4895085) q[2];
rz(0.70513519) q[3];
sx q[3];
rz(-1.4950246) q[3];
sx q[3];
rz(-0.42858538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9014277) q[0];
sx q[0];
rz(-0.21900284) q[0];
sx q[0];
rz(1.7967615) q[0];
rz(-0.52974686) q[1];
sx q[1];
rz(-1.2836645) q[1];
sx q[1];
rz(1.8240428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.844125) q[0];
sx q[0];
rz(-0.58124761) q[0];
sx q[0];
rz(1.4232924) q[0];
rz(-pi) q[1];
rz(1.3438247) q[2];
sx q[2];
rz(-1.3161873) q[2];
sx q[2];
rz(-1.8396149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81281137) q[1];
sx q[1];
rz(-1.9095632) q[1];
sx q[1];
rz(0.12943204) q[1];
x q[2];
rz(2.9530355) q[3];
sx q[3];
rz(-1.0692763) q[3];
sx q[3];
rz(-0.88392183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.162447) q[2];
sx q[2];
rz(-1.2787168) q[2];
sx q[2];
rz(1.0865967) q[2];
rz(-3.047591) q[3];
sx q[3];
rz(-1.1007997) q[3];
sx q[3];
rz(-0.47959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78733665) q[0];
sx q[0];
rz(-0.84699637) q[0];
sx q[0];
rz(0.081548668) q[0];
rz(1.5230007) q[1];
sx q[1];
rz(-1.0146419) q[1];
sx q[1];
rz(-0.60752216) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089585) q[0];
sx q[0];
rz(-2.3130401) q[0];
sx q[0];
rz(0.25018042) q[0];
rz(1.9532246) q[2];
sx q[2];
rz(-2.0086346) q[2];
sx q[2];
rz(3.0944648) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1178149) q[1];
sx q[1];
rz(-1.5842321) q[1];
sx q[1];
rz(1.4680844) q[1];
rz(-0.19316407) q[3];
sx q[3];
rz(-0.37988362) q[3];
sx q[3];
rz(-1.8233606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.088323204) q[2];
sx q[2];
rz(-0.20603267) q[2];
sx q[2];
rz(2.6982809) q[2];
rz(-3.0961224) q[3];
sx q[3];
rz(-1.9732405) q[3];
sx q[3];
rz(2.5209691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670253) q[0];
sx q[0];
rz(-0.98282951) q[0];
sx q[0];
rz(-0.56812754) q[0];
rz(1.2688295) q[1];
sx q[1];
rz(-1.1548235) q[1];
sx q[1];
rz(0.61187977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4418714) q[0];
sx q[0];
rz(-1.6590902) q[0];
sx q[0];
rz(-2.7631212) q[0];
rz(-2.955825) q[2];
sx q[2];
rz(-1.8427197) q[2];
sx q[2];
rz(2.853554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.438209) q[1];
sx q[1];
rz(-1.8132231) q[1];
sx q[1];
rz(-2.0722778) q[1];
rz(-pi) q[2];
rz(-0.17380865) q[3];
sx q[3];
rz(-1.7663071) q[3];
sx q[3];
rz(-1.3427092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.890471) q[2];
sx q[2];
rz(-1.2476363) q[2];
sx q[2];
rz(-2.9337511) q[2];
rz(-0.017092997) q[3];
sx q[3];
rz(-2.1199675) q[3];
sx q[3];
rz(-2.0489073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51181483) q[0];
sx q[0];
rz(-2.9371174) q[0];
sx q[0];
rz(1.7275607) q[0];
rz(3.0746225) q[1];
sx q[1];
rz(-1.8194865) q[1];
sx q[1];
rz(3.0019143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84565879) q[0];
sx q[0];
rz(-0.45546969) q[0];
sx q[0];
rz(2.3607872) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61105543) q[2];
sx q[2];
rz(-1.9381028) q[2];
sx q[2];
rz(-0.33314785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3744065) q[1];
sx q[1];
rz(-2.8303478) q[1];
sx q[1];
rz(2.3900476) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8627251) q[3];
sx q[3];
rz(-0.17660429) q[3];
sx q[3];
rz(2.7677329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4329873) q[2];
sx q[2];
rz(-1.2689509) q[2];
sx q[2];
rz(-0.73680669) q[2];
rz(-1.5271198) q[3];
sx q[3];
rz(-2.1325285) q[3];
sx q[3];
rz(-1.8710776) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0438185) q[0];
sx q[0];
rz(-2.173824) q[0];
sx q[0];
rz(-1.9737825) q[0];
rz(-2.4763926) q[1];
sx q[1];
rz(-1.2780814) q[1];
sx q[1];
rz(2.1591689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2806399) q[0];
sx q[0];
rz(-0.5361983) q[0];
sx q[0];
rz(1.1711189) q[0];
rz(-pi) q[1];
rz(2.5024611) q[2];
sx q[2];
rz(-1.35372) q[2];
sx q[2];
rz(-0.21737305) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.82653) q[1];
sx q[1];
rz(-0.86931397) q[1];
sx q[1];
rz(1.4884654) q[1];
rz(-pi) q[2];
x q[2];
rz(0.067422859) q[3];
sx q[3];
rz(-1.1734087) q[3];
sx q[3];
rz(-1.2681761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3490225) q[2];
sx q[2];
rz(-0.64705619) q[2];
sx q[2];
rz(1.6046074) q[2];
rz(-2.1215306) q[3];
sx q[3];
rz(-1.5217109) q[3];
sx q[3];
rz(2.9807239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0237324) q[0];
sx q[0];
rz(-1.5570138) q[0];
sx q[0];
rz(1.8099063) q[0];
rz(-2.3125519) q[1];
sx q[1];
rz(-1.4844898) q[1];
sx q[1];
rz(-0.97113562) q[1];
rz(0.42837684) q[2];
sx q[2];
rz(-2.4171792) q[2];
sx q[2];
rz(0.43424594) q[2];
rz(-2.7169124) q[3];
sx q[3];
rz(-1.1241354) q[3];
sx q[3];
rz(-0.507171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
