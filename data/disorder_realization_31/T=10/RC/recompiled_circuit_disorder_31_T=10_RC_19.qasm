OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52656093) q[0];
sx q[0];
rz(-2.5685413) q[0];
sx q[0];
rz(-0.84258643) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6095088) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(2.8340333) q[0];
rz(-pi) q[1];
rz(-0.61383944) q[2];
sx q[2];
rz(-1.5547353) q[2];
sx q[2];
rz(1.2889372) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9021437) q[1];
sx q[1];
rz(-1.7001171) q[1];
sx q[1];
rz(2.7620478) q[1];
x q[2];
rz(-2.5472766) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(-3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(1.1546086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-3.1161838) q[0];
sx q[0];
rz(2.4029762) q[0];
rz(-pi) q[1];
rz(-2.0915394) q[2];
sx q[2];
rz(-0.69486952) q[2];
sx q[2];
rz(-2.3827202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7696015) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(3.0197057) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9566544) q[3];
sx q[3];
rz(-1.3122845) q[3];
sx q[3];
rz(-1.0227433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(-2.2581805) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(0.60107636) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(1.0916969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1611623) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(-2.3479793) q[0];
rz(-0.91471471) q[2];
sx q[2];
rz(-1.2351742) q[2];
sx q[2];
rz(2.6732973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.162902) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(-2.6497926) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.320257) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(-1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-0.31035796) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46582857) q[0];
sx q[0];
rz(-0.40256631) q[0];
sx q[0];
rz(0.34253828) q[0];
rz(0.68508673) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(0.098066559) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.031361) q[1];
sx q[1];
rz(-2.7831315) q[1];
sx q[1];
rz(1.5178174) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1060018) q[3];
sx q[3];
rz(-1.4471874) q[3];
sx q[3];
rz(1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0115396) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(-1.594054) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(1.5769618) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2244959) q[0];
sx q[0];
rz(-1.3670237) q[0];
sx q[0];
rz(2.8208371) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7587897) q[2];
sx q[2];
rz(-0.88289875) q[2];
sx q[2];
rz(-2.5685513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2493077) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(0.73242964) q[1];
x q[2];
rz(1.8364041) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(0.8997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(1.4917096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314635) q[0];
sx q[0];
rz(-1.9995814) q[0];
sx q[0];
rz(1.462912) q[0];
rz(-pi) q[1];
rz(1.1509402) q[2];
sx q[2];
rz(-0.93517762) q[2];
sx q[2];
rz(2.2395526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3867492) q[1];
sx q[1];
rz(-2.1347087) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8317354) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(-0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(-3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(-0.89404026) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-3.1076028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695553) q[0];
sx q[0];
rz(-2.3453658) q[0];
sx q[0];
rz(-1.7982593) q[0];
x q[1];
rz(-0.88862822) q[2];
sx q[2];
rz(-2.3805328) q[2];
sx q[2];
rz(1.6737446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36044869) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(-3.0767246) q[1];
rz(-pi) q[2];
rz(-0.61492413) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(1.8254335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(0.87654385) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(2.7476655) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.4454909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60364265) q[0];
sx q[0];
rz(-2.0041487) q[0];
sx q[0];
rz(-3.135878) q[0];
rz(-pi) q[1];
rz(-0.94930737) q[2];
sx q[2];
rz(-2.0645803) q[2];
sx q[2];
rz(2.2311503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3280914) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(0.94533841) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-0.64627796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(-0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(2.8318185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4315902) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(2.8938328) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7477112) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(-2.6975346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8621091) q[1];
sx q[1];
rz(-0.21394193) q[1];
sx q[1];
rz(-1.2941542) q[1];
x q[2];
rz(-0.36848948) q[3];
sx q[3];
rz(-0.62928761) q[3];
sx q[3];
rz(-3.1143509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(-0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.9357095) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.6419798) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6584872) q[0];
sx q[0];
rz(-0.30880901) q[0];
sx q[0];
rz(1.9103861) q[0];
x q[1];
rz(0.76403107) q[2];
sx q[2];
rz(-1.6198297) q[2];
sx q[2];
rz(1.3438091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3824532) q[1];
sx q[1];
rz(-2.0331953) q[1];
sx q[1];
rz(0.23666246) q[1];
x q[2];
rz(-1.5564735) q[3];
sx q[3];
rz(-0.17281547) q[3];
sx q[3];
rz(0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-3.070667) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(-0.46438607) q[2];
sx q[2];
rz(-2.0225564) q[2];
sx q[2];
rz(0.49973942) q[2];
rz(-2.0675038) q[3];
sx q[3];
rz(-2.2278193) q[3];
sx q[3];
rz(-0.57886119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
