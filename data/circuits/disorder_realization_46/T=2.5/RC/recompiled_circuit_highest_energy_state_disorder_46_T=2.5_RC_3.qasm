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
rz(-2.6645339) q[0];
sx q[0];
rz(-2.8395489) q[0];
sx q[0];
rz(-1.8855236) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(-1.3060275) q[1];
sx q[1];
rz(2.0130872) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80917847) q[0];
sx q[0];
rz(-2.1398394) q[0];
sx q[0];
rz(2.8079985) q[0];
rz(-pi) q[1];
rz(-2.520325) q[2];
sx q[2];
rz(-1.8227907) q[2];
sx q[2];
rz(-0.14185837) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0417718) q[1];
sx q[1];
rz(-1.9536293) q[1];
sx q[1];
rz(-2.113093) q[1];
x q[2];
rz(0.63367356) q[3];
sx q[3];
rz(-1.9320935) q[3];
sx q[3];
rz(2.6665045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6065373) q[2];
sx q[2];
rz(-0.29416072) q[2];
sx q[2];
rz(1.199022) q[2];
rz(0.27896518) q[3];
sx q[3];
rz(-1.4594892) q[3];
sx q[3];
rz(3.103106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8798384) q[0];
sx q[0];
rz(-2.6848875) q[0];
sx q[0];
rz(-2.7754011) q[0];
rz(-1.2884033) q[1];
sx q[1];
rz(-0.90694702) q[1];
sx q[1];
rz(-2.8790976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3655626) q[0];
sx q[0];
rz(-1.7509116) q[0];
sx q[0];
rz(0.30091648) q[0];
x q[1];
rz(1.8196202) q[2];
sx q[2];
rz(-1.5142528) q[2];
sx q[2];
rz(-0.6565993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4387063) q[1];
sx q[1];
rz(-1.2186495) q[1];
sx q[1];
rz(-1.8789324) q[1];
rz(-pi) q[2];
rz(-0.025215312) q[3];
sx q[3];
rz(-1.8256911) q[3];
sx q[3];
rz(-2.1372634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3890106) q[2];
sx q[2];
rz(-1.2685308) q[2];
sx q[2];
rz(2.9627723) q[2];
rz(3.1255152) q[3];
sx q[3];
rz(-2.6060846) q[3];
sx q[3];
rz(0.9308365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43612424) q[0];
sx q[0];
rz(-3.0146764) q[0];
sx q[0];
rz(2.5674168) q[0];
rz(-0.15159675) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(-1.2281598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235665) q[0];
sx q[0];
rz(-1.0581338) q[0];
sx q[0];
rz(2.961757) q[0];
rz(-0.97790678) q[2];
sx q[2];
rz(-1.5139069) q[2];
sx q[2];
rz(-1.6757193) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0361402) q[1];
sx q[1];
rz(-2.5220118) q[1];
sx q[1];
rz(1.416227) q[1];
rz(-2.853573) q[3];
sx q[3];
rz(-1.1041118) q[3];
sx q[3];
rz(2.6096901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14153081) q[2];
sx q[2];
rz(-0.63566339) q[2];
sx q[2];
rz(-2.2504811) q[2];
rz(2.7530503) q[3];
sx q[3];
rz(-1.6308234) q[3];
sx q[3];
rz(-2.8019606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163088) q[0];
sx q[0];
rz(-0.044965222) q[0];
sx q[0];
rz(2.6594824) q[0];
rz(-0.66857839) q[1];
sx q[1];
rz(-1.5336288) q[1];
sx q[1];
rz(-2.5040928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8746045) q[0];
sx q[0];
rz(-0.98867304) q[0];
sx q[0];
rz(-0.28018392) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0257865) q[2];
sx q[2];
rz(-2.2409332) q[2];
sx q[2];
rz(1.1095123) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5920221) q[1];
sx q[1];
rz(-1.5369104) q[1];
sx q[1];
rz(-3.0646044) q[1];
x q[2];
rz(-2.1737764) q[3];
sx q[3];
rz(-0.42438904) q[3];
sx q[3];
rz(2.5924204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.94347) q[2];
sx q[2];
rz(-1.024647) q[2];
sx q[2];
rz(-3.0812954) q[2];
rz(-2.8382525) q[3];
sx q[3];
rz(-0.077849418) q[3];
sx q[3];
rz(-2.8764603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.1497134) q[0];
sx q[0];
rz(-0.35357058) q[0];
sx q[0];
rz(-2.7283332) q[0];
rz(-2.7464271) q[1];
sx q[1];
rz(-1.3576077) q[1];
sx q[1];
rz(-1.0023592) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43267879) q[0];
sx q[0];
rz(-2.1122547) q[0];
sx q[0];
rz(2.0869291) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4031314) q[2];
sx q[2];
rz(-2.9288026) q[2];
sx q[2];
rz(-2.8820702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2152948) q[1];
sx q[1];
rz(-1.7855254) q[1];
sx q[1];
rz(-2.2533827) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8931562) q[3];
sx q[3];
rz(-1.3731706) q[3];
sx q[3];
rz(-1.2569497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62785441) q[2];
sx q[2];
rz(-2.0444874) q[2];
sx q[2];
rz(-1.6045137) q[2];
rz(-1.9419935) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(2.3698923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62726778) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(-2.7728873) q[0];
rz(2.3983224) q[1];
sx q[1];
rz(-2.5786046) q[1];
sx q[1];
rz(-2.7875913) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10121541) q[0];
sx q[0];
rz(-2.4987552) q[0];
sx q[0];
rz(-1.4835711) q[0];
x q[1];
rz(2.6145489) q[2];
sx q[2];
rz(-1.8436925) q[2];
sx q[2];
rz(-0.54293406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56997715) q[1];
sx q[1];
rz(-1.0034402) q[1];
sx q[1];
rz(2.9302271) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50320585) q[3];
sx q[3];
rz(-0.70069289) q[3];
sx q[3];
rz(-2.1029187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.02504286) q[2];
sx q[2];
rz(-1.2094689) q[2];
sx q[2];
rz(-0.0054736007) q[2];
rz(-0.51760751) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(3.1143809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7356877) q[0];
sx q[0];
rz(-0.55957496) q[0];
sx q[0];
rz(-2.9285808) q[0];
rz(-1.6464015) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(-0.84943736) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6005858) q[0];
sx q[0];
rz(-1.5274421) q[0];
sx q[0];
rz(-2.5227273) q[0];
rz(-1.6140217) q[2];
sx q[2];
rz(-1.6014978) q[2];
sx q[2];
rz(-0.30393201) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7905459) q[1];
sx q[1];
rz(-0.67706579) q[1];
sx q[1];
rz(-0.76404567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3091503) q[3];
sx q[3];
rz(-2.5980686) q[3];
sx q[3];
rz(-1.3268918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9269632) q[2];
sx q[2];
rz(-1.233485) q[2];
sx q[2];
rz(-0.54369533) q[2];
rz(0.26116535) q[3];
sx q[3];
rz(-1.5800913) q[3];
sx q[3];
rz(0.15682596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549304) q[0];
sx q[0];
rz(-2.08827) q[0];
sx q[0];
rz(-2.9539811) q[0];
rz(2.018351) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(0.16256464) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190348) q[0];
sx q[0];
rz(-1.8936689) q[0];
sx q[0];
rz(-0.69299181) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0948903) q[2];
sx q[2];
rz(-1.419029) q[2];
sx q[2];
rz(-2.2693315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8423262) q[1];
sx q[1];
rz(-2.5262621) q[1];
sx q[1];
rz(-1.6493114) q[1];
rz(-pi) q[2];
rz(0.1127693) q[3];
sx q[3];
rz(-2.164132) q[3];
sx q[3];
rz(-1.6360374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8673914) q[2];
sx q[2];
rz(-0.31495366) q[2];
sx q[2];
rz(-0.20142041) q[2];
rz(2.6650186) q[3];
sx q[3];
rz(-1.3786992) q[3];
sx q[3];
rz(-0.99982888) q[3];
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
rz(1.1808712) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(-2.5647105) q[0];
rz(-2.8374529) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(-2.0374128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0074771) q[0];
sx q[0];
rz(-0.38950697) q[0];
sx q[0];
rz(1.3358119) q[0];
x q[1];
rz(-1.6921894) q[2];
sx q[2];
rz(-1.0052201) q[2];
sx q[2];
rz(2.7944416) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7420099) q[1];
sx q[1];
rz(-0.72232692) q[1];
sx q[1];
rz(-0.96256017) q[1];
x q[2];
rz(-1.9619759) q[3];
sx q[3];
rz(-1.5558387) q[3];
sx q[3];
rz(-2.7103031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96678174) q[2];
sx q[2];
rz(-2.8947783) q[2];
sx q[2];
rz(-2.4583869) q[2];
rz(1.5481868) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(-0.21827503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9504647) q[0];
sx q[0];
rz(-2.2948965) q[0];
sx q[0];
rz(2.6306613) q[0];
rz(-1.9966985) q[1];
sx q[1];
rz(-1.1868008) q[1];
sx q[1];
rz(-1.6330382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36105267) q[0];
sx q[0];
rz(-1.0298488) q[0];
sx q[0];
rz(-1.8622647) q[0];
rz(-pi) q[1];
rz(0.26692674) q[2];
sx q[2];
rz(-1.1014001) q[2];
sx q[2];
rz(1.5482582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.252987) q[1];
sx q[1];
rz(-1.3934416) q[1];
sx q[1];
rz(0.4528404) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4436117) q[3];
sx q[3];
rz(-0.38376402) q[3];
sx q[3];
rz(1.1205387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42023429) q[2];
sx q[2];
rz(-2.5098462) q[2];
sx q[2];
rz(-1.7101804) q[2];
rz(-3.1259649) q[3];
sx q[3];
rz(-1.8764007) q[3];
sx q[3];
rz(-2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90841993) q[0];
sx q[0];
rz(-1.40042) q[0];
sx q[0];
rz(-0.9216876) q[0];
rz(1.5064916) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(2.8534081) q[2];
sx q[2];
rz(-1.9389887) q[2];
sx q[2];
rz(-1.5252319) q[2];
rz(1.4931645) q[3];
sx q[3];
rz(-1.5356728) q[3];
sx q[3];
rz(0.079726797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
