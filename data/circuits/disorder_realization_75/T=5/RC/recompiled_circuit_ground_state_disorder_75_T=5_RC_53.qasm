OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(6.5919696) q[0];
sx q[0];
rz(6.0433521) q[0];
rz(-0.56107768) q[1];
sx q[1];
rz(-2.3993888) q[1];
sx q[1];
rz(-1.5129369) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1334907) q[0];
sx q[0];
rz(-1.2043796) q[0];
sx q[0];
rz(0.13704637) q[0];
rz(1.0885029) q[2];
sx q[2];
rz(-0.52327195) q[2];
sx q[2];
rz(0.54973093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4862772) q[1];
sx q[1];
rz(-1.4239171) q[1];
sx q[1];
rz(2.6492061) q[1];
rz(-3.0835864) q[3];
sx q[3];
rz(-1.7388302) q[3];
sx q[3];
rz(-1.1917654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1538887) q[2];
sx q[2];
rz(-0.99606267) q[2];
sx q[2];
rz(-1.055701) q[2];
rz(2.740247) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(1.5041941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317247) q[0];
sx q[0];
rz(-0.26917502) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(-0.74613219) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(3.0768118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822183) q[0];
sx q[0];
rz(-1.1752593) q[0];
sx q[0];
rz(0.48898029) q[0];
rz(1.4006813) q[2];
sx q[2];
rz(-0.89808849) q[2];
sx q[2];
rz(-1.1748888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8967817) q[1];
sx q[1];
rz(-2.9771388) q[1];
sx q[1];
rz(-2.0816148) q[1];
x q[2];
rz(1.719789) q[3];
sx q[3];
rz(-0.21315609) q[3];
sx q[3];
rz(0.085863559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88227162) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-3.1003013) q[2];
rz(2.0222372) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(-2.4422755) q[0];
rz(-0.62272561) q[1];
sx q[1];
rz(-0.8131578) q[1];
sx q[1];
rz(-2.8831388) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7454433) q[0];
sx q[0];
rz(-0.93702261) q[0];
sx q[0];
rz(1.5859912) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11451377) q[2];
sx q[2];
rz(-0.92161676) q[2];
sx q[2];
rz(-0.69930062) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.099951) q[1];
sx q[1];
rz(-2.0196242) q[1];
sx q[1];
rz(-1.5931409) q[1];
rz(-0.036470099) q[3];
sx q[3];
rz(-1.1410645) q[3];
sx q[3];
rz(0.018485395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93566018) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(0.70651954) q[2];
rz(-2.5126854) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(-2.0188873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(-2.1674147) q[0];
rz(-1.4840508) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(2.1407703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054996252) q[0];
sx q[0];
rz(-0.97584134) q[0];
sx q[0];
rz(-0.18878285) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77299574) q[2];
sx q[2];
rz(-2.2064035) q[2];
sx q[2];
rz(2.5781812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19480579) q[1];
sx q[1];
rz(-1.1071702) q[1];
sx q[1];
rz(-0.63321106) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9934512) q[3];
sx q[3];
rz(-0.32399789) q[3];
sx q[3];
rz(1.2202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4546844) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(-2.3243813) q[2];
rz(2.5456083) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(1.7433085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3572094) q[0];
sx q[0];
rz(-1.4056118) q[0];
sx q[0];
rz(1.0719517) q[0];
rz(-2.0042073) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(0.17328182) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77861518) q[0];
sx q[0];
rz(-0.6181051) q[0];
sx q[0];
rz(-2.419986) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91867515) q[2];
sx q[2];
rz(-2.5498418) q[2];
sx q[2];
rz(2.7482035) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4049591) q[1];
sx q[1];
rz(-1.4990605) q[1];
sx q[1];
rz(-0.81291764) q[1];
rz(-pi) q[2];
rz(-1.6638443) q[3];
sx q[3];
rz(-1.4085839) q[3];
sx q[3];
rz(-1.7013659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5002354) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(0.87453169) q[2];
rz(2.4747961) q[3];
sx q[3];
rz(-1.6801445) q[3];
sx q[3];
rz(0.82505208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2905529) q[0];
sx q[0];
rz(-2.9930826) q[0];
sx q[0];
rz(0.17459757) q[0];
rz(-0.54706508) q[1];
sx q[1];
rz(-2.6262296) q[1];
sx q[1];
rz(-1.2778767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9491315) q[0];
sx q[0];
rz(-1.8572079) q[0];
sx q[0];
rz(0.064563036) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75210394) q[2];
sx q[2];
rz(-1.9255203) q[2];
sx q[2];
rz(-2.6238476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6656832) q[1];
sx q[1];
rz(-1.922125) q[1];
sx q[1];
rz(-2.372588) q[1];
rz(-pi) q[2];
rz(-0.84482362) q[3];
sx q[3];
rz(-0.81171821) q[3];
sx q[3];
rz(-0.38540781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7586907) q[2];
sx q[2];
rz(-2.8768657) q[2];
sx q[2];
rz(-2.7721789) q[2];
rz(2.3139125) q[3];
sx q[3];
rz(-1.3134495) q[3];
sx q[3];
rz(2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740042) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(-0.23707238) q[0];
rz(-1.392662) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(-1.0038092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0802409) q[0];
sx q[0];
rz(-0.46886849) q[0];
sx q[0];
rz(2.7659225) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71356524) q[2];
sx q[2];
rz(-1.9305561) q[2];
sx q[2];
rz(0.52116115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29664111) q[1];
sx q[1];
rz(-1.2403204) q[1];
sx q[1];
rz(-2.5806694) q[1];
rz(-0.69606651) q[3];
sx q[3];
rz(-0.68361527) q[3];
sx q[3];
rz(-3.0974577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5281333) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(-2.6252739) q[2];
rz(-2.3033219) q[3];
sx q[3];
rz(-1.4811938) q[3];
sx q[3];
rz(0.18394884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2876005) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(0.91424346) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.4062358) q[1];
sx q[1];
rz(-0.90726888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7210684) q[0];
sx q[0];
rz(-2.4926909) q[0];
sx q[0];
rz(-1.6008928) q[0];
rz(-2.5261717) q[2];
sx q[2];
rz(-1.713101) q[2];
sx q[2];
rz(2.1880544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3104035) q[1];
sx q[1];
rz(-2.0210279) q[1];
sx q[1];
rz(-0.72527171) q[1];
rz(-pi) q[2];
rz(0.97213718) q[3];
sx q[3];
rz(-1.3479509) q[3];
sx q[3];
rz(0.11892648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3153136) q[2];
sx q[2];
rz(-1.4297012) q[2];
sx q[2];
rz(-0.86177525) q[2];
rz(-1.9050542) q[3];
sx q[3];
rz(-0.084241353) q[3];
sx q[3];
rz(-1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64060193) q[0];
sx q[0];
rz(-2.3589098) q[0];
sx q[0];
rz(1.030141) q[0];
rz(2.8252699) q[1];
sx q[1];
rz(-1.373469) q[1];
sx q[1];
rz(2.8533459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9194946) q[0];
sx q[0];
rz(-1.7159425) q[0];
sx q[0];
rz(-2.7371489) q[0];
rz(2.4924303) q[2];
sx q[2];
rz(-1.9948975) q[2];
sx q[2];
rz(2.8555388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0287778) q[1];
sx q[1];
rz(-0.76422404) q[1];
sx q[1];
rz(2.1995596) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.796631) q[3];
sx q[3];
rz(-0.91167456) q[3];
sx q[3];
rz(-0.82443217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2890702) q[2];
sx q[2];
rz(-1.728629) q[2];
sx q[2];
rz(0.22656013) q[2];
rz(-1.6864927) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(-0.69211012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15014547) q[0];
sx q[0];
rz(-1.2757855) q[0];
sx q[0];
rz(-2.5164497) q[0];
rz(-1.217968) q[1];
sx q[1];
rz(-0.70911276) q[1];
sx q[1];
rz(1.2120754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21975133) q[0];
sx q[0];
rz(-0.80972396) q[0];
sx q[0];
rz(2.452616) q[0];
rz(-0.67135129) q[2];
sx q[2];
rz(-1.3666479) q[2];
sx q[2];
rz(-2.8189557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21704504) q[1];
sx q[1];
rz(-0.44483652) q[1];
sx q[1];
rz(-1.7052827) q[1];
rz(-1.8623705) q[3];
sx q[3];
rz(-1.2816888) q[3];
sx q[3];
rz(1.5856727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43442279) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(1.0506857) q[2];
rz(0.55189842) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(-1.8570073) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(2.3552409) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(1.9342061) q[2];
sx q[2];
rz(-1.5354234) q[2];
sx q[2];
rz(-1.0101049) q[2];
rz(0.10765392) q[3];
sx q[3];
rz(-0.68895491) q[3];
sx q[3];
rz(0.72910492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
