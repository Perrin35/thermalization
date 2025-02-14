OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(-3.0734835) q[0];
sx q[0];
rz(2.3647302) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(-1.8196655) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711219) q[0];
sx q[0];
rz(-1.4323455) q[0];
sx q[0];
rz(-2.2884503) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.916531) q[2];
sx q[2];
rz(-2.0421197) q[2];
sx q[2];
rz(-0.87475785) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3275431) q[1];
sx q[1];
rz(-2.7316878) q[1];
sx q[1];
rz(-1.6442714) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1573651) q[3];
sx q[3];
rz(-2.7968458) q[3];
sx q[3];
rz(-2.1380827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.063735828) q[2];
sx q[2];
rz(-1.3602164) q[2];
sx q[2];
rz(-2.784101) q[2];
rz(-2.5743971) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(0.43958694) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95328632) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(-1.8081007) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(-0.83998799) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6463722) q[0];
sx q[0];
rz(-2.6073808) q[0];
sx q[0];
rz(1.87508) q[0];
x q[1];
rz(2.5469668) q[2];
sx q[2];
rz(-1.42808) q[2];
sx q[2];
rz(-0.12446257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6691554) q[1];
sx q[1];
rz(-1.0095114) q[1];
sx q[1];
rz(2.2538633) q[1];
rz(-pi) q[2];
rz(1.5072823) q[3];
sx q[3];
rz(-0.5507142) q[3];
sx q[3];
rz(-2.7515819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(-1.7355512) q[2];
rz(0.16215912) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(0.46419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(2.719847) q[0];
rz(-2.2987507) q[1];
sx q[1];
rz(-1.755736) q[1];
sx q[1];
rz(-2.8657894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3746412) q[0];
sx q[0];
rz(-1.8948484) q[0];
sx q[0];
rz(1.814591) q[0];
x q[1];
rz(-0.47068542) q[2];
sx q[2];
rz(-1.7690725) q[2];
sx q[2];
rz(-0.40775611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2308064) q[1];
sx q[1];
rz(-1.2503997) q[1];
sx q[1];
rz(-1.6256534) q[1];
rz(2.2725676) q[3];
sx q[3];
rz(-1.6459673) q[3];
sx q[3];
rz(1.6628671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80321035) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(-2.4513643) q[2];
rz(-0.35987443) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(-2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2279219) q[0];
sx q[0];
rz(-1.0424732) q[0];
sx q[0];
rz(-0.87164718) q[0];
rz(0.40920416) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-0.31235487) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.271628) q[0];
sx q[0];
rz(-1.4532491) q[0];
sx q[0];
rz(2.6948476) q[0];
rz(1.6705475) q[2];
sx q[2];
rz(-1.776223) q[2];
sx q[2];
rz(-0.86459407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.732199) q[1];
sx q[1];
rz(-2.7229561) q[1];
sx q[1];
rz(-0.51987176) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.254989) q[3];
sx q[3];
rz(-2.4772128) q[3];
sx q[3];
rz(-0.57306108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.006134) q[2];
sx q[2];
rz(-1.6929071) q[2];
sx q[2];
rz(-0.82497605) q[2];
rz(1.7168335) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(-0.43535522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(0.87442526) q[0];
rz(0.32886109) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(1.0985451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1060396) q[0];
sx q[0];
rz(-1.4604202) q[0];
sx q[0];
rz(-2.7644185) q[0];
rz(-pi) q[1];
rz(1.5308762) q[2];
sx q[2];
rz(-2.4062059) q[2];
sx q[2];
rz(-3.0288754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9096127) q[1];
sx q[1];
rz(-1.8706338) q[1];
sx q[1];
rz(-1.058387) q[1];
rz(-pi) q[2];
rz(-0.27733488) q[3];
sx q[3];
rz(-0.99858741) q[3];
sx q[3];
rz(-0.25857833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.141779) q[2];
sx q[2];
rz(-1.1835316) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(-1.6541803) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(0.28356799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312626) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(-1.3856101) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(-1.7562235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7840609) q[0];
sx q[0];
rz(-2.1936532) q[0];
sx q[0];
rz(2.8297367) q[0];
rz(1.6567066) q[2];
sx q[2];
rz(-1.0738147) q[2];
sx q[2];
rz(2.6045406) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5934903) q[1];
sx q[1];
rz(-1.8917276) q[1];
sx q[1];
rz(1.7744508) q[1];
x q[2];
rz(2.9117111) q[3];
sx q[3];
rz(-0.15860367) q[3];
sx q[3];
rz(1.6386709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.063227) q[2];
sx q[2];
rz(-0.93457064) q[2];
sx q[2];
rz(0.37609491) q[2];
rz(1.1345351) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-2.334972) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1217839) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(-0.10251775) q[0];
rz(2.5566697) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.9580511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4393504) q[0];
sx q[0];
rz(-2.9493414) q[0];
sx q[0];
rz(0.73701145) q[0];
rz(1.905422) q[2];
sx q[2];
rz(-2.2058704) q[2];
sx q[2];
rz(2.2872567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.592748) q[1];
sx q[1];
rz(-1.3164489) q[1];
sx q[1];
rz(-0.91491048) q[1];
x q[2];
rz(2.6747796) q[3];
sx q[3];
rz(-1.8100421) q[3];
sx q[3];
rz(-1.5249263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3465603) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(2.9417876) q[2];
rz(-1.0605158) q[3];
sx q[3];
rz(-2.6001055) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26315966) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(0.31295452) q[0];
rz(0.9494268) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.688028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1412555) q[0];
sx q[0];
rz(-0.20788684) q[0];
sx q[0];
rz(2.6882386) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4916999) q[2];
sx q[2];
rz(-0.95273113) q[2];
sx q[2];
rz(1.5072418) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9447228) q[1];
sx q[1];
rz(-2.2579282) q[1];
sx q[1];
rz(1.8623167) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62409921) q[3];
sx q[3];
rz(-2.5988376) q[3];
sx q[3];
rz(-1.8763157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7877385) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(-1.806951) q[2];
rz(-1.8169962) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(-1.2142396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1052833) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(2.4435254) q[0];
rz(-2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(0.36044136) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3981913) q[0];
sx q[0];
rz(-2.0071021) q[0];
sx q[0];
rz(-2.4453336) q[0];
rz(-pi) q[1];
rz(1.6351885) q[2];
sx q[2];
rz(-1.9528104) q[2];
sx q[2];
rz(-2.7329993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4354168) q[1];
sx q[1];
rz(-0.98505322) q[1];
sx q[1];
rz(-2.6800568) q[1];
x q[2];
rz(0.27020755) q[3];
sx q[3];
rz(-2.5761309) q[3];
sx q[3];
rz(2.5266441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91161072) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(2.3243375) q[2];
rz(-0.83003712) q[3];
sx q[3];
rz(-2.6006112) q[3];
sx q[3];
rz(-0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63477409) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-3.1095374) q[0];
rz(-1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(-2.4874036) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830156) q[0];
sx q[0];
rz(-2.9111283) q[0];
sx q[0];
rz(-2.5925267) q[0];
rz(2.6903321) q[2];
sx q[2];
rz(-0.84091821) q[2];
sx q[2];
rz(0.40715363) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5947551) q[1];
sx q[1];
rz(-1.0887301) q[1];
sx q[1];
rz(-1.0422816) q[1];
rz(-pi) q[2];
rz(-0.93864949) q[3];
sx q[3];
rz(-1.8018) q[3];
sx q[3];
rz(1.3579901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1303945) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(1.0065669) q[2];
rz(-2.448163) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(-1.2815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4257767) q[0];
sx q[0];
rz(-0.67831138) q[0];
sx q[0];
rz(-0.074445733) q[0];
rz(0.88059942) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(-0.975561) q[2];
sx q[2];
rz(-2.3481993) q[2];
sx q[2];
rz(1.3994155) q[2];
rz(-2.6941723) q[3];
sx q[3];
rz(-0.74759737) q[3];
sx q[3];
rz(2.916593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
