OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.10211927) q[0];
sx q[0];
rz(4.6146225) q[0];
sx q[0];
rz(6.1518402) q[0];
rz(-3.1047473) q[1];
sx q[1];
rz(-2.614202) q[1];
sx q[1];
rz(-1.8324469) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2764171) q[0];
sx q[0];
rz(-2.2593479) q[0];
sx q[0];
rz(-0.49746969) q[0];
rz(-pi) q[1];
rz(-0.45707656) q[2];
sx q[2];
rz(-2.0684048) q[2];
sx q[2];
rz(-0.08809419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6908474) q[1];
sx q[1];
rz(-1.2589129) q[1];
sx q[1];
rz(1.8212832) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74983024) q[3];
sx q[3];
rz(-2.8932126) q[3];
sx q[3];
rz(1.5313784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(2.479539) q[2];
rz(-2.3946297) q[3];
sx q[3];
rz(-0.91619879) q[3];
sx q[3];
rz(-1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5544283) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(1.1821049) q[0];
rz(-2.3960522) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(3.0381957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7000843) q[0];
sx q[0];
rz(-1.147195) q[0];
sx q[0];
rz(1.3157428) q[0];
x q[1];
rz(2.1076043) q[2];
sx q[2];
rz(-2.4721382) q[2];
sx q[2];
rz(0.34399271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89241806) q[1];
sx q[1];
rz(-1.9016445) q[1];
sx q[1];
rz(-2.3187917) q[1];
x q[2];
rz(-2.313226) q[3];
sx q[3];
rz(-1.6992555) q[3];
sx q[3];
rz(1.0297694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55249247) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(-0.44431552) q[2];
rz(1.0537423) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(-1.9452555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0788197) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(-2.479082) q[0];
rz(1.6569116) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(-2.9248765) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22610006) q[0];
sx q[0];
rz(-2.2558172) q[0];
sx q[0];
rz(-0.95695509) q[0];
rz(-1.1323053) q[2];
sx q[2];
rz(-2.2542033) q[2];
sx q[2];
rz(2.8029798) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3308476) q[1];
sx q[1];
rz(-2.243089) q[1];
sx q[1];
rz(0.051203392) q[1];
rz(1.7314579) q[3];
sx q[3];
rz(-1.0815797) q[3];
sx q[3];
rz(-0.9582203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4776939) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(0.23507512) q[2];
rz(2.7663686) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5594056) q[0];
sx q[0];
rz(-3.0538054) q[0];
sx q[0];
rz(-1.0002332) q[0];
rz(-1.2031215) q[1];
sx q[1];
rz(-1.6327881) q[1];
sx q[1];
rz(0.54471725) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31123589) q[0];
sx q[0];
rz(-2.2831627) q[0];
sx q[0];
rz(1.5115378) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80766957) q[2];
sx q[2];
rz(-2.2414506) q[2];
sx q[2];
rz(-2.9602221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9904321) q[1];
sx q[1];
rz(-2.9176788) q[1];
sx q[1];
rz(-1.0099221) q[1];
rz(-pi) q[2];
rz(1.9849378) q[3];
sx q[3];
rz(-2.4865287) q[3];
sx q[3];
rz(2.6998346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7228221) q[2];
sx q[2];
rz(-2.3493769) q[2];
sx q[2];
rz(2.344632) q[2];
rz(-0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(0.42118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5274984) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(-1.2689137) q[0];
rz(2.1753963) q[1];
sx q[1];
rz(-1.393001) q[1];
sx q[1];
rz(2.6752245) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.120935) q[0];
sx q[0];
rz(-1.0758721) q[0];
sx q[0];
rz(1.0878956) q[0];
rz(-1.9210468) q[2];
sx q[2];
rz(-1.1333457) q[2];
sx q[2];
rz(-1.3910363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9590145) q[1];
sx q[1];
rz(-1.2532638) q[1];
sx q[1];
rz(0.93341254) q[1];
x q[2];
rz(-2.0556695) q[3];
sx q[3];
rz(-2.4182662) q[3];
sx q[3];
rz(1.5736027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7352778) q[2];
sx q[2];
rz(-0.80253989) q[2];
sx q[2];
rz(-2.3373513) q[2];
rz(-0.4869701) q[3];
sx q[3];
rz(-1.0420957) q[3];
sx q[3];
rz(0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93953472) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(-0.95788389) q[0];
rz(1.5127888) q[1];
sx q[1];
rz(-1.0572546) q[1];
sx q[1];
rz(1.5527976) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68715019) q[0];
sx q[0];
rz(-1.5479364) q[0];
sx q[0];
rz(1.6506565) q[0];
rz(-pi) q[1];
rz(0.79926305) q[2];
sx q[2];
rz(-2.7529319) q[2];
sx q[2];
rz(0.51539076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1576288) q[1];
sx q[1];
rz(-2.2215054) q[1];
sx q[1];
rz(-1.7939066) q[1];
x q[2];
rz(-2.5143322) q[3];
sx q[3];
rz(-1.903844) q[3];
sx q[3];
rz(1.8740538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2019041) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(-0.73516694) q[2];
rz(-0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(-0.86400509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.76148024) q[0];
sx q[0];
rz(-2.1367456) q[0];
sx q[0];
rz(-2.5349706) q[0];
rz(0.78423777) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(-1.3444791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494251) q[0];
sx q[0];
rz(-1.0972404) q[0];
sx q[0];
rz(0.72832941) q[0];
rz(-pi) q[1];
rz(1.3188348) q[2];
sx q[2];
rz(-1.0875541) q[2];
sx q[2];
rz(1.9117219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5658907) q[1];
sx q[1];
rz(-1.1739578) q[1];
sx q[1];
rz(1.8116519) q[1];
rz(-1.3913888) q[3];
sx q[3];
rz(-0.47166892) q[3];
sx q[3];
rz(0.52667945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6091696) q[2];
sx q[2];
rz(-2.6340941) q[2];
sx q[2];
rz(1.3896821) q[2];
rz(2.9621647) q[3];
sx q[3];
rz(-2.1556985) q[3];
sx q[3];
rz(-0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0591902) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(0.75505906) q[0];
rz(-0.45626196) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(2.0645352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2307325) q[0];
sx q[0];
rz(-0.60711475) q[0];
sx q[0];
rz(-0.40672983) q[0];
rz(0.78387129) q[2];
sx q[2];
rz(-1.6567363) q[2];
sx q[2];
rz(-0.29168561) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3216599) q[1];
sx q[1];
rz(-2.0256307) q[1];
sx q[1];
rz(1.2534035) q[1];
rz(-pi) q[2];
rz(0.15865008) q[3];
sx q[3];
rz(-2.1057743) q[3];
sx q[3];
rz(-1.3033997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-1.1001526) q[2];
sx q[2];
rz(2.4033578) q[2];
rz(1.632656) q[3];
sx q[3];
rz(-1.8176707) q[3];
sx q[3];
rz(-1.9807321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410925) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(-0.68967462) q[1];
sx q[1];
rz(-2.4318305) q[1];
sx q[1];
rz(0.41742691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8807424) q[0];
sx q[0];
rz(-2.1584903) q[0];
sx q[0];
rz(-1.2686612) q[0];
rz(-0.29977389) q[2];
sx q[2];
rz(-1.6437141) q[2];
sx q[2];
rz(-2.656183) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3950353) q[1];
sx q[1];
rz(-0.42320028) q[1];
sx q[1];
rz(2.8854135) q[1];
rz(-1.1754009) q[3];
sx q[3];
rz(-1.1987276) q[3];
sx q[3];
rz(1.3517861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4735585) q[2];
sx q[2];
rz(-1.7201951) q[2];
sx q[2];
rz(0.045698015) q[2];
rz(0.33347305) q[3];
sx q[3];
rz(-0.57531753) q[3];
sx q[3];
rz(-3.0157109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
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
rz(-2.4561653) q[0];
sx q[0];
rz(-0.5394772) q[0];
sx q[0];
rz(1.9239377) q[0];
rz(2.894891) q[1];
sx q[1];
rz(-1.291899) q[1];
sx q[1];
rz(-2.6620679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8586774) q[0];
sx q[0];
rz(-1.1395795) q[0];
sx q[0];
rz(0.89003508) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6423655) q[2];
sx q[2];
rz(-2.2031234) q[2];
sx q[2];
rz(1.3030753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7361141) q[1];
sx q[1];
rz(-0.60156721) q[1];
sx q[1];
rz(-2.6094237) q[1];
rz(-2.6761618) q[3];
sx q[3];
rz(-1.3579391) q[3];
sx q[3];
rz(1.928291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55629998) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(-1.7362107) q[2];
rz(-2.4893238) q[3];
sx q[3];
rz(-2.8758719) q[3];
sx q[3];
rz(0.71776596) q[3];
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
rz(-pi/2) q[0];
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
rz(1.2300867) q[0];
sx q[0];
rz(-1.390504) q[0];
sx q[0];
rz(2.6934296) q[0];
rz(-0.74408342) q[1];
sx q[1];
rz(-0.71744812) q[1];
sx q[1];
rz(1.7070028) q[1];
rz(0.23575959) q[2];
sx q[2];
rz(-2.5452062) q[2];
sx q[2];
rz(1.7231981) q[2];
rz(-2.7884095) q[3];
sx q[3];
rz(-1.0241057) q[3];
sx q[3];
rz(3.1288341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
