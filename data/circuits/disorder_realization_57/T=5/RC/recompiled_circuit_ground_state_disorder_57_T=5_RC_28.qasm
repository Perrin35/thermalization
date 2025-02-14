OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(-2.804629) q[0];
rz(-2.872074) q[1];
sx q[1];
rz(-0.94269204) q[1];
sx q[1];
rz(1.8820794) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2452606) q[0];
sx q[0];
rz(-1.0335644) q[0];
sx q[0];
rz(2.2396829) q[0];
rz(-0.23064166) q[2];
sx q[2];
rz(-0.75048026) q[2];
sx q[2];
rz(-0.41507687) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6549398) q[1];
sx q[1];
rz(-1.6702685) q[1];
sx q[1];
rz(-0.54648262) q[1];
rz(-pi) q[2];
rz(0.97136949) q[3];
sx q[3];
rz(-1.2232336) q[3];
sx q[3];
rz(0.83272782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17244615) q[2];
sx q[2];
rz(-1.5643876) q[2];
sx q[2];
rz(-1.9077612) q[2];
rz(-1.5305758) q[3];
sx q[3];
rz(-1.7350035) q[3];
sx q[3];
rz(2.5764901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3259786) q[0];
sx q[0];
rz(-2.1765206) q[0];
sx q[0];
rz(-0.24398971) q[0];
rz(-1.9832393) q[1];
sx q[1];
rz(-2.127425) q[1];
sx q[1];
rz(0.44833952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2650484) q[0];
sx q[0];
rz(-0.52406132) q[0];
sx q[0];
rz(2.7830809) q[0];
x q[1];
rz(-0.82996394) q[2];
sx q[2];
rz(-1.1051851) q[2];
sx q[2];
rz(0.59620406) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.073001677) q[1];
sx q[1];
rz(-1.0319971) q[1];
sx q[1];
rz(-0.52937845) q[1];
rz(0.86037029) q[3];
sx q[3];
rz(-0.67796889) q[3];
sx q[3];
rz(1.9063661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9490732) q[2];
sx q[2];
rz(-0.039704617) q[2];
sx q[2];
rz(2.330759) q[2];
rz(-0.0090946322) q[3];
sx q[3];
rz(-1.5261212) q[3];
sx q[3];
rz(-2.8240589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4119165) q[0];
sx q[0];
rz(-1.7971973) q[0];
sx q[0];
rz(-0.94773951) q[0];
rz(2.2432227) q[1];
sx q[1];
rz(-0.71304524) q[1];
sx q[1];
rz(2.9352442) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8977933) q[0];
sx q[0];
rz(-0.81598213) q[0];
sx q[0];
rz(0.010494626) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3487665) q[2];
sx q[2];
rz(-1.2812876) q[2];
sx q[2];
rz(0.0057366554) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5134597) q[1];
sx q[1];
rz(-1.1185137) q[1];
sx q[1];
rz(-2.9123405) q[1];
rz(-pi) q[2];
rz(3.0037155) q[3];
sx q[3];
rz(-2.8123224) q[3];
sx q[3];
rz(0.38789685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4121033) q[2];
sx q[2];
rz(-2.1495843) q[2];
sx q[2];
rz(-2.1573529) q[2];
rz(-0.49736831) q[3];
sx q[3];
rz(-1.8822742) q[3];
sx q[3];
rz(0.74500144) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618149) q[0];
sx q[0];
rz(-0.091876939) q[0];
sx q[0];
rz(0.20633695) q[0];
rz(1.4641209) q[1];
sx q[1];
rz(-0.76281491) q[1];
sx q[1];
rz(-0.62613553) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38202823) q[0];
sx q[0];
rz(-0.81597486) q[0];
sx q[0];
rz(1.0797281) q[0];
rz(-pi) q[1];
rz(0.017109326) q[2];
sx q[2];
rz(-3.0387911) q[2];
sx q[2];
rz(0.67827618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4724222) q[1];
sx q[1];
rz(-2.0991994) q[1];
sx q[1];
rz(1.0086568) q[1];
x q[2];
rz(-1.3317165) q[3];
sx q[3];
rz(-2.7446973) q[3];
sx q[3];
rz(-1.3439076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76554406) q[2];
sx q[2];
rz(-2.3638793) q[2];
sx q[2];
rz(-2.1518339) q[2];
rz(2.1073714) q[3];
sx q[3];
rz(-1.1690305) q[3];
sx q[3];
rz(-2.6160252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.702221) q[0];
sx q[0];
rz(-1.419743) q[0];
sx q[0];
rz(1.1628344) q[0];
rz(0.16712664) q[1];
sx q[1];
rz(-1.5252599) q[1];
sx q[1];
rz(0.67536813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7279926) q[0];
sx q[0];
rz(-2.4420847) q[0];
sx q[0];
rz(-0.34297717) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33950342) q[2];
sx q[2];
rz(-1.9304747) q[2];
sx q[2];
rz(-2.8946257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7168419) q[1];
sx q[1];
rz(-1.2256241) q[1];
sx q[1];
rz(2.0276245) q[1];
rz(-pi) q[2];
rz(-2.1618103) q[3];
sx q[3];
rz(-1.3434935) q[3];
sx q[3];
rz(-2.1066372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9250138) q[2];
sx q[2];
rz(-1.719097) q[2];
sx q[2];
rz(0.50755802) q[2];
rz(-3.1223068) q[3];
sx q[3];
rz(-0.79500335) q[3];
sx q[3];
rz(-1.219334) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60292202) q[0];
sx q[0];
rz(-2.5141073) q[0];
sx q[0];
rz(2.4399309) q[0];
rz(-0.64274669) q[1];
sx q[1];
rz(-1.9822491) q[1];
sx q[1];
rz(-1.3210375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3468332) q[0];
sx q[0];
rz(-0.94330314) q[0];
sx q[0];
rz(-0.015608257) q[0];
rz(0.3568383) q[2];
sx q[2];
rz(-0.90274397) q[2];
sx q[2];
rz(-0.84260637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5112662) q[1];
sx q[1];
rz(-2.0980586) q[1];
sx q[1];
rz(-2.7106337) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81084588) q[3];
sx q[3];
rz(-1.7827991) q[3];
sx q[3];
rz(0.9190587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8464437) q[2];
sx q[2];
rz(-2.6602793) q[2];
sx q[2];
rz(-0.63977891) q[2];
rz(2.4844737) q[3];
sx q[3];
rz(-1.6491456) q[3];
sx q[3];
rz(-0.28321987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039728634) q[0];
sx q[0];
rz(-1.855408) q[0];
sx q[0];
rz(-2.3611327) q[0];
rz(-1.2377493) q[1];
sx q[1];
rz(-1.8433808) q[1];
sx q[1];
rz(-2.9775528) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6934768) q[0];
sx q[0];
rz(-2.0882029) q[0];
sx q[0];
rz(0.2736926) q[0];
rz(-2.3140644) q[2];
sx q[2];
rz(-2.568576) q[2];
sx q[2];
rz(0.94486754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61172813) q[1];
sx q[1];
rz(-2.0862038) q[1];
sx q[1];
rz(0.5926253) q[1];
rz(-pi) q[2];
rz(-2.3497054) q[3];
sx q[3];
rz(-1.2372176) q[3];
sx q[3];
rz(-1.6543948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8728472) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(-2.2609113) q[2];
rz(1.4745447) q[3];
sx q[3];
rz(-2.0737952) q[3];
sx q[3];
rz(-1.4145981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47377652) q[0];
sx q[0];
rz(-2.9528463) q[0];
sx q[0];
rz(1.1291946) q[0];
rz(-2.1881564) q[1];
sx q[1];
rz(-0.9402746) q[1];
sx q[1];
rz(-0.58852351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58266334) q[0];
sx q[0];
rz(-1.869794) q[0];
sx q[0];
rz(1.0010946) q[0];
rz(-pi) q[1];
rz(0.97739525) q[2];
sx q[2];
rz(-2.7468312) q[2];
sx q[2];
rz(-0.66937689) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5626763) q[1];
sx q[1];
rz(-2.2136218) q[1];
sx q[1];
rz(1.1688656) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7993808) q[3];
sx q[3];
rz(-2.9116063) q[3];
sx q[3];
rz(-0.019345779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15527655) q[2];
sx q[2];
rz(-1.0078112) q[2];
sx q[2];
rz(-1.0922095) q[2];
rz(-2.3327667) q[3];
sx q[3];
rz(-1.0816962) q[3];
sx q[3];
rz(-1.4441747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5539733) q[0];
sx q[0];
rz(-2.0646136) q[0];
sx q[0];
rz(-0.67163604) q[0];
rz(2.6486168) q[1];
sx q[1];
rz(-2.2608345) q[1];
sx q[1];
rz(1.2395476) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7412146) q[0];
sx q[0];
rz(-1.5682055) q[0];
sx q[0];
rz(1.8080312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.049767869) q[2];
sx q[2];
rz(-2.9441212) q[2];
sx q[2];
rz(1.2021499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2375657) q[1];
sx q[1];
rz(-0.75566245) q[1];
sx q[1];
rz(0.66046884) q[1];
x q[2];
rz(1.7335692) q[3];
sx q[3];
rz(-2.3331169) q[3];
sx q[3];
rz(2.2089434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12941831) q[2];
sx q[2];
rz(-1.6678026) q[2];
sx q[2];
rz(2.5950281) q[2];
rz(-2.2187345) q[3];
sx q[3];
rz(-0.62896362) q[3];
sx q[3];
rz(-1.9950689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.4839812) q[0];
sx q[0];
rz(-2.681499) q[0];
sx q[0];
rz(-0.42903236) q[0];
rz(-1.2826762) q[1];
sx q[1];
rz(-1.3653283) q[1];
sx q[1];
rz(-3.0603337) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6264296) q[0];
sx q[0];
rz(-1.7395955) q[0];
sx q[0];
rz(-1.6506881) q[0];
rz(1.5052028) q[2];
sx q[2];
rz(-2.6806405) q[2];
sx q[2];
rz(-2.2923802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9489021) q[1];
sx q[1];
rz(-1.16445) q[1];
sx q[1];
rz(2.3188616) q[1];
rz(-2.6762415) q[3];
sx q[3];
rz(-2.029226) q[3];
sx q[3];
rz(-1.1155333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.405534) q[2];
sx q[2];
rz(-2.2106876) q[2];
sx q[2];
rz(-1.1116568) q[2];
rz(1.9723655) q[3];
sx q[3];
rz(-2.4284095) q[3];
sx q[3];
rz(-0.14298239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0040759857) q[0];
sx q[0];
rz(-0.62340323) q[0];
sx q[0];
rz(2.2829983) q[0];
rz(-2.2463592) q[1];
sx q[1];
rz(-1.3139071) q[1];
sx q[1];
rz(0.039176686) q[1];
rz(1.5262114) q[2];
sx q[2];
rz(-2.8845061) q[2];
sx q[2];
rz(1.3731879) q[2];
rz(-0.64057401) q[3];
sx q[3];
rz(-1.6770029) q[3];
sx q[3];
rz(-2.3322662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
