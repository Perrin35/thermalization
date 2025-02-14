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
rz(1.9795228) q[0];
sx q[0];
rz(3.3764163) q[0];
sx q[0];
rz(11.233815) q[0];
rz(1.3136343) q[1];
sx q[1];
rz(-1.87275) q[1];
sx q[1];
rz(-2.3743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.965897) q[0];
sx q[0];
rz(-1.5724369) q[0];
sx q[0];
rz(2.7465435) q[0];
rz(-pi) q[1];
rz(1.1512773) q[2];
sx q[2];
rz(-1.1184208) q[2];
sx q[2];
rz(0.6587761) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0815107) q[1];
sx q[1];
rz(-2.7742371) q[1];
sx q[1];
rz(-1.0698331) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6758322) q[3];
sx q[3];
rz(-0.6094555) q[3];
sx q[3];
rz(-2.6026562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2904498) q[2];
sx q[2];
rz(-0.81633154) q[2];
sx q[2];
rz(1.2510703) q[2];
rz(-2.3276681) q[3];
sx q[3];
rz(-1.1000752) q[3];
sx q[3];
rz(-1.4936911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55452764) q[0];
sx q[0];
rz(-1.4545414) q[0];
sx q[0];
rz(2.5982507) q[0];
rz(-2.4045565) q[1];
sx q[1];
rz(-0.68865028) q[1];
sx q[1];
rz(0.063311689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30697) q[0];
sx q[0];
rz(-2.1181014) q[0];
sx q[0];
rz(0.91816781) q[0];
rz(-pi) q[1];
rz(2.43119) q[2];
sx q[2];
rz(-0.36493976) q[2];
sx q[2];
rz(-2.0055564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1406537) q[1];
sx q[1];
rz(-1.88511) q[1];
sx q[1];
rz(-0.94625603) q[1];
rz(-pi) q[2];
rz(1.3593986) q[3];
sx q[3];
rz(-1.711004) q[3];
sx q[3];
rz(-1.141215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43725499) q[2];
sx q[2];
rz(-2.0180118) q[2];
sx q[2];
rz(0.91427461) q[2];
rz(-2.8546913) q[3];
sx q[3];
rz(-0.56921452) q[3];
sx q[3];
rz(0.60681075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27980003) q[0];
sx q[0];
rz(-1.2759811) q[0];
sx q[0];
rz(1.3943425) q[0];
rz(-2.2718248) q[1];
sx q[1];
rz(-0.52647796) q[1];
sx q[1];
rz(-0.24551749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.190827) q[0];
sx q[0];
rz(-0.47266911) q[0];
sx q[0];
rz(-0.60421519) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1046529) q[2];
sx q[2];
rz(-2.1203342) q[2];
sx q[2];
rz(0.20333268) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88855381) q[1];
sx q[1];
rz(-1.2978857) q[1];
sx q[1];
rz(2.7103579) q[1];
rz(-pi) q[2];
rz(0.80733733) q[3];
sx q[3];
rz(-1.5745947) q[3];
sx q[3];
rz(-1.2646874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3123582) q[2];
sx q[2];
rz(-1.1623397) q[2];
sx q[2];
rz(-1.8355628) q[2];
rz(-2.8570789) q[3];
sx q[3];
rz(-1.6407216) q[3];
sx q[3];
rz(1.3792926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156085) q[0];
sx q[0];
rz(-0.59399501) q[0];
sx q[0];
rz(2.3748412) q[0];
rz(-1.7150257) q[1];
sx q[1];
rz(-1.7568935) q[1];
sx q[1];
rz(-1.0624622) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7257197) q[0];
sx q[0];
rz(-1.5860646) q[0];
sx q[0];
rz(-1.589561) q[0];
rz(1.3704186) q[2];
sx q[2];
rz(-0.60094072) q[2];
sx q[2];
rz(-2.9487425) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.26902449) q[1];
sx q[1];
rz(-1.8360129) q[1];
sx q[1];
rz(-2.0716297) q[1];
rz(-pi) q[2];
rz(-2.7492529) q[3];
sx q[3];
rz(-0.15695394) q[3];
sx q[3];
rz(-0.63877393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1876424) q[2];
sx q[2];
rz(-0.01454777) q[2];
sx q[2];
rz(-3.0237831) q[2];
rz(2.5976962) q[3];
sx q[3];
rz(-2.1102648) q[3];
sx q[3];
rz(-0.73345524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2619005) q[0];
sx q[0];
rz(-1.6629135) q[0];
sx q[0];
rz(-2.6035768) q[0];
rz(-0.063848786) q[1];
sx q[1];
rz(-1.0212746) q[1];
sx q[1];
rz(-1.341238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8842476) q[0];
sx q[0];
rz(-1.0228511) q[0];
sx q[0];
rz(2.1764285) q[0];
x q[1];
rz(-0.55300216) q[2];
sx q[2];
rz(-1.2931839) q[2];
sx q[2];
rz(1.0652519) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76829925) q[1];
sx q[1];
rz(-1.6316669) q[1];
sx q[1];
rz(2.5817691) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4405386) q[3];
sx q[3];
rz(-1.1410115) q[3];
sx q[3];
rz(-0.48359475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6824048) q[2];
sx q[2];
rz(-2.8509199) q[2];
sx q[2];
rz(0.9064557) q[2];
rz(1.9393548) q[3];
sx q[3];
rz(-1.5463444) q[3];
sx q[3];
rz(-1.3902339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9850013) q[0];
sx q[0];
rz(-2.5522794) q[0];
sx q[0];
rz(-0.44664788) q[0];
rz(-0.43529549) q[1];
sx q[1];
rz(-0.63658249) q[1];
sx q[1];
rz(0.84904233) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34693107) q[0];
sx q[0];
rz(-1.4022572) q[0];
sx q[0];
rz(-3.0656673) q[0];
rz(-pi) q[1];
rz(1.4150494) q[2];
sx q[2];
rz(-2.097599) q[2];
sx q[2];
rz(2.0447363) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0739797) q[1];
sx q[1];
rz(-0.22691209) q[1];
sx q[1];
rz(3.0935442) q[1];
rz(-1.3983512) q[3];
sx q[3];
rz(-1.9151805) q[3];
sx q[3];
rz(-1.9454959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0173831) q[2];
sx q[2];
rz(-0.97627348) q[2];
sx q[2];
rz(-1.3503831) q[2];
rz(-0.24556686) q[3];
sx q[3];
rz(-2.1731302) q[3];
sx q[3];
rz(2.811331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6733112) q[0];
sx q[0];
rz(-2.249233) q[0];
sx q[0];
rz(0.44912502) q[0];
rz(-1.4324073) q[1];
sx q[1];
rz(-2.1604373) q[1];
sx q[1];
rz(-1.2264576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0104727) q[0];
sx q[0];
rz(-0.95632271) q[0];
sx q[0];
rz(3.0486186) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6388387) q[2];
sx q[2];
rz(-1.9079676) q[2];
sx q[2];
rz(1.7751116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82348541) q[1];
sx q[1];
rz(-0.13581443) q[1];
sx q[1];
rz(-2.4450355) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38568078) q[3];
sx q[3];
rz(-1.2880518) q[3];
sx q[3];
rz(-1.677142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2655098) q[2];
sx q[2];
rz(-0.86388695) q[2];
sx q[2];
rz(0.49501219) q[2];
rz(1.0692976) q[3];
sx q[3];
rz(-2.4102231) q[3];
sx q[3];
rz(0.67720145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7345562) q[0];
sx q[0];
rz(-2.7249536) q[0];
sx q[0];
rz(-2.6805342) q[0];
rz(3.0126493) q[1];
sx q[1];
rz(-0.72595969) q[1];
sx q[1];
rz(2.2393548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0490992) q[0];
sx q[0];
rz(-1.5295795) q[0];
sx q[0];
rz(-1.5638086) q[0];
x q[1];
rz(-0.058470825) q[2];
sx q[2];
rz(-1.2035511) q[2];
sx q[2];
rz(0.98352393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.167978) q[1];
sx q[1];
rz(-1.2149286) q[1];
sx q[1];
rz(-2.0116429) q[1];
rz(-pi) q[2];
rz(-1.447313) q[3];
sx q[3];
rz(-2.4816645) q[3];
sx q[3];
rz(-2.8483158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43877131) q[2];
sx q[2];
rz(-0.92239014) q[2];
sx q[2];
rz(1.2141466) q[2];
rz(2.7625648) q[3];
sx q[3];
rz(-1.6774079) q[3];
sx q[3];
rz(-1.3488784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2269065) q[0];
sx q[0];
rz(-1.5808957) q[0];
sx q[0];
rz(-3.1353986) q[0];
rz(-0.50503039) q[1];
sx q[1];
rz(-2.6513702) q[1];
sx q[1];
rz(-0.44713155) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3792586) q[0];
sx q[0];
rz(-2.0885944) q[0];
sx q[0];
rz(0.8014265) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6880264) q[2];
sx q[2];
rz(-2.0868868) q[2];
sx q[2];
rz(2.1150401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8033235) q[1];
sx q[1];
rz(-1.4191966) q[1];
sx q[1];
rz(3.1123509) q[1];
x q[2];
rz(-0.38180967) q[3];
sx q[3];
rz(-0.47857943) q[3];
sx q[3];
rz(3.0604387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3588244) q[2];
sx q[2];
rz(-2.0147822) q[2];
sx q[2];
rz(0.44754851) q[2];
rz(-2.0284082) q[3];
sx q[3];
rz(-2.1035078) q[3];
sx q[3];
rz(0.32570496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96414763) q[0];
sx q[0];
rz(-0.064346813) q[0];
sx q[0];
rz(2.294975) q[0];
rz(2.8586491) q[1];
sx q[1];
rz(-1.6694262) q[1];
sx q[1];
rz(2.4683594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9603014) q[0];
sx q[0];
rz(-0.01894572) q[0];
sx q[0];
rz(0.51769729) q[0];
rz(-2.368957) q[2];
sx q[2];
rz(-1.4062721) q[2];
sx q[2];
rz(-2.5009837) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9615508) q[1];
sx q[1];
rz(-0.93856877) q[1];
sx q[1];
rz(-0.2747196) q[1];
rz(-0.85906927) q[3];
sx q[3];
rz(-0.43032703) q[3];
sx q[3];
rz(-0.31498614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8474569) q[2];
sx q[2];
rz(-2.2723276) q[2];
sx q[2];
rz(0.67226234) q[2];
rz(1.8654035) q[3];
sx q[3];
rz(-1.0381235) q[3];
sx q[3];
rz(-0.41260317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857674) q[0];
sx q[0];
rz(-2.0073931) q[0];
sx q[0];
rz(1.8452992) q[0];
rz(3.1271707) q[1];
sx q[1];
rz(-1.2851234) q[1];
sx q[1];
rz(-1.8869225) q[1];
rz(0.6966656) q[2];
sx q[2];
rz(-1.2045384) q[2];
sx q[2];
rz(-0.62919553) q[2];
rz(1.779535) q[3];
sx q[3];
rz(-1.345428) q[3];
sx q[3];
rz(-0.71161436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
