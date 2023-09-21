OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(2.8300571) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(-2.5826366) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5200978) q[0];
sx q[0];
rz(-1.9858452) q[0];
sx q[0];
rz(0.15226224) q[0];
x q[1];
rz(-1.5833601) q[2];
sx q[2];
rz(-1.0135092) q[2];
sx q[2];
rz(-1.81665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6701339) q[1];
sx q[1];
rz(-1.0308627) q[1];
sx q[1];
rz(2.259841) q[1];
x q[2];
rz(2.5823703) q[3];
sx q[3];
rz(-2.6359574) q[3];
sx q[3];
rz(-1.7749744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(2.2926245) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-2.9887181) q[0];
rz(2.3846467) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(2.1551932) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2043641) q[0];
sx q[0];
rz(-1.6305271) q[0];
sx q[0];
rz(-1.5304969) q[0];
rz(-pi) q[1];
rz(1.0436922) q[2];
sx q[2];
rz(-2.2615848) q[2];
sx q[2];
rz(-1.8156798) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4437372) q[1];
sx q[1];
rz(-1.194423) q[1];
sx q[1];
rz(-1.8259551) q[1];
x q[2];
rz(-1.2419224) q[3];
sx q[3];
rz(-1.4423014) q[3];
sx q[3];
rz(-2.0957029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5622921) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(0.78312773) q[2];
rz(3.1230208) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(-2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(-0.96167481) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(3.0128984) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6260687) q[0];
sx q[0];
rz(-1.4883853) q[0];
sx q[0];
rz(3.0936196) q[0];
rz(-pi) q[1];
rz(1.0110858) q[2];
sx q[2];
rz(-0.84583827) q[2];
sx q[2];
rz(-0.17979187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5161908) q[1];
sx q[1];
rz(-2.0149724) q[1];
sx q[1];
rz(1.8810012) q[1];
x q[2];
rz(0.80273654) q[3];
sx q[3];
rz(-1.4736796) q[3];
sx q[3];
rz(3.0320398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.06015691) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.8998247) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-0.97755066) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63242763) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-1.084491) q[0];
rz(1.658461) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(0.09253563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1922069) q[0];
sx q[0];
rz(-1.3017676) q[0];
sx q[0];
rz(1.2775248) q[0];
x q[1];
rz(-2.0501577) q[2];
sx q[2];
rz(-1.2617246) q[2];
sx q[2];
rz(2.022559) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42379984) q[1];
sx q[1];
rz(-1.6289627) q[1];
sx q[1];
rz(0.34323005) q[1];
rz(0.16102287) q[3];
sx q[3];
rz(-2.2194214) q[3];
sx q[3];
rz(2.8341858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83355054) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(-2.8095424) q[2];
rz(1.0559233) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(2.5312996) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8191391) q[0];
sx q[0];
rz(-1.8503014) q[0];
sx q[0];
rz(0.21155587) q[0];
rz(-1.3062723) q[1];
sx q[1];
rz(-1.897656) q[1];
sx q[1];
rz(-0.64770118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0928597) q[0];
sx q[0];
rz(-1.6824241) q[0];
sx q[0];
rz(0.1603006) q[0];
x q[1];
rz(-1.4164045) q[2];
sx q[2];
rz(-1.2263745) q[2];
sx q[2];
rz(1.4167348) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60407818) q[1];
sx q[1];
rz(-1.068183) q[1];
sx q[1];
rz(2.9167049) q[1];
rz(-1.5557489) q[3];
sx q[3];
rz(-1.5516557) q[3];
sx q[3];
rz(2.0862938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56090474) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(-0.69331759) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(3.1366689) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-2.9142002) q[0];
sx q[0];
rz(-1.90907) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(0.17428621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586853) q[0];
sx q[0];
rz(-1.8795965) q[0];
sx q[0];
rz(-0.82682825) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3107489) q[2];
sx q[2];
rz(-0.85692642) q[2];
sx q[2];
rz(-1.4884399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7277158) q[1];
sx q[1];
rz(-1.9503647) q[1];
sx q[1];
rz(0.80645251) q[1];
rz(-pi) q[2];
x q[2];
rz(1.918119) q[3];
sx q[3];
rz(-0.33063774) q[3];
sx q[3];
rz(0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(-1.4060085) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.063868) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(0.96616191) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.5055515) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063232139) q[0];
sx q[0];
rz(-0.97023836) q[0];
sx q[0];
rz(-2.2230704) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0622919) q[2];
sx q[2];
rz(-2.7100025) q[2];
sx q[2];
rz(2.0183795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70430763) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(-2.1561161) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3397293) q[3];
sx q[3];
rz(-2.8390084) q[3];
sx q[3];
rz(-1.8031977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.85764) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(-1.3939259) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.39032787) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(-0.51399291) q[0];
rz(3.0184074) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(-2.2095912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31873736) q[0];
sx q[0];
rz(-1.1294951) q[0];
sx q[0];
rz(-2.9916828) q[0];
rz(-pi) q[1];
rz(-0.93512647) q[2];
sx q[2];
rz(-1.3526275) q[2];
sx q[2];
rz(-0.93014923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.060170505) q[1];
sx q[1];
rz(-2.0875071) q[1];
sx q[1];
rz(2.4992649) q[1];
rz(-1.3436414) q[3];
sx q[3];
rz(-2.4121768) q[3];
sx q[3];
rz(1.9279355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1121858) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(2.712148) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(-0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747303) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(-1.8027579) q[0];
rz(-2.4354637) q[1];
sx q[1];
rz(-1.8920205) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6317752) q[0];
sx q[0];
rz(-1.063856) q[0];
sx q[0];
rz(-2.1419873) q[0];
rz(-pi) q[1];
x q[1];
rz(0.075245113) q[2];
sx q[2];
rz(-1.1560455) q[2];
sx q[2];
rz(-3.0441949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49950019) q[1];
sx q[1];
rz(-1.7576808) q[1];
sx q[1];
rz(2.1128113) q[1];
rz(-2.8870228) q[3];
sx q[3];
rz(-1.1708784) q[3];
sx q[3];
rz(3.1018156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(0.48428145) q[2];
rz(-2.2144923) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(-2.5203729) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973688) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.9186585) q[1];
sx q[1];
rz(-0.71892175) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842429) q[0];
sx q[0];
rz(-1.6650668) q[0];
sx q[0];
rz(2.1988792) q[0];
x q[1];
rz(1.747379) q[2];
sx q[2];
rz(-2.546616) q[2];
sx q[2];
rz(2.8240311) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5303516) q[1];
sx q[1];
rz(-1.9396922) q[1];
sx q[1];
rz(1.7289274) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94902456) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(-2.4728647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.2021844) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(-2.2693999) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9185716) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(-0.37721286) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(-1.8833075) q[2];
sx q[2];
rz(-2.0682813) q[2];
sx q[2];
rz(-2.993519) q[2];
rz(-2.4242998) q[3];
sx q[3];
rz(-1.4581231) q[3];
sx q[3];
rz(0.067618528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];