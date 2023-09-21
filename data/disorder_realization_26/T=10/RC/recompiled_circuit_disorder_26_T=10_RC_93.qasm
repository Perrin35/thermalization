OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(4.8298782) q[0];
sx q[0];
rz(9.7363135) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5200978) q[0];
sx q[0];
rz(-1.9858452) q[0];
sx q[0];
rz(0.15226224) q[0];
rz(-pi) q[1];
rz(0.55732255) q[2];
sx q[2];
rz(-1.5601336) q[2];
sx q[2];
rz(-0.23920857) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6701339) q[1];
sx q[1];
rz(-1.0308627) q[1];
sx q[1];
rz(2.259841) q[1];
rz(-pi) q[2];
rz(0.55922237) q[3];
sx q[3];
rz(-0.50563522) q[3];
sx q[3];
rz(1.3666183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(2.5048845) q[2];
rz(-0.84896815) q[3];
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
x q[2];
rz(-pi/2) q[3];
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
rz(0.37671509) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(0.15287457) q[0];
rz(2.3846467) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(-0.98639948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7727535) q[0];
sx q[0];
rz(-1.5305688) q[0];
sx q[0];
rz(0.059779151) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54665357) q[2];
sx q[2];
rz(-0.84178998) q[2];
sx q[2];
rz(0.58568776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6978554) q[1];
sx q[1];
rz(-1.9471696) q[1];
sx q[1];
rz(-1.8259551) q[1];
x q[2];
rz(-0.13568474) q[3];
sx q[3];
rz(-1.8968582) q[3];
sx q[3];
rz(0.56860926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5793005) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(-0.78312773) q[2];
rz(-3.1230208) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0884393) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(-2.1799178) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(-3.0128984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0823682) q[0];
sx q[0];
rz(-1.6186065) q[0];
sx q[0];
rz(-1.6533018) q[0];
x q[1];
rz(0.80758904) q[2];
sx q[2];
rz(-1.979504) q[2];
sx q[2];
rz(1.3568211) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3330831) q[1];
sx q[1];
rz(-1.8500449) q[1];
sx q[1];
rz(0.46344325) q[1];
rz(-pi) q[2];
rz(1.4314753) q[3];
sx q[3];
rz(-2.3686667) q[3];
sx q[3];
rz(-1.3611925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.06015691) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(1.241768) q[2];
rz(-2.5545819) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.63242763) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(2.0571016) q[0];
rz(-1.4831316) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(-3.049057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10034427) q[0];
sx q[0];
rz(-2.7462602) q[0];
sx q[0];
rz(2.3325217) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1763457) q[2];
sx q[2];
rz(-2.5778228) q[2];
sx q[2];
rz(3.0639067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1558518) q[1];
sx q[1];
rz(-0.34793138) q[1];
sx q[1];
rz(-0.17133979) q[1];
rz(1.7792286) q[3];
sx q[3];
rz(-2.4760893) q[3];
sx q[3];
rz(-3.0968551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83355054) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(-0.33205024) q[2];
rz(-1.0559233) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-0.21155587) q[0];
rz(1.3062723) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(2.4938915) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081213148) q[0];
sx q[0];
rz(-0.19506422) q[0];
sx q[0];
rz(-2.5293406) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7933502) q[2];
sx q[2];
rz(-1.4255382) q[2];
sx q[2];
rz(-2.9350304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60407818) q[1];
sx q[1];
rz(-1.068183) q[1];
sx q[1];
rz(0.22488774) q[1];
rz(-1.5557489) q[3];
sx q[3];
rz(-1.5516557) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5806879) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(2.4482751) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.90907) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-2.0697846) q[1];
sx q[1];
rz(-0.17428621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5151383) q[0];
sx q[0];
rz(-2.2720552) q[0];
sx q[0];
rz(-0.40909543) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2767378) q[2];
sx q[2];
rz(-2.1055429) q[2];
sx q[2];
rz(0.45644444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41387687) q[1];
sx q[1];
rz(-1.1912279) q[1];
sx q[1];
rz(-0.80645251) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2586081) q[3];
sx q[3];
rz(-1.4600666) q[3];
sx q[3];
rz(2.1544416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(-0.19763395) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(-1.4060085) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063868) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(1.6360412) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2265598) q[0];
sx q[0];
rz(-1.0462927) q[0];
sx q[0];
rz(-2.4302308) q[0];
rz(1.534329) q[2];
sx q[2];
rz(-1.140653) q[2];
sx q[2];
rz(1.9311116) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4475496) q[1];
sx q[1];
rz(-0.72206891) q[1];
sx q[1];
rz(-0.85192792) q[1];
rz(-pi) q[2];
rz(-2.8553477) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(-3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.85764) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(0.097578438) q[2];
rz(1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7512648) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(-0.51399291) q[0];
rz(0.12318525) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(0.93200144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65864627) q[0];
sx q[0];
rz(-0.46447771) q[0];
sx q[0];
rz(-1.8770201) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2134238) q[2];
sx q[2];
rz(-0.66713152) q[2];
sx q[2];
rz(0.92598976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2145558) q[1];
sx q[1];
rz(-2.3408457) q[1];
sx q[1];
rz(-2.3826249) q[1];
rz(-pi) q[2];
rz(1.7979513) q[3];
sx q[3];
rz(-0.72941581) q[3];
sx q[3];
rz(1.2136572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1121858) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(2.712148) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(-0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(1.3388348) q[0];
rz(-0.70612899) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-2.1910117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50981748) q[0];
sx q[0];
rz(-1.063856) q[0];
sx q[0];
rz(-2.1419873) q[0];
rz(1.739903) q[2];
sx q[2];
rz(-0.42113129) q[2];
sx q[2];
rz(0.087547628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7710167) q[1];
sx q[1];
rz(-0.5702714) q[1];
sx q[1];
rz(-1.9221406) q[1];
rz(-pi) q[2];
rz(2.1080984) q[3];
sx q[3];
rz(-0.47035445) q[3];
sx q[3];
rz(0.62894097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4799698) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(-0.22928672) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(-2.4226709) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15764788) q[0];
sx q[0];
rz(-2.5074208) q[0];
sx q[0];
rz(-1.7303403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11833338) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(0.10533939) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.194866) q[1];
sx q[1];
rz(-2.7416639) q[1];
sx q[1];
rz(2.7547794) q[1];
rz(2.1925681) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(2.4728647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(-0.87219277) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9185716) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(-2.7643798) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-1.2582851) q[2];
sx q[2];
rz(-1.0733114) q[2];
sx q[2];
rz(0.14807362) q[2];
rz(0.71729284) q[3];
sx q[3];
rz(-1.4581231) q[3];
sx q[3];
rz(0.067618528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];