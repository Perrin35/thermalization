OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0225585) q[0];
sx q[0];
rz(2.0329539) q[0];
sx q[0];
rz(11.995765) q[0];
rz(1.624149) q[1];
sx q[1];
rz(-0.94944209) q[1];
sx q[1];
rz(-1.5631262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0235897) q[0];
sx q[0];
rz(-2.6896416) q[0];
sx q[0];
rz(-1.8382545) q[0];
x q[1];
rz(1.203166) q[2];
sx q[2];
rz(-0.7945294) q[2];
sx q[2];
rz(2.6472435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7540463) q[1];
sx q[1];
rz(-1.3943732) q[1];
sx q[1];
rz(0.37732084) q[1];
x q[2];
rz(1.9195032) q[3];
sx q[3];
rz(-2.152386) q[3];
sx q[3];
rz(-2.3720456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40330517) q[2];
sx q[2];
rz(-2.3052576) q[2];
sx q[2];
rz(-2.9264012) q[2];
rz(-0.0067986851) q[3];
sx q[3];
rz(-0.88383043) q[3];
sx q[3];
rz(2.1782037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5408408) q[0];
sx q[0];
rz(-2.231926) q[0];
sx q[0];
rz(1.2700861) q[0];
rz(2.0974244) q[1];
sx q[1];
rz(-2.7296598) q[1];
sx q[1];
rz(-0.56234908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5836307) q[0];
sx q[0];
rz(-1.4578125) q[0];
sx q[0];
rz(1.8336748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1435881) q[2];
sx q[2];
rz(-1.0173804) q[2];
sx q[2];
rz(1.9678022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.96404379) q[1];
sx q[1];
rz(-0.73638396) q[1];
sx q[1];
rz(0.089929637) q[1];
rz(2.8329912) q[3];
sx q[3];
rz(-1.8695334) q[3];
sx q[3];
rz(1.5200391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.0079785138) q[2];
sx q[2];
rz(-2.9972711) q[2];
sx q[2];
rz(-0.14145429) q[2];
rz(-1.0072297) q[3];
sx q[3];
rz(-1.6601345) q[3];
sx q[3];
rz(3.0032515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7815642) q[0];
sx q[0];
rz(-1.5046295) q[0];
sx q[0];
rz(3.0887399) q[0];
rz(-1.7618893) q[1];
sx q[1];
rz(-2.3638937) q[1];
sx q[1];
rz(0.27390972) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5316469) q[0];
sx q[0];
rz(-1.7478806) q[0];
sx q[0];
rz(0.037399423) q[0];
rz(-pi) q[1];
rz(0.81676151) q[2];
sx q[2];
rz(-1.5026983) q[2];
sx q[2];
rz(0.20055972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1812836) q[1];
sx q[1];
rz(-0.77086721) q[1];
sx q[1];
rz(0.038357448) q[1];
rz(0.11577167) q[3];
sx q[3];
rz(-1.8325509) q[3];
sx q[3];
rz(3.0512471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.86458796) q[2];
sx q[2];
rz(-1.6463582) q[2];
sx q[2];
rz(0.60025364) q[2];
rz(-0.65435919) q[3];
sx q[3];
rz(-0.31939605) q[3];
sx q[3];
rz(1.6681685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3510503) q[0];
sx q[0];
rz(-1.6449787) q[0];
sx q[0];
rz(0.45781621) q[0];
rz(-2.5425725) q[1];
sx q[1];
rz(-0.51282239) q[1];
sx q[1];
rz(-0.27214989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7632338) q[0];
sx q[0];
rz(-0.97017479) q[0];
sx q[0];
rz(-1.8778223) q[0];
rz(-pi) q[1];
rz(0.53109269) q[2];
sx q[2];
rz(-2.0390687) q[2];
sx q[2];
rz(-2.7837769) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.052201) q[1];
sx q[1];
rz(-1.3864582) q[1];
sx q[1];
rz(1.4505397) q[1];
rz(-pi) q[2];
rz(2.5013345) q[3];
sx q[3];
rz(-1.7096412) q[3];
sx q[3];
rz(1.1306896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6279471) q[2];
sx q[2];
rz(-0.38868913) q[2];
sx q[2];
rz(1.909931) q[2];
rz(1.7197459) q[3];
sx q[3];
rz(-1.809241) q[3];
sx q[3];
rz(-2.1739912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49895367) q[0];
sx q[0];
rz(-0.068878219) q[0];
sx q[0];
rz(-2.4054085) q[0];
rz(-0.63756293) q[1];
sx q[1];
rz(-1.9393189) q[1];
sx q[1];
rz(0.41086248) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4291358) q[0];
sx q[0];
rz(-2.4596058) q[0];
sx q[0];
rz(2.3739967) q[0];
rz(-pi) q[1];
rz(-2.1396808) q[2];
sx q[2];
rz(-1.7497154) q[2];
sx q[2];
rz(-0.093947336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4102829) q[1];
sx q[1];
rz(-1.5775133) q[1];
sx q[1];
rz(-2.5120458) q[1];
rz(-pi) q[2];
rz(2.210064) q[3];
sx q[3];
rz(-2.1327104) q[3];
sx q[3];
rz(-0.84980295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32834443) q[2];
sx q[2];
rz(-1.594955) q[2];
sx q[2];
rz(-0.52097875) q[2];
rz(0.10884918) q[3];
sx q[3];
rz(-1.0131001) q[3];
sx q[3];
rz(-2.525135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2206889) q[0];
sx q[0];
rz(-1.2257129) q[0];
sx q[0];
rz(1.2589681) q[0];
rz(1.0733696) q[1];
sx q[1];
rz(-1.6683038) q[1];
sx q[1];
rz(2.4980256) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94144017) q[0];
sx q[0];
rz(-2.6905031) q[0];
sx q[0];
rz(3.1300225) q[0];
x q[1];
rz(0.37677209) q[2];
sx q[2];
rz(-1.7724282) q[2];
sx q[2];
rz(-1.0472991) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38222183) q[1];
sx q[1];
rz(-2.367451) q[1];
sx q[1];
rz(1.3355119) q[1];
x q[2];
rz(-2.4197641) q[3];
sx q[3];
rz(-2.4284644) q[3];
sx q[3];
rz(-1.4827181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7966938) q[2];
sx q[2];
rz(-2.0774136) q[2];
sx q[2];
rz(1.3382443) q[2];
rz(-1.6898588) q[3];
sx q[3];
rz(-1.9386049) q[3];
sx q[3];
rz(3.0411327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79784262) q[0];
sx q[0];
rz(-1.2748953) q[0];
sx q[0];
rz(0.95329681) q[0];
rz(-2.9816309) q[1];
sx q[1];
rz(-2.5123031) q[1];
sx q[1];
rz(-0.86872855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9729444) q[0];
sx q[0];
rz(-2.4356844) q[0];
sx q[0];
rz(1.9078518) q[0];
rz(-pi) q[1];
rz(-0.8192058) q[2];
sx q[2];
rz(-0.98759606) q[2];
sx q[2];
rz(1.6811586) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9142469) q[1];
sx q[1];
rz(-0.29623756) q[1];
sx q[1];
rz(3.0692718) q[1];
rz(-2.8817435) q[3];
sx q[3];
rz(-0.9161549) q[3];
sx q[3];
rz(2.8878867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6877044) q[2];
sx q[2];
rz(-1.1331465) q[2];
sx q[2];
rz(0.49393168) q[2];
rz(-1.6535053) q[3];
sx q[3];
rz(-0.86797124) q[3];
sx q[3];
rz(-0.021763703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9138551) q[0];
sx q[0];
rz(-1.7182925) q[0];
sx q[0];
rz(-0.29105759) q[0];
rz(3.011138) q[1];
sx q[1];
rz(-0.37965241) q[1];
sx q[1];
rz(1.5721273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5792517) q[0];
sx q[0];
rz(-1.6214193) q[0];
sx q[0];
rz(0.013112092) q[0];
x q[1];
rz(2.2443971) q[2];
sx q[2];
rz(-1.7056124) q[2];
sx q[2];
rz(-1.7133985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80940217) q[1];
sx q[1];
rz(-2.5173264) q[1];
sx q[1];
rz(2.9430976) q[1];
rz(2.7588284) q[3];
sx q[3];
rz(-2.8105002) q[3];
sx q[3];
rz(2.4352026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36720413) q[2];
sx q[2];
rz(-2.9640894) q[2];
sx q[2];
rz(-0.11167488) q[2];
rz(2.3292134) q[3];
sx q[3];
rz(-1.2923391) q[3];
sx q[3];
rz(-1.5264548) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77686253) q[0];
sx q[0];
rz(-0.22219816) q[0];
sx q[0];
rz(0.36866933) q[0];
rz(-2.5313976) q[1];
sx q[1];
rz(-1.5564324) q[1];
sx q[1];
rz(-1.0367702) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2431738) q[0];
sx q[0];
rz(-2.5372525) q[0];
sx q[0];
rz(2.739796) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5620932) q[2];
sx q[2];
rz(-0.51299757) q[2];
sx q[2];
rz(-0.41288105) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2543176) q[1];
sx q[1];
rz(-2.3000082) q[1];
sx q[1];
rz(-0.18730728) q[1];
rz(-1.0843363) q[3];
sx q[3];
rz(-1.7170625) q[3];
sx q[3];
rz(1.676596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6250352) q[2];
sx q[2];
rz(-0.80356193) q[2];
sx q[2];
rz(-0.45300031) q[2];
rz(-2.1895444) q[3];
sx q[3];
rz(-1.1166078) q[3];
sx q[3];
rz(-3.0208352) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9348738) q[0];
sx q[0];
rz(-3.0937338) q[0];
sx q[0];
rz(-0.55024838) q[0];
rz(-2.0568559) q[1];
sx q[1];
rz(-0.58894283) q[1];
sx q[1];
rz(-1.750607) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6577639) q[0];
sx q[0];
rz(-0.50301174) q[0];
sx q[0];
rz(-1.6526954) q[0];
x q[1];
rz(-2.1065305) q[2];
sx q[2];
rz(-1.7155871) q[2];
sx q[2];
rz(2.5404251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3251614) q[1];
sx q[1];
rz(-1.0712938) q[1];
sx q[1];
rz(0.71814594) q[1];
rz(-pi) q[2];
rz(1.6272561) q[3];
sx q[3];
rz(-1.5096231) q[3];
sx q[3];
rz(-0.45555761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37798205) q[2];
sx q[2];
rz(-2.4201329) q[2];
sx q[2];
rz(-1.9788431) q[2];
rz(2.5871596) q[3];
sx q[3];
rz(-0.88397637) q[3];
sx q[3];
rz(-2.673431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32600317) q[0];
sx q[0];
rz(-1.3023888) q[0];
sx q[0];
rz(1.6525255) q[0];
rz(2.9266657) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(1.3431637) q[2];
sx q[2];
rz(-0.91522436) q[2];
sx q[2];
rz(2.1686423) q[2];
rz(-0.074474143) q[3];
sx q[3];
rz(-2.4182033) q[3];
sx q[3];
rz(-0.084666336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
