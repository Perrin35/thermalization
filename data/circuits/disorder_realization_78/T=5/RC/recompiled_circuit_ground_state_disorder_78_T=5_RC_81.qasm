OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(-2.3999441) q[1];
sx q[1];
rz(2.9902966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1681874) q[0];
sx q[0];
rz(-3.0136913) q[0];
sx q[0];
rz(-1.3692877) q[0];
x q[1];
rz(0.1163255) q[2];
sx q[2];
rz(-1.0243729) q[2];
sx q[2];
rz(0.31508581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9090167) q[1];
sx q[1];
rz(-1.0006418) q[1];
sx q[1];
rz(2.2297165) q[1];
rz(0.80161867) q[3];
sx q[3];
rz(-1.773917) q[3];
sx q[3];
rz(2.5020848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0155045) q[2];
sx q[2];
rz(-1.2903004) q[2];
sx q[2];
rz(0.31910953) q[2];
rz(-1.1110405) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(1.8266953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023271712) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(2.5530489) q[0];
rz(-0.59666807) q[1];
sx q[1];
rz(-1.3304973) q[1];
sx q[1];
rz(0.24135022) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7241192) q[0];
sx q[0];
rz(-0.80097526) q[0];
sx q[0];
rz(2.1885314) q[0];
x q[1];
rz(1.0701837) q[2];
sx q[2];
rz(-1.6334849) q[2];
sx q[2];
rz(2.3586065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4171364) q[1];
sx q[1];
rz(-1.1358741) q[1];
sx q[1];
rz(3.0318854) q[1];
rz(2.1961658) q[3];
sx q[3];
rz(-1.6540048) q[3];
sx q[3];
rz(-1.041846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0236464) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(3.0493951) q[2];
rz(0.89208952) q[3];
sx q[3];
rz(-2.2690319) q[3];
sx q[3];
rz(1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43111619) q[0];
sx q[0];
rz(-1.6350063) q[0];
sx q[0];
rz(-3.0431252) q[0];
rz(-0.59421986) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(2.1509511) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4020549) q[0];
sx q[0];
rz(-1.6574114) q[0];
sx q[0];
rz(1.195407) q[0];
rz(-0.57035302) q[2];
sx q[2];
rz(-2.2945171) q[2];
sx q[2];
rz(3.0511659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.65730598) q[1];
sx q[1];
rz(-2.4542744) q[1];
sx q[1];
rz(1.9621852) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22728592) q[3];
sx q[3];
rz(-2.0281938) q[3];
sx q[3];
rz(-1.6136839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46382612) q[2];
sx q[2];
rz(-0.78654424) q[2];
sx q[2];
rz(-2.3393935) q[2];
rz(-1.5471316) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.39740729) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(-1.8205951) q[0];
rz(-1.6162704) q[1];
sx q[1];
rz(-0.95322144) q[1];
sx q[1];
rz(2.9811409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4686543) q[0];
sx q[0];
rz(-2.3529691) q[0];
sx q[0];
rz(-1.9140052) q[0];
x q[1];
rz(2.552659) q[2];
sx q[2];
rz(-2.9351165) q[2];
sx q[2];
rz(1.7084509) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.077549283) q[1];
sx q[1];
rz(-2.4740015) q[1];
sx q[1];
rz(-1.8197219) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53628199) q[3];
sx q[3];
rz(-1.9347526) q[3];
sx q[3];
rz(-0.82897794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3204331) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(-1.9197561) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4959167) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(-2.3107279) q[0];
rz(-1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(-2.1499706) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9308335) q[0];
sx q[0];
rz(-2.2607339) q[0];
sx q[0];
rz(-2.7372975) q[0];
rz(-pi) q[1];
rz(2.2575602) q[2];
sx q[2];
rz(-1.7900538) q[2];
sx q[2];
rz(0.21603661) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8960921) q[1];
sx q[1];
rz(-2.292521) q[1];
sx q[1];
rz(0.20770276) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2701459) q[3];
sx q[3];
rz(-2.0021523) q[3];
sx q[3];
rz(0.68701216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7603989) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(-2.7823616) q[2];
rz(-2.4257816) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(1.1904967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.0642218) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(-2.6203058) q[0];
rz(-2.298666) q[1];
sx q[1];
rz(-1.1784252) q[1];
sx q[1];
rz(-0.11046031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960707) q[0];
sx q[0];
rz(-1.8233436) q[0];
sx q[0];
rz(0.70742328) q[0];
x q[1];
rz(-0.02968024) q[2];
sx q[2];
rz(-0.92845193) q[2];
sx q[2];
rz(-2.6386767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4723052) q[1];
sx q[1];
rz(-0.51176039) q[1];
sx q[1];
rz(1.9838347) q[1];
rz(-pi) q[2];
rz(0.25909781) q[3];
sx q[3];
rz(-2.0784335) q[3];
sx q[3];
rz(2.6258084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20737401) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(0.99384394) q[2];
rz(0.55073109) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(1.6186835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9024502) q[0];
sx q[0];
rz(-2.7681594) q[0];
sx q[0];
rz(-0.19110876) q[0];
rz(2.7725819) q[1];
sx q[1];
rz(-1.399682) q[1];
sx q[1];
rz(0.63327995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2804945) q[0];
sx q[0];
rz(-1.8095353) q[0];
sx q[0];
rz(-1.6709836) q[0];
rz(-0.72004135) q[2];
sx q[2];
rz(-1.9332464) q[2];
sx q[2];
rz(-1.8099648) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25200462) q[1];
sx q[1];
rz(-1.7637296) q[1];
sx q[1];
rz(1.361752) q[1];
x q[2];
rz(0.69890396) q[3];
sx q[3];
rz(-1.6972542) q[3];
sx q[3];
rz(-1.5918283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9621027) q[2];
sx q[2];
rz(-1.3733764) q[2];
sx q[2];
rz(2.7607259) q[2];
rz(1.7535836) q[3];
sx q[3];
rz(-2.1151147) q[3];
sx q[3];
rz(-1.4306205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66985828) q[0];
sx q[0];
rz(-2.7405881) q[0];
sx q[0];
rz(-1.5420472) q[0];
rz(1.764864) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(0.9309887) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.087983) q[0];
sx q[0];
rz(-1.6585104) q[0];
sx q[0];
rz(-2.5045519) q[0];
rz(-2.6561894) q[2];
sx q[2];
rz(-0.50386643) q[2];
sx q[2];
rz(0.47270838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.37759763) q[1];
sx q[1];
rz(-0.74255172) q[1];
sx q[1];
rz(-1.5729891) q[1];
rz(-pi) q[2];
rz(2.2263636) q[3];
sx q[3];
rz(-1.9977565) q[3];
sx q[3];
rz(-2.2674436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57637438) q[2];
sx q[2];
rz(-2.5810676) q[2];
sx q[2];
rz(-3.0726688) q[2];
rz(1.8779514) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3987592) q[0];
sx q[0];
rz(-0.95785207) q[0];
sx q[0];
rz(-0.051890705) q[0];
rz(0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(0.35194078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0028241) q[0];
sx q[0];
rz(-1.616298) q[0];
sx q[0];
rz(2.7753434) q[0];
rz(-2.1743618) q[2];
sx q[2];
rz(-2.1082507) q[2];
sx q[2];
rz(3.1307901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5173147) q[1];
sx q[1];
rz(-1.0464962) q[1];
sx q[1];
rz(-0.45714) q[1];
rz(1.3424344) q[3];
sx q[3];
rz(-2.8274837) q[3];
sx q[3];
rz(-0.64947646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(2.7049098) q[2];
rz(2.7501578) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(0.88633886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201685) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(-0.11908764) q[0];
rz(1.8424312) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(1.7838759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0562818) q[0];
sx q[0];
rz(-1.1650024) q[0];
sx q[0];
rz(0.88078518) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7804271) q[2];
sx q[2];
rz(-2.2884011) q[2];
sx q[2];
rz(-1.7055562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75666243) q[1];
sx q[1];
rz(-2.1637056) q[1];
sx q[1];
rz(1.1641527) q[1];
rz(-1.0064899) q[3];
sx q[3];
rz(-2.3052633) q[3];
sx q[3];
rz(1.2655972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9327717) q[2];
sx q[2];
rz(-1.5522771) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(1.6299853) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83196249) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(0.070214179) q[2];
sx q[2];
rz(-1.6754985) q[2];
sx q[2];
rz(-1.4509261) q[2];
rz(2.3336505) q[3];
sx q[3];
rz(-1.0026889) q[3];
sx q[3];
rz(-1.9849594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
