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
rz(2.1863565) q[1];
sx q[1];
rz(-0.74164852) q[1];
sx q[1];
rz(0.15129605) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3713105) q[0];
sx q[0];
rz(-1.6960959) q[0];
sx q[0];
rz(0.025733982) q[0];
rz(-pi) q[1];
rz(1.7593616) q[2];
sx q[2];
rz(-2.5841568) q[2];
sx q[2];
rz(-3.0477026) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94231168) q[1];
sx q[1];
rz(-1.0293055) q[1];
sx q[1];
rz(2.460206) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27917464) q[3];
sx q[3];
rz(-2.320259) q[3];
sx q[3];
rz(-0.73842919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0155045) q[2];
sx q[2];
rz(-1.2903004) q[2];
sx q[2];
rz(-2.8224831) q[2];
rz(1.1110405) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023271712) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(-2.5530489) q[0];
rz(-2.5449246) q[1];
sx q[1];
rz(-1.3304973) q[1];
sx q[1];
rz(-0.24135022) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4174735) q[0];
sx q[0];
rz(-0.80097526) q[0];
sx q[0];
rz(-0.95306122) q[0];
x q[1];
rz(2.0714089) q[2];
sx q[2];
rz(-1.5081078) q[2];
sx q[2];
rz(2.3586065) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4687522) q[1];
sx q[1];
rz(-2.6939055) q[1];
sx q[1];
rz(-1.3393558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7123132) q[3];
sx q[3];
rz(-2.5114473) q[3];
sx q[3];
rz(2.7272448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0236464) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(-0.092197593) q[2];
rz(-0.89208952) q[3];
sx q[3];
rz(-2.2690319) q[3];
sx q[3];
rz(-1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104765) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-0.098467501) q[0];
rz(-0.59421986) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(2.1509511) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38483175) q[0];
sx q[0];
rz(-2.7568026) q[0];
sx q[0];
rz(1.8033474) q[0];
x q[1];
rz(-2.3805101) q[2];
sx q[2];
rz(-1.154261) q[2];
sx q[2];
rz(1.0786982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2222683) q[1];
sx q[1];
rz(-1.8152555) q[1];
sx q[1];
rz(-2.2198841) q[1];
rz(0.22728592) q[3];
sx q[3];
rz(-2.0281938) q[3];
sx q[3];
rz(-1.5279087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46382612) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(2.3393935) q[2];
rz(1.5471316) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39740729) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(-1.3209976) q[0];
rz(1.6162704) q[1];
sx q[1];
rz(-0.95322144) q[1];
sx q[1];
rz(-2.9811409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65112075) q[0];
sx q[0];
rz(-1.8118389) q[0];
sx q[0];
rz(-2.3293877) q[0];
rz(-pi) q[1];
rz(-1.686626) q[2];
sx q[2];
rz(-1.3994851) q[2];
sx q[2];
rz(-2.0320323) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0640434) q[1];
sx q[1];
rz(-0.66759118) q[1];
sx q[1];
rz(1.8197219) q[1];
x q[2];
rz(2.5009551) q[3];
sx q[3];
rz(-0.63797073) q[3];
sx q[3];
rz(0.20221329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3204331) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(1.2218366) q[3];
sx q[3];
rz(-0.33908436) q[3];
sx q[3];
rz(0.1327742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64567599) q[0];
sx q[0];
rz(-0.78302947) q[0];
sx q[0];
rz(2.3107279) q[0];
rz(-1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(0.99162203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9308335) q[0];
sx q[0];
rz(-0.88085876) q[0];
sx q[0];
rz(0.4042952) q[0];
x q[1];
rz(2.8610271) q[2];
sx q[2];
rz(-2.2380851) q[2];
sx q[2];
rz(-1.1782805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93672559) q[1];
sx q[1];
rz(-0.74581742) q[1];
sx q[1];
rz(1.3406483) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6001106) q[3];
sx q[3];
rz(-0.94621822) q[3];
sx q[3];
rz(-1.9196212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7603989) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(2.7823616) q[2];
rz(2.4257816) q[3];
sx q[3];
rz(-2.3575213) q[3];
sx q[3];
rz(1.1904967) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0642218) q[0];
sx q[0];
rz(-1.8697898) q[0];
sx q[0];
rz(2.6203058) q[0];
rz(0.84292665) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(-3.0311323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9006289) q[0];
sx q[0];
rz(-0.74375737) q[0];
sx q[0];
rz(2.7636011) q[0];
rz(2.213352) q[2];
sx q[2];
rz(-1.5945598) q[2];
sx q[2];
rz(1.0500963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1350491) q[1];
sx q[1];
rz(-1.1056756) q[1];
sx q[1];
rz(0.22174447) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0486097) q[3];
sx q[3];
rz(-1.7966101) q[3];
sx q[3];
rz(-1.9584394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20737401) q[2];
sx q[2];
rz(-1.5997581) q[2];
sx q[2];
rz(2.1477487) q[2];
rz(-2.5908616) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(-1.5229092) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2391424) q[0];
sx q[0];
rz(-2.7681594) q[0];
sx q[0];
rz(-0.19110876) q[0];
rz(-0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-0.63327995) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12145081) q[0];
sx q[0];
rz(-2.883054) q[0];
sx q[0];
rz(-0.38991897) q[0];
rz(-2.0379637) q[2];
sx q[2];
rz(-0.90641253) q[2];
sx q[2];
rz(2.6009022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25200462) q[1];
sx q[1];
rz(-1.377863) q[1];
sx q[1];
rz(-1.7798406) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19509372) q[3];
sx q[3];
rz(-2.4332402) q[3];
sx q[3];
rz(0.1699902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17948991) q[2];
sx q[2];
rz(-1.3733764) q[2];
sx q[2];
rz(-2.7607259) q[2];
rz(-1.7535836) q[3];
sx q[3];
rz(-1.0264779) q[3];
sx q[3];
rz(1.7109722) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66985828) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(-1.5420472) q[0];
rz(-1.3767287) q[1];
sx q[1];
rz(-1.3841261) q[1];
sx q[1];
rz(-0.9309887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.087983) q[0];
sx q[0];
rz(-1.6585104) q[0];
sx q[0];
rz(-0.63704078) q[0];
rz(-1.8225708) q[2];
sx q[2];
rz(-2.012017) q[2];
sx q[2];
rz(-3.0721498) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1915832) q[1];
sx q[1];
rz(-1.5693136) q[1];
sx q[1];
rz(-2.3133469) q[1];
x q[2];
rz(0.91522907) q[3];
sx q[3];
rz(-1.1438362) q[3];
sx q[3];
rz(-2.2674436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57637438) q[2];
sx q[2];
rz(-2.5810676) q[2];
sx q[2];
rz(0.068923846) q[2];
rz(1.8779514) q[3];
sx q[3];
rz(-2.7694323) q[3];
sx q[3];
rz(0.69839683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3987592) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(3.0897019) q[0];
rz(0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-2.7896519) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.591065) q[0];
sx q[0];
rz(-1.9366486) q[0];
sx q[0];
rz(-1.5220674) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62656709) q[2];
sx q[2];
rz(-2.0800903) q[2];
sx q[2];
rz(-1.9208822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4365153) q[1];
sx q[1];
rz(-1.1787346) q[1];
sx q[1];
rz(0.99832557) q[1];
x q[2];
rz(0.073411302) q[3];
sx q[3];
rz(-1.8764827) q[3];
sx q[3];
rz(-0.88915196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45741442) q[2];
sx q[2];
rz(-2.4922721) q[2];
sx q[2];
rz(2.7049098) q[2];
rz(-0.39143482) q[3];
sx q[3];
rz(-1.6792363) q[3];
sx q[3];
rz(-0.88633886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142414) q[0];
sx q[0];
rz(-0.7069718) q[0];
sx q[0];
rz(0.11908764) q[0];
rz(-1.2991615) q[1];
sx q[1];
rz(-1.3928587) q[1];
sx q[1];
rz(1.3577168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0393667) q[0];
sx q[0];
rz(-0.78327562) q[0];
sx q[0];
rz(0.97707763) q[0];
rz(-pi) q[1];
rz(-2.3215649) q[2];
sx q[2];
rz(-1.301328) q[2];
sx q[2];
rz(2.7633689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4132781) q[1];
sx q[1];
rz(-2.4367146) q[1];
sx q[1];
rz(0.5308297) q[1];
rz(-pi) q[2];
rz(2.3230053) q[3];
sx q[3];
rz(-1.162611) q[3];
sx q[3];
rz(3.0455923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.208821) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(1.5116073) q[3];
sx q[3];
rz(-0.67174086) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83196249) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(-2.108719) q[1];
sx q[1];
rz(-1.1943457) q[1];
sx q[1];
rz(1.7051382) q[1];
rz(0.070214179) q[2];
sx q[2];
rz(-1.6754985) q[2];
sx q[2];
rz(-1.4509261) q[2];
rz(-0.72335106) q[3];
sx q[3];
rz(-0.94905973) q[3];
sx q[3];
rz(2.2524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
