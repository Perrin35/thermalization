OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3044337) q[0];
sx q[0];
rz(-1.8456012) q[0];
sx q[0];
rz(-2.0438099) q[0];
rz(2.1863565) q[1];
sx q[1];
rz(-0.74164852) q[1];
sx q[1];
rz(0.15129605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97340527) q[0];
sx q[0];
rz(-3.0136913) q[0];
sx q[0];
rz(-1.3692877) q[0];
x q[1];
rz(-3.0252671) q[2];
sx q[2];
rz(-1.0243729) q[2];
sx q[2];
rz(0.31508581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1946434) q[1];
sx q[1];
rz(-2.2990755) q[1];
sx q[1];
rz(2.3792653) q[1];
rz(0.80161867) q[3];
sx q[3];
rz(-1.3676757) q[3];
sx q[3];
rz(-2.5020848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1260881) q[2];
sx q[2];
rz(-1.2903004) q[2];
sx q[2];
rz(-0.31910953) q[2];
rz(1.1110405) q[3];
sx q[3];
rz(-0.55912656) q[3];
sx q[3];
rz(-1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023271712) q[0];
sx q[0];
rz(-0.51902223) q[0];
sx q[0];
rz(0.58854377) q[0];
rz(0.59666807) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(0.24135022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7241192) q[0];
sx q[0];
rz(-0.80097526) q[0];
sx q[0];
rz(-2.1885314) q[0];
rz(-1.4407519) q[2];
sx q[2];
rz(-2.6374014) q[2];
sx q[2];
rz(2.2397704) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0343098) q[1];
sx q[1];
rz(-1.4713381) q[1];
sx q[1];
rz(-1.1335655) q[1];
x q[2];
rz(2.1961658) q[3];
sx q[3];
rz(-1.4875879) q[3];
sx q[3];
rz(-2.0997467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1179463) q[2];
sx q[2];
rz(-1.6829374) q[2];
sx q[2];
rz(-0.092197593) q[2];
rz(-0.89208952) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43111619) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-3.0431252) q[0];
rz(-2.5473728) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(-2.1509511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7567609) q[0];
sx q[0];
rz(-0.38479003) q[0];
sx q[0];
rz(-1.8033474) q[0];
rz(-pi) q[1];
rz(-0.76108257) q[2];
sx q[2];
rz(-1.154261) q[2];
sx q[2];
rz(2.0628945) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16690635) q[1];
sx q[1];
rz(-2.1975127) q[1];
sx q[1];
rz(2.8381366) q[1];
rz(-pi) q[2];
rz(0.22728592) q[3];
sx q[3];
rz(-1.1133988) q[3];
sx q[3];
rz(-1.6136839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6777665) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(-0.80219913) q[2];
rz(-1.5471316) q[3];
sx q[3];
rz(-1.0651257) q[3];
sx q[3];
rz(0.86435634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-0.39740729) q[0];
sx q[0];
rz(-2.7825401) q[0];
sx q[0];
rz(-1.3209976) q[0];
rz(-1.6162704) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(-2.9811409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67293834) q[0];
sx q[0];
rz(-2.3529691) q[0];
sx q[0];
rz(-1.9140052) q[0];
rz(-0.58893369) q[2];
sx q[2];
rz(-2.9351165) q[2];
sx q[2];
rz(-1.4331417) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7509979) q[1];
sx q[1];
rz(-0.92727755) q[1];
sx q[1];
rz(2.949763) q[1];
x q[2];
rz(-1.9879278) q[3];
sx q[3];
rz(-2.0685745) q[3];
sx q[3];
rz(-2.1912632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82115951) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(-0.56905812) q[2];
rz(-1.9197561) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(-0.1327742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64567599) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(-0.83086479) q[0];
rz(1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(-0.99162203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2107591) q[0];
sx q[0];
rz(-0.88085876) q[0];
sx q[0];
rz(2.7372975) q[0];
rz(-pi) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(0.93672559) q[1];
sx q[1];
rz(-0.74581742) q[1];
sx q[1];
rz(-1.8009443) q[1];
x q[2];
rz(2.1915221) q[3];
sx q[3];
rz(-2.3394428) q[3];
sx q[3];
rz(0.42250326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3811938) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(2.7823616) q[2];
rz(-0.71581101) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(-1.1904967) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0642218) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(-2.6203058) q[0];
rz(-0.84292665) q[1];
sx q[1];
rz(-1.1784252) q[1];
sx q[1];
rz(-3.0311323) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045521988) q[0];
sx q[0];
rz(-1.8233436) q[0];
sx q[0];
rz(2.4341694) q[0];
rz(-pi) q[1];
rz(-3.1119124) q[2];
sx q[2];
rz(-0.92845193) q[2];
sx q[2];
rz(2.6386767) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6692874) q[1];
sx q[1];
rz(-2.6298323) q[1];
sx q[1];
rz(1.9838347) q[1];
x q[2];
rz(2.0024226) q[3];
sx q[3];
rz(-0.56474287) q[3];
sx q[3];
rz(3.125001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9342186) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(0.99384394) q[2];
rz(-2.5908616) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(1.6186835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2391424) q[0];
sx q[0];
rz(-2.7681594) q[0];
sx q[0];
rz(0.19110876) q[0];
rz(0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-2.5083127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8610982) q[0];
sx q[0];
rz(-1.3320574) q[0];
sx q[0];
rz(1.4706091) q[0];
rz(-pi) q[1];
rz(2.6197144) q[2];
sx q[2];
rz(-2.3503135) q[2];
sx q[2];
rz(0.14497862) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7821473) q[1];
sx q[1];
rz(-1.3656865) q[1];
sx q[1];
rz(-0.19711776) q[1];
x q[2];
rz(1.406226) q[3];
sx q[3];
rz(-0.87858445) q[3];
sx q[3];
rz(-0.084567955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4717344) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(1.764864) q[1];
sx q[1];
rz(-1.3841261) q[1];
sx q[1];
rz(2.210604) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.72351) q[0];
sx q[0];
rz(-2.2049954) q[0];
sx q[0];
rz(-1.67976) q[0];
rz(2.6561894) q[2];
sx q[2];
rz(-2.6377262) q[2];
sx q[2];
rz(-2.6688843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37759763) q[1];
sx q[1];
rz(-0.74255172) q[1];
sx q[1];
rz(1.5729891) q[1];
rz(-0.91522907) q[3];
sx q[3];
rz(-1.9977565) q[3];
sx q[3];
rz(-2.2674436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5652183) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(3.0726688) q[2];
rz(1.8779514) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(-0.69839683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7428335) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(-3.0897019) q[0];
rz(-0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(2.7896519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.591065) q[0];
sx q[0];
rz(-1.204944) q[0];
sx q[0];
rz(1.5220674) q[0];
rz(-pi) q[1];
rz(-2.1743618) q[2];
sx q[2];
rz(-1.033342) q[2];
sx q[2];
rz(-3.1307901) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4365153) q[1];
sx q[1];
rz(-1.962858) q[1];
sx q[1];
rz(0.99832557) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8772576) q[3];
sx q[3];
rz(-1.6407986) q[3];
sx q[3];
rz(2.4378192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(0.43668288) q[2];
rz(0.39143482) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6201685) q[0];
sx q[0];
rz(-0.7069718) q[0];
sx q[0];
rz(-0.11908764) q[0];
rz(-1.8424312) q[1];
sx q[1];
rz(-1.3928587) q[1];
sx q[1];
rz(1.7838759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8004476) q[0];
sx q[0];
rz(-0.94609944) q[0];
sx q[0];
rz(0.50826061) q[0];
rz(-1.18612) q[2];
sx q[2];
rz(-2.3529077) q[2];
sx q[2];
rz(-2.2269611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57798856) q[1];
sx q[1];
rz(-1.23659) q[1];
sx q[1];
rz(-2.5086705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6068654) q[3];
sx q[3];
rz(-2.2486454) q[3];
sx q[3];
rz(-1.1191561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9327717) q[2];
sx q[2];
rz(-1.5522771) q[2];
sx q[2];
rz(2.5271752) q[2];
rz(-1.5116073) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.83196249) q[0];
sx q[0];
rz(-0.93292581) q[0];
sx q[0];
rz(0.25136872) q[0];
rz(1.0328737) q[1];
sx q[1];
rz(-1.1943457) q[1];
sx q[1];
rz(1.7051382) q[1];
rz(3.0713785) q[2];
sx q[2];
rz(-1.4660942) q[2];
sx q[2];
rz(1.6906665) q[2];
rz(-2.3165807) q[3];
sx q[3];
rz(-2.2259983) q[3];
sx q[3];
rz(-3.0430924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
