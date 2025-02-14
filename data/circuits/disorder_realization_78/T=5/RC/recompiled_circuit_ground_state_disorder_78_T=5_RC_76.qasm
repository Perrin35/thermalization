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
rz(1.0977828) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(-2.3999441) q[1];
sx q[1];
rz(2.9902966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1681874) q[0];
sx q[0];
rz(-0.12790132) q[0];
sx q[0];
rz(1.772305) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0252671) q[2];
sx q[2];
rz(-1.0243729) q[2];
sx q[2];
rz(-2.8265068) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.199281) q[1];
sx q[1];
rz(-1.0293055) q[1];
sx q[1];
rz(-0.6813867) q[1];
x q[2];
rz(-0.27917464) q[3];
sx q[3];
rz(-2.320259) q[3];
sx q[3];
rz(-0.73842919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1260881) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(2.8224831) q[2];
rz(-2.0305521) q[3];
sx q[3];
rz(-0.55912656) q[3];
sx q[3];
rz(1.8266953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023271712) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(-2.5530489) q[0];
rz(2.5449246) q[1];
sx q[1];
rz(-1.3304973) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38743153) q[0];
sx q[0];
rz(-1.9997055) q[0];
sx q[0];
rz(2.2700381) q[0];
rz(-pi) q[1];
rz(1.7008408) q[2];
sx q[2];
rz(-0.50419129) q[2];
sx q[2];
rz(0.90182226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4171364) q[1];
sx q[1];
rz(-1.1358741) q[1];
sx q[1];
rz(0.10970727) q[1];
rz(1.7123132) q[3];
sx q[3];
rz(-2.5114473) q[3];
sx q[3];
rz(-2.7272448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0236464) q[2];
sx q[2];
rz(-1.6829374) q[2];
sx q[2];
rz(0.092197593) q[2];
rz(0.89208952) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(-1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43111619) q[0];
sx q[0];
rz(-1.6350063) q[0];
sx q[0];
rz(-0.098467501) q[0];
rz(0.59421986) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(0.99064151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7567609) q[0];
sx q[0];
rz(-0.38479003) q[0];
sx q[0];
rz(1.8033474) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57035302) q[2];
sx q[2];
rz(-2.2945171) q[2];
sx q[2];
rz(0.090426771) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65730598) q[1];
sx q[1];
rz(-2.4542744) q[1];
sx q[1];
rz(1.1794075) q[1];
rz(-0.22728592) q[3];
sx q[3];
rz(-2.0281938) q[3];
sx q[3];
rz(-1.6136839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6777665) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(-0.80219913) q[2];
rz(1.5944611) q[3];
sx q[3];
rz(-1.0651257) q[3];
sx q[3];
rz(0.86435634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39740729) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(1.8205951) q[0];
rz(1.6162704) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(-0.1604518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4686543) q[0];
sx q[0];
rz(-2.3529691) q[0];
sx q[0];
rz(-1.2275874) q[0];
rz(2.9691485) q[2];
sx q[2];
rz(-1.4566696) q[2];
sx q[2];
rz(0.44140377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2962118) q[1];
sx q[1];
rz(-1.7239162) q[1];
sx q[1];
rz(0.91836849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53628199) q[3];
sx q[3];
rz(-1.9347526) q[3];
sx q[3];
rz(0.82897794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82115951) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(-0.56905812) q[2];
rz(1.9197561) q[3];
sx q[3];
rz(-0.33908436) q[3];
sx q[3];
rz(-0.1327742) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959167) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(-0.83086479) q[0];
rz(-1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(0.99162203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9308335) q[0];
sx q[0];
rz(-0.88085876) q[0];
sx q[0];
rz(0.4042952) q[0];
rz(-1.2328202) q[2];
sx q[2];
rz(-2.4261279) q[2];
sx q[2];
rz(-1.6140661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2048671) q[1];
sx q[1];
rz(-0.74581742) q[1];
sx q[1];
rz(-1.3406483) q[1];
rz(-pi) q[2];
rz(-0.95007054) q[3];
sx q[3];
rz(-0.80214989) q[3];
sx q[3];
rz(2.7190894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3811938) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(-0.35923108) q[2];
rz(-2.4257816) q[3];
sx q[3];
rz(-2.3575213) q[3];
sx q[3];
rz(1.951096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0773709) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(0.52128681) q[0];
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
rz(-0.045521988) q[0];
sx q[0];
rz(-1.8233436) q[0];
sx q[0];
rz(0.70742328) q[0];
rz(-pi) q[1];
rz(-0.02968024) q[2];
sx q[2];
rz(-2.2131407) q[2];
sx q[2];
rz(2.6386767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6781194) q[1];
sx q[1];
rz(-1.7686525) q[1];
sx q[1];
rz(-2.045928) q[1];
rz(2.8824948) q[3];
sx q[3];
rz(-1.0631592) q[3];
sx q[3];
rz(-0.51578427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20737401) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(-0.99384394) q[2];
rz(0.55073109) q[3];
sx q[3];
rz(-2.6271074) q[3];
sx q[3];
rz(-1.6186835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2391424) q[0];
sx q[0];
rz(-2.7681594) q[0];
sx q[0];
rz(0.19110876) q[0];
rz(0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(0.63327995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2804945) q[0];
sx q[0];
rz(-1.8095353) q[0];
sx q[0];
rz(-1.4706091) q[0];
x q[1];
rz(-2.6197144) q[2];
sx q[2];
rz(-0.79127914) q[2];
sx q[2];
rz(-2.996614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25200462) q[1];
sx q[1];
rz(-1.377863) q[1];
sx q[1];
rz(-1.7798406) q[1];
x q[2];
rz(-1.406226) q[3];
sx q[3];
rz(-0.87858445) q[3];
sx q[3];
rz(0.084567955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9621027) q[2];
sx q[2];
rz(-1.3733764) q[2];
sx q[2];
rz(-0.38086677) q[2];
rz(-1.3880091) q[3];
sx q[3];
rz(-2.1151147) q[3];
sx q[3];
rz(1.7109722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66985828) q[0];
sx q[0];
rz(-2.7405881) q[0];
sx q[0];
rz(-1.5420472) q[0];
rz(-1.764864) q[1];
sx q[1];
rz(-1.3841261) q[1];
sx q[1];
rz(-2.210604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0536097) q[0];
sx q[0];
rz(-1.4830822) q[0];
sx q[0];
rz(-0.63704078) q[0];
rz(-0.48540326) q[2];
sx q[2];
rz(-0.50386643) q[2];
sx q[2];
rz(2.6688843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37759763) q[1];
sx q[1];
rz(-2.3990409) q[1];
sx q[1];
rz(1.5686036) q[1];
rz(-pi) q[2];
rz(2.6205711) q[3];
sx q[3];
rz(-0.98257321) q[3];
sx q[3];
rz(-2.7532492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57637438) q[2];
sx q[2];
rz(-2.5810676) q[2];
sx q[2];
rz(-3.0726688) q[2];
rz(-1.8779514) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(0.69839683) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7428335) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(0.051890705) q[0];
rz(0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-2.7896519) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55052763) q[0];
sx q[0];
rz(-1.9366486) q[0];
sx q[0];
rz(1.5220674) q[0];
rz(-pi) q[1];
rz(-0.96723084) q[2];
sx q[2];
rz(-1.033342) q[2];
sx q[2];
rz(3.1307901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7406977) q[1];
sx q[1];
rz(-2.4603421) q[1];
sx q[1];
rz(-2.2227312) q[1];
rz(-pi) q[2];
rz(-0.073411302) q[3];
sx q[3];
rz(-1.8764827) q[3];
sx q[3];
rz(-2.2524407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6841782) q[2];
sx q[2];
rz(-2.4922721) q[2];
sx q[2];
rz(-0.43668288) q[2];
rz(0.39143482) q[3];
sx q[3];
rz(-1.6792363) q[3];
sx q[3];
rz(-2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142414) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(-3.022505) q[0];
rz(1.2991615) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(1.3577168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3411451) q[0];
sx q[0];
rz(-0.94609944) q[0];
sx q[0];
rz(0.50826061) q[0];
rz(-2.3215649) q[2];
sx q[2];
rz(-1.301328) q[2];
sx q[2];
rz(2.7633689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7283145) q[1];
sx q[1];
rz(-0.70487805) q[1];
sx q[1];
rz(2.610763) q[1];
x q[2];
rz(2.3230053) q[3];
sx q[3];
rz(-1.9789816) q[3];
sx q[3];
rz(0.096000324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9327717) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(-1.5116073) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3096302) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(-1.0328737) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(0.9821427) q[2];
sx q[2];
rz(-3.0155984) q[2];
sx q[2];
rz(2.2830008) q[2];
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
