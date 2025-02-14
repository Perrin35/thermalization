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
rz(-0.9552362) q[1];
sx q[1];
rz(7.0248338) q[1];
sx q[1];
rz(6.4344814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79729743) q[0];
sx q[0];
rz(-1.5963285) q[0];
sx q[0];
rz(-1.4454557) q[0];
x q[1];
rz(3.0252671) q[2];
sx q[2];
rz(-1.0243729) q[2];
sx q[2];
rz(2.8265068) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23257593) q[1];
sx q[1];
rz(-1.0006418) q[1];
sx q[1];
rz(2.2297165) q[1];
rz(-2.339974) q[3];
sx q[3];
rz(-1.3676757) q[3];
sx q[3];
rz(-2.5020848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0155045) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(2.8224831) q[2];
rz(1.1110405) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023271712) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(0.58854377) q[0];
rz(2.5449246) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(2.9002424) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198181) q[0];
sx q[0];
rz(-2.1960917) q[0];
sx q[0];
rz(0.53859512) q[0];
x q[1];
rz(-2.0714089) q[2];
sx q[2];
rz(-1.6334849) q[2];
sx q[2];
rz(-0.78298616) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4171364) q[1];
sx q[1];
rz(-1.1358741) q[1];
sx q[1];
rz(0.10970727) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1961658) q[3];
sx q[3];
rz(-1.4875879) q[3];
sx q[3];
rz(2.0997467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1179463) q[2];
sx q[2];
rz(-1.6829374) q[2];
sx q[2];
rz(3.0493951) q[2];
rz(0.89208952) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(-1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104765) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-0.098467501) q[0];
rz(-2.5473728) q[1];
sx q[1];
rz(-1.0585982) q[1];
sx q[1];
rz(2.1509511) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0069283) q[0];
sx q[0];
rz(-1.9447088) q[0];
sx q[0];
rz(-0.093061826) q[0];
rz(-pi) q[1];
rz(-1.0223234) q[2];
sx q[2];
rz(-0.88829852) q[2];
sx q[2];
rz(-2.2816531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65730598) q[1];
sx q[1];
rz(-2.4542744) q[1];
sx q[1];
rz(-1.9621852) q[1];
rz(2.9143067) q[3];
sx q[3];
rz(-2.0281938) q[3];
sx q[3];
rz(-1.6136839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46382612) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(-0.80219913) q[2];
rz(1.5944611) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.35905251) q[0];
sx q[0];
rz(1.3209976) q[0];
rz(-1.6162704) q[1];
sx q[1];
rz(-0.95322144) q[1];
sx q[1];
rz(2.9811409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67293834) q[0];
sx q[0];
rz(-2.3529691) q[0];
sx q[0];
rz(-1.9140052) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9691485) q[2];
sx q[2];
rz(-1.4566696) q[2];
sx q[2];
rz(0.44140377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8453809) q[1];
sx q[1];
rz(-1.4176765) q[1];
sx q[1];
rz(0.91836849) q[1];
x q[2];
rz(-1.1536648) q[3];
sx q[3];
rz(-1.0730181) q[3];
sx q[3];
rz(-2.1912632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82115951) q[2];
sx q[2];
rz(-1.4904138) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(1.2218366) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(-0.1327742) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959167) q[0];
sx q[0];
rz(-0.78302947) q[0];
sx q[0];
rz(-2.3107279) q[0];
rz(1.9310541) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(0.99162203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0474393) q[0];
sx q[0];
rz(-1.8790566) q[0];
sx q[0];
rz(0.83931132) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9087725) q[2];
sx q[2];
rz(-0.71546474) q[2];
sx q[2];
rz(-1.6140661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2455006) q[1];
sx q[1];
rz(-2.292521) q[1];
sx q[1];
rz(-0.20770276) q[1];
x q[2];
rz(-2.1915221) q[3];
sx q[3];
rz(-2.3394428) q[3];
sx q[3];
rz(2.7190894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3811938) q[2];
sx q[2];
rz(-0.92279592) q[2];
sx q[2];
rz(-0.35923108) q[2];
rz(-2.4257816) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(-1.951096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.0642218) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(2.6203058) q[0];
rz(2.298666) q[1];
sx q[1];
rz(-1.1784252) q[1];
sx q[1];
rz(0.11046031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960707) q[0];
sx q[0];
rz(-1.318249) q[0];
sx q[0];
rz(2.4341694) q[0];
rz(-pi) q[1];
rz(-3.1119124) q[2];
sx q[2];
rz(-0.92845193) q[2];
sx q[2];
rz(-0.50291598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6692874) q[1];
sx q[1];
rz(-2.6298323) q[1];
sx q[1];
rz(-1.157758) q[1];
x q[2];
rz(-2.0024226) q[3];
sx q[3];
rz(-0.56474287) q[3];
sx q[3];
rz(0.016591681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9342186) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(-0.99384394) q[2];
rz(-0.55073109) q[3];
sx q[3];
rz(-2.6271074) q[3];
sx q[3];
rz(1.6186835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2391424) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(-0.19110876) q[0];
rz(-0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(2.5083127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8610982) q[0];
sx q[0];
rz(-1.3320574) q[0];
sx q[0];
rz(-1.6709836) q[0];
rz(-0.52187829) q[2];
sx q[2];
rz(-2.3503135) q[2];
sx q[2];
rz(0.14497862) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5577199) q[1];
sx q[1];
rz(-0.28350949) q[1];
sx q[1];
rz(-2.3260172) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19509372) q[3];
sx q[3];
rz(-2.4332402) q[3];
sx q[3];
rz(2.9716024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9621027) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(-2.7607259) q[2];
rz(-1.3880091) q[3];
sx q[3];
rz(-2.1151147) q[3];
sx q[3];
rz(1.7109722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.4717344) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(1.764864) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(0.9309887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0536097) q[0];
sx q[0];
rz(-1.6585104) q[0];
sx q[0];
rz(-0.63704078) q[0];
rz(1.3190218) q[2];
sx q[2];
rz(-2.012017) q[2];
sx q[2];
rz(0.069442858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38057391) q[1];
sx q[1];
rz(-0.82824681) q[1];
sx q[1];
rz(3.1395802) q[1];
rz(-pi) q[2];
rz(2.2119207) q[3];
sx q[3];
rz(-2.3768209) q[3];
sx q[3];
rz(1.9509893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5652183) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(0.068923846) q[2];
rz(-1.2636412) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3987592) q[0];
sx q[0];
rz(-0.95785207) q[0];
sx q[0];
rz(-0.051890705) q[0];
rz(-2.9715111) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(0.35194078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55052763) q[0];
sx q[0];
rz(-1.204944) q[0];
sx q[0];
rz(-1.6195253) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76100112) q[2];
sx q[2];
rz(-2.3563851) q[2];
sx q[2];
rz(2.1987555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62427795) q[1];
sx q[1];
rz(-2.0950965) q[1];
sx q[1];
rz(-0.45714) q[1];
rz(-pi) q[2];
rz(-1.8772576) q[3];
sx q[3];
rz(-1.500794) q[3];
sx q[3];
rz(-0.7037735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.45741442) q[2];
sx q[2];
rz(-2.4922721) q[2];
sx q[2];
rz(-0.43668288) q[2];
rz(0.39143482) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(-0.88633886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142414) q[0];
sx q[0];
rz(-0.7069718) q[0];
sx q[0];
rz(-0.11908764) q[0];
rz(-1.8424312) q[1];
sx q[1];
rz(-1.3928587) q[1];
sx q[1];
rz(-1.3577168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0562818) q[0];
sx q[0];
rz(-1.9765903) q[0];
sx q[0];
rz(-0.88078518) q[0];
rz(-pi) q[1];
rz(-1.18612) q[2];
sx q[2];
rz(-0.78868491) q[2];
sx q[2];
rz(2.2269611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5636041) q[1];
sx q[1];
rz(-1.23659) q[1];
sx q[1];
rz(2.5086705) q[1];
rz(-pi) q[2];
rz(2.6068654) q[3];
sx q[3];
rz(-0.89294723) q[3];
sx q[3];
rz(1.1191561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9327717) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(-0.61441747) q[2];
rz(1.5116073) q[3];
sx q[3];
rz(-0.67174086) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83196249) q[0];
sx q[0];
rz(-2.2086668) q[0];
sx q[0];
rz(-2.8902239) q[0];
rz(-1.0328737) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(-1.6757552) q[2];
sx q[2];
rz(-1.6406254) q[2];
sx q[2];
rz(0.1272203) q[2];
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
