OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7991601) q[0];
sx q[0];
rz(-2.7437796) q[0];
sx q[0];
rz(1.2678658) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(1.9222577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.08377) q[0];
sx q[0];
rz(-1.738213) q[0];
sx q[0];
rz(2.1517702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5930685) q[2];
sx q[2];
rz(-1.4269575) q[2];
sx q[2];
rz(-1.5576236) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0525236) q[1];
sx q[1];
rz(-2.7929651) q[1];
sx q[1];
rz(-2.2050956) q[1];
x q[2];
rz(2.2279927) q[3];
sx q[3];
rz(-1.7570279) q[3];
sx q[3];
rz(2.5600764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5775602) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.3684028) q[2];
rz(-0.91156256) q[3];
sx q[3];
rz(-1.3699968) q[3];
sx q[3];
rz(-0.024356775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1034705) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(-0.11904112) q[0];
rz(-2.0114404) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(-1.7134604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16258612) q[0];
sx q[0];
rz(-1.1631199) q[0];
sx q[0];
rz(-2.5752221) q[0];
rz(0.4006673) q[2];
sx q[2];
rz(-1.180449) q[2];
sx q[2];
rz(2.0290749) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5806313) q[1];
sx q[1];
rz(-2.033618) q[1];
sx q[1];
rz(0.14974071) q[1];
x q[2];
rz(1.2327594) q[3];
sx q[3];
rz(-1.2814953) q[3];
sx q[3];
rz(-2.3551331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.52593645) q[2];
sx q[2];
rz(-1.4639414) q[2];
sx q[2];
rz(2.3770135) q[2];
rz(0.48314759) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(-2.29276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.1044384) q[0];
sx q[0];
rz(-2.8762682) q[0];
sx q[0];
rz(-1.4720488) q[0];
rz(-2.5813591) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(1.4405506) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9870696) q[0];
sx q[0];
rz(-1.3366342) q[0];
sx q[0];
rz(-2.9375141) q[0];
x q[1];
rz(0.0071390634) q[2];
sx q[2];
rz(-2.4666365) q[2];
sx q[2];
rz(2.1330041) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8941108) q[1];
sx q[1];
rz(-1.9097345) q[1];
sx q[1];
rz(2.1268658) q[1];
rz(-0.23594956) q[3];
sx q[3];
rz(-1.7914158) q[3];
sx q[3];
rz(1.556213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6560087) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(2.7991925) q[2];
rz(1.7229236) q[3];
sx q[3];
rz(-0.72903052) q[3];
sx q[3];
rz(-2.5820144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90958428) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.5268071) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5680983) q[1];
sx q[1];
rz(-1.9400914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9634092) q[0];
sx q[0];
rz(-2.3187175) q[0];
sx q[0];
rz(-1.9848518) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3551201) q[2];
sx q[2];
rz(-2.3209089) q[2];
sx q[2];
rz(-2.414626) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2953342) q[1];
sx q[1];
rz(-0.82678079) q[1];
sx q[1];
rz(1.0215525) q[1];
rz(-0.19530794) q[3];
sx q[3];
rz(-2.8805519) q[3];
sx q[3];
rz(1.9633499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3379007) q[2];
sx q[2];
rz(-0.98946977) q[2];
sx q[2];
rz(1.2376415) q[2];
rz(1.8303653) q[3];
sx q[3];
rz(-1.6236191) q[3];
sx q[3];
rz(-1.1782014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(2.232724) q[0];
rz(0.88242775) q[1];
sx q[1];
rz(-2.4274223) q[1];
sx q[1];
rz(2.4095101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9764938) q[0];
sx q[0];
rz(-1.1015176) q[0];
sx q[0];
rz(2.3114572) q[0];
rz(-1.8514439) q[2];
sx q[2];
rz(-0.82260231) q[2];
sx q[2];
rz(-1.4952687) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.656658) q[1];
sx q[1];
rz(-2.0044494) q[1];
sx q[1];
rz(-1.1364494) q[1];
rz(-pi) q[2];
rz(-0.83314912) q[3];
sx q[3];
rz(-2.1705496) q[3];
sx q[3];
rz(0.38755998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0714134) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(-1.8522235) q[2];
rz(-1.4687126) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(-0.41951352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5759204) q[0];
sx q[0];
rz(-1.1812295) q[0];
sx q[0];
rz(-2.0400203) q[0];
rz(0.22483243) q[1];
sx q[1];
rz(-1.3460527) q[1];
sx q[1];
rz(-0.98519957) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48822847) q[0];
sx q[0];
rz(-0.66500992) q[0];
sx q[0];
rz(-1.8652161) q[0];
rz(-pi) q[1];
rz(1.2217789) q[2];
sx q[2];
rz(-0.98528457) q[2];
sx q[2];
rz(1.7988009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0877062) q[1];
sx q[1];
rz(-2.3017677) q[1];
sx q[1];
rz(2.5165416) q[1];
rz(-pi) q[2];
rz(0.043126338) q[3];
sx q[3];
rz(-0.52845983) q[3];
sx q[3];
rz(0.915574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3605986) q[2];
sx q[2];
rz(-2.4807319) q[2];
sx q[2];
rz(1.9471656) q[2];
rz(-3.0418975) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7128971) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(2.0275443) q[0];
rz(1.3975337) q[1];
sx q[1];
rz(-1.3797398) q[1];
sx q[1];
rz(-2.8278415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6614051) q[0];
sx q[0];
rz(-1.620694) q[0];
sx q[0];
rz(1.2779092) q[0];
rz(-pi) q[1];
rz(1.4391104) q[2];
sx q[2];
rz(-1.734231) q[2];
sx q[2];
rz(1.2800467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75668979) q[1];
sx q[1];
rz(-0.80763615) q[1];
sx q[1];
rz(2.7212423) q[1];
rz(-pi) q[2];
rz(1.3082771) q[3];
sx q[3];
rz(-0.90210017) q[3];
sx q[3];
rz(1.4895542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.067001192) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(-2.5355549) q[2];
rz(1.0780942) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(1.7830474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17230497) q[0];
sx q[0];
rz(-2.2242039) q[0];
sx q[0];
rz(-1.1267927) q[0];
rz(-0.91398319) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(-0.21805683) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1153961) q[0];
sx q[0];
rz(-0.65803981) q[0];
sx q[0];
rz(0.81320073) q[0];
rz(-2.3534691) q[2];
sx q[2];
rz(-1.0255073) q[2];
sx q[2];
rz(0.1696378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8149643) q[1];
sx q[1];
rz(-1.6964629) q[1];
sx q[1];
rz(-0.92355048) q[1];
rz(-pi) q[2];
rz(1.7084684) q[3];
sx q[3];
rz(-1.4542504) q[3];
sx q[3];
rz(-2.2046409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6286991) q[2];
sx q[2];
rz(-2.4662374) q[2];
sx q[2];
rz(-2.8524354) q[2];
rz(2.1469927) q[3];
sx q[3];
rz(-1.9528439) q[3];
sx q[3];
rz(0.4579671) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7116123) q[0];
sx q[0];
rz(-2.0531605) q[0];
sx q[0];
rz(-2.3751538) q[0];
rz(-1.537716) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(-2.2235353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1963475) q[0];
sx q[0];
rz(-0.96789384) q[0];
sx q[0];
rz(-2.5640998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7154022) q[2];
sx q[2];
rz(-1.40398) q[2];
sx q[2];
rz(1.6289935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54967105) q[1];
sx q[1];
rz(-1.3240314) q[1];
sx q[1];
rz(-0.73149577) q[1];
rz(-pi) q[2];
rz(1.3206743) q[3];
sx q[3];
rz(-1.4682574) q[3];
sx q[3];
rz(-0.96270442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90091577) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(0.051699836) q[2];
rz(-0.85203552) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(0.26028546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8155415) q[0];
sx q[0];
rz(-1.6764078) q[0];
sx q[0];
rz(-3.0503804) q[0];
rz(0.7971898) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(-2.7104654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.618396) q[0];
sx q[0];
rz(3.0831227) q[0];
rz(-pi) q[1];
rz(2.2195039) q[2];
sx q[2];
rz(-1.6493634) q[2];
sx q[2];
rz(0.16115133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23049244) q[1];
sx q[1];
rz(-2.6395032) q[1];
sx q[1];
rz(2.6761495) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43041269) q[3];
sx q[3];
rz(-1.0135883) q[3];
sx q[3];
rz(2.1490974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.572523) q[2];
sx q[2];
rz(-1.5745682) q[2];
sx q[2];
rz(1.265906) q[2];
rz(1.2931394) q[3];
sx q[3];
rz(-2.0796516) q[3];
sx q[3];
rz(-0.40614793) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1311998) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(-2.5776183) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(-1.3971267) q[2];
sx q[2];
rz(-2.7026855) q[2];
sx q[2];
rz(2.4365946) q[2];
rz(-2.1544477) q[3];
sx q[3];
rz(-0.55716438) q[3];
sx q[3];
rz(-1.1732994) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
