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
rz(-2.2351216) q[0];
sx q[0];
rz(-1.6989166) q[0];
sx q[0];
rz(-0.28191167) q[0];
rz(0.52892041) q[1];
sx q[1];
rz(4.79098) q[1];
sx q[1];
rz(11.00287) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53858763) q[0];
sx q[0];
rz(-1.2055802) q[0];
sx q[0];
rz(0.19293789) q[0];
rz(-pi) q[1];
rz(-0.018997832) q[2];
sx q[2];
rz(-2.8667031) q[2];
sx q[2];
rz(-1.887448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2927478) q[1];
sx q[1];
rz(-2.1969921) q[1];
sx q[1];
rz(1.7262154) q[1];
rz(-0.36564499) q[3];
sx q[3];
rz(-2.4198101) q[3];
sx q[3];
rz(-2.4810098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4549183) q[2];
sx q[2];
rz(-0.29759559) q[2];
sx q[2];
rz(2.5285524) q[2];
rz(2.6679299) q[3];
sx q[3];
rz(-1.1981755) q[3];
sx q[3];
rz(-1.7090428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050201) q[0];
sx q[0];
rz(-2.9614145) q[0];
sx q[0];
rz(-2.3905684) q[0];
rz(-2.660102) q[1];
sx q[1];
rz(-2.057169) q[1];
sx q[1];
rz(2.1717333) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4879726) q[0];
sx q[0];
rz(-0.75838415) q[0];
sx q[0];
rz(-2.5874596) q[0];
rz(1.5300691) q[2];
sx q[2];
rz(-1.7864405) q[2];
sx q[2];
rz(3.0037896) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10895059) q[1];
sx q[1];
rz(-1.8284997) q[1];
sx q[1];
rz(0.9089246) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5447846) q[3];
sx q[3];
rz(-1.1973901) q[3];
sx q[3];
rz(1.663409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91886175) q[2];
sx q[2];
rz(-1.4758045) q[2];
sx q[2];
rz(0.53885031) q[2];
rz(-0.017596267) q[3];
sx q[3];
rz(-2.9774057) q[3];
sx q[3];
rz(-2.2331451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410412) q[0];
sx q[0];
rz(-0.38527641) q[0];
sx q[0];
rz(3.0803296) q[0];
rz(1.6368658) q[1];
sx q[1];
rz(-0.69925362) q[1];
sx q[1];
rz(1.9452852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17249566) q[0];
sx q[0];
rz(-2.3845256) q[0];
sx q[0];
rz(1.2189381) q[0];
rz(-pi) q[1];
rz(-1.5012791) q[2];
sx q[2];
rz(-1.3090773) q[2];
sx q[2];
rz(1.4756605) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4704901) q[1];
sx q[1];
rz(-1.2467978) q[1];
sx q[1];
rz(-0.98172202) q[1];
rz(-pi) q[2];
rz(-1.1617834) q[3];
sx q[3];
rz(-2.4065285) q[3];
sx q[3];
rz(0.72003698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4521744) q[2];
sx q[2];
rz(-1.3070561) q[2];
sx q[2];
rz(0.43251953) q[2];
rz(1.0906667) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(1.0961016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5525621) q[0];
sx q[0];
rz(-2.2046389) q[0];
sx q[0];
rz(0.20137782) q[0];
rz(-2.6248113) q[1];
sx q[1];
rz(-2.7613381) q[1];
sx q[1];
rz(2.1929599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.082371) q[0];
sx q[0];
rz(-1.9192524) q[0];
sx q[0];
rz(-2.2414464) q[0];
x q[1];
rz(-1.2327551) q[2];
sx q[2];
rz(-2.4043407) q[2];
sx q[2];
rz(0.33589943) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3310254) q[1];
sx q[1];
rz(-2.0750891) q[1];
sx q[1];
rz(1.585998) q[1];
rz(-1.6679156) q[3];
sx q[3];
rz(-2.0896974) q[3];
sx q[3];
rz(3.1176381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8676694) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(-0.55740994) q[2];
rz(-3.0607767) q[3];
sx q[3];
rz(-1.2405453) q[3];
sx q[3];
rz(-1.935299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7252561) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(2.985785) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(-2.9373998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9568528) q[0];
sx q[0];
rz(-2.889688) q[0];
sx q[0];
rz(2.8624318) q[0];
rz(2.5956018) q[2];
sx q[2];
rz(-0.24790774) q[2];
sx q[2];
rz(-1.966983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62120092) q[1];
sx q[1];
rz(-1.0869496) q[1];
sx q[1];
rz(-1.8041759) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3315807) q[3];
sx q[3];
rz(-0.81221928) q[3];
sx q[3];
rz(-1.1438469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2121409) q[2];
sx q[2];
rz(-1.4931623) q[2];
sx q[2];
rz(0.41352752) q[2];
rz(0.84189576) q[3];
sx q[3];
rz(-1.7588408) q[3];
sx q[3];
rz(2.1690185) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095116422) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(0.0055775642) q[0];
rz(-1.7976409) q[1];
sx q[1];
rz(-0.19897142) q[1];
sx q[1];
rz(-2.4406348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0049135) q[0];
sx q[0];
rz(-0.71109301) q[0];
sx q[0];
rz(-0.85791608) q[0];
x q[1];
rz(0.92971071) q[2];
sx q[2];
rz(-1.2107009) q[2];
sx q[2];
rz(1.6508697) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.13853874) q[1];
sx q[1];
rz(-1.8459311) q[1];
sx q[1];
rz(-1.5281926) q[1];
rz(-2.088955) q[3];
sx q[3];
rz(-2.0817882) q[3];
sx q[3];
rz(2.4571153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31412101) q[2];
sx q[2];
rz(-0.33385971) q[2];
sx q[2];
rz(-1.1673048) q[2];
rz(2.2589034) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(2.816443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90671396) q[0];
sx q[0];
rz(-2.0218847) q[0];
sx q[0];
rz(0.085302189) q[0];
rz(2.3872497) q[1];
sx q[1];
rz(-0.5032379) q[1];
sx q[1];
rz(2.1139961) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37759532) q[0];
sx q[0];
rz(-2.9288284) q[0];
sx q[0];
rz(0.048325267) q[0];
rz(0.59664388) q[2];
sx q[2];
rz(-1.5663356) q[2];
sx q[2];
rz(0.73847929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.69313777) q[1];
sx q[1];
rz(-1.4488646) q[1];
sx q[1];
rz(-1.8961402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1291162) q[3];
sx q[3];
rz(-1.9838247) q[3];
sx q[3];
rz(-0.2193887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9895642) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(-0.41934553) q[2];
rz(-2.5069405) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(-2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80970508) q[0];
sx q[0];
rz(-0.87600791) q[0];
sx q[0];
rz(-2.4494655) q[0];
rz(0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(1.7587761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.064474) q[0];
sx q[0];
rz(-3.0713284) q[0];
sx q[0];
rz(1.8514567) q[0];
x q[1];
rz(2.369129) q[2];
sx q[2];
rz(-2.303745) q[2];
sx q[2];
rz(0.39847429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18305138) q[1];
sx q[1];
rz(-2.1223096) q[1];
sx q[1];
rz(0.27490487) q[1];
rz(-pi) q[2];
rz(-1.0442078) q[3];
sx q[3];
rz(-0.79681891) q[3];
sx q[3];
rz(-2.6897893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77070037) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(-2.6254081) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(-1.9019351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6462964) q[0];
sx q[0];
rz(-1.7072059) q[0];
sx q[0];
rz(2.7700951) q[0];
rz(-0.83483541) q[1];
sx q[1];
rz(-2.3284262) q[1];
sx q[1];
rz(3.1239948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478701) q[0];
sx q[0];
rz(-1.2945064) q[0];
sx q[0];
rz(1.0280861) q[0];
rz(-pi) q[1];
rz(2.405898) q[2];
sx q[2];
rz(-1.383505) q[2];
sx q[2];
rz(-0.012030727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4092651) q[1];
sx q[1];
rz(-1.0068934) q[1];
sx q[1];
rz(1.1281518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.892421) q[3];
sx q[3];
rz(-2.1545046) q[3];
sx q[3];
rz(-2.5584588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1060433) q[2];
sx q[2];
rz(-2.273166) q[2];
sx q[2];
rz(3.0412728) q[2];
rz(2.912168) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(-2.6917698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.68596524) q[0];
sx q[0];
rz(-2.8534511) q[0];
sx q[0];
rz(2.567754) q[0];
rz(2.0876743) q[1];
sx q[1];
rz(-1.046448) q[1];
sx q[1];
rz(-2.4651249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8343) q[0];
sx q[0];
rz(-0.3737475) q[0];
sx q[0];
rz(2.4227002) q[0];
rz(-0.25679882) q[2];
sx q[2];
rz(-1.9635824) q[2];
sx q[2];
rz(1.2744255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65634644) q[1];
sx q[1];
rz(-2.606483) q[1];
sx q[1];
rz(-0.17406221) q[1];
rz(-pi) q[2];
rz(1.2419534) q[3];
sx q[3];
rz(-2.5848205) q[3];
sx q[3];
rz(-2.013735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68710589) q[2];
sx q[2];
rz(-0.7578308) q[2];
sx q[2];
rz(1.494361) q[2];
rz(2.2729661) q[3];
sx q[3];
rz(-2.4354911) q[3];
sx q[3];
rz(2.6715265) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710707) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(1.3337878) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(1.6619353) q[2];
sx q[2];
rz(-2.3341134) q[2];
sx q[2];
rz(-2.3011617) q[2];
rz(1.7464433) q[3];
sx q[3];
rz(-1.5271389) q[3];
sx q[3];
rz(2.3300119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
