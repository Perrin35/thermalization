OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34243256) q[0];
sx q[0];
rz(-0.39781308) q[0];
sx q[0];
rz(-1.2678658) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(-1.2193349) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3800664) q[0];
sx q[0];
rz(-0.60193578) q[0];
sx q[0];
rz(-1.2720889) q[0];
rz(-pi) q[1];
rz(-2.8706495) q[2];
sx q[2];
rz(-2.5763955) q[2];
sx q[2];
rz(-0.2172367) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87675512) q[1];
sx q[1];
rz(-1.7746468) q[1];
sx q[1];
rz(1.8556103) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9079535) q[3];
sx q[3];
rz(-0.92687449) q[3];
sx q[3];
rz(-1.131191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56403247) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(0.91156256) q[3];
sx q[3];
rz(-1.3699968) q[3];
sx q[3];
rz(-3.1172359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1034705) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(-3.0225515) q[0];
rz(2.0114404) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(1.7134604) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9790065) q[0];
sx q[0];
rz(-1.9784728) q[0];
sx q[0];
rz(2.5752221) q[0];
rz(0.4006673) q[2];
sx q[2];
rz(-1.9611437) q[2];
sx q[2];
rz(1.1125178) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5806313) q[1];
sx q[1];
rz(-1.1079746) q[1];
sx q[1];
rz(-0.14974071) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9088332) q[3];
sx q[3];
rz(-1.2814953) q[3];
sx q[3];
rz(0.78645951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52593645) q[2];
sx q[2];
rz(-1.4639414) q[2];
sx q[2];
rz(-2.3770135) q[2];
rz(0.48314759) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(0.84883261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0371542) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(1.6695439) q[0];
rz(-2.5813591) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(1.4405506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2587145) q[0];
sx q[0];
rz(-2.8322161) q[0];
sx q[0];
rz(2.2750654) q[0];
rz(-pi) q[1];
rz(2.4666489) q[2];
sx q[2];
rz(-1.5663354) q[2];
sx q[2];
rz(0.5677815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24748188) q[1];
sx q[1];
rz(-1.2318582) q[1];
sx q[1];
rz(2.1268658) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3441003) q[3];
sx q[3];
rz(-1.8009225) q[3];
sx q[3];
rz(-0.067148681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(0.34240016) q[2];
rz(-1.418669) q[3];
sx q[3];
rz(-2.4125621) q[3];
sx q[3];
rz(2.5820144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320084) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.5268071) q[0];
rz(0.80742637) q[1];
sx q[1];
rz(-1.5680983) q[1];
sx q[1];
rz(1.2015013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9634092) q[0];
sx q[0];
rz(-2.3187175) q[0];
sx q[0];
rz(1.1567409) q[0];
rz(-2.2204578) q[2];
sx q[2];
rz(-2.1138722) q[2];
sx q[2];
rz(1.7004101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0240942) q[1];
sx q[1];
rz(-1.9650241) q[1];
sx q[1];
rz(-2.3181163) q[1];
rz(-pi) q[2];
rz(1.5189999) q[3];
sx q[3];
rz(-1.8267617) q[3];
sx q[3];
rz(-1.761375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3379007) q[2];
sx q[2];
rz(-0.98946977) q[2];
sx q[2];
rz(-1.2376415) q[2];
rz(-1.8303653) q[3];
sx q[3];
rz(-1.6236191) q[3];
sx q[3];
rz(-1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118498) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(-2.232724) q[0];
rz(-0.88242775) q[1];
sx q[1];
rz(-2.4274223) q[1];
sx q[1];
rz(-2.4095101) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9764938) q[0];
sx q[0];
rz(-2.0400751) q[0];
sx q[0];
rz(-2.3114572) q[0];
rz(-pi) q[1];
rz(0.76811789) q[2];
sx q[2];
rz(-1.3663732) q[2];
sx q[2];
rz(0.26917514) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27837023) q[1];
sx q[1];
rz(-1.178982) q[1];
sx q[1];
rz(2.6696221) q[1];
x q[2];
rz(2.3084435) q[3];
sx q[3];
rz(-0.97104302) q[3];
sx q[3];
rz(-0.38755998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.070179209) q[2];
sx q[2];
rz(-2.3290403) q[2];
sx q[2];
rz(1.8522235) q[2];
rz(1.4687126) q[3];
sx q[3];
rz(-1.2082992) q[3];
sx q[3];
rz(-0.41951352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656723) q[0];
sx q[0];
rz(-1.1812295) q[0];
sx q[0];
rz(2.0400203) q[0];
rz(0.22483243) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(-2.1563931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0211739) q[0];
sx q[0];
rz(-2.2024931) q[0];
sx q[0];
rz(2.9178502) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9198138) q[2];
sx q[2];
rz(-0.98528457) q[2];
sx q[2];
rz(1.7988009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0538865) q[1];
sx q[1];
rz(-0.83982491) q[1];
sx q[1];
rz(0.62505109) q[1];
x q[2];
rz(1.5456301) q[3];
sx q[3];
rz(-1.0428793) q[3];
sx q[3];
rz(-2.2759469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3605986) q[2];
sx q[2];
rz(-2.4807319) q[2];
sx q[2];
rz(-1.9471656) q[2];
rz(0.099695168) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(2.1626507) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128971) q[0];
sx q[0];
rz(-2.6850057) q[0];
sx q[0];
rz(2.0275443) q[0];
rz(-1.7440589) q[1];
sx q[1];
rz(-1.7618529) q[1];
sx q[1];
rz(-0.31375113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6614051) q[0];
sx q[0];
rz(-1.620694) q[0];
sx q[0];
rz(1.8636835) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9767562) q[2];
sx q[2];
rz(-1.7007174) q[2];
sx q[2];
rz(-0.31229737) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3849029) q[1];
sx q[1];
rz(-2.3339565) q[1];
sx q[1];
rz(0.4203504) q[1];
rz(-pi) q[2];
rz(-2.8242495) q[3];
sx q[3];
rz(-2.4306524) q[3];
sx q[3];
rz(1.2430354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0745915) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(-0.60603777) q[2];
rz(-1.0780942) q[3];
sx q[3];
rz(-0.86229101) q[3];
sx q[3];
rz(-1.3585453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17230497) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(-1.1267927) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(-0.21805683) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.954531) q[0];
sx q[0];
rz(-2.0045223) q[0];
sx q[0];
rz(-1.0591255) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78812353) q[2];
sx q[2];
rz(-1.0255073) q[2];
sx q[2];
rz(-0.1696378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3266284) q[1];
sx q[1];
rz(-1.4451298) q[1];
sx q[1];
rz(0.92355048) q[1];
x q[2];
rz(0.86450926) q[3];
sx q[3];
rz(-0.1801404) q[3];
sx q[3];
rz(-0.064398191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-2.4662374) q[2];
sx q[2];
rz(-0.2891573) q[2];
rz(-2.1469927) q[3];
sx q[3];
rz(-1.9528439) q[3];
sx q[3];
rz(2.6836256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.91805735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.799918) q[0];
sx q[0];
rz(-0.8091439) q[0];
sx q[0];
rz(-2.2412712) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16854281) q[2];
sx q[2];
rz(-1.7133811) q[2];
sx q[2];
rz(-0.082372168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80528211) q[1];
sx q[1];
rz(-2.2754221) q[1];
sx q[1];
rz(1.8971838) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0357846) q[3];
sx q[3];
rz(-1.322016) q[3];
sx q[3];
rz(-2.5073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2406769) q[2];
sx q[2];
rz(-0.46611163) q[2];
sx q[2];
rz(-0.051699836) q[2];
rz(2.2895571) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(0.26028546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8155415) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(0.091212243) q[0];
rz(-2.3444029) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(0.43112722) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.618396) q[0];
sx q[0];
rz(3.0831227) q[0];
rz(2.2195039) q[2];
sx q[2];
rz(-1.6493634) q[2];
sx q[2];
rz(-2.9804413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9111002) q[1];
sx q[1];
rz(-0.5020895) q[1];
sx q[1];
rz(-2.6761495) q[1];
rz(-pi) q[2];
rz(0.43041269) q[3];
sx q[3];
rz(-2.1280043) q[3];
sx q[3];
rz(-2.1490974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.572523) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(1.8756867) q[2];
rz(-1.2931394) q[3];
sx q[3];
rz(-2.0796516) q[3];
sx q[3];
rz(-2.7354447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1311998) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(0.56397437) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(1.1376913) q[2];
sx q[2];
rz(-1.4972996) q[2];
sx q[2];
rz(1.0232915) q[2];
rz(1.0913783) q[3];
sx q[3];
rz(-1.2751083) q[3];
sx q[3];
rz(0.90838065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
