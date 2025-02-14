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
rz(0.46718207) q[0];
sx q[0];
rz(-1.2527569) q[0];
sx q[0];
rz(7.0897515) q[0];
rz(2.2944577) q[1];
sx q[1];
rz(-2.6557014) q[1];
sx q[1];
rz(-0.86056217) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33968494) q[0];
sx q[0];
rz(-2.1640477) q[0];
sx q[0];
rz(-2.9723333) q[0];
x q[1];
rz(-2.2769558) q[2];
sx q[2];
rz(-1.3792017) q[2];
sx q[2];
rz(0.77048555) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5570506) q[1];
sx q[1];
rz(-1.8452541) q[1];
sx q[1];
rz(-1.6226136) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2673554) q[3];
sx q[3];
rz(-2.4421066) q[3];
sx q[3];
rz(-1.6523838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90031558) q[2];
sx q[2];
rz(-2.6130455) q[2];
sx q[2];
rz(2.6701374) q[2];
rz(-0.49247646) q[3];
sx q[3];
rz(-1.9086647) q[3];
sx q[3];
rz(-1.8878149) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660698) q[0];
sx q[0];
rz(-2.7560784) q[0];
sx q[0];
rz(-2.8634014) q[0];
rz(1.9506075) q[1];
sx q[1];
rz(-1.8089801) q[1];
sx q[1];
rz(-1.2954378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2076072) q[0];
sx q[0];
rz(-1.1046358) q[0];
sx q[0];
rz(-0.68434478) q[0];
rz(-pi) q[1];
rz(2.5668199) q[2];
sx q[2];
rz(-3.0468371) q[2];
sx q[2];
rz(3.0916328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1300704) q[1];
sx q[1];
rz(-1.5685625) q[1];
sx q[1];
rz(1.190459) q[1];
x q[2];
rz(-2.2601028) q[3];
sx q[3];
rz(-0.75493287) q[3];
sx q[3];
rz(0.024294446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94464716) q[2];
sx q[2];
rz(-1.4713919) q[2];
sx q[2];
rz(-0.50986457) q[2];
rz(-1.6950131) q[3];
sx q[3];
rz(-0.46049419) q[3];
sx q[3];
rz(2.6367326) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.533621) q[0];
sx q[0];
rz(-1.6225659) q[0];
sx q[0];
rz(-0.98440379) q[0];
rz(-3.1254752) q[1];
sx q[1];
rz(-1.1382444) q[1];
sx q[1];
rz(1.791753) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4181217) q[0];
sx q[0];
rz(-1.4051361) q[0];
sx q[0];
rz(0.70933527) q[0];
x q[1];
rz(-1.0764112) q[2];
sx q[2];
rz(-0.62389031) q[2];
sx q[2];
rz(-2.9817493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.35914) q[1];
sx q[1];
rz(-1.5305007) q[1];
sx q[1];
rz(-0.69922478) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0257019) q[3];
sx q[3];
rz(-2.3724764) q[3];
sx q[3];
rz(-2.9562841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71780378) q[2];
sx q[2];
rz(-2.4393647) q[2];
sx q[2];
rz(2.2595937) q[2];
rz(1.9574022) q[3];
sx q[3];
rz(-0.94112527) q[3];
sx q[3];
rz(2.4052896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475567) q[0];
sx q[0];
rz(-1.2268257) q[0];
sx q[0];
rz(0.077022821) q[0];
rz(2.9432964) q[1];
sx q[1];
rz(-1.501333) q[1];
sx q[1];
rz(2.95453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88019218) q[0];
sx q[0];
rz(-0.37366707) q[0];
sx q[0];
rz(-2.6970661) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9603029) q[2];
sx q[2];
rz(-1.9167056) q[2];
sx q[2];
rz(-0.51900253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7834082) q[1];
sx q[1];
rz(-0.71137127) q[1];
sx q[1];
rz(3.0510805) q[1];
rz(-pi) q[2];
rz(2.8253205) q[3];
sx q[3];
rz(-1.894884) q[3];
sx q[3];
rz(-0.23095953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5932811) q[2];
sx q[2];
rz(-2.686794) q[2];
sx q[2];
rz(-1.6853257) q[2];
rz(-2.4233387) q[3];
sx q[3];
rz(-1.7536438) q[3];
sx q[3];
rz(-2.3072402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.89183557) q[0];
sx q[0];
rz(-1.8064073) q[0];
sx q[0];
rz(2.7030113) q[0];
rz(-2.7843685) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(1.8908148) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6092458) q[0];
sx q[0];
rz(-0.43918434) q[0];
sx q[0];
rz(-1.3339551) q[0];
x q[1];
rz(1.3322796) q[2];
sx q[2];
rz(-0.3198238) q[2];
sx q[2];
rz(2.8826098) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4080271) q[1];
sx q[1];
rz(-1.9434612) q[1];
sx q[1];
rz(-1.9994451) q[1];
x q[2];
rz(-2.881348) q[3];
sx q[3];
rz(-1.4422999) q[3];
sx q[3];
rz(-2.5082626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9095416) q[2];
sx q[2];
rz(-0.48018685) q[2];
sx q[2];
rz(1.6298182) q[2];
rz(-0.016228598) q[3];
sx q[3];
rz(-2.0461693) q[3];
sx q[3];
rz(-2.3487263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4474354) q[0];
sx q[0];
rz(-1.2801535) q[0];
sx q[0];
rz(1.9919027) q[0];
rz(2.7206874) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(-0.044205753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71647787) q[0];
sx q[0];
rz(-1.2477861) q[0];
sx q[0];
rz(-0.72854211) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46385997) q[2];
sx q[2];
rz(-2.0630699) q[2];
sx q[2];
rz(-2.4459723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.468535) q[1];
sx q[1];
rz(-1.6516764) q[1];
sx q[1];
rz(-1.9616425) q[1];
rz(-pi) q[2];
rz(1.8474691) q[3];
sx q[3];
rz(-1.0545316) q[3];
sx q[3];
rz(0.54126213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1786903) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(0.94513354) q[2];
rz(1.6804228) q[3];
sx q[3];
rz(-1.9943359) q[3];
sx q[3];
rz(-0.51819658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1951676) q[0];
sx q[0];
rz(-0.31038809) q[0];
sx q[0];
rz(-0.84939605) q[0];
rz(-1.3769582) q[1];
sx q[1];
rz(-2.5701249) q[1];
sx q[1];
rz(1.8570541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97235128) q[0];
sx q[0];
rz(-1.961876) q[0];
sx q[0];
rz(-0.53000662) q[0];
x q[1];
rz(2.8148267) q[2];
sx q[2];
rz(-0.32354718) q[2];
sx q[2];
rz(1.7015333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4139288) q[1];
sx q[1];
rz(-1.5352121) q[1];
sx q[1];
rz(0.52558454) q[1];
rz(-pi) q[2];
rz(-2.2734145) q[3];
sx q[3];
rz(-0.68553151) q[3];
sx q[3];
rz(-3.1044416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66594243) q[2];
sx q[2];
rz(-1.7721756) q[2];
sx q[2];
rz(-1.5671889) q[2];
rz(2.808908) q[3];
sx q[3];
rz(-1.0654457) q[3];
sx q[3];
rz(0.98193297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5653) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(-1.0885619) q[0];
rz(-2.2125878) q[1];
sx q[1];
rz(-1.6477511) q[1];
sx q[1];
rz(-2.8378024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4171241) q[0];
sx q[0];
rz(-1.1982293) q[0];
sx q[0];
rz(2.5346816) q[0];
rz(-pi) q[1];
rz(3.0225176) q[2];
sx q[2];
rz(-1.4085253) q[2];
sx q[2];
rz(-1.281126) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15081295) q[1];
sx q[1];
rz(-0.44027281) q[1];
sx q[1];
rz(1.3516462) q[1];
rz(-pi) q[2];
rz(-0.63173826) q[3];
sx q[3];
rz(-2.4186385) q[3];
sx q[3];
rz(0.50448862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2527689) q[2];
sx q[2];
rz(-1.6656275) q[2];
sx q[2];
rz(-0.26407537) q[2];
rz(2.0090328) q[3];
sx q[3];
rz(-0.3370291) q[3];
sx q[3];
rz(-2.0866709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7842512) q[0];
sx q[0];
rz(-0.34264523) q[0];
sx q[0];
rz(-1.0145048) q[0];
rz(2.2134589) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(-0.49044213) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883915) q[0];
sx q[0];
rz(-1.7176796) q[0];
sx q[0];
rz(-3.0753972) q[0];
rz(-2.2955206) q[2];
sx q[2];
rz(-2.6451689) q[2];
sx q[2];
rz(-1.7299394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0874493) q[1];
sx q[1];
rz(-1.5296525) q[1];
sx q[1];
rz(-1.270134) q[1];
rz(-pi) q[2];
rz(2.1098299) q[3];
sx q[3];
rz(-1.7787053) q[3];
sx q[3];
rz(-0.62139213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32533112) q[2];
sx q[2];
rz(-1.0518495) q[2];
sx q[2];
rz(-3.0391147) q[2];
rz(1.8898194) q[3];
sx q[3];
rz(-1.3636369) q[3];
sx q[3];
rz(-1.4832835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935299) q[0];
sx q[0];
rz(-0.04638014) q[0];
sx q[0];
rz(-1.2903794) q[0];
rz(2.6875467) q[1];
sx q[1];
rz(-1.8171277) q[1];
sx q[1];
rz(-0.18347278) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.030122) q[0];
sx q[0];
rz(-0.18387499) q[0];
sx q[0];
rz(-1.6746111) q[0];
rz(-pi) q[1];
rz(-1.5857592) q[2];
sx q[2];
rz(-0.97285473) q[2];
sx q[2];
rz(0.19549616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.390229) q[1];
sx q[1];
rz(-1.5946663) q[1];
sx q[1];
rz(-2.6493401) q[1];
x q[2];
rz(1.4278891) q[3];
sx q[3];
rz(-0.57051728) q[3];
sx q[3];
rz(2.4874668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4219249) q[2];
sx q[2];
rz(-2.0823961) q[2];
sx q[2];
rz(-0.025040778) q[2];
rz(0.54343623) q[3];
sx q[3];
rz(-2.4380324) q[3];
sx q[3];
rz(-2.8652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.339879) q[0];
sx q[0];
rz(-2.1262953) q[0];
sx q[0];
rz(-2.3332818) q[0];
rz(-2.8080151) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(1.2110151) q[2];
sx q[2];
rz(-0.3508437) q[2];
sx q[2];
rz(2.99101) q[2];
rz(2.3357441) q[3];
sx q[3];
rz(-0.83559201) q[3];
sx q[3];
rz(1.5960485) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
