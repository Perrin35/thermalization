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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.08377) q[0];
sx q[0];
rz(-1.738213) q[0];
sx q[0];
rz(2.1517702) q[0];
rz(-pi) q[1];
rz(1.7389347) q[2];
sx q[2];
rz(-1.0285707) q[2];
sx q[2];
rz(-3.0410492) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3883512) q[1];
sx q[1];
rz(-1.8495535) q[1];
sx q[1];
rz(-0.21215794) q[1];
rz(-pi) q[2];
rz(-0.23363913) q[3];
sx q[3];
rz(-2.2147182) q[3];
sx q[3];
rz(-2.0104017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56403247) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(2.2300301) q[3];
sx q[3];
rz(-1.7715958) q[3];
sx q[3];
rz(0.024356775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0381222) q[0];
sx q[0];
rz(-2.7439674) q[0];
sx q[0];
rz(-3.0225515) q[0];
rz(1.1301522) q[1];
sx q[1];
rz(-0.96769133) q[1];
sx q[1];
rz(-1.4281323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16258612) q[0];
sx q[0];
rz(-1.9784728) q[0];
sx q[0];
rz(2.5752221) q[0];
x q[1];
rz(-0.4006673) q[2];
sx q[2];
rz(-1.9611437) q[2];
sx q[2];
rz(-1.1125178) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5609614) q[1];
sx q[1];
rz(-1.1079746) q[1];
sx q[1];
rz(2.9918519) q[1];
rz(-pi) q[2];
rz(0.30562206) q[3];
sx q[3];
rz(-1.247331) q[3];
sx q[3];
rz(-0.68438578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52593645) q[2];
sx q[2];
rz(-1.4639414) q[2];
sx q[2];
rz(-2.3770135) q[2];
rz(0.48314759) q[3];
sx q[3];
rz(-1.316148) q[3];
sx q[3];
rz(-0.84883261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1044384) q[0];
sx q[0];
rz(-2.8762682) q[0];
sx q[0];
rz(1.4720488) q[0];
rz(2.5813591) q[1];
sx q[1];
rz(-1.7543703) q[1];
sx q[1];
rz(-1.701042) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2587145) q[0];
sx q[0];
rz(-2.8322161) q[0];
sx q[0];
rz(0.86652722) q[0];
rz(-0.0071390634) q[2];
sx q[2];
rz(-2.4666365) q[2];
sx q[2];
rz(1.0085886) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8941108) q[1];
sx q[1];
rz(-1.2318582) q[1];
sx q[1];
rz(1.0147269) q[1];
rz(-1.3441003) q[3];
sx q[3];
rz(-1.3406702) q[3];
sx q[3];
rz(-3.074444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(-0.34240016) q[2];
rz(1.418669) q[3];
sx q[3];
rz(-0.72903052) q[3];
sx q[3];
rz(-0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2320084) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.6147856) q[0];
rz(2.3341663) q[1];
sx q[1];
rz(-1.5734943) q[1];
sx q[1];
rz(-1.9400914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9634092) q[0];
sx q[0];
rz(-0.82287517) q[0];
sx q[0];
rz(1.1567409) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2204578) q[2];
sx q[2];
rz(-1.0277205) q[2];
sx q[2];
rz(1.4411826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2953342) q[1];
sx q[1];
rz(-0.82678079) q[1];
sx q[1];
rz(-1.0215525) q[1];
rz(-pi) q[2];
rz(1.6225927) q[3];
sx q[3];
rz(-1.3148309) q[3];
sx q[3];
rz(1.3802176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3379007) q[2];
sx q[2];
rz(-2.1521229) q[2];
sx q[2];
rz(-1.2376415) q[2];
rz(1.3112274) q[3];
sx q[3];
rz(-1.5179736) q[3];
sx q[3];
rz(-1.1782014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.7685522) q[0];
sx q[0];
rz(0.9088687) q[0];
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
rz(3.1279739) q[0];
sx q[0];
rz(-2.2165856) q[0];
sx q[0];
rz(-0.60198204) q[0];
rz(0.2899828) q[2];
sx q[2];
rz(-2.3521543) q[2];
sx q[2];
rz(-2.0470999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4914507) q[1];
sx q[1];
rz(-0.60381266) q[1];
sx q[1];
rz(-0.73765386) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3084435) q[3];
sx q[3];
rz(-2.1705496) q[3];
sx q[3];
rz(-2.7540327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.070179209) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(1.2893691) q[2];
rz(-1.4687126) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656723) q[0];
sx q[0];
rz(-1.1812295) q[0];
sx q[0];
rz(-2.0400203) q[0];
rz(-0.22483243) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(-0.98519957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48822847) q[0];
sx q[0];
rz(-0.66500992) q[0];
sx q[0];
rz(1.2763766) q[0];
rz(1.9198138) q[2];
sx q[2];
rz(-0.98528457) q[2];
sx q[2];
rz(-1.7988009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0877062) q[1];
sx q[1];
rz(-0.83982491) q[1];
sx q[1];
rz(-0.62505109) q[1];
rz(-pi) q[2];
rz(3.0984663) q[3];
sx q[3];
rz(-0.52845983) q[3];
sx q[3];
rz(-0.915574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3605986) q[2];
sx q[2];
rz(-0.66086078) q[2];
sx q[2];
rz(-1.9471656) q[2];
rz(-0.099695168) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(0.97894198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.926696) q[0];
sx q[0];
rz(-2.844606) q[0];
sx q[0];
rz(1.7420705) q[0];
rz(2.4691441) q[2];
sx q[2];
rz(-2.932076) q[2];
sx q[2];
rz(2.5447951) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75668979) q[1];
sx q[1];
rz(-2.3339565) q[1];
sx q[1];
rz(-2.7212423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31734316) q[3];
sx q[3];
rz(-0.71094027) q[3];
sx q[3];
rz(1.8985572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.067001192) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(0.60603777) q[2];
rz(-2.0634985) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(1.7830474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9692877) q[0];
sx q[0];
rz(-2.2242039) q[0];
sx q[0];
rz(2.0148) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-2.7028449) q[1];
sx q[1];
rz(0.21805683) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15197922) q[0];
sx q[0];
rz(-1.1104062) q[0];
sx q[0];
rz(-0.48826852) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3534691) q[2];
sx q[2];
rz(-1.0255073) q[2];
sx q[2];
rz(2.9719549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3266284) q[1];
sx q[1];
rz(-1.6964629) q[1];
sx q[1];
rz(-2.2180422) q[1];
rz(-pi) q[2];
rz(1.7084684) q[3];
sx q[3];
rz(-1.4542504) q[3];
sx q[3];
rz(0.93695177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6286991) q[2];
sx q[2];
rz(-0.67535526) q[2];
sx q[2];
rz(-2.8524354) q[2];
rz(2.1469927) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(-0.4579671) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4299803) q[0];
sx q[0];
rz(-2.0531605) q[0];
sx q[0];
rz(2.3751538) q[0];
rz(1.6038766) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(-2.2235353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1963475) q[0];
sx q[0];
rz(-0.96789384) q[0];
sx q[0];
rz(2.5640998) q[0];
x q[1];
rz(-0.70785825) q[2];
sx q[2];
rz(-2.9212657) q[2];
sx q[2];
rz(-0.79263055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80528211) q[1];
sx q[1];
rz(-2.2754221) q[1];
sx q[1];
rz(1.2444088) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(0.90091577) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(-0.051699836) q[2];
rz(2.2895571) q[3];
sx q[3];
rz(-1.4308735) q[3];
sx q[3];
rz(2.8813072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.8155415) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(0.091212243) q[0];
rz(-2.3444029) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(-2.7104654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.618396) q[0];
sx q[0];
rz(3.0831227) q[0];
rz(-pi) q[1];
x q[1];
rz(0.098478949) q[2];
sx q[2];
rz(-0.92442552) q[2];
sx q[2];
rz(1.7913833) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7550274) q[1];
sx q[1];
rz(-1.3530827) q[1];
sx q[1];
rz(2.6854807) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1608716) q[3];
sx q[3];
rz(-2.4517165) q[3];
sx q[3];
rz(1.4341314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.572523) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(-1.265906) q[2];
rz(1.2931394) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(0.40614793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010392808) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(0.56397437) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(-3.0606506) q[2];
sx q[2];
rz(-1.1389393) q[2];
sx q[2];
rz(2.6280279) q[2];
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
