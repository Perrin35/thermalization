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
rz(1.8737268) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(2.3550912) q[1];
sx q[1];
rz(13.785706) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76152626) q[0];
sx q[0];
rz(-2.5396569) q[0];
sx q[0];
rz(-1.2720889) q[0];
rz(1.4026579) q[2];
sx q[2];
rz(-2.113022) q[2];
sx q[2];
rz(0.10054345) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75324149) q[1];
sx q[1];
rz(-1.8495535) q[1];
sx q[1];
rz(-2.9294347) q[1];
x q[2];
rz(1.8699617) q[3];
sx q[3];
rz(-2.4623021) q[3];
sx q[3];
rz(-2.3878179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5775602) q[2];
sx q[2];
rz(-1.5215678) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(2.2300301) q[3];
sx q[3];
rz(-1.7715958) q[3];
sx q[3];
rz(-3.1172359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9790065) q[0];
sx q[0];
rz(-1.1631199) q[0];
sx q[0];
rz(-2.5752221) q[0];
rz(-0.4006673) q[2];
sx q[2];
rz(-1.9611437) q[2];
sx q[2];
rz(2.0290749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88682879) q[1];
sx q[1];
rz(-0.48476754) q[1];
sx q[1];
rz(1.280275) q[1];
rz(-2.8359706) q[3];
sx q[3];
rz(-1.8942617) q[3];
sx q[3];
rz(0.68438578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6156562) q[2];
sx q[2];
rz(-1.4639414) q[2];
sx q[2];
rz(2.3770135) q[2];
rz(2.6584451) q[3];
sx q[3];
rz(-1.316148) q[3];
sx q[3];
rz(-2.29276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0371542) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(-1.4720488) q[0];
rz(2.5813591) q[1];
sx q[1];
rz(-1.7543703) q[1];
sx q[1];
rz(1.4405506) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.773303) q[0];
sx q[0];
rz(-1.3723626) q[0];
sx q[0];
rz(-1.8097359) q[0];
rz(-pi) q[1];
rz(-0.67494377) q[2];
sx q[2];
rz(-1.5663354) q[2];
sx q[2];
rz(0.5677815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3095197) q[1];
sx q[1];
rz(-0.64180556) q[1];
sx q[1];
rz(2.1596396) q[1];
x q[2];
rz(-2.3769242) q[3];
sx q[3];
rz(-2.8199785) q[3];
sx q[3];
rz(-0.72383108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(-0.34240016) q[2];
rz(-1.7229236) q[3];
sx q[3];
rz(-0.72903052) q[3];
sx q[3];
rz(-0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320084) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.5268071) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5680983) q[1];
sx q[1];
rz(-1.9400914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9634092) q[0];
sx q[0];
rz(-2.3187175) q[0];
sx q[0];
rz(1.1567409) q[0];
x q[1];
rz(-2.3551201) q[2];
sx q[2];
rz(-0.82068372) q[2];
sx q[2];
rz(2.414626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11749841) q[1];
sx q[1];
rz(-1.9650241) q[1];
sx q[1];
rz(0.82347639) q[1];
rz(1.5189999) q[3];
sx q[3];
rz(-1.8267617) q[3];
sx q[3];
rz(-1.761375) q[3];
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
rz(1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0230947) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(2.232724) q[0];
rz(2.2591649) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(2.4095101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2765512) q[0];
sx q[0];
rz(-0.85231977) q[0];
sx q[0];
rz(-0.9263692) q[0];
x q[1];
rz(2.3734748) q[2];
sx q[2];
rz(-1.7752194) q[2];
sx q[2];
rz(0.26917514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8632224) q[1];
sx q[1];
rz(-1.9626106) q[1];
sx q[1];
rz(-2.6696221) q[1];
x q[2];
rz(2.3957) q[3];
sx q[3];
rz(-2.159366) q[3];
sx q[3];
rz(2.4323127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.070179209) q[2];
sx q[2];
rz(-2.3290403) q[2];
sx q[2];
rz(-1.8522235) q[2];
rz(-1.4687126) q[3];
sx q[3];
rz(-1.2082992) q[3];
sx q[3];
rz(-2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5759204) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(-1.1015724) q[0];
rz(-2.9167602) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(0.98519957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12041872) q[0];
sx q[0];
rz(-0.93909953) q[0];
sx q[0];
rz(-2.9178502) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9198138) q[2];
sx q[2];
rz(-2.1563081) q[2];
sx q[2];
rz(1.7988009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0538865) q[1];
sx q[1];
rz(-2.3017677) q[1];
sx q[1];
rz(2.5165416) q[1];
rz(-1.5959625) q[3];
sx q[3];
rz(-2.0987134) q[3];
sx q[3];
rz(2.2759469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3605986) q[2];
sx q[2];
rz(-2.4807319) q[2];
sx q[2];
rz(1.194427) q[2];
rz(3.0418975) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(-2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128971) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(-2.0275443) q[0];
rz(-1.3975337) q[1];
sx q[1];
rz(-1.3797398) q[1];
sx q[1];
rz(2.8278415) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.926696) q[0];
sx q[0];
rz(-0.29698661) q[0];
sx q[0];
rz(1.7420705) q[0];
rz(1.4391104) q[2];
sx q[2];
rz(-1.4073616) q[2];
sx q[2];
rz(-1.2800467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9589428) q[1];
sx q[1];
rz(-2.291276) q[1];
sx q[1];
rz(1.1675323) q[1];
rz(-pi) q[2];
rz(0.68571949) q[3];
sx q[3];
rz(-1.3657394) q[3];
sx q[3];
rz(-3.0577537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0745915) q[2];
sx q[2];
rz(-0.40412298) q[2];
sx q[2];
rz(-2.5355549) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9692877) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(1.1267927) q[0];
rz(-2.2276095) q[1];
sx q[1];
rz(-2.7028449) q[1];
sx q[1];
rz(2.9235358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954531) q[0];
sx q[0];
rz(-2.0045223) q[0];
sx q[0];
rz(1.0591255) q[0];
x q[1];
rz(2.3534691) q[2];
sx q[2];
rz(-2.1160853) q[2];
sx q[2];
rz(0.1696378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7331024) q[1];
sx q[1];
rz(-0.6576076) q[1];
sx q[1];
rz(1.7773184) q[1];
x q[2];
rz(-1.7084684) q[3];
sx q[3];
rz(-1.6873423) q[3];
sx q[3];
rz(-2.2046409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-2.4662374) q[2];
sx q[2];
rz(2.8524354) q[2];
rz(-0.9946) q[3];
sx q[3];
rz(-1.9528439) q[3];
sx q[3];
rz(-2.6836256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7116123) q[0];
sx q[0];
rz(-1.0884322) q[0];
sx q[0];
rz(-0.76643884) q[0];
rz(-1.537716) q[1];
sx q[1];
rz(-1.3130554) q[1];
sx q[1];
rz(-0.91805735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.799918) q[0];
sx q[0];
rz(-2.3324488) q[0];
sx q[0];
rz(2.2412712) q[0];
x q[1];
rz(2.4337344) q[2];
sx q[2];
rz(-0.22032693) q[2];
sx q[2];
rz(0.79263055) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80528211) q[1];
sx q[1];
rz(-0.86617058) q[1];
sx q[1];
rz(1.8971838) q[1];
rz(-1.3206743) q[3];
sx q[3];
rz(-1.4682574) q[3];
sx q[3];
rz(0.96270442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2406769) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(0.051699836) q[2];
rz(-2.2895571) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(2.8813072) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32605115) q[0];
sx q[0];
rz(-1.6764078) q[0];
sx q[0];
rz(3.0503804) q[0];
rz(-0.7971898) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(-0.43112722) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9190311) q[0];
sx q[0];
rz(-3.0662144) q[0];
sx q[0];
rz(-0.68392174) q[0];
rz(-1.7003785) q[2];
sx q[2];
rz(-0.65276546) q[2];
sx q[2];
rz(-1.5127986) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3865652) q[1];
sx q[1];
rz(-1.7885099) q[1];
sx q[1];
rz(-2.6854807) q[1];
x q[2];
rz(2.1608716) q[3];
sx q[3];
rz(-0.68987615) q[3];
sx q[3];
rz(-1.4341314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.572523) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(-1.8756867) q[2];
rz(1.8484533) q[3];
sx q[3];
rz(-2.0796516) q[3];
sx q[3];
rz(0.40614793) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1311998) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(-0.56397437) q[1];
sx q[1];
rz(-1.0697983) q[1];
sx q[1];
rz(-1.0615798) q[1];
rz(-0.080942091) q[2];
sx q[2];
rz(-2.0026534) q[2];
sx q[2];
rz(-0.51356471) q[2];
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
