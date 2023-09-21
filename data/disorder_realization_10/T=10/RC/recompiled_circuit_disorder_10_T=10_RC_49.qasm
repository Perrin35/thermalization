OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(-1.9300652) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053781833) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(1.5009297) q[0];
x q[1];
rz(1.3053427) q[2];
sx q[2];
rz(-0.34398088) q[2];
sx q[2];
rz(1.6814107) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28042291) q[1];
sx q[1];
rz(-2.0457595) q[1];
sx q[1];
rz(-2.1850987) q[1];
rz(-0.88793036) q[3];
sx q[3];
rz(-3.008932) q[3];
sx q[3];
rz(-2.5761029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(2.1477264) q[2];
rz(2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(-2.6699064) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90549201) q[0];
sx q[0];
rz(-0.24370757) q[0];
sx q[0];
rz(-0.36867504) q[0];
rz(-pi) q[1];
rz(-1.291044) q[2];
sx q[2];
rz(-1.0824167) q[2];
sx q[2];
rz(0.56088698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7447409) q[1];
sx q[1];
rz(-0.62780118) q[1];
sx q[1];
rz(0.46698924) q[1];
rz(-pi) q[2];
rz(0.15700335) q[3];
sx q[3];
rz(-0.96884851) q[3];
sx q[3];
rz(2.4412145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5919684) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(0.96898752) q[2];
rz(2.5668868) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(-2.95978) q[0];
rz(-1.1026985) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(-1.4556494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239379) q[0];
sx q[0];
rz(-1.0929937) q[0];
sx q[0];
rz(2.2295203) q[0];
rz(-0.74080148) q[2];
sx q[2];
rz(-1.891279) q[2];
sx q[2];
rz(-1.4893116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1118483) q[1];
sx q[1];
rz(-1.6964456) q[1];
sx q[1];
rz(1.5602342) q[1];
rz(-0.62526838) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(2.3141935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-0.96763119) q[2];
rz(-2.4140221) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(2.0137285) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.2329873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013279182) q[0];
sx q[0];
rz(-2.0648801) q[0];
sx q[0];
rz(-1.4923151) q[0];
rz(-pi) q[1];
rz(-1.0225251) q[2];
sx q[2];
rz(-2.0892482) q[2];
sx q[2];
rz(0.099629121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7113263) q[1];
sx q[1];
rz(-1.749199) q[1];
sx q[1];
rz(3.1049411) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2635157) q[3];
sx q[3];
rz(-1.746776) q[3];
sx q[3];
rz(0.39719492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91810742) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-2.510575) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994613) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(-1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(0.2125425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4310303) q[0];
sx q[0];
rz(-3.115603) q[0];
sx q[0];
rz(-0.89528577) q[0];
x q[1];
rz(-1.5151305) q[2];
sx q[2];
rz(-0.60534436) q[2];
sx q[2];
rz(2.2948613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5571641) q[1];
sx q[1];
rz(-0.85313988) q[1];
sx q[1];
rz(2.3889956) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84380031) q[3];
sx q[3];
rz(-1.2146597) q[3];
sx q[3];
rz(-0.86405495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(-0.5212211) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(1.8732171) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.3549995) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(2.7640142) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8416653) q[0];
sx q[0];
rz(-2.8746434) q[0];
sx q[0];
rz(1.6542692) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2874561) q[2];
sx q[2];
rz(-1.524458) q[2];
sx q[2];
rz(3.1099144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0203591) q[1];
sx q[1];
rz(-1.6755783) q[1];
sx q[1];
rz(2.0038414) q[1];
rz(0.029929786) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(-0.078775725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(0.48103508) q[2];
rz(2.7379819) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(-0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(0.21952195) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(0.25442466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1628111) q[0];
sx q[0];
rz(-0.33919507) q[0];
sx q[0];
rz(-1.2374872) q[0];
rz(2.1830325) q[2];
sx q[2];
rz(-1.4441635) q[2];
sx q[2];
rz(-2.3938092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2563045) q[1];
sx q[1];
rz(-1.2997775) q[1];
sx q[1];
rz(0.46896743) q[1];
rz(-pi) q[2];
rz(-0.45215996) q[3];
sx q[3];
rz(-1.1434571) q[3];
sx q[3];
rz(-0.51963453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.8257726) q[2];
rz(2.2655462) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(2.1379437) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5664068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67714171) q[0];
sx q[0];
rz(-1.8832708) q[0];
sx q[0];
rz(2.421445) q[0];
x q[1];
rz(1.2775531) q[2];
sx q[2];
rz(-0.55985057) q[2];
sx q[2];
rz(0.44550371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6187001) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(-0.79843847) q[1];
rz(-pi) q[2];
rz(1.1493707) q[3];
sx q[3];
rz(-2.0391658) q[3];
sx q[3];
rz(-2.3295662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(-2.8640462) q[2];
rz(1.4510441) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(-2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16185109) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(1.45654) q[0];
rz(-0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(-1.1368407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4515848) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(2.8103229) q[0];
rz(-0.85176977) q[2];
sx q[2];
rz(-2.6034144) q[2];
sx q[2];
rz(-2.1449508) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1946823) q[1];
sx q[1];
rz(-1.2543251) q[1];
sx q[1];
rz(3.0419993) q[1];
rz(-pi) q[2];
rz(-1.198248) q[3];
sx q[3];
rz(-0.87009831) q[3];
sx q[3];
rz(-2.8230599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4920766) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(0.31203976) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-0.94605207) q[0];
rz(-2.229915) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8004868) q[0];
sx q[0];
rz(-2.2414811) q[0];
sx q[0];
rz(0.97408803) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.811118) q[2];
sx q[2];
rz(-1.9069306) q[2];
sx q[2];
rz(-1.0898332) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8814197) q[1];
sx q[1];
rz(-0.91671413) q[1];
sx q[1];
rz(1.4152848) q[1];
x q[2];
rz(-1.9480115) q[3];
sx q[3];
rz(-1.0672788) q[3];
sx q[3];
rz(2.6138888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0570021) q[2];
sx q[2];
rz(-2.4995063) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(-0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6012797) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(-2.3868949) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(2.1379708) q[2];
sx q[2];
rz(-2.1688609) q[2];
sx q[2];
rz(-1.4458956) q[2];
rz(-1.8176953) q[3];
sx q[3];
rz(-0.3331475) q[3];
sx q[3];
rz(-3.0039136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
