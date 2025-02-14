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
rz(-0.15260881) q[0];
sx q[0];
rz(5.0835291) q[0];
sx q[0];
rz(9.1191376) q[0];
rz(-2.8513554) q[1];
sx q[1];
rz(-1.4479535) q[1];
sx q[1];
rz(-1.0040959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0073402) q[0];
sx q[0];
rz(-1.5135845) q[0];
sx q[0];
rz(-1.5438616) q[0];
x q[1];
rz(2.9485116) q[2];
sx q[2];
rz(-1.3261137) q[2];
sx q[2];
rz(-2.8249403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6184147) q[1];
sx q[1];
rz(-2.6175099) q[1];
sx q[1];
rz(2.2487703) q[1];
rz(-0.070843971) q[3];
sx q[3];
rz(-0.53991441) q[3];
sx q[3];
rz(1.8749464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64757887) q[2];
sx q[2];
rz(-0.61654377) q[2];
sx q[2];
rz(2.5971863) q[2];
rz(-0.38862774) q[3];
sx q[3];
rz(-1.2231239) q[3];
sx q[3];
rz(0.94038928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87043864) q[0];
sx q[0];
rz(-2.1156023) q[0];
sx q[0];
rz(-2.0197268) q[0];
rz(2.0815966) q[1];
sx q[1];
rz(-2.6385939) q[1];
sx q[1];
rz(-0.88567919) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6933357) q[0];
sx q[0];
rz(-1.8739032) q[0];
sx q[0];
rz(1.2616377) q[0];
x q[1];
rz(-2.3212385) q[2];
sx q[2];
rz(-2.1848896) q[2];
sx q[2];
rz(0.49763735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52351511) q[1];
sx q[1];
rz(-1.4898058) q[1];
sx q[1];
rz(0.063132719) q[1];
rz(2.2996385) q[3];
sx q[3];
rz(-2.9837787) q[3];
sx q[3];
rz(-2.006881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9649967) q[2];
sx q[2];
rz(-1.6792038) q[2];
sx q[2];
rz(0.029335955) q[2];
rz(2.7756179) q[3];
sx q[3];
rz(-0.47593203) q[3];
sx q[3];
rz(0.16987814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40816864) q[0];
sx q[0];
rz(-1.3329196) q[0];
sx q[0];
rz(0.15342203) q[0];
rz(-0.19639213) q[1];
sx q[1];
rz(-1.6729665) q[1];
sx q[1];
rz(-0.82082716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9388537) q[0];
sx q[0];
rz(-1.0717055) q[0];
sx q[0];
rz(-2.3853154) q[0];
rz(-2.2024764) q[2];
sx q[2];
rz(-1.16776) q[2];
sx q[2];
rz(0.86954651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7559636) q[1];
sx q[1];
rz(-1.6106933) q[1];
sx q[1];
rz(-2.7261158) q[1];
x q[2];
rz(-0.12755022) q[3];
sx q[3];
rz(-0.3693499) q[3];
sx q[3];
rz(0.41778676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11744943) q[2];
sx q[2];
rz(-1.8154181) q[2];
sx q[2];
rz(-0.93309012) q[2];
rz(2.7945331) q[3];
sx q[3];
rz(-2.0324028) q[3];
sx q[3];
rz(-0.6412653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1953122) q[0];
sx q[0];
rz(-2.0144137) q[0];
sx q[0];
rz(-1.0585002) q[0];
rz(-0.095254101) q[1];
sx q[1];
rz(-1.9922099) q[1];
sx q[1];
rz(-0.20406318) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51523717) q[0];
sx q[0];
rz(-1.7254819) q[0];
sx q[0];
rz(1.3503287) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45943164) q[2];
sx q[2];
rz(-1.1387912) q[2];
sx q[2];
rz(2.9615796) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5285572) q[1];
sx q[1];
rz(-1.370346) q[1];
sx q[1];
rz(0.58581183) q[1];
rz(-pi) q[2];
rz(-2.8974523) q[3];
sx q[3];
rz(-2.4286037) q[3];
sx q[3];
rz(-2.0363852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.05693398) q[2];
sx q[2];
rz(-1.441842) q[2];
sx q[2];
rz(-1.8208246) q[2];
rz(-1.002958) q[3];
sx q[3];
rz(-1.2489677) q[3];
sx q[3];
rz(-2.5761719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4249307) q[0];
sx q[0];
rz(-1.2053763) q[0];
sx q[0];
rz(-2.5004814) q[0];
rz(-2.5504316) q[1];
sx q[1];
rz(-1.4407651) q[1];
sx q[1];
rz(1.444918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49316367) q[0];
sx q[0];
rz(-0.38721353) q[0];
sx q[0];
rz(3.0680543) q[0];
x q[1];
rz(0.88052184) q[2];
sx q[2];
rz(-1.7362927) q[2];
sx q[2];
rz(1.5269765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7406297) q[1];
sx q[1];
rz(-1.1006121) q[1];
sx q[1];
rz(-1.8484673) q[1];
x q[2];
rz(1.0743477) q[3];
sx q[3];
rz(-0.58234057) q[3];
sx q[3];
rz(2.989188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5833907) q[2];
sx q[2];
rz(-2.148197) q[2];
sx q[2];
rz(2.2885585) q[2];
rz(-2.9663626) q[3];
sx q[3];
rz(-1.3220738) q[3];
sx q[3];
rz(2.5186553) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22195062) q[0];
sx q[0];
rz(-2.3355244) q[0];
sx q[0];
rz(0.44575086) q[0];
rz(2.1469965) q[1];
sx q[1];
rz(-1.0044731) q[1];
sx q[1];
rz(1.6533143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22791187) q[0];
sx q[0];
rz(-0.74800038) q[0];
sx q[0];
rz(-2.8196003) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6449261) q[2];
sx q[2];
rz(-1.6152116) q[2];
sx q[2];
rz(-1.6629459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46785746) q[1];
sx q[1];
rz(-2.7731491) q[1];
sx q[1];
rz(2.175827) q[1];
rz(1.8976977) q[3];
sx q[3];
rz(-1.4964678) q[3];
sx q[3];
rz(-1.9063708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17066869) q[2];
sx q[2];
rz(-1.2430151) q[2];
sx q[2];
rz(-2.1155913) q[2];
rz(2.3435727) q[3];
sx q[3];
rz(-0.84037104) q[3];
sx q[3];
rz(0.26871267) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84665027) q[0];
sx q[0];
rz(-1.3947399) q[0];
sx q[0];
rz(0.41644874) q[0];
rz(-1.7448447) q[1];
sx q[1];
rz(-1.6114707) q[1];
sx q[1];
rz(1.4769311) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.653841) q[0];
sx q[0];
rz(-1.6076822) q[0];
sx q[0];
rz(2.1247618) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8139122) q[2];
sx q[2];
rz(-1.1736571) q[2];
sx q[2];
rz(-2.5458592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3116377) q[1];
sx q[1];
rz(-0.88407239) q[1];
sx q[1];
rz(-2.555189) q[1];
rz(-pi) q[2];
rz(1.6536413) q[3];
sx q[3];
rz(-1.9085064) q[3];
sx q[3];
rz(0.089525539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2485409) q[2];
sx q[2];
rz(-0.38022843) q[2];
sx q[2];
rz(0.60379544) q[2];
rz(0.90886146) q[3];
sx q[3];
rz(-1.7085608) q[3];
sx q[3];
rz(-1.0994256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73868442) q[0];
sx q[0];
rz(-1.383902) q[0];
sx q[0];
rz(0.72147328) q[0];
rz(1.3353222) q[1];
sx q[1];
rz(-2.0029009) q[1];
sx q[1];
rz(1.8849751) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2892814) q[0];
sx q[0];
rz(-0.86332488) q[0];
sx q[0];
rz(1.1723164) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0851791) q[2];
sx q[2];
rz(-2.0943421) q[2];
sx q[2];
rz(-0.26748891) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2366888) q[1];
sx q[1];
rz(-1.8042548) q[1];
sx q[1];
rz(-1.1462565) q[1];
x q[2];
rz(0.04560915) q[3];
sx q[3];
rz(-1.6001587) q[3];
sx q[3];
rz(1.9085397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98715544) q[2];
sx q[2];
rz(-0.64266959) q[2];
sx q[2];
rz(-1.8670234) q[2];
rz(0.68862033) q[3];
sx q[3];
rz(-1.6504811) q[3];
sx q[3];
rz(3.137818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0140728) q[0];
sx q[0];
rz(-1.2496244) q[0];
sx q[0];
rz(0.67181146) q[0];
rz(1.6067243) q[1];
sx q[1];
rz(-0.22767362) q[1];
sx q[1];
rz(0.52275503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76394698) q[0];
sx q[0];
rz(-1.3113759) q[0];
sx q[0];
rz(2.460145) q[0];
rz(-pi) q[1];
rz(-1.2272862) q[2];
sx q[2];
rz(-1.6186423) q[2];
sx q[2];
rz(-0.44105083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8485496) q[1];
sx q[1];
rz(-0.58401744) q[1];
sx q[1];
rz(-2.4021781) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32355752) q[3];
sx q[3];
rz(-1.3731628) q[3];
sx q[3];
rz(1.4909286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25040024) q[2];
sx q[2];
rz(-0.6468536) q[2];
sx q[2];
rz(0.35987535) q[2];
rz(-1.322809) q[3];
sx q[3];
rz(-1.6836932) q[3];
sx q[3];
rz(0.047164269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5322402) q[0];
sx q[0];
rz(-2.1873964) q[0];
sx q[0];
rz(0.94616079) q[0];
rz(3.1006475) q[1];
sx q[1];
rz(-1.2766726) q[1];
sx q[1];
rz(1.2650222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.301279) q[0];
sx q[0];
rz(-1.095533) q[0];
sx q[0];
rz(1.6874357) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6640915) q[2];
sx q[2];
rz(-0.5029808) q[2];
sx q[2];
rz(1.2625842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62502669) q[1];
sx q[1];
rz(-2.0989162) q[1];
sx q[1];
rz(-0.81104802) q[1];
x q[2];
rz(1.446883) q[3];
sx q[3];
rz(-2.765345) q[3];
sx q[3];
rz(0.28722426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1653183) q[2];
sx q[2];
rz(-1.9860622) q[2];
sx q[2];
rz(-2.9602642) q[2];
rz(-2.2881962) q[3];
sx q[3];
rz(-2.3386164) q[3];
sx q[3];
rz(1.8027423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637909) q[0];
sx q[0];
rz(-1.3297357) q[0];
sx q[0];
rz(1.7672675) q[0];
rz(1.455066) q[1];
sx q[1];
rz(-1.4604026) q[1];
sx q[1];
rz(-0.30053465) q[1];
rz(-2.5004432) q[2];
sx q[2];
rz(-1.1417626) q[2];
sx q[2];
rz(-0.7456197) q[2];
rz(-1.3131014) q[3];
sx q[3];
rz(-1.8492167) q[3];
sx q[3];
rz(1.3067018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
