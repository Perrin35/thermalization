OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.328182) q[0];
sx q[0];
rz(5.2433104) q[0];
sx q[0];
rz(5.1007895) q[0];
rz(0.27529588) q[1];
sx q[1];
rz(7.2841865) q[1];
sx q[1];
rz(10.871973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49791928) q[0];
sx q[0];
rz(-0.6090954) q[0];
sx q[0];
rz(-0.66225556) q[0];
rz(-pi) q[1];
rz(-2.2859226) q[2];
sx q[2];
rz(-1.8026094) q[2];
sx q[2];
rz(2.135853) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1684678) q[1];
sx q[1];
rz(-0.68328062) q[1];
sx q[1];
rz(-0.82078187) q[1];
x q[2];
rz(-2.1547444) q[3];
sx q[3];
rz(-0.60923558) q[3];
sx q[3];
rz(1.9690994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64624661) q[2];
sx q[2];
rz(-2.4301811) q[2];
sx q[2];
rz(-1.6515674) q[2];
rz(-1.9735769) q[3];
sx q[3];
rz(-1.370859) q[3];
sx q[3];
rz(1.5994387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9769576) q[0];
sx q[0];
rz(-0.63666207) q[0];
sx q[0];
rz(-2.013999) q[0];
rz(-0.98059869) q[1];
sx q[1];
rz(-0.53739986) q[1];
sx q[1];
rz(-1.604039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1297596) q[0];
sx q[0];
rz(-0.3922222) q[0];
sx q[0];
rz(1.5480199) q[0];
x q[1];
rz(-1.3959742) q[2];
sx q[2];
rz(-3.0478301) q[2];
sx q[2];
rz(-1.7435058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.6863197) q[1];
sx q[1];
rz(-1.4497633) q[1];
sx q[1];
rz(-1.0334121) q[1];
rz(-pi) q[2];
rz(1.2885154) q[3];
sx q[3];
rz(-1.936463) q[3];
sx q[3];
rz(2.7136841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.299122) q[2];
sx q[2];
rz(-2.292558) q[2];
sx q[2];
rz(-0.15057286) q[2];
rz(-1.8922837) q[3];
sx q[3];
rz(-2.1842897) q[3];
sx q[3];
rz(-1.5875491) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4318749) q[0];
sx q[0];
rz(-1.4282325) q[0];
sx q[0];
rz(3.0336483) q[0];
rz(0.81186324) q[1];
sx q[1];
rz(-0.64473647) q[1];
sx q[1];
rz(2.937584) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0013179) q[0];
sx q[0];
rz(-1.5840127) q[0];
sx q[0];
rz(0.93854882) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7008867) q[2];
sx q[2];
rz(-1.9153898) q[2];
sx q[2];
rz(-1.454892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94212342) q[1];
sx q[1];
rz(-0.55049039) q[1];
sx q[1];
rz(-1.8041736) q[1];
rz(-pi) q[2];
rz(2.3164204) q[3];
sx q[3];
rz(-1.3047855) q[3];
sx q[3];
rz(1.1893742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72309368) q[2];
sx q[2];
rz(-2.8321224) q[2];
sx q[2];
rz(1.2014368) q[2];
rz(-0.57032436) q[3];
sx q[3];
rz(-1.0446965) q[3];
sx q[3];
rz(0.23049155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0975321) q[0];
sx q[0];
rz(-0.965913) q[0];
sx q[0];
rz(3.0934546) q[0];
rz(-0.33787456) q[1];
sx q[1];
rz(-0.66403762) q[1];
sx q[1];
rz(-2.7545676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95985896) q[0];
sx q[0];
rz(-0.67200586) q[0];
sx q[0];
rz(1.9197587) q[0];
rz(2.0104644) q[2];
sx q[2];
rz(-2.2626468) q[2];
sx q[2];
rz(0.95702167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6070575) q[1];
sx q[1];
rz(-1.9448154) q[1];
sx q[1];
rz(0.97403557) q[1];
rz(-pi) q[2];
rz(-0.13366661) q[3];
sx q[3];
rz(-1.2959305) q[3];
sx q[3];
rz(0.2256338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7016697) q[2];
sx q[2];
rz(-0.54653007) q[2];
sx q[2];
rz(0.92865357) q[2];
rz(-0.52740151) q[3];
sx q[3];
rz(-1.6583859) q[3];
sx q[3];
rz(-1.5206913) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034390282) q[0];
sx q[0];
rz(-0.90029383) q[0];
sx q[0];
rz(-0.1874371) q[0];
rz(1.772359) q[1];
sx q[1];
rz(-1.0405468) q[1];
sx q[1];
rz(-1.7231411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56031449) q[0];
sx q[0];
rz(-2.1933103) q[0];
sx q[0];
rz(0.050943663) q[0];
rz(0.64308738) q[2];
sx q[2];
rz(-1.4355112) q[2];
sx q[2];
rz(0.39150086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.975357) q[1];
sx q[1];
rz(-1.72962) q[1];
sx q[1];
rz(1.3162196) q[1];
rz(0.62331919) q[3];
sx q[3];
rz(-1.5006362) q[3];
sx q[3];
rz(2.1753785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9886542) q[2];
sx q[2];
rz(-0.31854892) q[2];
sx q[2];
rz(0.62206507) q[2];
rz(-2.7234744) q[3];
sx q[3];
rz(-1.6303635) q[3];
sx q[3];
rz(-2.6869584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.76573265) q[0];
sx q[0];
rz(-1.8143761) q[0];
sx q[0];
rz(2.3151929) q[0];
rz(-1.3580258) q[1];
sx q[1];
rz(-2.2260901) q[1];
sx q[1];
rz(0.079039097) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86980692) q[0];
sx q[0];
rz(-0.1087063) q[0];
sx q[0];
rz(-0.73819091) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68361552) q[2];
sx q[2];
rz(-1.9208462) q[2];
sx q[2];
rz(-1.8253714) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0519281) q[1];
sx q[1];
rz(-0.66510495) q[1];
sx q[1];
rz(-1.0699349) q[1];
rz(-pi) q[2];
rz(-2.6812234) q[3];
sx q[3];
rz(-1.7475024) q[3];
sx q[3];
rz(-1.1713771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0283811) q[2];
sx q[2];
rz(-1.6190642) q[2];
sx q[2];
rz(-1.9436504) q[2];
rz(-1.0269264) q[3];
sx q[3];
rz(-2.9668861) q[3];
sx q[3];
rz(-0.13897707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748264) q[0];
sx q[0];
rz(-2.3737895) q[0];
sx q[0];
rz(0.36198947) q[0];
rz(-2.844574) q[1];
sx q[1];
rz(-1.2535008) q[1];
sx q[1];
rz(-2.6716935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7720761) q[0];
sx q[0];
rz(-2.2945991) q[0];
sx q[0];
rz(-0.93627255) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1761492) q[2];
sx q[2];
rz(-1.5684853) q[2];
sx q[2];
rz(2.1049079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5928157) q[1];
sx q[1];
rz(-1.5783184) q[1];
sx q[1];
rz(-1.8846017) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.04977241) q[3];
sx q[3];
rz(-0.94954606) q[3];
sx q[3];
rz(1.5240108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.517259) q[2];
sx q[2];
rz(-2.2225311) q[2];
sx q[2];
rz(-2.2993235) q[2];
rz(1.9159348) q[3];
sx q[3];
rz(-1.3282447) q[3];
sx q[3];
rz(1.8153502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4351742) q[0];
sx q[0];
rz(-1.5442727) q[0];
sx q[0];
rz(0.73961863) q[0];
rz(0.30532062) q[1];
sx q[1];
rz(-1.0605597) q[1];
sx q[1];
rz(-1.7230497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.212767) q[0];
sx q[0];
rz(-0.97090844) q[0];
sx q[0];
rz(1.3296574) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3408634) q[2];
sx q[2];
rz(-0.62494102) q[2];
sx q[2];
rz(2.2228624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3736014) q[1];
sx q[1];
rz(-1.8680675) q[1];
sx q[1];
rz(0.41042787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3548325) q[3];
sx q[3];
rz(-1.4359754) q[3];
sx q[3];
rz(2.1444837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6366987) q[2];
sx q[2];
rz(-0.86504522) q[2];
sx q[2];
rz(-1.0587586) q[2];
rz(-0.59075683) q[3];
sx q[3];
rz(-2.1754913) q[3];
sx q[3];
rz(3.0920715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010794086) q[0];
sx q[0];
rz(-0.32242355) q[0];
sx q[0];
rz(1.7266493) q[0];
rz(2.2825799) q[1];
sx q[1];
rz(-2.450727) q[1];
sx q[1];
rz(-2.1376999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0738132) q[0];
sx q[0];
rz(-1.8260024) q[0];
sx q[0];
rz(-1.4007481) q[0];
rz(-pi) q[1];
rz(-2.48433) q[2];
sx q[2];
rz(-1.8397775) q[2];
sx q[2];
rz(1.2439234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2768663) q[1];
sx q[1];
rz(-1.4819615) q[1];
sx q[1];
rz(2.4125742) q[1];
rz(-pi) q[2];
rz(1.4850158) q[3];
sx q[3];
rz(-0.82786938) q[3];
sx q[3];
rz(-1.08094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88177195) q[2];
sx q[2];
rz(-1.5871781) q[2];
sx q[2];
rz(0.47015321) q[2];
rz(-1.638089) q[3];
sx q[3];
rz(-0.94109002) q[3];
sx q[3];
rz(2.9269384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.692602) q[0];
sx q[0];
rz(-0.34994352) q[0];
sx q[0];
rz(1.0119447) q[0];
rz(0.28255209) q[1];
sx q[1];
rz(-2.2239182) q[1];
sx q[1];
rz(0.73175398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35039038) q[0];
sx q[0];
rz(-1.6178432) q[0];
sx q[0];
rz(-1.5806517) q[0];
x q[1];
rz(-0.75586478) q[2];
sx q[2];
rz(-0.24768344) q[2];
sx q[2];
rz(0.29049722) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2866757) q[1];
sx q[1];
rz(-2.1647506) q[1];
sx q[1];
rz(-2.0718365) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.326272) q[3];
sx q[3];
rz(-0.36416194) q[3];
sx q[3];
rz(0.18362602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4628576) q[2];
sx q[2];
rz(-1.5827468) q[2];
sx q[2];
rz(-2.8850214) q[2];
rz(-2.3738142) q[3];
sx q[3];
rz(-1.0662181) q[3];
sx q[3];
rz(-2.0496875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0302122) q[0];
sx q[0];
rz(-0.82652265) q[0];
sx q[0];
rz(-0.65435456) q[0];
rz(-0.70814537) q[1];
sx q[1];
rz(-1.000052) q[1];
sx q[1];
rz(2.4194385) q[1];
rz(-1.234645) q[2];
sx q[2];
rz(-0.79094255) q[2];
sx q[2];
rz(-3.0002181) q[2];
rz(1.5453399) q[3];
sx q[3];
rz(-0.40332356) q[3];
sx q[3];
rz(-1.0298503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
