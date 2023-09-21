OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(-0.95681325) q[0];
sx q[0];
rz(1.4332888) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(1.1448316) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0639609) q[0];
sx q[0];
rz(-1.6368757) q[0];
sx q[0];
rz(-1.8128916) q[0];
rz(-1.1510017) q[2];
sx q[2];
rz(-1.8325873) q[2];
sx q[2];
rz(0.37444886) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.634234) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(-1.6139612) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47031109) q[3];
sx q[3];
rz(-2.5458126) q[3];
sx q[3];
rz(2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-2.3117075) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(0.84709644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704437) q[0];
sx q[0];
rz(-2.0111472) q[0];
sx q[0];
rz(-2.0106993) q[0];
rz(-pi) q[1];
rz(0.77913021) q[2];
sx q[2];
rz(-1.7841633) q[2];
sx q[2];
rz(-1.2029635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21287316) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(0.54547711) q[1];
rz(0.91244016) q[3];
sx q[3];
rz(-0.53386253) q[3];
sx q[3];
rz(-2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
rz(1.8711202) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8910599) q[0];
sx q[0];
rz(-0.30936229) q[0];
sx q[0];
rz(-3.1191349) q[0];
rz(2.4685523) q[2];
sx q[2];
rz(-1.8290142) q[2];
sx q[2];
rz(-2.8845127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61422435) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(0.36348344) q[1];
x q[2];
rz(-1.8970044) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(-2.0446442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(2.4148338) q[2];
rz(2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.5184901) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-0.8262659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5465281) q[0];
sx q[0];
rz(-2.1273158) q[0];
sx q[0];
rz(2.3331649) q[0];
rz(-3.0777061) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(2.137616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(-0.080288447) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9598947) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(-2.1967595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-1.1337093) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782398) q[0];
sx q[0];
rz(-3.087128) q[0];
sx q[0];
rz(0.99476238) q[0];
x q[1];
rz(-1.1281625) q[2];
sx q[2];
rz(-2.1652522) q[2];
sx q[2];
rz(0.15304676) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87318201) q[1];
sx q[1];
rz(-2.1408484) q[1];
sx q[1];
rz(-0.33318712) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9660452) q[3];
sx q[3];
rz(-1.5834337) q[3];
sx q[3];
rz(-0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0952547) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(-2.6020223) q[0];
rz(0.18896582) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(-1.8408066) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.99299586) q[1];
sx q[1];
rz(-0.33278782) q[1];
sx q[1];
rz(-2.6246043) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6052386) q[3];
sx q[3];
rz(-2.648571) q[3];
sx q[3];
rz(2.8480414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95773762) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.6523223) q[0];
sx q[0];
rz(1.3315014) q[0];
x q[1];
rz(-0.046182403) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(-1.6398167) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0446339) q[1];
sx q[1];
rz(-0.63475906) q[1];
sx q[1];
rz(-1.8083014) q[1];
rz(-1.8764898) q[3];
sx q[3];
rz(-1.8743268) q[3];
sx q[3];
rz(0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0030901) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(1.2069758) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.5054024) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36639402) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(1.9154857) q[0];
x q[1];
rz(-0.82377388) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.6549695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-0.69680981) q[1];
sx q[1];
rz(-0.050683024) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1721341) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(-2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(1.6020417) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(0.64569965) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18683534) q[0];
sx q[0];
rz(-1.1823913) q[0];
sx q[0];
rz(0.79083058) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5876797) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(1.6162789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7449194) q[1];
sx q[1];
rz(-1.5349689) q[1];
sx q[1];
rz(-2.9958535) q[1];
rz(0.42616578) q[3];
sx q[3];
rz(-2.6821972) q[3];
sx q[3];
rz(0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.6773178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7510371) q[0];
sx q[0];
rz(-1.3532191) q[0];
sx q[0];
rz(1.6839954) q[0];
rz(-pi) q[1];
rz(-0.22428959) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(-1.6155417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36832419) q[1];
sx q[1];
rz(-1.3001633) q[1];
sx q[1];
rz(-2.2776105) q[1];
x q[2];
rz(3.0562923) q[3];
sx q[3];
rz(-1.395426) q[3];
sx q[3];
rz(0.264197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-2.0132813) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-0.55143572) q[2];
sx q[2];
rz(-0.79179344) q[2];
sx q[2];
rz(1.0168016) q[2];
rz(2.5682156) q[3];
sx q[3];
rz(-1.2497414) q[3];
sx q[3];
rz(-2.9336815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
