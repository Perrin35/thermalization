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
rz(1.9749405) q[0];
sx q[0];
rz(-2.6986172) q[0];
sx q[0];
rz(0.38323453) q[0];
rz(-1.2869599) q[1];
sx q[1];
rz(-1.1416963) q[1];
sx q[1];
rz(-0.72007522) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543862) q[0];
sx q[0];
rz(-1.1875679) q[0];
sx q[0];
rz(-2.711074) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7699446) q[2];
sx q[2];
rz(-1.3825644) q[2];
sx q[2];
rz(0.07478274) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9125058) q[1];
sx q[1];
rz(-0.47059637) q[1];
sx q[1];
rz(2.5060593) q[1];
x q[2];
rz(2.1813857) q[3];
sx q[3];
rz(-1.7060602) q[3];
sx q[3];
rz(-1.3162515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8522475) q[2];
sx q[2];
rz(-0.87163681) q[2];
sx q[2];
rz(2.3343425) q[2];
rz(-1.6926258) q[3];
sx q[3];
rz(-2.2117386) q[3];
sx q[3];
rz(0.71678954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221231) q[0];
sx q[0];
rz(-1.2998281) q[0];
sx q[0];
rz(1.3053373) q[0];
rz(-2.3985825) q[1];
sx q[1];
rz(-2.4842255) q[1];
sx q[1];
rz(1.8441127) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8332307) q[0];
sx q[0];
rz(-1.9675281) q[0];
sx q[0];
rz(-0.19773592) q[0];
rz(-pi) q[1];
rz(-0.99886151) q[2];
sx q[2];
rz(-1.1872269) q[2];
sx q[2];
rz(0.072173031) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0540794) q[1];
sx q[1];
rz(-0.98615188) q[1];
sx q[1];
rz(-2.784009) q[1];
rz(-pi) q[2];
rz(-2.7801974) q[3];
sx q[3];
rz(-1.4618498) q[3];
sx q[3];
rz(1.4055085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3798736) q[2];
sx q[2];
rz(-2.1642809) q[2];
sx q[2];
rz(-1.0404111) q[2];
rz(0.31446332) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(1.2578806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8085025) q[0];
sx q[0];
rz(-2.1603778) q[0];
sx q[0];
rz(1.3370978) q[0];
rz(-0.62726504) q[1];
sx q[1];
rz(-1.4074872) q[1];
sx q[1];
rz(-1.6076535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1264135) q[0];
sx q[0];
rz(-0.94278383) q[0];
sx q[0];
rz(-2.6005484) q[0];
rz(-pi) q[1];
rz(-2.7159458) q[2];
sx q[2];
rz(-1.8432277) q[2];
sx q[2];
rz(2.6657588) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2308637) q[1];
sx q[1];
rz(-1.7553991) q[1];
sx q[1];
rz(2.6631825) q[1];
rz(1.1008784) q[3];
sx q[3];
rz(-0.4199314) q[3];
sx q[3];
rz(0.41929276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4072121) q[2];
sx q[2];
rz(-0.89443365) q[2];
sx q[2];
rz(-1.2073995) q[2];
rz(0.95737988) q[3];
sx q[3];
rz(-1.9060241) q[3];
sx q[3];
rz(-0.68937075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6974238) q[0];
sx q[0];
rz(-0.95132315) q[0];
sx q[0];
rz(-3.0498258) q[0];
rz(0.16608876) q[1];
sx q[1];
rz(-1.2867915) q[1];
sx q[1];
rz(2.1002358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9516141) q[0];
sx q[0];
rz(-2.3654177) q[0];
sx q[0];
rz(1.5086195) q[0];
rz(-2.3538386) q[2];
sx q[2];
rz(-1.8623036) q[2];
sx q[2];
rz(-1.1109835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77063054) q[1];
sx q[1];
rz(-0.47415552) q[1];
sx q[1];
rz(1.4956135) q[1];
rz(0.54817537) q[3];
sx q[3];
rz(-2.976307) q[3];
sx q[3];
rz(-2.838244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.912821) q[2];
sx q[2];
rz(-2.3232338) q[2];
sx q[2];
rz(1.6010326) q[2];
rz(-2.1238964) q[3];
sx q[3];
rz(-1.8120268) q[3];
sx q[3];
rz(1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-1.3200662) q[0];
sx q[0];
rz(-1.7463266) q[0];
sx q[0];
rz(-2.9030002) q[0];
rz(0.78856167) q[1];
sx q[1];
rz(-2.8193654) q[1];
sx q[1];
rz(-2.2408748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64553211) q[0];
sx q[0];
rz(-1.1405611) q[0];
sx q[0];
rz(-1.2853464) q[0];
rz(0.094036799) q[2];
sx q[2];
rz(-1.9584736) q[2];
sx q[2];
rz(2.2212818) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2701931) q[1];
sx q[1];
rz(-2.7084623) q[1];
sx q[1];
rz(0.30850839) q[1];
x q[2];
rz(2.2969518) q[3];
sx q[3];
rz(-0.9880522) q[3];
sx q[3];
rz(-0.83707367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0939402) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(1.368604) q[2];
rz(-0.37745825) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(-1.1641758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0193598) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(2.9121616) q[0];
rz(0.58652985) q[1];
sx q[1];
rz(-0.49085453) q[1];
sx q[1];
rz(1.6612926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585715) q[0];
sx q[0];
rz(-0.93556306) q[0];
sx q[0];
rz(-0.44376377) q[0];
rz(-2.7469559) q[2];
sx q[2];
rz(-2.5355314) q[2];
sx q[2];
rz(0.74575033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52119285) q[1];
sx q[1];
rz(-1.7971276) q[1];
sx q[1];
rz(-1.5160376) q[1];
x q[2];
rz(2.924891) q[3];
sx q[3];
rz(-2.2335972) q[3];
sx q[3];
rz(-0.13706438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71969879) q[2];
sx q[2];
rz(-0.59013683) q[2];
sx q[2];
rz(1.3479007) q[2];
rz(-0.19139309) q[3];
sx q[3];
rz(-2.8189711) q[3];
sx q[3];
rz(-1.9026683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58452463) q[0];
sx q[0];
rz(-1.6275591) q[0];
sx q[0];
rz(2.6051482) q[0];
rz(0.084511936) q[1];
sx q[1];
rz(-0.5785431) q[1];
sx q[1];
rz(-1.1892148) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093350323) q[0];
sx q[0];
rz(-0.95365253) q[0];
sx q[0];
rz(2.9523152) q[0];
rz(1.1841449) q[2];
sx q[2];
rz(-2.2154567) q[2];
sx q[2];
rz(1.8895011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3757194) q[1];
sx q[1];
rz(-0.60085591) q[1];
sx q[1];
rz(2.145776) q[1];
x q[2];
rz(-2.0246451) q[3];
sx q[3];
rz(-0.50411187) q[3];
sx q[3];
rz(-0.6716118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9540003) q[2];
sx q[2];
rz(-1.8686998) q[2];
sx q[2];
rz(-0.083258955) q[2];
rz(1.2620874) q[3];
sx q[3];
rz(-0.83674651) q[3];
sx q[3];
rz(-2.4038103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5872203) q[0];
sx q[0];
rz(-1.502259) q[0];
sx q[0];
rz(-2.3578405) q[0];
rz(-1.4844249) q[1];
sx q[1];
rz(-2.617651) q[1];
sx q[1];
rz(-2.8481683) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83111887) q[0];
sx q[0];
rz(-1.9378512) q[0];
sx q[0];
rz(-2.0622753) q[0];
rz(1.2905471) q[2];
sx q[2];
rz(-2.4538605) q[2];
sx q[2];
rz(2.6838145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30006584) q[1];
sx q[1];
rz(-1.7151388) q[1];
sx q[1];
rz(1.1529902) q[1];
x q[2];
rz(-1.8752304) q[3];
sx q[3];
rz(-0.67825156) q[3];
sx q[3];
rz(1.1224318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36711127) q[2];
sx q[2];
rz(-1.7256871) q[2];
sx q[2];
rz(1.492738) q[2];
rz(1.9914918) q[3];
sx q[3];
rz(-0.29905683) q[3];
sx q[3];
rz(0.90112954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25006008) q[0];
sx q[0];
rz(-1.0279259) q[0];
sx q[0];
rz(-3.133339) q[0];
rz(-3.1390269) q[1];
sx q[1];
rz(-0.51104128) q[1];
sx q[1];
rz(2.0770238) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8966347) q[0];
sx q[0];
rz(-1.9799398) q[0];
sx q[0];
rz(-2.4372502) q[0];
rz(-1.8202174) q[2];
sx q[2];
rz(-1.5787243) q[2];
sx q[2];
rz(-0.2570105) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9656532) q[1];
sx q[1];
rz(-2.4373131) q[1];
sx q[1];
rz(-2.2910816) q[1];
rz(-pi) q[2];
rz(1.4225619) q[3];
sx q[3];
rz(-0.98662107) q[3];
sx q[3];
rz(-2.119182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0801487) q[2];
sx q[2];
rz(-1.6950357) q[2];
sx q[2];
rz(0.20305571) q[2];
rz(1.2114581) q[3];
sx q[3];
rz(-1.0695499) q[3];
sx q[3];
rz(-3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9515297) q[0];
sx q[0];
rz(-2.0385346) q[0];
sx q[0];
rz(-0.67361012) q[0];
rz(2.8296962) q[1];
sx q[1];
rz(-1.936828) q[1];
sx q[1];
rz(1.9685251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2791084) q[0];
sx q[0];
rz(-2.0983834) q[0];
sx q[0];
rz(2.8010445) q[0];
rz(-pi) q[1];
rz(-2.8946213) q[2];
sx q[2];
rz(-2.3701043) q[2];
sx q[2];
rz(-2.889973) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38531479) q[1];
sx q[1];
rz(-0.42149252) q[1];
sx q[1];
rz(-1.7569067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42822631) q[3];
sx q[3];
rz(-0.90598327) q[3];
sx q[3];
rz(-1.6237824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46549002) q[2];
sx q[2];
rz(-0.61773053) q[2];
sx q[2];
rz(0.22107302) q[2];
rz(0.57340932) q[3];
sx q[3];
rz(-0.87380496) q[3];
sx q[3];
rz(-2.0360086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40972805) q[0];
sx q[0];
rz(-1.2936214) q[0];
sx q[0];
rz(-3.1142942) q[0];
rz(1.191054) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(1.2192192) q[2];
sx q[2];
rz(-0.58813358) q[2];
sx q[2];
rz(1.0434601) q[2];
rz(2.1656169) q[3];
sx q[3];
rz(-1.2479758) q[3];
sx q[3];
rz(1.6231619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
