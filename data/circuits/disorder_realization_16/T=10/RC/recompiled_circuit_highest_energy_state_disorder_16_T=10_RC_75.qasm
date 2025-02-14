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
rz(1.8546328) q[1];
sx q[1];
rz(-1.9998963) q[1];
sx q[1];
rz(0.72007522) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135277) q[0];
sx q[0];
rz(-1.9682471) q[0];
sx q[0];
rz(1.1532408) q[0];
rz(-1.371648) q[2];
sx q[2];
rz(-1.3825644) q[2];
sx q[2];
rz(-3.0668099) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.076700828) q[1];
sx q[1];
rz(-1.2982839) q[1];
sx q[1];
rz(2.7530159) q[1];
x q[2];
rz(-0.96020697) q[3];
sx q[3];
rz(-1.7060602) q[3];
sx q[3];
rz(-1.3162515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2893452) q[2];
sx q[2];
rz(-2.2699558) q[2];
sx q[2];
rz(-2.3343425) q[2];
rz(-1.6926258) q[3];
sx q[3];
rz(-2.2117386) q[3];
sx q[3];
rz(-2.4248031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221231) q[0];
sx q[0];
rz(-1.8417646) q[0];
sx q[0];
rz(1.8362554) q[0];
rz(-0.7430101) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(1.8441127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7867048) q[0];
sx q[0];
rz(-2.700665) q[0];
sx q[0];
rz(-2.0092677) q[0];
rz(0.93012419) q[2];
sx q[2];
rz(-0.67652297) q[2];
sx q[2];
rz(-2.1695824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6492122) q[1];
sx q[1];
rz(-2.4673176) q[1];
sx q[1];
rz(1.0843305) q[1];
x q[2];
rz(1.4543919) q[3];
sx q[3];
rz(-1.9299514) q[3];
sx q[3];
rz(-2.9352279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7617191) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(2.1011815) q[2];
rz(2.8271293) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(1.8837121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8085025) q[0];
sx q[0];
rz(-0.98121488) q[0];
sx q[0];
rz(-1.3370978) q[0];
rz(0.62726504) q[1];
sx q[1];
rz(-1.4074872) q[1];
sx q[1];
rz(-1.5339392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.92534) q[0];
sx q[0];
rz(-2.0006764) q[0];
sx q[0];
rz(-0.86801189) q[0];
x q[1];
rz(2.5467403) q[2];
sx q[2];
rz(-0.50083465) q[2];
sx q[2];
rz(-2.5822617) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2308637) q[1];
sx q[1];
rz(-1.3861935) q[1];
sx q[1];
rz(0.47841013) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.942105) q[3];
sx q[3];
rz(-1.1987682) q[3];
sx q[3];
rz(-2.2147199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4072121) q[2];
sx q[2];
rz(-2.247159) q[2];
sx q[2];
rz(1.9341932) q[2];
rz(-2.1842128) q[3];
sx q[3];
rz(-1.2355685) q[3];
sx q[3];
rz(0.68937075) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44416881) q[0];
sx q[0];
rz(-2.1902695) q[0];
sx q[0];
rz(-3.0498258) q[0];
rz(0.16608876) q[1];
sx q[1];
rz(-1.8548012) q[1];
sx q[1];
rz(1.0413569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18997858) q[0];
sx q[0];
rz(-0.77617499) q[0];
sx q[0];
rz(-1.6329732) q[0];
x q[1];
rz(-1.9729593) q[2];
sx q[2];
rz(-2.3170174) q[2];
sx q[2];
rz(-2.4006725) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2745121) q[1];
sx q[1];
rz(-1.6050982) q[1];
sx q[1];
rz(2.0438036) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6575069) q[3];
sx q[3];
rz(-1.4299046) q[3];
sx q[3];
rz(-0.85764641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2287717) q[2];
sx q[2];
rz(-0.81835881) q[2];
sx q[2];
rz(-1.6010326) q[2];
rz(2.1238964) q[3];
sx q[3];
rz(-1.3295659) q[3];
sx q[3];
rz(1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8215264) q[0];
sx q[0];
rz(-1.7463266) q[0];
sx q[0];
rz(-2.9030002) q[0];
rz(2.353031) q[1];
sx q[1];
rz(-2.8193654) q[1];
sx q[1];
rz(2.2408748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4960605) q[0];
sx q[0];
rz(-1.1405611) q[0];
sx q[0];
rz(1.2853464) q[0];
rz(-0.094036799) q[2];
sx q[2];
rz(-1.9584736) q[2];
sx q[2];
rz(-2.2212818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2701931) q[1];
sx q[1];
rz(-0.43313036) q[1];
sx q[1];
rz(-0.30850839) q[1];
x q[2];
rz(0.84464083) q[3];
sx q[3];
rz(-0.9880522) q[3];
sx q[3];
rz(0.83707367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0939402) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(1.368604) q[2];
rz(0.37745825) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(-1.9774168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222328) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(2.9121616) q[0];
rz(-2.5550628) q[1];
sx q[1];
rz(-0.49085453) q[1];
sx q[1];
rz(-1.4803001) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585715) q[0];
sx q[0];
rz(-0.93556306) q[0];
sx q[0];
rz(-0.44376377) q[0];
x q[1];
rz(-1.3103799) q[2];
sx q[2];
rz(-2.1245109) q[2];
sx q[2];
rz(-1.2147724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6203998) q[1];
sx q[1];
rz(-1.7971276) q[1];
sx q[1];
rz(1.625555) q[1];
rz(2.2451083) q[3];
sx q[3];
rz(-1.7411045) q[3];
sx q[3];
rz(-1.8425106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4218939) q[2];
sx q[2];
rz(-0.59013683) q[2];
sx q[2];
rz(-1.7936919) q[2];
rz(2.9501996) q[3];
sx q[3];
rz(-2.8189711) q[3];
sx q[3];
rz(1.2389244) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58452463) q[0];
sx q[0];
rz(-1.5140336) q[0];
sx q[0];
rz(2.6051482) q[0];
rz(3.0570807) q[1];
sx q[1];
rz(-0.5785431) q[1];
sx q[1];
rz(-1.9523778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093350323) q[0];
sx q[0];
rz(-0.95365253) q[0];
sx q[0];
rz(0.18927745) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68185735) q[2];
sx q[2];
rz(-1.8769662) q[2];
sx q[2];
rz(-0.078753565) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3757194) q[1];
sx q[1];
rz(-2.5407367) q[1];
sx q[1];
rz(-0.99581666) q[1];
rz(-pi) q[2];
rz(0.23730324) q[3];
sx q[3];
rz(-2.0198698) q[3];
sx q[3];
rz(0.16332353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1875923) q[2];
sx q[2];
rz(-1.8686998) q[2];
sx q[2];
rz(0.083258955) q[2];
rz(-1.2620874) q[3];
sx q[3];
rz(-0.83674651) q[3];
sx q[3];
rz(2.4038103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5872203) q[0];
sx q[0];
rz(-1.502259) q[0];
sx q[0];
rz(2.3578405) q[0];
rz(1.4844249) q[1];
sx q[1];
rz(-0.52394167) q[1];
sx q[1];
rz(0.2934244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14908168) q[0];
sx q[0];
rz(-0.60428491) q[0];
sx q[0];
rz(0.88715951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2905471) q[2];
sx q[2];
rz(-2.4538605) q[2];
sx q[2];
rz(2.6838145) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8070911) q[1];
sx q[1];
rz(-1.9839904) q[1];
sx q[1];
rz(2.9838802) q[1];
x q[2];
rz(-2.226172) q[3];
sx q[3];
rz(-1.3815945) q[3];
sx q[3];
rz(-0.20841852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7744814) q[2];
sx q[2];
rz(-1.4159055) q[2];
sx q[2];
rz(-1.492738) q[2];
rz(1.9914918) q[3];
sx q[3];
rz(-2.8425358) q[3];
sx q[3];
rz(2.2404631) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25006008) q[0];
sx q[0];
rz(-2.1136668) q[0];
sx q[0];
rz(3.133339) q[0];
rz(3.1390269) q[1];
sx q[1];
rz(-2.6305514) q[1];
sx q[1];
rz(-1.0645688) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417404) q[0];
sx q[0];
rz(-2.2069262) q[0];
sx q[0];
rz(-1.0534691) q[0];
rz(1.3213753) q[2];
sx q[2];
rz(-1.5787243) q[2];
sx q[2];
rz(2.8845822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67978102) q[1];
sx q[1];
rz(-2.0790599) q[1];
sx q[1];
rz(0.51080461) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4225619) q[3];
sx q[3];
rz(-0.98662107) q[3];
sx q[3];
rz(1.0224107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0801487) q[2];
sx q[2];
rz(-1.6950357) q[2];
sx q[2];
rz(2.9385369) q[2];
rz(1.2114581) q[3];
sx q[3];
rz(-1.0695499) q[3];
sx q[3];
rz(-3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-0.31189648) q[1];
sx q[1];
rz(-1.2047647) q[1];
sx q[1];
rz(-1.9685251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4757901) q[0];
sx q[0];
rz(-0.6190933) q[0];
sx q[0];
rz(2.0912916) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3854957) q[2];
sx q[2];
rz(-1.399516) q[2];
sx q[2];
rz(-1.1403699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0153326) q[1];
sx q[1];
rz(-1.646572) q[1];
sx q[1];
rz(1.9858236) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0836175) q[3];
sx q[3];
rz(-2.3688032) q[3];
sx q[3];
rz(-2.1548398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6761026) q[2];
sx q[2];
rz(-2.5238621) q[2];
sx q[2];
rz(-2.9205196) q[2];
rz(-0.57340932) q[3];
sx q[3];
rz(-0.87380496) q[3];
sx q[3];
rz(2.0360086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40972805) q[0];
sx q[0];
rz(-1.8479713) q[0];
sx q[0];
rz(0.027298409) q[0];
rz(-1.9505386) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(-2.9158557) q[2];
sx q[2];
rz(-1.0229243) q[2];
sx q[2];
rz(-1.682874) q[2];
rz(-1.0326019) q[3];
sx q[3];
rz(-0.66734869) q[3];
sx q[3];
rz(-2.6507631) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
