OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(2.7288781) q[0];
sx q[0];
rz(5.4469845) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(-2.1332512) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9158186) q[0];
sx q[0];
rz(-1.649722) q[0];
sx q[0];
rz(1.1962587) q[0];
x q[1];
rz(-0.37990976) q[2];
sx q[2];
rz(-1.3515944) q[2];
sx q[2];
rz(0.70218147) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28532449) q[1];
sx q[1];
rz(-2.6372077) q[1];
sx q[1];
rz(-2.8125202) q[1];
rz(-1.3775695) q[3];
sx q[3];
rz(-1.591103) q[3];
sx q[3];
rz(1.3042252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3081554) q[2];
sx q[2];
rz(-0.94352949) q[2];
sx q[2];
rz(-2.4885139) q[2];
rz(0.11450442) q[3];
sx q[3];
rz(-1.3936035) q[3];
sx q[3];
rz(2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36650518) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(1.3436226) q[0];
rz(-1.1603181) q[1];
sx q[1];
rz(-2.7313045) q[1];
sx q[1];
rz(0.62058273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4575093) q[0];
sx q[0];
rz(-0.39827049) q[0];
sx q[0];
rz(-1.0578367) q[0];
x q[1];
rz(-2.7730915) q[2];
sx q[2];
rz(-0.93612367) q[2];
sx q[2];
rz(2.2919185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3239336) q[1];
sx q[1];
rz(-0.91485564) q[1];
sx q[1];
rz(2.0327912) q[1];
rz(-pi) q[2];
rz(-1.1813404) q[3];
sx q[3];
rz(-2.4886294) q[3];
sx q[3];
rz(-0.12784004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4429984) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(-2.911574) q[2];
rz(-1.5864141) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(-0.27568278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25552937) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(2.0607167) q[0];
rz(0.64322645) q[1];
sx q[1];
rz(-2.3422362) q[1];
sx q[1];
rz(0.08531514) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236508) q[0];
sx q[0];
rz(-2.1408014) q[0];
sx q[0];
rz(-0.78582363) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5719937) q[2];
sx q[2];
rz(-1.7503009) q[2];
sx q[2];
rz(0.010590503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6416723) q[1];
sx q[1];
rz(-2.9292078) q[1];
sx q[1];
rz(-2.4714064) q[1];
x q[2];
rz(-1.5253228) q[3];
sx q[3];
rz(-2.1915428) q[3];
sx q[3];
rz(2.5705702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2049415) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(0.19634518) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(2.096368) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333756) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(2.7561482) q[0];
rz(-1.7281945) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(2.8916496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88943931) q[0];
sx q[0];
rz(-1.1595386) q[0];
sx q[0];
rz(-2.2338339) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5416597) q[2];
sx q[2];
rz(-3.0351808) q[2];
sx q[2];
rz(1.0628884) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.476772) q[1];
sx q[1];
rz(-0.52174924) q[1];
sx q[1];
rz(-2.797965) q[1];
rz(0.19574768) q[3];
sx q[3];
rz(-2.3382814) q[3];
sx q[3];
rz(2.3565528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46127737) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(2.9052022) q[2];
rz(-1.2605028) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(2.4642956) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054955) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(2.6564823) q[0];
rz(1.9612034) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(-0.83650437) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9125875) q[0];
sx q[0];
rz(-1.2558172) q[0];
sx q[0];
rz(-0.37512414) q[0];
x q[1];
rz(2.8233158) q[2];
sx q[2];
rz(-1.1851382) q[2];
sx q[2];
rz(-1.1952656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5768535) q[1];
sx q[1];
rz(-0.67138636) q[1];
sx q[1];
rz(-1.4935843) q[1];
rz(-0.14797609) q[3];
sx q[3];
rz(-0.40468369) q[3];
sx q[3];
rz(-1.2860822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(0.76954976) q[2];
rz(0.92929333) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(-0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1997851) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(-2.3727681) q[0];
rz(-1.3517514) q[1];
sx q[1];
rz(-2.6685346) q[1];
sx q[1];
rz(-2.2519055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6883967) q[0];
sx q[0];
rz(-1.9510117) q[0];
sx q[0];
rz(2.9158457) q[0];
x q[1];
rz(-3.1402428) q[2];
sx q[2];
rz(-1.5264411) q[2];
sx q[2];
rz(2.7173017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5701712) q[1];
sx q[1];
rz(-1.4246877) q[1];
sx q[1];
rz(0.04695462) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6906312) q[3];
sx q[3];
rz(-1.6743252) q[3];
sx q[3];
rz(0.24117392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6512904) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(-1.5865405) q[2];
rz(-1.2706612) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(-0.13203013) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1726058) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(1.9975115) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(2.3544618) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89881508) q[0];
sx q[0];
rz(-1.4393596) q[0];
sx q[0];
rz(1.4351033) q[0];
rz(-pi) q[1];
rz(-0.88655858) q[2];
sx q[2];
rz(-1.1349003) q[2];
sx q[2];
rz(0.8800216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5435782) q[1];
sx q[1];
rz(-2.0435963) q[1];
sx q[1];
rz(-0.531457) q[1];
rz(1.9612585) q[3];
sx q[3];
rz(-2.0956846) q[3];
sx q[3];
rz(-0.75492862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0234915) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(3.1084295) q[2];
rz(-1.5432594) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(-1.5020812) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892266) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(1.7484885) q[0];
rz(-2.6847367) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(-1.1005864) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7129214) q[0];
sx q[0];
rz(-1.572519) q[0];
sx q[0];
rz(3.0655906) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4068836) q[2];
sx q[2];
rz(-0.8438973) q[2];
sx q[2];
rz(-1.5555842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3708027) q[1];
sx q[1];
rz(-0.90401559) q[1];
sx q[1];
rz(3.0984224) q[1];
rz(-1.4310775) q[3];
sx q[3];
rz(-1.1828239) q[3];
sx q[3];
rz(-3.1041077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.042772375) q[2];
sx q[2];
rz(-2.6440812) q[2];
sx q[2];
rz(-1.5808606) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(2.1055351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40912691) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(0.29148802) q[0];
rz(-1.2358707) q[1];
sx q[1];
rz(-1.196685) q[1];
sx q[1];
rz(-3.018697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0451806) q[0];
sx q[0];
rz(-1.2244321) q[0];
sx q[0];
rz(1.5161432) q[0];
rz(-pi) q[1];
rz(-0.22237402) q[2];
sx q[2];
rz(-1.8463328) q[2];
sx q[2];
rz(1.0433407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6195685) q[1];
sx q[1];
rz(-1.3641197) q[1];
sx q[1];
rz(0.11264888) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7154152) q[3];
sx q[3];
rz(-1.5767617) q[3];
sx q[3];
rz(-0.98814135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62873944) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(0.17744803) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(-1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43750957) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(1.7105239) q[0];
rz(-0.60072947) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(0.61753714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3445471) q[0];
sx q[0];
rz(-2.8918242) q[0];
sx q[0];
rz(-2.0956371) q[0];
x q[1];
rz(0.52458101) q[2];
sx q[2];
rz(-2.2426105) q[2];
sx q[2];
rz(-1.8730522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70676196) q[1];
sx q[1];
rz(-2.1036356) q[1];
sx q[1];
rz(3.117242) q[1];
rz(-pi) q[2];
rz(-0.57281877) q[3];
sx q[3];
rz(-1.3098992) q[3];
sx q[3];
rz(-3.0134921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2403229) q[2];
sx q[2];
rz(-2.1412886) q[2];
sx q[2];
rz(-0.97736248) q[2];
rz(1.0824925) q[3];
sx q[3];
rz(-0.87477028) q[3];
sx q[3];
rz(3.0860743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8025773) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(-3.0958685) q[1];
sx q[1];
rz(-0.8174236) q[1];
sx q[1];
rz(1.5541706) q[1];
rz(2.634766) q[2];
sx q[2];
rz(-1.097737) q[2];
sx q[2];
rz(-3.1180827) q[2];
rz(2.9228899) q[3];
sx q[3];
rz(-0.66853157) q[3];
sx q[3];
rz(2.4562277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
