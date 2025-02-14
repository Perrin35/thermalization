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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(0.32455197) q[0];
rz(2.1895154) q[1];
sx q[1];
rz(-2.9157186) q[1];
sx q[1];
rz(-1.3083375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65303923) q[0];
sx q[0];
rz(-1.7069611) q[0];
sx q[0];
rz(-0.67739886) q[0];
x q[1];
rz(-1.4412854) q[2];
sx q[2];
rz(-1.0900094) q[2];
sx q[2];
rz(-3.0134984) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5950299) q[1];
sx q[1];
rz(-1.1076704) q[1];
sx q[1];
rz(2.8519277) q[1];
x q[2];
rz(1.7115373) q[3];
sx q[3];
rz(-0.87647299) q[3];
sx q[3];
rz(-1.7870045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7294881) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(-0.35120249) q[2];
rz(1.2900194) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7555162) q[0];
sx q[0];
rz(-1.9367243) q[0];
sx q[0];
rz(2.2112041) q[0];
rz(1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-2.5321541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.45935) q[0];
sx q[0];
rz(-1.1455904) q[0];
sx q[0];
rz(-3.0283079) q[0];
rz(-pi) q[1];
rz(-1.8687042) q[2];
sx q[2];
rz(-1.6570083) q[2];
sx q[2];
rz(2.6484368) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6516782) q[1];
sx q[1];
rz(-2.6082572) q[1];
sx q[1];
rz(2.3467605) q[1];
rz(-1.0004382) q[3];
sx q[3];
rz(-1.7118239) q[3];
sx q[3];
rz(1.0424249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(-2.6212202) q[3];
sx q[3];
rz(-0.64295355) q[3];
sx q[3];
rz(-0.61409942) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60270131) q[0];
sx q[0];
rz(-2.4245872) q[0];
sx q[0];
rz(0.30211788) q[0];
rz(0.27613861) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(1.709323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75452226) q[0];
sx q[0];
rz(-1.1898479) q[0];
sx q[0];
rz(-0.98636287) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8319857) q[2];
sx q[2];
rz(-0.99203324) q[2];
sx q[2];
rz(-0.28489339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70967445) q[1];
sx q[1];
rz(-2.1248105) q[1];
sx q[1];
rz(2.9565503) q[1];
rz(-3.0705419) q[3];
sx q[3];
rz(-2.0356352) q[3];
sx q[3];
rz(-2.0185061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8351195) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(1.845537) q[2];
rz(-2.6895788) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(-0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52962676) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(1.8091328) q[0];
rz(2.0924856) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(-1.5528991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4175966) q[0];
sx q[0];
rz(-1.231047) q[0];
sx q[0];
rz(-1.2509565) q[0];
rz(-pi) q[1];
rz(2.0715782) q[2];
sx q[2];
rz(-1.5588648) q[2];
sx q[2];
rz(-0.8296166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3192352) q[1];
sx q[1];
rz(-2.2062782) q[1];
sx q[1];
rz(0.24680071) q[1];
x q[2];
rz(-2.2838628) q[3];
sx q[3];
rz(-1.1489043) q[3];
sx q[3];
rz(0.5798012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4730452) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(0.2505396) q[2];
rz(0.9969095) q[3];
sx q[3];
rz(-2.0018115) q[3];
sx q[3];
rz(-0.38058773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32960358) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(-1.3478152) q[0];
rz(-2.7373121) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(-2.0124729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383352) q[0];
sx q[0];
rz(-1.1764488) q[0];
sx q[0];
rz(-2.3189937) q[0];
x q[1];
rz(0.26910946) q[2];
sx q[2];
rz(-0.92293149) q[2];
sx q[2];
rz(-1.8108167) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4562005) q[1];
sx q[1];
rz(-2.2397352) q[1];
sx q[1];
rz(3.0078956) q[1];
x q[2];
rz(0.31480883) q[3];
sx q[3];
rz(-1.6757351) q[3];
sx q[3];
rz(0.85280692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79444844) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(-1.8974737) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619549) q[0];
sx q[0];
rz(-0.10529101) q[0];
sx q[0];
rz(0.42718497) q[0];
rz(-0.034491388) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(-0.010206612) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50960582) q[0];
sx q[0];
rz(-1.5684396) q[0];
sx q[0];
rz(1.5690593) q[0];
x q[1];
rz(0.9385061) q[2];
sx q[2];
rz(-1.7006497) q[2];
sx q[2];
rz(2.1406895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8275244) q[1];
sx q[1];
rz(-1.6351103) q[1];
sx q[1];
rz(-2.2576648) q[1];
rz(-pi) q[2];
rz(-0.7179435) q[3];
sx q[3];
rz(-2.4417851) q[3];
sx q[3];
rz(-2.1199302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5222142) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(2.348032) q[2];
rz(-0.86269745) q[3];
sx q[3];
rz(-1.5746652) q[3];
sx q[3];
rz(-0.70393744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70235395) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(2.0080361) q[0];
rz(-1.0776862) q[1];
sx q[1];
rz(-0.51529854) q[1];
sx q[1];
rz(2.8394707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0691119) q[0];
sx q[0];
rz(-1.2337419) q[0];
sx q[0];
rz(-1.2500136) q[0];
rz(-2.0321192) q[2];
sx q[2];
rz(-1.0240842) q[2];
sx q[2];
rz(0.49671587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1715917) q[1];
sx q[1];
rz(-1.0495032) q[1];
sx q[1];
rz(-1.820351) q[1];
rz(-0.98432912) q[3];
sx q[3];
rz(-2.1696089) q[3];
sx q[3];
rz(-0.58517712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9083531) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(1.3471777) q[2];
rz(1.2498648) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(-3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5938479) q[0];
sx q[0];
rz(-0.43314728) q[0];
sx q[0];
rz(-0.10840848) q[0];
rz(2.0938865) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(2.122706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.326556) q[0];
sx q[0];
rz(-0.78062764) q[0];
sx q[0];
rz(2.3579979) q[0];
x q[1];
rz(-0.76129976) q[2];
sx q[2];
rz(-1.5660962) q[2];
sx q[2];
rz(-2.0705303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8727303) q[1];
sx q[1];
rz(-2.514317) q[1];
sx q[1];
rz(-2.4226047) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4298444) q[3];
sx q[3];
rz(-2.3059855) q[3];
sx q[3];
rz(-2.6648389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63742796) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(0.49312433) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(1.5639308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2757932) q[0];
sx q[0];
rz(-1.9075305) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(-1.5419143) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(2.7154198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3329553) q[0];
sx q[0];
rz(-1.5530968) q[0];
sx q[0];
rz(1.5541535) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8502281) q[2];
sx q[2];
rz(-0.93932187) q[2];
sx q[2];
rz(-0.15635083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8746169) q[1];
sx q[1];
rz(-0.24493453) q[1];
sx q[1];
rz(2.0744978) q[1];
rz(-pi) q[2];
rz(2.3087193) q[3];
sx q[3];
rz(-0.98372059) q[3];
sx q[3];
rz(0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6466732) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(-1.0529998) q[2];
rz(2.7827175) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7786355) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(-2.2208075) q[0];
rz(2.6329363) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(-0.42053929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3497999) q[0];
sx q[0];
rz(-1.0726439) q[0];
sx q[0];
rz(-0.80043261) q[0];
x q[1];
rz(-2.2013046) q[2];
sx q[2];
rz(-1.7675401) q[2];
sx q[2];
rz(-1.8041174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.072445446) q[1];
sx q[1];
rz(-2.0305567) q[1];
sx q[1];
rz(-2.2531177) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4484506) q[3];
sx q[3];
rz(-2.4469355) q[3];
sx q[3];
rz(2.8784424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4113808) q[2];
sx q[2];
rz(-0.79019848) q[2];
sx q[2];
rz(-0.31663695) q[2];
rz(0.5591048) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(2.3234698) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.485514) q[0];
sx q[0];
rz(-1.0777127) q[0];
sx q[0];
rz(1.0832473) q[0];
rz(2.6651233) q[1];
sx q[1];
rz(-2.10119) q[1];
sx q[1];
rz(1.4269921) q[1];
rz(1.477735) q[2];
sx q[2];
rz(-2.3882967) q[2];
sx q[2];
rz(-0.49327539) q[2];
rz(3.0477897) q[3];
sx q[3];
rz(-1.9524912) q[3];
sx q[3];
rz(-2.5759202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
