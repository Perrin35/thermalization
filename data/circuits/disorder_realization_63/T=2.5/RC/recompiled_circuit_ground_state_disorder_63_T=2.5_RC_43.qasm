OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0286921) q[0];
sx q[0];
rz(-2.8960462) q[0];
sx q[0];
rz(-3.088933) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(-0.3408365) q[1];
sx q[1];
rz(-1.0885106) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3707293) q[0];
sx q[0];
rz(-1.5077871) q[0];
sx q[0];
rz(-1.1746688) q[0];
x q[1];
rz(-1.4289189) q[2];
sx q[2];
rz(-1.6798225) q[2];
sx q[2];
rz(-3.1224868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8177942) q[1];
sx q[1];
rz(-1.7355738) q[1];
sx q[1];
rz(0.64357693) q[1];
x q[2];
rz(-0.0040715374) q[3];
sx q[3];
rz(-1.1511251) q[3];
sx q[3];
rz(1.3908902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6344305) q[2];
sx q[2];
rz(-1.7796702) q[2];
sx q[2];
rz(-2.5334899) q[2];
rz(-2.0173343) q[3];
sx q[3];
rz(-1.9367155) q[3];
sx q[3];
rz(-2.2500136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41481498) q[0];
sx q[0];
rz(-1.5685273) q[0];
sx q[0];
rz(-0.60348764) q[0];
rz(-0.51015774) q[1];
sx q[1];
rz(-2.3699103) q[1];
sx q[1];
rz(1.0268432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.802216) q[0];
sx q[0];
rz(-0.23288865) q[0];
sx q[0];
rz(0.73407895) q[0];
rz(-pi) q[1];
rz(0.85058327) q[2];
sx q[2];
rz(-2.2424978) q[2];
sx q[2];
rz(1.6784061) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34705838) q[1];
sx q[1];
rz(-1.0640941) q[1];
sx q[1];
rz(-2.3757412) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2193416) q[3];
sx q[3];
rz(-2.3664858) q[3];
sx q[3];
rz(2.4668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5828731) q[2];
sx q[2];
rz(-1.1483973) q[2];
sx q[2];
rz(0.64644512) q[2];
rz(1.8630113) q[3];
sx q[3];
rz(-0.94978142) q[3];
sx q[3];
rz(1.6897374) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2090787) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(0.34700829) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(-2.2775547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0873252) q[0];
sx q[0];
rz(-1.010276) q[0];
sx q[0];
rz(-1.3232687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50484294) q[2];
sx q[2];
rz(-2.6528203) q[2];
sx q[2];
rz(-0.29112838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62926312) q[1];
sx q[1];
rz(-2.1844668) q[1];
sx q[1];
rz(2.3113109) q[1];
rz(1.044831) q[3];
sx q[3];
rz(-1.4792031) q[3];
sx q[3];
rz(-2.4909702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0916834) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(-1.4570215) q[2];
rz(2.125804) q[3];
sx q[3];
rz(-0.53272811) q[3];
sx q[3];
rz(-0.31639019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006573) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(1.578791) q[0];
rz(-1.0904301) q[1];
sx q[1];
rz(-2.5999887) q[1];
sx q[1];
rz(-1.1515559) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59127677) q[0];
sx q[0];
rz(-1.7999066) q[0];
sx q[0];
rz(-2.5510049) q[0];
rz(-pi) q[1];
rz(-0.76049034) q[2];
sx q[2];
rz(-1.9273238) q[2];
sx q[2];
rz(-0.52887756) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.71660605) q[1];
sx q[1];
rz(-2.0946436) q[1];
sx q[1];
rz(0.040158466) q[1];
x q[2];
rz(1.5308446) q[3];
sx q[3];
rz(-2.5343347) q[3];
sx q[3];
rz(2.0366397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84732032) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(1.4997743) q[2];
rz(-0.13449399) q[3];
sx q[3];
rz(-1.2765063) q[3];
sx q[3];
rz(-0.78181481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.0791557) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(-0.4656747) q[0];
rz(-1.8648719) q[1];
sx q[1];
rz(-2.1379505) q[1];
sx q[1];
rz(2.1086955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55556923) q[0];
sx q[0];
rz(-2.9831671) q[0];
sx q[0];
rz(1.5044295) q[0];
x q[1];
rz(-0.6200142) q[2];
sx q[2];
rz(-2.9699932) q[2];
sx q[2];
rz(0.5172587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8414265) q[1];
sx q[1];
rz(-2.4304217) q[1];
sx q[1];
rz(-0.39022846) q[1];
x q[2];
rz(0.9744076) q[3];
sx q[3];
rz(-0.63035175) q[3];
sx q[3];
rz(-0.74385616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5501962) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-0.93185321) q[2];
rz(-3.0469117) q[3];
sx q[3];
rz(-0.90016142) q[3];
sx q[3];
rz(1.3100821) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4513627) q[0];
sx q[0];
rz(-0.19915038) q[0];
sx q[0];
rz(-0.1061826) q[0];
rz(0.55446082) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(-1.6000481) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743177) q[0];
sx q[0];
rz(-1.5467001) q[0];
sx q[0];
rz(2.7510608) q[0];
rz(-pi) q[1];
rz(1.0511398) q[2];
sx q[2];
rz(-2.6259632) q[2];
sx q[2];
rz(-1.3517429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9408343) q[1];
sx q[1];
rz(-1.6929827) q[1];
sx q[1];
rz(-0.69735511) q[1];
rz(2.9227436) q[3];
sx q[3];
rz(-0.25222029) q[3];
sx q[3];
rz(0.24009934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9083378) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(-1.6439269) q[2];
rz(2.6805367) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(-1.3155931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6719565) q[0];
sx q[0];
rz(-2.3905583) q[0];
sx q[0];
rz(2.0627956) q[0];
rz(0.59393334) q[1];
sx q[1];
rz(-0.97336951) q[1];
sx q[1];
rz(2.914391) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5626635) q[0];
sx q[0];
rz(-1.4592071) q[0];
sx q[0];
rz(1.7203549) q[0];
rz(-pi) q[1];
rz(0.14102139) q[2];
sx q[2];
rz(-1.2407836) q[2];
sx q[2];
rz(-2.7591456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0409483) q[1];
sx q[1];
rz(-2.4923263) q[1];
sx q[1];
rz(-1.0430286) q[1];
rz(-pi) q[2];
rz(1.3183837) q[3];
sx q[3];
rz(-0.38312437) q[3];
sx q[3];
rz(-2.8069404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6195153) q[2];
sx q[2];
rz(-2.7232813) q[2];
sx q[2];
rz(0.72959161) q[2];
rz(-1.4531762) q[3];
sx q[3];
rz(-1.4540693) q[3];
sx q[3];
rz(-0.61354536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5486117) q[0];
sx q[0];
rz(-0.91658968) q[0];
sx q[0];
rz(0.35162893) q[0];
rz(-2.8278606) q[1];
sx q[1];
rz(-1.9351363) q[1];
sx q[1];
rz(-1.7650013) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976395) q[0];
sx q[0];
rz(-0.40297976) q[0];
sx q[0];
rz(-0.12059327) q[0];
rz(-pi) q[1];
rz(2.1765472) q[2];
sx q[2];
rz(-0.20119431) q[2];
sx q[2];
rz(-2.4810042) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2089858) q[1];
sx q[1];
rz(-2.5521005) q[1];
sx q[1];
rz(-0.28372753) q[1];
x q[2];
rz(-1.3198765) q[3];
sx q[3];
rz(-2.7610215) q[3];
sx q[3];
rz(-1.5151092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5967963) q[2];
sx q[2];
rz(-2.4185541) q[2];
sx q[2];
rz(1.6205622) q[2];
rz(-2.9936786) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(1.773905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77383298) q[0];
sx q[0];
rz(-1.8748883) q[0];
sx q[0];
rz(1.2441147) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(-0.22020766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51802902) q[0];
sx q[0];
rz(-1.9401778) q[0];
sx q[0];
rz(-2.2165038) q[0];
rz(-pi) q[1];
rz(2.979109) q[2];
sx q[2];
rz(-2.0318609) q[2];
sx q[2];
rz(-1.7562255) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0427093) q[1];
sx q[1];
rz(-2.115887) q[1];
sx q[1];
rz(-0.49867757) q[1];
rz(-pi) q[2];
rz(-1.4768019) q[3];
sx q[3];
rz(-1.4729744) q[3];
sx q[3];
rz(1.0536915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9571017) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(-1.5000878) q[2];
rz(-3.0575276) q[3];
sx q[3];
rz(-1.8819239) q[3];
sx q[3];
rz(0.070092289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67749196) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(0.41859928) q[0];
rz(0.94340008) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(-2.8387866) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3269791) q[0];
sx q[0];
rz(-1.7075065) q[0];
sx q[0];
rz(-1.4025406) q[0];
rz(-pi) q[1];
rz(-0.19089107) q[2];
sx q[2];
rz(-1.7451236) q[2];
sx q[2];
rz(-0.088100351) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3258241) q[1];
sx q[1];
rz(-3.0481955) q[1];
sx q[1];
rz(-2.6572822) q[1];
rz(-pi) q[2];
rz(-1.8859768) q[3];
sx q[3];
rz(-0.36096301) q[3];
sx q[3];
rz(2.7415467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6454978) q[2];
sx q[2];
rz(-2.2690513) q[2];
sx q[2];
rz(2.3910451) q[2];
rz(1.6802855) q[3];
sx q[3];
rz(-2.3716726) q[3];
sx q[3];
rz(-0.51699483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207598) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(-1.2278521) q[1];
sx q[1];
rz(-0.49929437) q[1];
sx q[1];
rz(1.1566537) q[1];
rz(-2.5748809) q[2];
sx q[2];
rz(-1.7253582) q[2];
sx q[2];
rz(2.2127989) q[2];
rz(-0.1838818) q[3];
sx q[3];
rz(-1.2533698) q[3];
sx q[3];
rz(-1.7141242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
