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
rz(0.25843698) q[0];
sx q[0];
rz(2.8539477) q[0];
sx q[0];
rz(10.532425) q[0];
rz(1.3232752) q[1];
sx q[1];
rz(3.4634436) q[1];
sx q[1];
rz(8.7965214) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1028672) q[0];
sx q[0];
rz(-2.9280781) q[0];
sx q[0];
rz(2.1597014) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4187367) q[2];
sx q[2];
rz(-1.8713163) q[2];
sx q[2];
rz(-0.84096694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19753708) q[1];
sx q[1];
rz(-1.5345795) q[1];
sx q[1];
rz(0.22878583) q[1];
x q[2];
rz(1.0367582) q[3];
sx q[3];
rz(-1.1849355) q[3];
sx q[3];
rz(2.2552668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85311741) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(2.8004069) q[2];
rz(0.8902542) q[3];
sx q[3];
rz(-1.8743926) q[3];
sx q[3];
rz(-0.53513479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.61966908) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(-2.4271915) q[0];
rz(1.5419386) q[1];
sx q[1];
rz(-0.49595141) q[1];
sx q[1];
rz(-3.0468859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3886914) q[0];
sx q[0];
rz(-1.6651648) q[0];
sx q[0];
rz(1.8917985) q[0];
x q[1];
rz(0.96927283) q[2];
sx q[2];
rz(-1.8516632) q[2];
sx q[2];
rz(1.2796677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7357121) q[1];
sx q[1];
rz(-0.42197126) q[1];
sx q[1];
rz(2.5338245) q[1];
x q[2];
rz(-1.7294782) q[3];
sx q[3];
rz(-1.0668716) q[3];
sx q[3];
rz(-2.9014587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9793205) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(-2.6675513) q[2];
rz(0.76215172) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(-2.6704085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2786461) q[0];
sx q[0];
rz(-3.0248248) q[0];
sx q[0];
rz(0.85292029) q[0];
rz(2.9184753) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(0.51582897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11027656) q[0];
sx q[0];
rz(-3.0678684) q[0];
sx q[0];
rz(-2.323193) q[0];
rz(-pi) q[1];
rz(0.25571172) q[2];
sx q[2];
rz(-0.94921321) q[2];
sx q[2];
rz(0.58178025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3577376) q[1];
sx q[1];
rz(-0.92354316) q[1];
sx q[1];
rz(0.81070645) q[1];
x q[2];
rz(-1.9733377) q[3];
sx q[3];
rz(-2.4022837) q[3];
sx q[3];
rz(-2.4420247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8082661) q[2];
sx q[2];
rz(-0.9767248) q[2];
sx q[2];
rz(-3.0466363) q[2];
rz(2.700108) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(-1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3917291) q[0];
sx q[0];
rz(-2.5666105) q[0];
sx q[0];
rz(1.0561426) q[0];
rz(0.9437584) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(1.6897078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34925845) q[0];
sx q[0];
rz(-2.1448488) q[0];
sx q[0];
rz(-0.7416772) q[0];
rz(-pi) q[1];
x q[1];
rz(2.679616) q[2];
sx q[2];
rz(-1.9236132) q[2];
sx q[2];
rz(-1.3188254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31409594) q[1];
sx q[1];
rz(-0.24330595) q[1];
sx q[1];
rz(2.9929586) q[1];
rz(-pi) q[2];
rz(-1.5632252) q[3];
sx q[3];
rz(-1.4601344) q[3];
sx q[3];
rz(2.8286434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4132495) q[2];
sx q[2];
rz(-0.80060935) q[2];
sx q[2];
rz(-2.441414) q[2];
rz(2.2705196) q[3];
sx q[3];
rz(-1.0659404) q[3];
sx q[3];
rz(1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313318) q[0];
sx q[0];
rz(-0.49322525) q[0];
sx q[0];
rz(-2.6569271) q[0];
rz(2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(-2.7427618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054886) q[0];
sx q[0];
rz(-1.7057422) q[0];
sx q[0];
rz(0.55079726) q[0];
rz(-pi) q[1];
rz(-0.70339922) q[2];
sx q[2];
rz(-1.8336772) q[2];
sx q[2];
rz(0.92353283) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3276766) q[1];
sx q[1];
rz(-2.3688201) q[1];
sx q[1];
rz(2.1111958) q[1];
rz(-pi) q[2];
rz(2.5957803) q[3];
sx q[3];
rz(-0.94792507) q[3];
sx q[3];
rz(-0.92178492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6681119) q[2];
sx q[2];
rz(-1.5405737) q[2];
sx q[2];
rz(-0.87781805) q[2];
rz(0.38464883) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(-1.1788684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8574852) q[0];
sx q[0];
rz(-2.6326038) q[0];
sx q[0];
rz(2.8711163) q[0];
rz(-0.30666223) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(-0.57672966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79918843) q[0];
sx q[0];
rz(-1.6584089) q[0];
sx q[0];
rz(1.7566998) q[0];
rz(-pi) q[1];
rz(-2.2606196) q[2];
sx q[2];
rz(-1.2861797) q[2];
sx q[2];
rz(-1.6135482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19837436) q[1];
sx q[1];
rz(-2.4369708) q[1];
sx q[1];
rz(2.2072029) q[1];
rz(-pi) q[2];
rz(-0.0402952) q[3];
sx q[3];
rz(-0.3350122) q[3];
sx q[3];
rz(-2.1744436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7755255) q[2];
sx q[2];
rz(-2.0719353) q[2];
sx q[2];
rz(0.97037399) q[2];
rz(-2.7905285) q[3];
sx q[3];
rz(-2.909436) q[3];
sx q[3];
rz(-2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.905726) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(-2.6406777) q[0];
rz(-2.4210335) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(0.93773425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2116579) q[0];
sx q[0];
rz(-1.5296361) q[0];
sx q[0];
rz(1.7743054) q[0];
rz(-pi) q[1];
rz(2.4831122) q[2];
sx q[2];
rz(-0.57504932) q[2];
sx q[2];
rz(-2.5031896) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.42635065) q[1];
sx q[1];
rz(-0.75255005) q[1];
sx q[1];
rz(2.1133711) q[1];
rz(-pi) q[2];
rz(0.012091919) q[3];
sx q[3];
rz(-0.68688697) q[3];
sx q[3];
rz(1.505132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0096036) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(2.2930938) q[2];
rz(2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(1.2257303) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7312412) q[0];
sx q[0];
rz(-2.1827965) q[0];
sx q[0];
rz(0.89286667) q[0];
rz(0.061575312) q[1];
sx q[1];
rz(-2.4696746) q[1];
sx q[1];
rz(2.9796013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91770691) q[0];
sx q[0];
rz(-0.10334238) q[0];
sx q[0];
rz(-1.7211821) q[0];
rz(-pi) q[1];
rz(0.70190491) q[2];
sx q[2];
rz(-1.4830952) q[2];
sx q[2];
rz(3.1021995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3969035) q[1];
sx q[1];
rz(-1.2619881) q[1];
sx q[1];
rz(1.4267322) q[1];
x q[2];
rz(1.6747159) q[3];
sx q[3];
rz(-1.6641518) q[3];
sx q[3];
rz(-0.68478497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.050345) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(-0.72489911) q[2];
rz(0.27724087) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(-3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8898833) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(3.1157893) q[0];
rz(1.7426527) q[1];
sx q[1];
rz(-1.6394697) q[1];
sx q[1];
rz(-0.10765156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336595) q[0];
sx q[0];
rz(-1.7042394) q[0];
sx q[0];
rz(2.1428277) q[0];
rz(-1.1321105) q[2];
sx q[2];
rz(-2.0033547) q[2];
sx q[2];
rz(-2.3985942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.746826) q[1];
sx q[1];
rz(-1.0977543) q[1];
sx q[1];
rz(1.6113144) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4156451) q[3];
sx q[3];
rz(-0.11501139) q[3];
sx q[3];
rz(1.0105159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.671635) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(2.0476445) q[2];
rz(-2.5661902) q[3];
sx q[3];
rz(-0.83991528) q[3];
sx q[3];
rz(-1.1168787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-0.58654439) q[0];
sx q[0];
rz(-2.4285512) q[0];
rz(2.9167922) q[1];
sx q[1];
rz(-2.5495922) q[1];
sx q[1];
rz(-2.0483268) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7894692) q[0];
sx q[0];
rz(-2.2575245) q[0];
sx q[0];
rz(2.5144469) q[0];
x q[1];
rz(1.8453127) q[2];
sx q[2];
rz(-1.3764672) q[2];
sx q[2];
rz(-2.7467665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.374129) q[1];
sx q[1];
rz(-1.1157562) q[1];
sx q[1];
rz(-1.6947075) q[1];
rz(-pi) q[2];
rz(-1.9841927) q[3];
sx q[3];
rz(-1.6781312) q[3];
sx q[3];
rz(-1.3232376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1405868) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(0.33983964) q[2];
rz(0.042424399) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(-2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18702678) q[0];
sx q[0];
rz(-1.621959) q[0];
sx q[0];
rz(1.8931615) q[0];
rz(-0.78202248) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(2.7693979) q[2];
sx q[2];
rz(-1.2584465) q[2];
sx q[2];
rz(-1.0312205) q[2];
rz(-2.8911968) q[3];
sx q[3];
rz(-1.9848817) q[3];
sx q[3];
rz(3.1105697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
