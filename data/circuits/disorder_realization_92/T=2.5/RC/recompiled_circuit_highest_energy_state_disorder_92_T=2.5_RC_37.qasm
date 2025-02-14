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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1028672) q[0];
sx q[0];
rz(-2.9280781) q[0];
sx q[0];
rz(-0.98189129) q[0];
rz(-0.43810644) q[2];
sx q[2];
rz(-0.77229653) q[2];
sx q[2];
rz(0.4060678) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3648277) q[1];
sx q[1];
rz(-1.7994295) q[1];
sx q[1];
rz(-1.6079812) q[1];
rz(-pi) q[2];
rz(2.1048344) q[3];
sx q[3];
rz(-1.9566571) q[3];
sx q[3];
rz(2.2552668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2884752) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(0.34118578) q[2];
rz(0.8902542) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(-2.6064579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61966908) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(-0.71440119) q[0];
rz(1.5419386) q[1];
sx q[1];
rz(-0.49595141) q[1];
sx q[1];
rz(0.094706789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7529012) q[0];
sx q[0];
rz(-1.6651648) q[0];
sx q[0];
rz(-1.8917985) q[0];
rz(-pi) q[1];
rz(0.33659597) q[2];
sx q[2];
rz(-0.99592745) q[2];
sx q[2];
rz(-2.6624555) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4112139) q[1];
sx q[1];
rz(-1.3347365) q[1];
sx q[1];
rz(2.7884931) q[1];
x q[2];
rz(-1.4121145) q[3];
sx q[3];
rz(-2.074721) q[3];
sx q[3];
rz(-2.9014587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16227214) q[2];
sx q[2];
rz(-1.2459735) q[2];
sx q[2];
rz(0.47404131) q[2];
rz(0.76215172) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(-2.6704085) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-0.85292029) q[0];
rz(-0.22311738) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(0.51582897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2775622) q[0];
sx q[0];
rz(-1.5169965) q[0];
sx q[0];
rz(-3.0911616) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8858809) q[2];
sx q[2];
rz(-0.94921321) q[2];
sx q[2];
rz(-0.58178025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8341537) q[1];
sx q[1];
rz(-0.98886988) q[1];
sx q[1];
rz(0.80640275) q[1];
x q[2];
rz(2.7985057) q[3];
sx q[3];
rz(-2.2394387) q[3];
sx q[3];
rz(-0.17681387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8082661) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(-0.094956368) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74986356) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(-2.0854501) q[0];
rz(2.1978343) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(1.4518849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34925845) q[0];
sx q[0];
rz(-0.99674388) q[0];
sx q[0];
rz(2.3999155) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4511172) q[2];
sx q[2];
rz(-0.5734517) q[2];
sx q[2];
rz(0.8586463) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8274967) q[1];
sx q[1];
rz(-0.24330595) q[1];
sx q[1];
rz(-2.9929586) q[1];
x q[2];
rz(1.5783674) q[3];
sx q[3];
rz(-1.4601344) q[3];
sx q[3];
rz(-0.31294926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(0.70017868) q[2];
rz(2.2705196) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(1.357249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313318) q[0];
sx q[0];
rz(-0.49322525) q[0];
sx q[0];
rz(0.48466551) q[0];
rz(2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(-2.7427618) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1171378) q[0];
sx q[0];
rz(-2.1160191) q[0];
sx q[0];
rz(-1.7288037) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4381934) q[2];
sx q[2];
rz(-1.3079155) q[2];
sx q[2];
rz(2.2180598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3276766) q[1];
sx q[1];
rz(-0.77277257) q[1];
sx q[1];
rz(-2.1111958) q[1];
x q[2];
rz(0.871931) q[3];
sx q[3];
rz(-1.1355577) q[3];
sx q[3];
rz(-0.98952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47348076) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(-0.87781805) q[2];
rz(2.7569438) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(-1.9627242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8574852) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(0.2704764) q[0];
rz(-0.30666223) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(2.564863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79918843) q[0];
sx q[0];
rz(-1.6584089) q[0];
sx q[0];
rz(-1.7566998) q[0];
rz(-pi) q[1];
rz(-0.36251601) q[2];
sx q[2];
rz(-0.91372817) q[2];
sx q[2];
rz(-0.18489472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8851848) q[1];
sx q[1];
rz(-1.1756304) q[1];
sx q[1];
rz(-0.97104071) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8068325) q[3];
sx q[3];
rz(-1.5575512) q[3];
sx q[3];
rz(2.5760026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3660672) q[2];
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
sx q[3];
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
rz(1.2358667) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(-2.6406777) q[0];
rz(-0.72055912) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(2.2038584) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9792694) q[0];
sx q[0];
rz(-2.9340194) q[0];
sx q[0];
rz(-1.3697778) q[0];
rz(2.6679084) q[2];
sx q[2];
rz(-1.2315183) q[2];
sx q[2];
rz(2.7850399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.42635065) q[1];
sx q[1];
rz(-0.75255005) q[1];
sx q[1];
rz(-2.1133711) q[1];
rz(1.5807125) q[3];
sx q[3];
rz(-2.2576233) q[3];
sx q[3];
rz(1.4894942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.131989) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(2.2930938) q[2];
rz(2.2961473) q[3];
sx q[3];
rz(-2.4852821) q[3];
sx q[3];
rz(1.9158624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41035143) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(0.89286667) q[0];
rz(-0.061575312) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(2.9796013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50349405) q[0];
sx q[0];
rz(-1.5553405) q[0];
sx q[0];
rz(-1.6729805) q[0];
rz(0.70190491) q[2];
sx q[2];
rz(-1.4830952) q[2];
sx q[2];
rz(-0.039393124) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74468915) q[1];
sx q[1];
rz(-1.2619881) q[1];
sx q[1];
rz(1.4267322) q[1];
rz(-pi) q[2];
rz(0.093858899) q[3];
sx q[3];
rz(-1.6742618) q[3];
sx q[3];
rz(-0.89573345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.050345) q[2];
sx q[2];
rz(-2.7140736) q[2];
sx q[2];
rz(-2.4166935) q[2];
rz(2.8643518) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8898833) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(3.1157893) q[0];
rz(1.39894) q[1];
sx q[1];
rz(-1.6394697) q[1];
sx q[1];
rz(0.10765156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9190034) q[0];
sx q[0];
rz(-2.1371142) q[0];
sx q[0];
rz(-2.983271) q[0];
rz(-2.669966) q[2];
sx q[2];
rz(-1.9666858) q[2];
sx q[2];
rz(-2.5079923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3947666) q[1];
sx q[1];
rz(-1.0977543) q[1];
sx q[1];
rz(1.6113144) q[1];
rz(-pi) q[2];
x q[2];
rz(0.017849542) q[3];
sx q[3];
rz(-1.4571725) q[3];
sx q[3];
rz(-0.85434948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.671635) q[2];
sx q[2];
rz(-1.8822972) q[2];
sx q[2];
rz(-2.0476445) q[2];
rz(2.5661902) q[3];
sx q[3];
rz(-0.83991528) q[3];
sx q[3];
rz(1.1168787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(-0.71304148) q[0];
rz(-0.22480045) q[1];
sx q[1];
rz(-2.5495922) q[1];
sx q[1];
rz(-2.0483268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9295392) q[0];
sx q[0];
rz(-1.0997547) q[0];
sx q[0];
rz(-0.77917288) q[0];
rz(-0.20168882) q[2];
sx q[2];
rz(-1.301577) q[2];
sx q[2];
rz(-1.2303011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8902407) q[1];
sx q[1];
rz(-1.6820434) q[1];
sx q[1];
rz(0.45808582) q[1];
x q[2];
rz(-1.3087464) q[3];
sx q[3];
rz(-0.42632494) q[3];
sx q[3];
rz(-2.654512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1405868) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(2.801753) q[2];
rz(-3.0991683) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(0.62780082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.4154631) q[2];
sx q[2];
rz(-2.6604322) q[2];
sx q[2];
rz(-0.12728035) q[2];
rz(1.9966765) q[3];
sx q[3];
rz(-1.3419587) q[3];
sx q[3];
rz(1.6423196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
