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
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(-2.5133361) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4391651) q[0];
sx q[0];
rz(-1.3936696) q[0];
sx q[0];
rz(-0.11985531) q[0];
x q[1];
rz(-1.9626748) q[2];
sx q[2];
rz(-0.88681839) q[2];
sx q[2];
rz(-0.98525233) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3648277) q[1];
sx q[1];
rz(-1.3421632) q[1];
sx q[1];
rz(1.5336114) q[1];
rz(-1.0367582) q[3];
sx q[3];
rz(-1.9566571) q[3];
sx q[3];
rz(-0.88632583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85311741) q[2];
sx q[2];
rz(-1.709781) q[2];
sx q[2];
rz(0.34118578) q[2];
rz(-2.2513385) q[3];
sx q[3];
rz(-1.8743926) q[3];
sx q[3];
rz(-0.53513479) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5219236) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(0.71440119) q[0];
rz(-1.5419386) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(0.094706789) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15078397) q[0];
sx q[0];
rz(-1.8903194) q[0];
sx q[0];
rz(-0.099415778) q[0];
x q[1];
rz(2.8049967) q[2];
sx q[2];
rz(-2.1456652) q[2];
sx q[2];
rz(0.47913715) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.73037877) q[1];
sx q[1];
rz(-1.8068562) q[1];
sx q[1];
rz(-0.35309955) q[1];
rz(-pi) q[2];
rz(2.6323039) q[3];
sx q[3];
rz(-1.7096161) q[3];
sx q[3];
rz(1.8880488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16227214) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(0.47404131) q[2];
rz(0.76215172) q[3];
sx q[3];
rz(-0.63665974) q[3];
sx q[3];
rz(2.6704085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2786461) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-2.2886724) q[0];
rz(0.22311738) q[1];
sx q[1];
rz(-0.47839636) q[1];
sx q[1];
rz(-2.6257637) q[1];
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
rz(2.8858809) q[2];
sx q[2];
rz(-2.1923794) q[2];
sx q[2];
rz(-2.5598124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3577376) q[1];
sx q[1];
rz(-0.92354316) q[1];
sx q[1];
rz(-2.3308862) q[1];
rz(-pi) q[2];
rz(-2.7985057) q[3];
sx q[3];
rz(-2.2394387) q[3];
sx q[3];
rz(-2.9647788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8082661) q[2];
sx q[2];
rz(-0.9767248) q[2];
sx q[2];
rz(-0.094956368) q[2];
rz(-0.44148463) q[3];
sx q[3];
rz(-1.7850826) q[3];
sx q[3];
rz(-1.1536185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3917291) q[0];
sx q[0];
rz(-2.5666105) q[0];
sx q[0];
rz(-2.0854501) q[0];
rz(-2.1978343) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(-1.4518849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7565931) q[0];
sx q[0];
rz(-0.90314048) q[0];
sx q[0];
rz(0.7636015) q[0];
rz(-2.4511172) q[2];
sx q[2];
rz(-0.5734517) q[2];
sx q[2];
rz(0.8586463) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6744212) q[1];
sx q[1];
rz(-1.3302263) q[1];
sx q[1];
rz(-1.6075385) q[1];
x q[2];
rz(-0.068030997) q[3];
sx q[3];
rz(-3.0306731) q[3];
sx q[3];
rz(-2.7601931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(2.441414) q[2];
rz(0.87107301) q[3];
sx q[3];
rz(-1.0659404) q[3];
sx q[3];
rz(-1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313318) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(-0.48466551) q[0];
rz(2.2635745) q[1];
sx q[1];
rz(-1.9819825) q[1];
sx q[1];
rz(-0.39883089) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8190552) q[0];
sx q[0];
rz(-0.56542552) q[0];
sx q[0];
rz(-0.25382332) q[0];
x q[1];
rz(0.39426784) q[2];
sx q[2];
rz(-2.3986001) q[2];
sx q[2];
rz(-2.791648) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4927401) q[1];
sx q[1];
rz(-1.2034186) q[1];
sx q[1];
rz(2.2671764) q[1];
x q[2];
rz(-0.871931) q[3];
sx q[3];
rz(-2.0060349) q[3];
sx q[3];
rz(2.1520681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6681119) q[2];
sx q[2];
rz(-1.5405737) q[2];
sx q[2];
rz(0.87781805) q[2];
rz(2.7569438) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(1.1788684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28410742) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(0.2704764) q[0];
rz(-2.8349304) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(-2.564863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79918843) q[0];
sx q[0];
rz(-1.4831838) q[0];
sx q[0];
rz(1.3848928) q[0];
x q[1];
rz(0.36251601) q[2];
sx q[2];
rz(-0.91372817) q[2];
sx q[2];
rz(0.18489472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57177793) q[1];
sx q[1];
rz(-1.0228436) q[1];
sx q[1];
rz(0.46787365) q[1];
rz(-1.5567729) q[3];
sx q[3];
rz(-1.2360667) q[3];
sx q[3];
rz(-2.1317792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7755255) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(2.1712187) q[2];
rz(-2.7905285) q[3];
sx q[3];
rz(-2.909436) q[3];
sx q[3];
rz(0.27950132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
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
rz(-0.84545285) q[1];
sx q[1];
rz(-2.2038584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739627) q[0];
sx q[0];
rz(-1.7741307) q[0];
sx q[0];
rz(-3.0995661) q[0];
x q[1];
rz(0.47368424) q[2];
sx q[2];
rz(-1.2315183) q[2];
sx q[2];
rz(-2.7850399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8775505) q[1];
sx q[1];
rz(-2.1960947) q[1];
sx q[1];
rz(0.45035119) q[1];
rz(0.012091919) q[3];
sx q[3];
rz(-2.4547057) q[3];
sx q[3];
rz(1.6364607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.131989) q[2];
sx q[2];
rz(-0.14543532) q[2];
sx q[2];
rz(-2.2930938) q[2];
rz(-2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(1.9158624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41035143) q[0];
sx q[0];
rz(-2.1827965) q[0];
sx q[0];
rz(-0.89286667) q[0];
rz(-0.061575312) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(2.9796013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50349405) q[0];
sx q[0];
rz(-1.5862521) q[0];
sx q[0];
rz(-1.6729805) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0062469) q[2];
sx q[2];
rz(-2.4351546) q[2];
sx q[2];
rz(-1.634623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74468915) q[1];
sx q[1];
rz(-1.8796046) q[1];
sx q[1];
rz(-1.4267322) q[1];
x q[2];
rz(-1.6747159) q[3];
sx q[3];
rz(-1.6641518) q[3];
sx q[3];
rz(0.68478497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25170931) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(-0.025803331) q[0];
rz(-1.7426527) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(-0.10765156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336595) q[0];
sx q[0];
rz(-1.4373533) q[0];
sx q[0];
rz(-2.1428277) q[0];
rz(-pi) q[1];
rz(-2.669966) q[2];
sx q[2];
rz(-1.1749069) q[2];
sx q[2];
rz(2.5079923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15756179) q[1];
sx q[1];
rz(-1.6068629) q[1];
sx q[1];
rz(0.47337516) q[1];
rz(-pi) q[2];
rz(0.017849542) q[3];
sx q[3];
rz(-1.4571725) q[3];
sx q[3];
rz(-0.85434948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.671635) q[2];
sx q[2];
rz(-1.8822972) q[2];
sx q[2];
rz(1.0939481) q[2];
rz(2.5661902) q[3];
sx q[3];
rz(-0.83991528) q[3];
sx q[3];
rz(-2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(2.4285512) q[0];
rz(0.22480045) q[1];
sx q[1];
rz(-2.5495922) q[1];
sx q[1];
rz(2.0483268) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2042052) q[0];
sx q[0];
rz(-2.2474216) q[0];
sx q[0];
rz(-0.94955523) q[0];
rz(-1.29628) q[2];
sx q[2];
rz(-1.3764672) q[2];
sx q[2];
rz(0.39482612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0436127) q[1];
sx q[1];
rz(-2.671125) q[1];
sx q[1];
rz(-0.2474252) q[1];
rz(-1.8328463) q[3];
sx q[3];
rz(-2.7152677) q[3];
sx q[3];
rz(0.48708068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0010058086) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(-2.801753) q[2];
rz(0.042424399) q[3];
sx q[3];
rz(-1.8570615) q[3];
sx q[3];
rz(-0.62780082) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18702678) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(0.78202248) q[1];
sx q[1];
rz(-2.2047058) q[1];
sx q[1];
rz(2.4155736) q[1];
rz(-0.72612957) q[2];
sx q[2];
rz(-0.48116044) q[2];
sx q[2];
rz(3.0143123) q[2];
rz(-0.25039582) q[3];
sx q[3];
rz(-1.156711) q[3];
sx q[3];
rz(-0.031022978) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
