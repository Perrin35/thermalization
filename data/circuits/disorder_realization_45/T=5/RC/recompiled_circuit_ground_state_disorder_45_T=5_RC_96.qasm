OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0319808) q[0];
sx q[0];
rz(4.1469753) q[0];
sx q[0];
rz(8.4561705) q[0];
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(3.1117575) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.077622) q[0];
sx q[0];
rz(-1.0322303) q[0];
sx q[0];
rz(1.2993815) q[0];
x q[1];
rz(1.4953243) q[2];
sx q[2];
rz(-2.3879037) q[2];
sx q[2];
rz(-1.8729295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7368841) q[1];
sx q[1];
rz(-0.88062693) q[1];
sx q[1];
rz(2.0815013) q[1];
rz(-pi) q[2];
rz(0.13600341) q[3];
sx q[3];
rz(-1.4823071) q[3];
sx q[3];
rz(-1.1565151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6385103) q[2];
sx q[2];
rz(-0.57380399) q[2];
sx q[2];
rz(2.8199675) q[2];
rz(2.5173748) q[3];
sx q[3];
rz(-1.7523242) q[3];
sx q[3];
rz(-1.483884) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6931848) q[0];
sx q[0];
rz(-2.8852561) q[0];
sx q[0];
rz(0.93628991) q[0];
rz(1.2721277) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(1.3880656) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427338) q[0];
sx q[0];
rz(-1.6279477) q[0];
sx q[0];
rz(1.3512035) q[0];
rz(-0.77735591) q[2];
sx q[2];
rz(-1.1749975) q[2];
sx q[2];
rz(2.0344025) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1202045) q[1];
sx q[1];
rz(-0.98155456) q[1];
sx q[1];
rz(-0.91597775) q[1];
x q[2];
rz(-0.8901435) q[3];
sx q[3];
rz(-0.33484866) q[3];
sx q[3];
rz(-2.072842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0973413) q[2];
sx q[2];
rz(-0.96434957) q[2];
sx q[2];
rz(0.67330366) q[2];
rz(-1.3168969) q[3];
sx q[3];
rz(-1.293106) q[3];
sx q[3];
rz(-2.3015658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8116542) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(3.0918308) q[0];
rz(-2.9352169) q[1];
sx q[1];
rz(-1.3866837) q[1];
sx q[1];
rz(1.4899563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0127209) q[0];
sx q[0];
rz(-1.6908592) q[0];
sx q[0];
rz(-2.524656) q[0];
x q[1];
rz(0.55402886) q[2];
sx q[2];
rz(-1.1893335) q[2];
sx q[2];
rz(1.2020401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4088285) q[1];
sx q[1];
rz(-1.8547426) q[1];
sx q[1];
rz(2.8676621) q[1];
x q[2];
rz(0.95155119) q[3];
sx q[3];
rz(-0.54446044) q[3];
sx q[3];
rz(-1.8306554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18996198) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(-3.087431) q[2];
rz(-1.1545898) q[3];
sx q[3];
rz(-1.8535987) q[3];
sx q[3];
rz(-2.0258928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.7390249) q[0];
sx q[0];
rz(-1.4250647) q[0];
sx q[0];
rz(2.2836852) q[0];
rz(-0.67445451) q[1];
sx q[1];
rz(-2.1916788) q[1];
sx q[1];
rz(-1.3309006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895578) q[0];
sx q[0];
rz(-1.7072132) q[0];
sx q[0];
rz(2.9035764) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9711645) q[2];
sx q[2];
rz(-1.870451) q[2];
sx q[2];
rz(-1.48502) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.808271) q[1];
sx q[1];
rz(-2.1262319) q[1];
sx q[1];
rz(2.9395648) q[1];
rz(-0.0033802948) q[3];
sx q[3];
rz(-1.1515318) q[3];
sx q[3];
rz(1.5730891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7669507) q[2];
sx q[2];
rz(-0.16051126) q[2];
sx q[2];
rz(-0.10425076) q[2];
rz(-2.5637964) q[3];
sx q[3];
rz(-2.0640852) q[3];
sx q[3];
rz(0.92208636) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0084956) q[0];
sx q[0];
rz(-1.47559) q[0];
sx q[0];
rz(2.1971028) q[0];
rz(-1.7780444) q[1];
sx q[1];
rz(-1.8133769) q[1];
sx q[1];
rz(-1.4469226) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5529842) q[0];
sx q[0];
rz(-2.406334) q[0];
sx q[0];
rz(-1.8553493) q[0];
x q[1];
rz(3.1051272) q[2];
sx q[2];
rz(-1.8210095) q[2];
sx q[2];
rz(1.6584876) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80396485) q[1];
sx q[1];
rz(-1.9586438) q[1];
sx q[1];
rz(-0.66168534) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78131494) q[3];
sx q[3];
rz(-2.4014946) q[3];
sx q[3];
rz(-0.61352713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1030582) q[2];
sx q[2];
rz(-2.6699622) q[2];
sx q[2];
rz(-0.72059694) q[2];
rz(-0.12901148) q[3];
sx q[3];
rz(-1.7804264) q[3];
sx q[3];
rz(1.8425997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98220372) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(-2.5250315) q[0];
rz(-2.189134) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(-1.8985101) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4774507) q[0];
sx q[0];
rz(-1.7074025) q[0];
sx q[0];
rz(-1.3090429) q[0];
x q[1];
rz(2.277563) q[2];
sx q[2];
rz(-1.2507157) q[2];
sx q[2];
rz(0.79369407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4259498) q[1];
sx q[1];
rz(-2.2897807) q[1];
sx q[1];
rz(1.4904316) q[1];
rz(1.467435) q[3];
sx q[3];
rz(-0.62400104) q[3];
sx q[3];
rz(2.0614908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9055966) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(-0.074404152) q[2];
rz(1.0807886) q[3];
sx q[3];
rz(-1.6483043) q[3];
sx q[3];
rz(1.0409482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8393132) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(0.062051274) q[0];
rz(1.4150367) q[1];
sx q[1];
rz(-0.77295417) q[1];
sx q[1];
rz(-0.089990377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22253972) q[0];
sx q[0];
rz(-0.74680674) q[0];
sx q[0];
rz(-1.9487621) q[0];
x q[1];
rz(-2.9758866) q[2];
sx q[2];
rz(-2.2519654) q[2];
sx q[2];
rz(-0.15128862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9518826) q[1];
sx q[1];
rz(-1.2219349) q[1];
sx q[1];
rz(2.2079218) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2136647) q[3];
sx q[3];
rz(-1.0786329) q[3];
sx q[3];
rz(-1.6097691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0061079582) q[2];
sx q[2];
rz(-1.2697479) q[2];
sx q[2];
rz(0.6832068) q[2];
rz(-1.3225383) q[3];
sx q[3];
rz(-2.1893978) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6739864) q[0];
sx q[0];
rz(-1.645393) q[0];
sx q[0];
rz(2.6737387) q[0];
rz(2.4521008) q[1];
sx q[1];
rz(-0.71527022) q[1];
sx q[1];
rz(0.47324866) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9625585) q[0];
sx q[0];
rz(-2.9518572) q[0];
sx q[0];
rz(1.9535854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32903657) q[2];
sx q[2];
rz(-0.46571443) q[2];
sx q[2];
rz(0.28458111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93750226) q[1];
sx q[1];
rz(-3.1181304) q[1];
sx q[1];
rz(-0.84659441) q[1];
rz(-1.3159479) q[3];
sx q[3];
rz(-1.8425966) q[3];
sx q[3];
rz(2.013828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1271992) q[2];
sx q[2];
rz(-1.3213804) q[2];
sx q[2];
rz(-0.42352208) q[2];
rz(1.8262156) q[3];
sx q[3];
rz(-0.49153057) q[3];
sx q[3];
rz(-1.2929085) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92726436) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(-1.3182053) q[1];
sx q[1];
rz(-2.2315836) q[1];
sx q[1];
rz(0.013890161) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0158169) q[0];
sx q[0];
rz(-2.7770574) q[0];
sx q[0];
rz(1.6423854) q[0];
rz(-pi) q[1];
rz(-0.35160323) q[2];
sx q[2];
rz(-1.7187864) q[2];
sx q[2];
rz(-2.6949712) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74456066) q[1];
sx q[1];
rz(-0.63493997) q[1];
sx q[1];
rz(-1.145361) q[1];
rz(-3.0180279) q[3];
sx q[3];
rz(-0.84513226) q[3];
sx q[3];
rz(2.7853876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58005971) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(-1.053099) q[2];
rz(-0.92755353) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(0.48228669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841566) q[0];
sx q[0];
rz(-0.3967537) q[0];
sx q[0];
rz(-2.5411153) q[0];
rz(-0.26957574) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(1.762766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81299128) q[0];
sx q[0];
rz(-2.6909851) q[0];
sx q[0];
rz(-2.2297321) q[0];
rz(-pi) q[1];
rz(2.7711033) q[2];
sx q[2];
rz(-1.8571808) q[2];
sx q[2];
rz(3.1305673) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8104756) q[1];
sx q[1];
rz(-1.1441191) q[1];
sx q[1];
rz(1.7716549) q[1];
rz(-2.5559101) q[3];
sx q[3];
rz(-2.437311) q[3];
sx q[3];
rz(-1.3552988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1158585) q[2];
sx q[2];
rz(-0.33418843) q[2];
sx q[2];
rz(0.96246976) q[2];
rz(0.95364237) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(1.3941221) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3772603) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(-1.3537021) q[1];
sx q[1];
rz(-1.1693015) q[1];
sx q[1];
rz(2.404626) q[1];
rz(-2.216966) q[2];
sx q[2];
rz(-1.007111) q[2];
sx q[2];
rz(2.5180857) q[2];
rz(-0.24772027) q[3];
sx q[3];
rz(-0.52023028) q[3];
sx q[3];
rz(0.46432555) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
