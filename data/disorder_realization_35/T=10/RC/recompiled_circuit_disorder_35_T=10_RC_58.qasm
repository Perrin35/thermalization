OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(-2.536474) q[1];
sx q[1];
rz(-2.6095698) q[1];
sx q[1];
rz(-1.1693118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7919851) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(2.4277359) q[0];
rz(-pi) q[1];
rz(0.77180441) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(-2.8231205) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77820233) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(1.1671288) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23006769) q[3];
sx q[3];
rz(-2.821273) q[3];
sx q[3];
rz(0.39428082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-2.1038726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1514725) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-2.5734076) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77387626) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(0.80640031) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91451007) q[1];
sx q[1];
rz(-2.7343379) q[1];
sx q[1];
rz(-1.6362908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5492937) q[3];
sx q[3];
rz(-2.2364738) q[3];
sx q[3];
rz(-0.75331068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304853) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7645435) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(-0.069805108) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77261749) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(2.2881743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55858559) q[1];
sx q[1];
rz(-1.188394) q[1];
sx q[1];
rz(1.26182) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40804789) q[3];
sx q[3];
rz(-1.1432262) q[3];
sx q[3];
rz(-1.5023295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-2.5909246) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6837316) q[0];
sx q[0];
rz(-2.3944003) q[0];
sx q[0];
rz(-1.9283717) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0998459) q[2];
sx q[2];
rz(-3.1051271) q[2];
sx q[2];
rz(1.0066102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1945222) q[1];
sx q[1];
rz(-2.2010989) q[1];
sx q[1];
rz(1.6987726) q[1];
rz(3.0691931) q[3];
sx q[3];
rz(-0.94745938) q[3];
sx q[3];
rz(2.6423531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.682122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5495816) q[0];
sx q[0];
rz(-1.3322543) q[0];
sx q[0];
rz(-2.9907988) q[0];
x q[1];
rz(0.51909165) q[2];
sx q[2];
rz(-0.95538288) q[2];
sx q[2];
rz(-1.5450684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1294714) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(2.4232037) q[1];
x q[2];
rz(-2.8206283) q[3];
sx q[3];
rz(-1.1043613) q[3];
sx q[3];
rz(2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(-0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(3.0029283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667144) q[0];
sx q[0];
rz(-0.020665558) q[0];
sx q[0];
rz(2.979435) q[0];
rz(-pi) q[1];
rz(1.4939098) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(2.4424057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2648592) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(-0.77244669) q[1];
rz(-pi) q[2];
rz(1.9964553) q[3];
sx q[3];
rz(-2.2666551) q[3];
sx q[3];
rz(1.3130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3665294) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(1.1602317) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(0.53652525) q[0];
rz(0.58553186) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-0.62430635) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3793959) q[0];
sx q[0];
rz(-1.3610098) q[0];
sx q[0];
rz(2.4136153) q[0];
rz(-pi) q[1];
rz(0.46632669) q[2];
sx q[2];
rz(-1.000324) q[2];
sx q[2];
rz(-1.8898659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1075322) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(0.88453102) q[1];
x q[2];
rz(-0.47847139) q[3];
sx q[3];
rz(-1.0020743) q[3];
sx q[3];
rz(2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6162993) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-2.7590511) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0807954) q[0];
sx q[0];
rz(-1.1757438) q[0];
sx q[0];
rz(0.11443826) q[0];
rz(2.244197) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(-0.71133864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8184549) q[1];
sx q[1];
rz(-1.8433246) q[1];
sx q[1];
rz(-2.6998991) q[1];
x q[2];
rz(0.73980476) q[3];
sx q[3];
rz(-1.5590258) q[3];
sx q[3];
rz(-1.7108325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-2.074266) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.8918442) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352167) q[0];
sx q[0];
rz(-1.750964) q[0];
sx q[0];
rz(0.10794497) q[0];
rz(1.1662912) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(1.5664958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2495888) q[1];
sx q[1];
rz(-0.87000404) q[1];
sx q[1];
rz(2.4904817) q[1];
x q[2];
rz(-2.8899367) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(0.8151527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0525557) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(0.54668033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0944259) q[0];
sx q[0];
rz(-1.7921653) q[0];
sx q[0];
rz(0.84035994) q[0];
rz(0.28995138) q[2];
sx q[2];
rz(-0.49294127) q[2];
sx q[2];
rz(1.8268367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(-1.5323557) q[1];
rz(1.7095079) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(2.9885938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13835779) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.2658723) q[2];
sx q[2];
rz(-3.0921428) q[2];
sx q[2];
rz(-2.5947239) q[2];
rz(2.1605282) q[3];
sx q[3];
rz(-2.5909501) q[3];
sx q[3];
rz(0.29028374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];