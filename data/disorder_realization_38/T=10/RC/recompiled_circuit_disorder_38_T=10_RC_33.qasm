OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(0.24917319) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(-2.1698706) q[1];
sx q[1];
rz(0.5501774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30480706) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(0.027713393) q[0];
rz(-0.89994853) q[2];
sx q[2];
rz(-2.4953825) q[2];
sx q[2];
rz(2.3044569) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33221153) q[1];
sx q[1];
rz(-0.84377938) q[1];
sx q[1];
rz(-2.0903281) q[1];
rz(-pi) q[2];
x q[2];
rz(0.083505587) q[3];
sx q[3];
rz(-2.180035) q[3];
sx q[3];
rz(3.1374251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(-1.6289904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727901) q[0];
sx q[0];
rz(-1.8913942) q[0];
sx q[0];
rz(0.94706236) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(-0.46797215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.43803793) q[1];
sx q[1];
rz(-0.8982656) q[1];
sx q[1];
rz(0.1708252) q[1];
x q[2];
rz(-1.2849502) q[3];
sx q[3];
rz(-1.3405521) q[3];
sx q[3];
rz(-2.2765991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2597044) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-2.846068) q[0];
rz(0.23513901) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(0.74584109) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94507664) q[0];
sx q[0];
rz(-1.7917504) q[0];
sx q[0];
rz(-3.0483079) q[0];
x q[1];
rz(-0.52650555) q[2];
sx q[2];
rz(-1.7130934) q[2];
sx q[2];
rz(-1.7775747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44887603) q[1];
sx q[1];
rz(-1.2926896) q[1];
sx q[1];
rz(0.49726185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8196705) q[3];
sx q[3];
rz(-2.0119152) q[3];
sx q[3];
rz(-2.6114024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(3.0912494) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(-0.27329683) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647588) q[0];
sx q[0];
rz(-2.8059373) q[0];
sx q[0];
rz(1.2980952) q[0];
rz(2.4457473) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(2.3455182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4218688) q[1];
sx q[1];
rz(-2.8641652) q[1];
sx q[1];
rz(1.006702) q[1];
rz(-pi) q[2];
rz(-2.9076505) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(0.58594221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6325729) q[0];
sx q[0];
rz(-0.97826725) q[0];
sx q[0];
rz(2.65034) q[0];
rz(-1.5613902) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(0.26088342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8122711) q[1];
sx q[1];
rz(-1.8835888) q[1];
sx q[1];
rz(-2.4815464) q[1];
rz(1.6721252) q[3];
sx q[3];
rz(-0.39468995) q[3];
sx q[3];
rz(2.8916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(2.9303072) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(0.791839) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.2794367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6364481) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(1.5178773) q[0];
x q[1];
rz(-1.860414) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(-2.9786125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4931603) q[1];
sx q[1];
rz(-1.6006159) q[1];
sx q[1];
rz(-1.0489419) q[1];
rz(-pi) q[2];
rz(0.1412973) q[3];
sx q[3];
rz(-1.7501037) q[3];
sx q[3];
rz(-2.3395204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(2.9260013) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(3.1380222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647033) q[0];
sx q[0];
rz(-1.9386374) q[0];
sx q[0];
rz(-0.44643114) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9622757) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(-1.5228524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4184985) q[1];
sx q[1];
rz(-1.8685891) q[1];
sx q[1];
rz(0.94350068) q[1];
x q[2];
rz(0.076038578) q[3];
sx q[3];
rz(-1.7685316) q[3];
sx q[3];
rz(2.7969489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(0.92010951) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(2.4628941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7413901) q[0];
sx q[0];
rz(-0.42660248) q[0];
sx q[0];
rz(-2.461117) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0661725) q[2];
sx q[2];
rz(-0.80427158) q[2];
sx q[2];
rz(1.417516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2862257) q[1];
sx q[1];
rz(-1.507842) q[1];
sx q[1];
rz(1.0463868) q[1];
rz(-2.7160866) q[3];
sx q[3];
rz(-1.4262428) q[3];
sx q[3];
rz(1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(-1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(0.52694595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066756847) q[0];
sx q[0];
rz(-1.2196676) q[0];
sx q[0];
rz(1.3394651) q[0];
rz(1.2277463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(-2.3234141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2745797) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(0.32151476) q[1];
x q[2];
rz(-2.5832289) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(-0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(0.3607761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4496778) q[0];
sx q[0];
rz(-1.5458108) q[0];
sx q[0];
rz(1.4033532) q[0];
rz(-3.0214494) q[2];
sx q[2];
rz(-1.3573682) q[2];
sx q[2];
rz(-0.10211589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1211179) q[1];
sx q[1];
rz(-1.7133683) q[1];
sx q[1];
rz(-0.74417443) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4264614) q[3];
sx q[3];
rz(-2.9188041) q[3];
sx q[3];
rz(-2.1905394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(-1.3600596) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(1.3654937) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(-2.2684569) q[3];
sx q[3];
rz(-0.94948873) q[3];
sx q[3];
rz(-0.82617847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
