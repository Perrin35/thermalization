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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(-0.43098488) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(-2.7117742) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34944158) q[0];
sx q[0];
rz(-0.49956736) q[0];
sx q[0];
rz(-2.9527305) q[0];
x q[1];
rz(-0.016021803) q[2];
sx q[2];
rz(-1.1384769) q[2];
sx q[2];
rz(-0.95903083) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4861075) q[1];
sx q[1];
rz(-1.361651) q[1];
sx q[1];
rz(-2.7746592) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1898891) q[3];
sx q[3];
rz(-1.105068) q[3];
sx q[3];
rz(-1.4066175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0777883) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(2.1303614) q[2];
rz(-1.0412591) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083199) q[0];
sx q[0];
rz(-2.1450873) q[0];
sx q[0];
rz(0.74747768) q[0];
rz(0.29391089) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(-0.25516137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6601488) q[0];
sx q[0];
rz(-2.143476) q[0];
sx q[0];
rz(0.19788589) q[0];
x q[1];
rz(2.601971) q[2];
sx q[2];
rz(-0.35225062) q[2];
sx q[2];
rz(0.43708153) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0415062) q[1];
sx q[1];
rz(-1.1250682) q[1];
sx q[1];
rz(0.84560945) q[1];
rz(-pi) q[2];
rz(-0.18682602) q[3];
sx q[3];
rz(-2.1869724) q[3];
sx q[3];
rz(-2.0957859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1110288) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(2.2500706) q[2];
rz(-2.8920178) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7340649) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(-3.0974645) q[0];
rz(-1.8924425) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(-2.6557907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6225548) q[0];
sx q[0];
rz(-1.7960894) q[0];
sx q[0];
rz(1.5507429) q[0];
rz(-0.93636958) q[2];
sx q[2];
rz(-2.9108725) q[2];
sx q[2];
rz(1.3619193) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0016077) q[1];
sx q[1];
rz(-2.1473653) q[1];
sx q[1];
rz(-0.92553161) q[1];
rz(-pi) q[2];
rz(2.3084675) q[3];
sx q[3];
rz(-0.55748788) q[3];
sx q[3];
rz(-2.0035898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(2.5658549) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(1.7339138) q[3];
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
rz(-2.1908252) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(0.29275352) q[0];
rz(0.5689019) q[1];
sx q[1];
rz(-2.3141373) q[1];
sx q[1];
rz(-2.2946766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305188) q[0];
sx q[0];
rz(-1.9785709) q[0];
sx q[0];
rz(0.084534377) q[0];
rz(-0.90130536) q[2];
sx q[2];
rz(-1.6274823) q[2];
sx q[2];
rz(-2.0079119) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0140527) q[1];
sx q[1];
rz(-1.3779281) q[1];
sx q[1];
rz(2.8991153) q[1];
x q[2];
rz(-0.067518721) q[3];
sx q[3];
rz(-1.3775) q[3];
sx q[3];
rz(0.75261469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34951052) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(0.2087896) q[2];
rz(-2.4236692) q[3];
sx q[3];
rz(-2.8585298) q[3];
sx q[3];
rz(-0.94055241) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402886) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(-0.33863246) q[0];
rz(0.70308095) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(2.3032761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49946872) q[0];
sx q[0];
rz(-1.1422061) q[0];
sx q[0];
rz(1.8502251) q[0];
rz(-pi) q[1];
rz(-0.19836004) q[2];
sx q[2];
rz(-0.83399978) q[2];
sx q[2];
rz(1.5892346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5456713) q[1];
sx q[1];
rz(-1.6041846) q[1];
sx q[1];
rz(-2.2708793) q[1];
rz(-pi) q[2];
rz(-0.60689599) q[3];
sx q[3];
rz(-1.2147153) q[3];
sx q[3];
rz(1.0676366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.68359739) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(1.0028769) q[2];
rz(-2.9966127) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(-0.063044757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.018464) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(-0.96328324) q[0];
rz(-2.9604498) q[1];
sx q[1];
rz(-2.2434442) q[1];
sx q[1];
rz(1.2394261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2116263) q[0];
sx q[0];
rz(-2.3986882) q[0];
sx q[0];
rz(2.5725098) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.507429) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(-1.5786174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2378578) q[1];
sx q[1];
rz(-1.9457726) q[1];
sx q[1];
rz(1.1695678) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7224947) q[3];
sx q[3];
rz(-1.6983508) q[3];
sx q[3];
rz(2.2383245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2565101) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(-2.2086842) q[2];
rz(-1.7289303) q[3];
sx q[3];
rz(-1.7191073) q[3];
sx q[3];
rz(-2.9322114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4466208) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(-2.3387261) q[0];
rz(-2.39957) q[1];
sx q[1];
rz(-1.6525533) q[1];
sx q[1];
rz(-0.41025695) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2712532) q[0];
sx q[0];
rz(-1.2294719) q[0];
sx q[0];
rz(0.56296157) q[0];
rz(0.76426498) q[2];
sx q[2];
rz(-1.4857123) q[2];
sx q[2];
rz(-2.6972023) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8150543) q[1];
sx q[1];
rz(-1.6449303) q[1];
sx q[1];
rz(1.4085517) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0954082) q[3];
sx q[3];
rz(-0.96964004) q[3];
sx q[3];
rz(0.22829311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2731169) q[2];
sx q[2];
rz(-1.5831524) q[2];
sx q[2];
rz(2.9203019) q[2];
rz(-0.41915974) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(2.6656849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047423) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(1.6118443) q[0];
rz(-1.8086241) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.453368) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7573563) q[0];
sx q[0];
rz(-2.1290551) q[0];
sx q[0];
rz(-2.4008958) q[0];
rz(-pi) q[1];
rz(0.55744268) q[2];
sx q[2];
rz(-2.9176783) q[2];
sx q[2];
rz(0.73665184) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.61078) q[1];
sx q[1];
rz(-1.7208092) q[1];
sx q[1];
rz(0.71041815) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0734208) q[3];
sx q[3];
rz(-2.8262071) q[3];
sx q[3];
rz(2.634545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2754485) q[2];
sx q[2];
rz(-0.25442213) q[2];
sx q[2];
rz(-0.010738372) q[2];
rz(-2.8209316) q[3];
sx q[3];
rz(-1.843822) q[3];
sx q[3];
rz(0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(-2.6668715) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-2.0359998) q[1];
sx q[1];
rz(1.4795823) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6410165) q[0];
sx q[0];
rz(-1.5307284) q[0];
sx q[0];
rz(0.03682076) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3298673) q[2];
sx q[2];
rz(-2.2325667) q[2];
sx q[2];
rz(-2.2994201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7254527) q[1];
sx q[1];
rz(-1.3403157) q[1];
sx q[1];
rz(2.6621755) q[1];
rz(-1.817486) q[3];
sx q[3];
rz(-1.7462915) q[3];
sx q[3];
rz(-1.7481723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.097229615) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(-2.7105159) q[2];
rz(-2.9512682) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(0.88461191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97656074) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(-2.9227559) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(-1.7494019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4138654) q[0];
sx q[0];
rz(-1.5194335) q[0];
sx q[0];
rz(-2.9482909) q[0];
x q[1];
rz(2.7213412) q[2];
sx q[2];
rz(-1.4388348) q[2];
sx q[2];
rz(-0.38334639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6619685) q[1];
sx q[1];
rz(-1.2733409) q[1];
sx q[1];
rz(-0.65315078) q[1];
rz(-pi) q[2];
rz(2.7110215) q[3];
sx q[3];
rz(-1.0007676) q[3];
sx q[3];
rz(-1.4206629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.86768156) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(1.8217746) q[2];
rz(-0.53226081) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(1.5925315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5821447) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(-2.9136912) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(1.8310905) q[2];
sx q[2];
rz(-1.9734662) q[2];
sx q[2];
rz(-3.037813) q[2];
rz(1.3263973) q[3];
sx q[3];
rz(-0.59323372) q[3];
sx q[3];
rz(-2.9842971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
