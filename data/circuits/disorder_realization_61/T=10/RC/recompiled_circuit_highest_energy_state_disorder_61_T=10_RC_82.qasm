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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(2.4866009) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(-1.8185599) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060299035) q[0];
sx q[0];
rz(-2.4954688) q[0];
sx q[0];
rz(2.6160014) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5147314) q[2];
sx q[2];
rz(-0.38598362) q[2];
sx q[2];
rz(0.71039334) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3016794) q[1];
sx q[1];
rz(-0.29407802) q[1];
sx q[1];
rz(-1.9211576) q[1];
rz(-pi) q[2];
rz(-0.58889525) q[3];
sx q[3];
rz(-2.3213904) q[3];
sx q[3];
rz(2.7429648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(-0.65832552) q[2];
rz(-0.63411921) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(-1.3707976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089791678) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(0.077089699) q[0];
rz(0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(1.9416521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.752992) q[0];
sx q[0];
rz(-1.3738828) q[0];
sx q[0];
rz(-1.1042751) q[0];
x q[1];
rz(1.7351772) q[2];
sx q[2];
rz(-1.5417678) q[2];
sx q[2];
rz(-0.2982225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5189831) q[1];
sx q[1];
rz(-0.72777343) q[1];
sx q[1];
rz(-2.1534797) q[1];
x q[2];
rz(-0.55910965) q[3];
sx q[3];
rz(-1.5190795) q[3];
sx q[3];
rz(2.8962108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2526907) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(0.67548951) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6460687) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(0.75463265) q[0];
rz(3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(-2.1307814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36463144) q[0];
sx q[0];
rz(-1.2662132) q[0];
sx q[0];
rz(1.0582032) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6057683) q[2];
sx q[2];
rz(-2.4596697) q[2];
sx q[2];
rz(-0.65820314) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52151187) q[1];
sx q[1];
rz(-2.3151868) q[1];
sx q[1];
rz(1.6537731) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95960404) q[3];
sx q[3];
rz(-2.2356847) q[3];
sx q[3];
rz(1.5184107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(-0.50673103) q[2];
rz(1.8690551) q[3];
sx q[3];
rz(-1.3576018) q[3];
sx q[3];
rz(-2.8633964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5469359) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(-2.0261672) q[0];
rz(1.2660654) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-0.87475264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0698358) q[0];
sx q[0];
rz(-1.4821782) q[0];
sx q[0];
rz(1.3184692) q[0];
rz(-pi) q[1];
x q[1];
rz(1.546193) q[2];
sx q[2];
rz(-2.4877791) q[2];
sx q[2];
rz(0.023651274) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1072448) q[1];
sx q[1];
rz(-0.49458359) q[1];
sx q[1];
rz(0.58652189) q[1];
rz(-pi) q[2];
rz(-0.86438365) q[3];
sx q[3];
rz(-1.4339851) q[3];
sx q[3];
rz(2.9446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11883417) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(-1.9258707) q[2];
rz(-1.7440965) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(-1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8656411) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.6708466) q[0];
rz(1.7631081) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.7313622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277692) q[0];
sx q[0];
rz(-0.54049379) q[0];
sx q[0];
rz(1.1590411) q[0];
x q[1];
rz(1.330014) q[2];
sx q[2];
rz(-2.2301939) q[2];
sx q[2];
rz(-1.3032152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8961337) q[1];
sx q[1];
rz(-1.2583548) q[1];
sx q[1];
rz(1.5691084) q[1];
x q[2];
rz(0.35378176) q[3];
sx q[3];
rz(-1.980034) q[3];
sx q[3];
rz(-1.080846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1269647) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(-1.8178168) q[2];
rz(1.7823559) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(-0.38715473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(3.1193745) q[0];
rz(-2.4413595) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(-0.95692316) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2882521) q[0];
sx q[0];
rz(-1.4256251) q[0];
sx q[0];
rz(0.70183825) q[0];
rz(-pi) q[1];
rz(0.21778743) q[2];
sx q[2];
rz(-1.7484669) q[2];
sx q[2];
rz(-1.9363994) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80920389) q[1];
sx q[1];
rz(-1.1706377) q[1];
sx q[1];
rz(0.2984115) q[1];
x q[2];
rz(2.2150008) q[3];
sx q[3];
rz(-2.6098688) q[3];
sx q[3];
rz(0.34429541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9243098) q[2];
sx q[2];
rz(-1.1223015) q[2];
sx q[2];
rz(-1.4420606) q[2];
rz(2.0349272) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(1.9895915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7149413) q[0];
sx q[0];
rz(-1.277667) q[0];
sx q[0];
rz(-2.6746124) q[0];
rz(-1.3106208) q[1];
sx q[1];
rz(-2.1881723) q[1];
sx q[1];
rz(-2.7511645) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6990427) q[0];
sx q[0];
rz(-0.14420284) q[0];
sx q[0];
rz(-0.22871916) q[0];
x q[1];
rz(0.10755928) q[2];
sx q[2];
rz(-1.2778408) q[2];
sx q[2];
rz(2.7358472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5891287) q[1];
sx q[1];
rz(-0.93163449) q[1];
sx q[1];
rz(2.5659849) q[1];
rz(-pi) q[2];
rz(-1.3620699) q[3];
sx q[3];
rz(-0.62297076) q[3];
sx q[3];
rz(-1.8546212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(-1.0170271) q[2];
rz(2.2022066) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(0.92923195) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4262714) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(-1.9287047) q[0];
rz(2.5384278) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0807665) q[0];
sx q[0];
rz(-2.2696324) q[0];
sx q[0];
rz(-3.0399975) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.064761) q[2];
sx q[2];
rz(-1.72121) q[2];
sx q[2];
rz(2.7422649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0871157) q[1];
sx q[1];
rz(-1.7810139) q[1];
sx q[1];
rz(1.3940548) q[1];
x q[2];
rz(-2.9717507) q[3];
sx q[3];
rz(-1.5780996) q[3];
sx q[3];
rz(0.23582349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6508871) q[2];
sx q[2];
rz(-2.8287973) q[2];
sx q[2];
rz(2.1631964) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(1.8470701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097718358) q[0];
sx q[0];
rz(-2.9053423) q[0];
sx q[0];
rz(-0.043721113) q[0];
rz(-1.1737191) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(2.3497605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39003644) q[0];
sx q[0];
rz(-1.1564213) q[0];
sx q[0];
rz(0.66585559) q[0];
x q[1];
rz(-0.20251198) q[2];
sx q[2];
rz(-1.2767451) q[2];
sx q[2];
rz(2.685315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89347311) q[1];
sx q[1];
rz(-1.9933102) q[1];
sx q[1];
rz(2.7759107) q[1];
rz(-pi) q[2];
rz(0.47747647) q[3];
sx q[3];
rz(-1.2315742) q[3];
sx q[3];
rz(0.71996688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23058471) q[2];
sx q[2];
rz(-0.30969301) q[2];
sx q[2];
rz(-1.2370375) q[2];
rz(2.6593995) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(-0.48875109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1173387) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(1.8035969) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(1.1504014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8088732) q[0];
sx q[0];
rz(-0.74902636) q[0];
sx q[0];
rz(2.3179884) q[0];
rz(-1.8878292) q[2];
sx q[2];
rz(-1.46508) q[2];
sx q[2];
rz(-2.6144165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.73455056) q[1];
sx q[1];
rz(-1.4311387) q[1];
sx q[1];
rz(1.1753083) q[1];
x q[2];
rz(-0.38693736) q[3];
sx q[3];
rz(-2.134857) q[3];
sx q[3];
rz(1.2880985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77356768) q[2];
sx q[2];
rz(-1.0057534) q[2];
sx q[2];
rz(-0.31039882) q[2];
rz(-2.5144905) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(-1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8695759) q[0];
sx q[0];
rz(-1.5740812) q[0];
sx q[0];
rz(1.5480702) q[0];
rz(1.9238453) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(-2.0704253) q[2];
sx q[2];
rz(-1.7555321) q[2];
sx q[2];
rz(2.1499014) q[2];
rz(-0.18219215) q[3];
sx q[3];
rz(-2.8145418) q[3];
sx q[3];
rz(-2.9180632) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
