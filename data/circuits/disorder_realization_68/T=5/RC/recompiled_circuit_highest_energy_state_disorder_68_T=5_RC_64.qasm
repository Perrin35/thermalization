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
rz(0.8402549) q[0];
sx q[0];
rz(4.1794887) q[0];
sx q[0];
rz(10.269796) q[0];
rz(0.71212274) q[1];
sx q[1];
rz(-2.1399838) q[1];
sx q[1];
rz(1.4955624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5069731) q[0];
sx q[0];
rz(-1.15043) q[0];
sx q[0];
rz(-2.2078321) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7579384) q[2];
sx q[2];
rz(-1.4315413) q[2];
sx q[2];
rz(-0.68472199) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6604267) q[1];
sx q[1];
rz(-1.1200953) q[1];
sx q[1];
rz(2.6066102) q[1];
x q[2];
rz(-0.59652416) q[3];
sx q[3];
rz(-2.6748195) q[3];
sx q[3];
rz(-1.0214361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40319765) q[2];
sx q[2];
rz(-1.1824111) q[2];
sx q[2];
rz(2.5851868) q[2];
rz(-0.49499908) q[3];
sx q[3];
rz(-0.35111108) q[3];
sx q[3];
rz(0.88465148) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86968386) q[0];
sx q[0];
rz(-0.10232919) q[0];
sx q[0];
rz(0.62865692) q[0];
rz(0.87720811) q[1];
sx q[1];
rz(-0.49867201) q[1];
sx q[1];
rz(-0.25310755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.077674) q[0];
sx q[0];
rz(-1.5627994) q[0];
sx q[0];
rz(1.5863938) q[0];
rz(-pi) q[1];
rz(-1.2727401) q[2];
sx q[2];
rz(-1.9949081) q[2];
sx q[2];
rz(-1.1828681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4483466) q[1];
sx q[1];
rz(-0.23138675) q[1];
sx q[1];
rz(2.6103002) q[1];
rz(1.289648) q[3];
sx q[3];
rz(-2.3901403) q[3];
sx q[3];
rz(-0.40159097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6796598) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(-2.0750462) q[2];
rz(-1.8343532) q[3];
sx q[3];
rz(-0.46590889) q[3];
sx q[3];
rz(-0.90827847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9268554) q[0];
sx q[0];
rz(-0.94811386) q[0];
sx q[0];
rz(0.98534775) q[0];
rz(-0.39380479) q[1];
sx q[1];
rz(-0.99291283) q[1];
sx q[1];
rz(0.52678144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7486289) q[0];
sx q[0];
rz(-0.53501463) q[0];
sx q[0];
rz(1.4070542) q[0];
rz(-pi) q[1];
rz(-2.5622732) q[2];
sx q[2];
rz(-1.123482) q[2];
sx q[2];
rz(-1.2944702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.078048) q[1];
sx q[1];
rz(-1.1495665) q[1];
sx q[1];
rz(-0.18555141) q[1];
x q[2];
rz(0.44528499) q[3];
sx q[3];
rz(-2.1213343) q[3];
sx q[3];
rz(2.8360644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7963316) q[2];
sx q[2];
rz(-2.3806206) q[2];
sx q[2];
rz(-2.995028) q[2];
rz(2.2740299) q[3];
sx q[3];
rz(-0.87351322) q[3];
sx q[3];
rz(0.7578907) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72916156) q[0];
sx q[0];
rz(-0.49598345) q[0];
sx q[0];
rz(-2.5878986) q[0];
rz(1.4750534) q[1];
sx q[1];
rz(-2.8429884) q[1];
sx q[1];
rz(-0.24436229) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15904832) q[0];
sx q[0];
rz(-1.9031018) q[0];
sx q[0];
rz(2.766702) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8342488) q[2];
sx q[2];
rz(-1.3400199) q[2];
sx q[2];
rz(2.1342056) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8148538) q[1];
sx q[1];
rz(-2.5757901) q[1];
sx q[1];
rz(2.4292355) q[1];
rz(-pi) q[2];
rz(-3.1126106) q[3];
sx q[3];
rz(-0.71235114) q[3];
sx q[3];
rz(2.1952352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0814521) q[2];
sx q[2];
rz(-1.1271366) q[2];
sx q[2];
rz(-2.3606908) q[2];
rz(2.7447356) q[3];
sx q[3];
rz(-0.01440993) q[3];
sx q[3];
rz(1.074033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9020554) q[0];
sx q[0];
rz(-2.9750415) q[0];
sx q[0];
rz(-0.2567513) q[0];
rz(-2.9638839) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(-2.1977052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0487918) q[0];
sx q[0];
rz(-1.208713) q[0];
sx q[0];
rz(1.0975773) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2805528) q[2];
sx q[2];
rz(-1.3849045) q[2];
sx q[2];
rz(-3.0772532) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1348995) q[1];
sx q[1];
rz(-2.1104782) q[1];
sx q[1];
rz(1.5758908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6750608) q[3];
sx q[3];
rz(-1.0634907) q[3];
sx q[3];
rz(0.093884163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.93861598) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(-0.20982783) q[2];
rz(-0.046791568) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(-2.7941424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059134722) q[0];
sx q[0];
rz(-0.81675285) q[0];
sx q[0];
rz(-1.2030075) q[0];
rz(-1.6251534) q[1];
sx q[1];
rz(-0.37767437) q[1];
sx q[1];
rz(3.1086521) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66844475) q[0];
sx q[0];
rz(-2.0347036) q[0];
sx q[0];
rz(2.8345246) q[0];
x q[1];
rz(0.85196544) q[2];
sx q[2];
rz(-2.4174066) q[2];
sx q[2];
rz(-2.9018096) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9666933) q[1];
sx q[1];
rz(-2.1191349) q[1];
sx q[1];
rz(1.95005) q[1];
rz(-0.88992702) q[3];
sx q[3];
rz(-2.0924727) q[3];
sx q[3];
rz(-1.8311178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5861627) q[2];
sx q[2];
rz(-1.0032434) q[2];
sx q[2];
rz(0.45247751) q[2];
rz(2.9943976) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(-1.8424621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066605695) q[0];
sx q[0];
rz(-0.43656483) q[0];
sx q[0];
rz(2.7845352) q[0];
rz(-1.0128516) q[1];
sx q[1];
rz(-1.169299) q[1];
sx q[1];
rz(-3.1049407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3507639) q[0];
sx q[0];
rz(-2.0022656) q[0];
sx q[0];
rz(-0.58833265) q[0];
rz(1.4674856) q[2];
sx q[2];
rz(-1.1432308) q[2];
sx q[2];
rz(-0.72139886) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.086163047) q[1];
sx q[1];
rz(-2.2698672) q[1];
sx q[1];
rz(0.12075874) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5753079) q[3];
sx q[3];
rz(-1.5399714) q[3];
sx q[3];
rz(2.2982321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4770294) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(-3.0518517) q[2];
rz(-1.2612032) q[3];
sx q[3];
rz(-1.3124876) q[3];
sx q[3];
rz(1.4242273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4795714) q[0];
sx q[0];
rz(-0.57858545) q[0];
sx q[0];
rz(-1.4805967) q[0];
rz(1.4501694) q[1];
sx q[1];
rz(-0.65933508) q[1];
sx q[1];
rz(0.75993842) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1131033) q[0];
sx q[0];
rz(-3.0249502) q[0];
sx q[0];
rz(1.6595429) q[0];
x q[1];
rz(-0.6282232) q[2];
sx q[2];
rz(-1.6464273) q[2];
sx q[2];
rz(0.13720195) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2719136) q[1];
sx q[1];
rz(-1.9909503) q[1];
sx q[1];
rz(-2.7019482) q[1];
x q[2];
rz(2.6989903) q[3];
sx q[3];
rz(-2.5902904) q[3];
sx q[3];
rz(-2.9568903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30457589) q[2];
sx q[2];
rz(-1.872007) q[2];
sx q[2];
rz(-0.059583511) q[2];
rz(-2.7130821) q[3];
sx q[3];
rz(-2.2969552) q[3];
sx q[3];
rz(0.60215157) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8932327) q[0];
sx q[0];
rz(-0.28286523) q[0];
sx q[0];
rz(0.34348139) q[0];
rz(0.19613014) q[1];
sx q[1];
rz(-0.88534147) q[1];
sx q[1];
rz(-2.3416065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.26639) q[0];
sx q[0];
rz(-0.81884844) q[0];
sx q[0];
rz(-2.1455392) q[0];
x q[1];
rz(0.38651379) q[2];
sx q[2];
rz(-2.5441351) q[2];
sx q[2];
rz(0.51624417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1175122) q[1];
sx q[1];
rz(-1.5264866) q[1];
sx q[1];
rz(-2.5818392) q[1];
x q[2];
rz(0.88801791) q[3];
sx q[3];
rz(-1.6245442) q[3];
sx q[3];
rz(-2.2537504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34667748) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(0.22870341) q[2];
rz(1.4165357) q[3];
sx q[3];
rz(-1.718947) q[3];
sx q[3];
rz(-2.4712839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5423841) q[0];
sx q[0];
rz(-2.6238361) q[0];
sx q[0];
rz(-1.9589348) q[0];
rz(0.54284894) q[1];
sx q[1];
rz(-3.019637) q[1];
sx q[1];
rz(-0.45546946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96129721) q[0];
sx q[0];
rz(-1.1396798) q[0];
sx q[0];
rz(2.6102553) q[0];
x q[1];
rz(-1.1588016) q[2];
sx q[2];
rz(-0.48589009) q[2];
sx q[2];
rz(-2.0249572) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23697904) q[1];
sx q[1];
rz(-2.0354668) q[1];
sx q[1];
rz(1.418271) q[1];
x q[2];
rz(-2.4157468) q[3];
sx q[3];
rz(-2.5739087) q[3];
sx q[3];
rz(2.3633702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21696422) q[2];
sx q[2];
rz(-2.2483726) q[2];
sx q[2];
rz(-0.42579892) q[2];
rz(2.9567772) q[3];
sx q[3];
rz(-2.5033689) q[3];
sx q[3];
rz(-0.027801175) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4375147) q[0];
sx q[0];
rz(-1.5298433) q[0];
sx q[0];
rz(-0.035506305) q[0];
rz(2.6620445) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(0.32301474) q[2];
sx q[2];
rz(-2.5217006) q[2];
sx q[2];
rz(2.9519242) q[2];
rz(0.95994031) q[3];
sx q[3];
rz(-1.452594) q[3];
sx q[3];
rz(0.95127524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
