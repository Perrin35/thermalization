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
rz(-2.4294699) q[1];
sx q[1];
rz(-1.0016088) q[1];
sx q[1];
rz(1.6460302) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43216773) q[0];
sx q[0];
rz(-2.3948095) q[0];
sx q[0];
rz(0.92632697) q[0];
x q[1];
rz(1.4207815) q[2];
sx q[2];
rz(-1.9505462) q[2];
sx q[2];
rz(0.94204547) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.54421669) q[1];
sx q[1];
rz(-2.4565182) q[1];
sx q[1];
rz(2.3822576) q[1];
rz(-pi) q[2];
rz(-2.5450685) q[3];
sx q[3];
rz(-2.6748195) q[3];
sx q[3];
rz(1.0214361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.738395) q[2];
sx q[2];
rz(-1.9591816) q[2];
sx q[2];
rz(-2.5851868) q[2];
rz(-2.6465936) q[3];
sx q[3];
rz(-2.7904816) q[3];
sx q[3];
rz(-2.2569412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.2719088) q[0];
sx q[0];
rz(-0.10232919) q[0];
sx q[0];
rz(-2.5129357) q[0];
rz(-0.87720811) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(2.8884851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0639187) q[0];
sx q[0];
rz(-1.5627994) q[0];
sx q[0];
rz(1.5863938) q[0];
rz(-1.8688525) q[2];
sx q[2];
rz(-1.1466845) q[2];
sx q[2];
rz(-1.1828681) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.991515) q[1];
sx q[1];
rz(-1.7698231) q[1];
sx q[1];
rz(1.4519889) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8879426) q[3];
sx q[3];
rz(-0.85542711) q[3];
sx q[3];
rz(2.3636092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6796598) q[2];
sx q[2];
rz(-1.9340645) q[2];
sx q[2];
rz(1.0665464) q[2];
rz(1.3072394) q[3];
sx q[3];
rz(-0.46590889) q[3];
sx q[3];
rz(-0.90827847) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9268554) q[0];
sx q[0];
rz(-0.94811386) q[0];
sx q[0];
rz(2.1562449) q[0];
rz(-2.7477879) q[1];
sx q[1];
rz(-2.1486798) q[1];
sx q[1];
rz(0.52678144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0366448) q[0];
sx q[0];
rz(-1.6540043) q[0];
sx q[0];
rz(1.0416688) q[0];
x q[1];
rz(0.57931945) q[2];
sx q[2];
rz(-1.123482) q[2];
sx q[2];
rz(1.8471225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0635447) q[1];
sx q[1];
rz(-1.1495665) q[1];
sx q[1];
rz(0.18555141) q[1];
rz(-pi) q[2];
rz(0.97350307) q[3];
sx q[3];
rz(-1.1949348) q[3];
sx q[3];
rz(1.6316538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3452611) q[2];
sx q[2];
rz(-0.76097208) q[2];
sx q[2];
rz(2.995028) q[2];
rz(-2.2740299) q[3];
sx q[3];
rz(-2.2680794) q[3];
sx q[3];
rz(-2.383702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72916156) q[0];
sx q[0];
rz(-2.6456092) q[0];
sx q[0];
rz(-2.5878986) q[0];
rz(1.4750534) q[1];
sx q[1];
rz(-0.29860425) q[1];
sx q[1];
rz(0.24436229) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15904832) q[0];
sx q[0];
rz(-1.9031018) q[0];
sx q[0];
rz(-2.766702) q[0];
rz(-1.2343132) q[2];
sx q[2];
rz(-0.76533106) q[2];
sx q[2];
rz(2.825277) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38589726) q[1];
sx q[1];
rz(-1.2127969) q[1];
sx q[1];
rz(2.6935607) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5958173) q[3];
sx q[3];
rz(-0.85880781) q[3];
sx q[3];
rz(-2.156949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.060140572) q[2];
sx q[2];
rz(-1.1271366) q[2];
sx q[2];
rz(2.3606908) q[2];
rz(0.39685708) q[3];
sx q[3];
rz(-0.01440993) q[3];
sx q[3];
rz(-1.074033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9020554) q[0];
sx q[0];
rz(-0.1665512) q[0];
sx q[0];
rz(0.2567513) q[0];
rz(-0.17770879) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(-0.94388747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0583874) q[0];
sx q[0];
rz(-0.58738601) q[0];
sx q[0];
rz(-2.2642231) q[0];
rz(-pi) q[1];
rz(-0.98942049) q[2];
sx q[2];
rz(-2.7983449) q[2];
sx q[2];
rz(-2.1892669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9967794) q[1];
sx q[1];
rz(-0.5397035) q[1];
sx q[1];
rz(3.1330879) q[1];
x q[2];
rz(-0.89045276) q[3];
sx q[3];
rz(-0.6750921) q[3];
sx q[3];
rz(-2.4321041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93861598) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(-0.20982783) q[2];
rz(-3.0948011) q[3];
sx q[3];
rz(-1.6822466) q[3];
sx q[3];
rz(0.34745026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059134722) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(1.2030075) q[0];
rz(-1.5164392) q[1];
sx q[1];
rz(-2.7639183) q[1];
sx q[1];
rz(3.1086521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66844475) q[0];
sx q[0];
rz(-1.1068891) q[0];
sx q[0];
rz(-0.3070681) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6141784) q[2];
sx q[2];
rz(-2.0927807) q[2];
sx q[2];
rz(-2.5185713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9666933) q[1];
sx q[1];
rz(-1.0224578) q[1];
sx q[1];
rz(1.95005) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63689639) q[3];
sx q[3];
rz(-0.99352443) q[3];
sx q[3];
rz(-3.0182216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.55543) q[2];
sx q[2];
rz(-1.0032434) q[2];
sx q[2];
rz(0.45247751) q[2];
rz(-0.14719506) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(1.2991306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.074987) q[0];
sx q[0];
rz(-0.43656483) q[0];
sx q[0];
rz(-2.7845352) q[0];
rz(-2.1287411) q[1];
sx q[1];
rz(-1.9722936) q[1];
sx q[1];
rz(-3.1049407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0520518) q[0];
sx q[0];
rz(-1.0424422) q[0];
sx q[0];
rz(2.0762877) q[0];
x q[1];
rz(0.42958625) q[2];
sx q[2];
rz(-1.4768147) q[2];
sx q[2];
rz(-2.2492301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7348902) q[1];
sx q[1];
rz(-1.478456) q[1];
sx q[1];
rz(-2.2734696) q[1];
rz(1.5662848) q[3];
sx q[3];
rz(-1.6016212) q[3];
sx q[3];
rz(-2.2982321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6645633) q[2];
sx q[2];
rz(-0.9534812) q[2];
sx q[2];
rz(3.0518517) q[2];
rz(1.2612032) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(-1.7173654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4795714) q[0];
sx q[0];
rz(-2.5630072) q[0];
sx q[0];
rz(1.660996) q[0];
rz(1.6914233) q[1];
sx q[1];
rz(-2.4822576) q[1];
sx q[1];
rz(-2.3816542) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284894) q[0];
sx q[0];
rz(-0.11664243) q[0];
sx q[0];
rz(1.4820497) q[0];
rz(2.5133695) q[2];
sx q[2];
rz(-1.6464273) q[2];
sx q[2];
rz(0.13720195) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48843417) q[1];
sx q[1];
rz(-1.1716845) q[1];
sx q[1];
rz(1.1121955) q[1];
x q[2];
rz(-1.8283056) q[3];
sx q[3];
rz(-1.0777359) q[3];
sx q[3];
rz(-0.3230394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8370168) q[2];
sx q[2];
rz(-1.2695856) q[2];
sx q[2];
rz(0.059583511) q[2];
rz(-2.7130821) q[3];
sx q[3];
rz(-0.84463745) q[3];
sx q[3];
rz(-0.60215157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2483599) q[0];
sx q[0];
rz(-0.28286523) q[0];
sx q[0];
rz(-0.34348139) q[0];
rz(-2.9454625) q[1];
sx q[1];
rz(-0.88534147) q[1];
sx q[1];
rz(-2.3416065) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752026) q[0];
sx q[0];
rz(-2.3227442) q[0];
sx q[0];
rz(0.99605346) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38651379) q[2];
sx q[2];
rz(-2.5441351) q[2];
sx q[2];
rz(-0.51624417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1175122) q[1];
sx q[1];
rz(-1.615106) q[1];
sx q[1];
rz(-2.5818392) q[1];
rz(-pi) q[2];
rz(-0.069234476) q[3];
sx q[3];
rz(-2.2524009) q[3];
sx q[3];
rz(0.72661663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7949152) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(-0.22870341) q[2];
rz(1.7250569) q[3];
sx q[3];
rz(-1.718947) q[3];
sx q[3];
rz(-0.67030877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.59920853) q[0];
sx q[0];
rz(-2.6238361) q[0];
sx q[0];
rz(1.1826578) q[0];
rz(-0.54284894) q[1];
sx q[1];
rz(-3.019637) q[1];
sx q[1];
rz(-2.6861232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1802954) q[0];
sx q[0];
rz(-2.0019128) q[0];
sx q[0];
rz(2.6102553) q[0];
rz(2.021505) q[2];
sx q[2];
rz(-1.7589065) q[2];
sx q[2];
rz(0.82291079) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7389981) q[1];
sx q[1];
rz(-1.4345503) q[1];
sx q[1];
rz(0.46936492) q[1];
rz(1.1703759) q[3];
sx q[3];
rz(-1.1569303) q[3];
sx q[3];
rz(3.1088157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21696422) q[2];
sx q[2];
rz(-0.89322007) q[2];
sx q[2];
rz(-0.42579892) q[2];
rz(0.18481542) q[3];
sx q[3];
rz(-0.63822377) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4375147) q[0];
sx q[0];
rz(-1.5298433) q[0];
sx q[0];
rz(-0.035506305) q[0];
rz(-0.47954814) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(-1.3479955) q[2];
sx q[2];
rz(-0.98735129) q[2];
sx q[2];
rz(-2.9410887) q[2];
rz(2.1816523) q[3];
sx q[3];
rz(-1.6889986) q[3];
sx q[3];
rz(-2.1903174) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
