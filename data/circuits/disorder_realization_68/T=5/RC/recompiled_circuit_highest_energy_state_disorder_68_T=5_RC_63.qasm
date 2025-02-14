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
rz(-2.1036966) q[0];
sx q[0];
rz(0.84501803) q[0];
rz(0.71212274) q[1];
sx q[1];
rz(4.1432015) q[1];
sx q[1];
rz(10.92034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43216773) q[0];
sx q[0];
rz(-2.3948095) q[0];
sx q[0];
rz(-0.92632697) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7579384) q[2];
sx q[2];
rz(-1.4315413) q[2];
sx q[2];
rz(-0.68472199) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.597376) q[1];
sx q[1];
rz(-2.4565182) q[1];
sx q[1];
rz(-2.3822576) q[1];
x q[2];
rz(-2.7466165) q[3];
sx q[3];
rz(-1.8263701) q[3];
sx q[3];
rz(3.1373484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40319765) q[2];
sx q[2];
rz(-1.9591816) q[2];
sx q[2];
rz(-0.5564059) q[2];
rz(-2.6465936) q[3];
sx q[3];
rz(-2.7904816) q[3];
sx q[3];
rz(-2.2569412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86968386) q[0];
sx q[0];
rz(-3.0392635) q[0];
sx q[0];
rz(2.5129357) q[0];
rz(2.2643845) q[1];
sx q[1];
rz(-0.49867201) q[1];
sx q[1];
rz(-2.8884851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49324711) q[0];
sx q[0];
rz(-1.5863933) q[0];
sx q[0];
rz(0.007997917) q[0];
rz(-pi) q[1];
rz(1.8688525) q[2];
sx q[2];
rz(-1.1466845) q[2];
sx q[2];
rz(-1.9587245) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15007764) q[1];
sx q[1];
rz(-1.7698231) q[1];
sx q[1];
rz(1.4519889) q[1];
x q[2];
rz(2.8879426) q[3];
sx q[3];
rz(-0.85542711) q[3];
sx q[3];
rz(-0.77798346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6796598) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(1.0665464) q[2];
rz(1.3072394) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(-2.2333142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2147373) q[0];
sx q[0];
rz(-2.1934788) q[0];
sx q[0];
rz(2.1562449) q[0];
rz(0.39380479) q[1];
sx q[1];
rz(-2.1486798) q[1];
sx q[1];
rz(0.52678144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1049479) q[0];
sx q[0];
rz(-1.6540043) q[0];
sx q[0];
rz(1.0416688) q[0];
rz(-2.4220139) q[2];
sx q[2];
rz(-2.4256896) q[2];
sx q[2];
rz(2.2810675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5722428) q[1];
sx q[1];
rz(-1.7399638) q[1];
sx q[1];
rz(1.9985241) q[1];
rz(2.1680896) q[3];
sx q[3];
rz(-1.1949348) q[3];
sx q[3];
rz(-1.6316538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7963316) q[2];
sx q[2];
rz(-2.3806206) q[2];
sx q[2];
rz(0.14656466) q[2];
rz(-0.86756271) q[3];
sx q[3];
rz(-0.87351322) q[3];
sx q[3];
rz(0.7578907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4124311) q[0];
sx q[0];
rz(-0.49598345) q[0];
sx q[0];
rz(0.55369401) q[0];
rz(1.4750534) q[1];
sx q[1];
rz(-2.8429884) q[1];
sx q[1];
rz(2.8972304) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.037652) q[0];
sx q[0];
rz(-2.6459011) q[0];
sx q[0];
rz(2.3858059) q[0];
rz(-0.8342488) q[2];
sx q[2];
rz(-1.8015727) q[2];
sx q[2];
rz(-1.007387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1235413) q[1];
sx q[1];
rz(-1.9885716) q[1];
sx q[1];
rz(-1.9642427) q[1];
x q[2];
rz(1.5958173) q[3];
sx q[3];
rz(-0.85880781) q[3];
sx q[3];
rz(0.98464363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0814521) q[2];
sx q[2];
rz(-1.1271366) q[2];
sx q[2];
rz(-2.3606908) q[2];
rz(0.39685708) q[3];
sx q[3];
rz(-3.1271827) q[3];
sx q[3];
rz(-2.0675596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395372) q[0];
sx q[0];
rz(-2.9750415) q[0];
sx q[0];
rz(2.8848414) q[0];
rz(-0.17770879) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(-0.94388747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2985757) q[0];
sx q[0];
rz(-1.1305048) q[0];
sx q[0];
rz(2.7392575) q[0];
rz(2.1521722) q[2];
sx q[2];
rz(-2.7983449) q[2];
sx q[2];
rz(-2.1892669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.56672102) q[1];
sx q[1];
rz(-1.5751667) q[1];
sx q[1];
rz(-2.6019051) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2511399) q[3];
sx q[3];
rz(-2.4665006) q[3];
sx q[3];
rz(-0.70948851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2029767) q[2];
sx q[2];
rz(-1.0304893) q[2];
sx q[2];
rz(2.9317648) q[2];
rz(3.0948011) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(0.34745026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.059134722) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(1.9385852) q[0];
rz(1.6251534) q[1];
sx q[1];
rz(-0.37767437) q[1];
sx q[1];
rz(0.032940544) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0982978) q[0];
sx q[0];
rz(-1.8445208) q[0];
sx q[0];
rz(-2.0541463) q[0];
x q[1];
rz(-2.2896272) q[2];
sx q[2];
rz(-2.4174066) q[2];
sx q[2];
rz(0.23978309) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1748993) q[1];
sx q[1];
rz(-2.1191349) q[1];
sx q[1];
rz(1.1915426) q[1];
x q[2];
rz(0.88992702) q[3];
sx q[3];
rz(-2.0924727) q[3];
sx q[3];
rz(-1.3104749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5861627) q[2];
sx q[2];
rz(-2.1383492) q[2];
sx q[2];
rz(2.6891151) q[2];
rz(-2.9943976) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(-1.2991306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.074987) q[0];
sx q[0];
rz(-2.7050278) q[0];
sx q[0];
rz(-2.7845352) q[0];
rz(-2.1287411) q[1];
sx q[1];
rz(-1.9722936) q[1];
sx q[1];
rz(-3.1049407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79082876) q[0];
sx q[0];
rz(-2.0022656) q[0];
sx q[0];
rz(-2.55326) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9190262) q[2];
sx q[2];
rz(-2.7024726) q[2];
sx q[2];
rz(2.665208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27254912) q[1];
sx q[1];
rz(-2.4339035) q[1];
sx q[1];
rz(1.7131192) q[1];
x q[2];
rz(-1.5753079) q[3];
sx q[3];
rz(-1.5399714) q[3];
sx q[3];
rz(-0.84336057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4770294) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(-3.0518517) q[2];
rz(-1.2612032) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(1.7173654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66202128) q[0];
sx q[0];
rz(-0.57858545) q[0];
sx q[0];
rz(-1.660996) q[0];
rz(-1.6914233) q[1];
sx q[1];
rz(-0.65933508) q[1];
sx q[1];
rz(0.75993842) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2024538) q[0];
sx q[0];
rz(-1.6869776) q[0];
sx q[0];
rz(-3.1312079) q[0];
x q[1];
rz(-1.6641812) q[2];
sx q[2];
rz(-2.1969446) q[2];
sx q[2];
rz(1.4884251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7260598) q[1];
sx q[1];
rz(-0.59846717) q[1];
sx q[1];
rz(-2.3319671) q[1];
rz(-2.6344243) q[3];
sx q[3];
rz(-1.7970603) q[3];
sx q[3];
rz(-1.3717701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8370168) q[2];
sx q[2];
rz(-1.2695856) q[2];
sx q[2];
rz(-0.059583511) q[2];
rz(-0.42851055) q[3];
sx q[3];
rz(-2.2969552) q[3];
sx q[3];
rz(2.5394411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.2483599) q[0];
sx q[0];
rz(-0.28286523) q[0];
sx q[0];
rz(0.34348139) q[0];
rz(-2.9454625) q[1];
sx q[1];
rz(-2.2562512) q[1];
sx q[1];
rz(-0.79998618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752026) q[0];
sx q[0];
rz(-0.81884844) q[0];
sx q[0];
rz(-0.99605346) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38651379) q[2];
sx q[2];
rz(-0.59745759) q[2];
sx q[2];
rz(-2.6253485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51896732) q[1];
sx q[1];
rz(-2.1299348) q[1];
sx q[1];
rz(-1.6230727) q[1];
rz(-pi) q[2];
rz(-1.6558582) q[3];
sx q[3];
rz(-0.68455212) q[3];
sx q[3];
rz(-2.5245997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34667748) q[2];
sx q[2];
rz(-2.3778043) q[2];
sx q[2];
rz(-0.22870341) q[2];
rz(1.7250569) q[3];
sx q[3];
rz(-1.718947) q[3];
sx q[3];
rz(2.4712839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5423841) q[0];
sx q[0];
rz(-2.6238361) q[0];
sx q[0];
rz(1.9589348) q[0];
rz(-0.54284894) q[1];
sx q[1];
rz(-3.019637) q[1];
sx q[1];
rz(-2.6861232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802954) q[0];
sx q[0];
rz(-1.1396798) q[0];
sx q[0];
rz(-0.53133734) q[0];
rz(-1.9827911) q[2];
sx q[2];
rz(-0.48589009) q[2];
sx q[2];
rz(2.0249572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9046136) q[1];
sx q[1];
rz(-1.1061258) q[1];
sx q[1];
rz(-1.7233217) q[1];
x q[2];
rz(1.1703759) q[3];
sx q[3];
rz(-1.1569303) q[3];
sx q[3];
rz(3.1088157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9246284) q[2];
sx q[2];
rz(-0.89322007) q[2];
sx q[2];
rz(0.42579892) q[2];
rz(0.18481542) q[3];
sx q[3];
rz(-0.63822377) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7040779) q[0];
sx q[0];
rz(-1.6117493) q[0];
sx q[0];
rz(3.1060863) q[0];
rz(2.6620445) q[1];
sx q[1];
rz(-1.5621114) q[1];
sx q[1];
rz(1.5893804) q[1];
rz(1.3479955) q[2];
sx q[2];
rz(-2.1542414) q[2];
sx q[2];
rz(0.20050394) q[2];
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
