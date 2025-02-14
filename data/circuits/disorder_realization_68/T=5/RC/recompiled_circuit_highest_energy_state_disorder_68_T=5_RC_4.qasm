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
rz(-2.3013378) q[0];
sx q[0];
rz(-1.037896) q[0];
sx q[0];
rz(2.2965746) q[0];
rz(-2.4294699) q[1];
sx q[1];
rz(-1.0016088) q[1];
sx q[1];
rz(1.6460302) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43216773) q[0];
sx q[0];
rz(-0.74678316) q[0];
sx q[0];
rz(0.92632697) q[0];
x q[1];
rz(2.7833013) q[2];
sx q[2];
rz(-0.40696496) q[2];
sx q[2];
rz(-2.586722) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.597376) q[1];
sx q[1];
rz(-0.68507441) q[1];
sx q[1];
rz(-2.3822576) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5450685) q[3];
sx q[3];
rz(-0.46677315) q[3];
sx q[3];
rz(1.0214361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.738395) q[2];
sx q[2];
rz(-1.9591816) q[2];
sx q[2];
rz(-0.5564059) q[2];
rz(0.49499908) q[3];
sx q[3];
rz(-0.35111108) q[3];
sx q[3];
rz(2.2569412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86968386) q[0];
sx q[0];
rz(-3.0392635) q[0];
sx q[0];
rz(-2.5129357) q[0];
rz(-2.2643845) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(0.25310755) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.077674) q[0];
sx q[0];
rz(-1.5787933) q[0];
sx q[0];
rz(1.5863938) q[0];
rz(-0.44127522) q[2];
sx q[2];
rz(-1.8417449) q[2];
sx q[2];
rz(0.51368062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15007764) q[1];
sx q[1];
rz(-1.3717695) q[1];
sx q[1];
rz(1.4519889) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8519446) q[3];
sx q[3];
rz(-2.3901403) q[3];
sx q[3];
rz(-0.40159097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46193281) q[2];
sx q[2];
rz(-1.2075281) q[2];
sx q[2];
rz(-1.0665464) q[2];
rz(1.3072394) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(-2.2333142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2147373) q[0];
sx q[0];
rz(-0.94811386) q[0];
sx q[0];
rz(0.98534775) q[0];
rz(-0.39380479) q[1];
sx q[1];
rz(-2.1486798) q[1];
sx q[1];
rz(2.6148112) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7486289) q[0];
sx q[0];
rz(-2.606578) q[0];
sx q[0];
rz(1.4070542) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57931945) q[2];
sx q[2];
rz(-2.0181106) q[2];
sx q[2];
rz(-1.8471225) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.078048) q[1];
sx q[1];
rz(-1.1495665) q[1];
sx q[1];
rz(0.18555141) q[1];
x q[2];
rz(-0.95895264) q[3];
sx q[3];
rz(-2.4482895) q[3];
sx q[3];
rz(-0.43402754) q[3];
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
rz(-0.87351322) q[3];
sx q[3];
rz(-0.7578907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4124311) q[0];
sx q[0];
rz(-2.6456092) q[0];
sx q[0];
rz(2.5878986) q[0];
rz(1.6665392) q[1];
sx q[1];
rz(-2.8429884) q[1];
sx q[1];
rz(0.24436229) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2840773) q[0];
sx q[0];
rz(-1.2173442) q[0];
sx q[0];
rz(-1.9259324) q[0];
x q[1];
rz(-0.8342488) q[2];
sx q[2];
rz(-1.3400199) q[2];
sx q[2];
rz(-2.1342056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0180514) q[1];
sx q[1];
rz(-1.1530211) q[1];
sx q[1];
rz(1.17735) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71214337) q[3];
sx q[3];
rz(-1.5518547) q[3];
sx q[3];
rz(-0.60250184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.060140572) q[2];
sx q[2];
rz(-2.0144561) q[2];
sx q[2];
rz(-2.3606908) q[2];
rz(-2.7447356) q[3];
sx q[3];
rz(-3.1271827) q[3];
sx q[3];
rz(1.074033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.9020554) q[0];
sx q[0];
rz(-0.1665512) q[0];
sx q[0];
rz(2.8848414) q[0];
rz(-0.17770879) q[1];
sx q[1];
rz(-2.5959028) q[1];
sx q[1];
rz(2.1977052) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8430169) q[0];
sx q[0];
rz(-1.1305048) q[0];
sx q[0];
rz(2.7392575) q[0];
x q[1];
rz(-0.19380865) q[2];
sx q[2];
rz(-1.8558981) q[2];
sx q[2];
rz(1.5616035) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.9967794) q[1];
sx q[1];
rz(-2.6018892) q[1];
sx q[1];
rz(3.1330879) q[1];
x q[2];
rz(1.0141171) q[3];
sx q[3];
rz(-1.9748439) q[3];
sx q[3];
rz(1.7168644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2029767) q[2];
sx q[2];
rz(-1.0304893) q[2];
sx q[2];
rz(-2.9317648) q[2];
rz(-3.0948011) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(-0.34745026) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059134722) q[0];
sx q[0];
rz(-0.81675285) q[0];
sx q[0];
rz(-1.2030075) q[0];
rz(-1.5164392) q[1];
sx q[1];
rz(-2.7639183) q[1];
sx q[1];
rz(3.1086521) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0982978) q[0];
sx q[0];
rz(-1.2970719) q[0];
sx q[0];
rz(1.0874463) q[0];
x q[1];
rz(0.85196544) q[2];
sx q[2];
rz(-2.4174066) q[2];
sx q[2];
rz(-2.9018096) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1748993) q[1];
sx q[1];
rz(-2.1191349) q[1];
sx q[1];
rz(1.95005) q[1];
x q[2];
rz(0.88992702) q[3];
sx q[3];
rz(-1.0491199) q[3];
sx q[3];
rz(-1.8311178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5861627) q[2];
sx q[2];
rz(-1.0032434) q[2];
sx q[2];
rz(2.6891151) q[2];
rz(2.9943976) q[3];
sx q[3];
rz(-0.23284027) q[3];
sx q[3];
rz(-1.8424621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066605695) q[0];
sx q[0];
rz(-0.43656483) q[0];
sx q[0];
rz(2.7845352) q[0];
rz(1.0128516) q[1];
sx q[1];
rz(-1.169299) q[1];
sx q[1];
rz(-0.036651932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9215611) q[0];
sx q[0];
rz(-2.4274094) q[0];
sx q[0];
rz(0.69252391) q[0];
rz(-pi) q[1];
rz(-0.22256644) q[2];
sx q[2];
rz(-0.43912008) q[2];
sx q[2];
rz(0.47638461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7348902) q[1];
sx q[1];
rz(-1.478456) q[1];
sx q[1];
rz(-0.86812302) q[1];
x q[2];
rz(0.030825214) q[3];
sx q[3];
rz(-1.5662869) q[3];
sx q[3];
rz(-2.4140178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6645633) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(3.0518517) q[2];
rz(-1.2612032) q[3];
sx q[3];
rz(-1.3124876) q[3];
sx q[3];
rz(-1.7173654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66202128) q[0];
sx q[0];
rz(-0.57858545) q[0];
sx q[0];
rz(1.660996) q[0];
rz(1.6914233) q[1];
sx q[1];
rz(-2.4822576) q[1];
sx q[1];
rz(0.75993842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.511139) q[0];
sx q[0];
rz(-1.5604815) q[0];
sx q[0];
rz(1.4546088) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6282232) q[2];
sx q[2];
rz(-1.4951653) q[2];
sx q[2];
rz(0.13720195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.41553283) q[1];
sx q[1];
rz(-2.5431255) q[1];
sx q[1];
rz(0.80962555) q[1];
rz(0.50716831) q[3];
sx q[3];
rz(-1.3445323) q[3];
sx q[3];
rz(1.3717701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8370168) q[2];
sx q[2];
rz(-1.872007) q[2];
sx q[2];
rz(-3.0820091) q[2];
rz(0.42851055) q[3];
sx q[3];
rz(-2.2969552) q[3];
sx q[3];
rz(-2.5394411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8932327) q[0];
sx q[0];
rz(-2.8587274) q[0];
sx q[0];
rz(0.34348139) q[0];
rz(-0.19613014) q[1];
sx q[1];
rz(-2.2562512) q[1];
sx q[1];
rz(0.79998618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8752026) q[0];
sx q[0];
rz(-0.81884844) q[0];
sx q[0];
rz(0.99605346) q[0];
rz(-2.5792511) q[2];
sx q[2];
rz(-1.7844756) q[2];
sx q[2];
rz(2.4116229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0240805) q[1];
sx q[1];
rz(-1.615106) q[1];
sx q[1];
rz(-2.5818392) q[1];
rz(-pi) q[2];
rz(-2.2535747) q[3];
sx q[3];
rz(-1.6245442) q[3];
sx q[3];
rz(-2.2537504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7949152) q[2];
sx q[2];
rz(-0.76378834) q[2];
sx q[2];
rz(-0.22870341) q[2];
rz(-1.7250569) q[3];
sx q[3];
rz(-1.4226457) q[3];
sx q[3];
rz(2.4712839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5423841) q[0];
sx q[0];
rz(-2.6238361) q[0];
sx q[0];
rz(-1.9589348) q[0];
rz(-2.5987437) q[1];
sx q[1];
rz(-3.019637) q[1];
sx q[1];
rz(2.6861232) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36866828) q[0];
sx q[0];
rz(-2.0491573) q[0];
sx q[0];
rz(-2.0609074) q[0];
x q[1];
rz(0.20840877) q[2];
sx q[2];
rz(-1.128607) q[2];
sx q[2];
rz(2.4839574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7389981) q[1];
sx q[1];
rz(-1.4345503) q[1];
sx q[1];
rz(-2.6722277) q[1];
rz(-1.1703759) q[3];
sx q[3];
rz(-1.1569303) q[3];
sx q[3];
rz(0.032776959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21696422) q[2];
sx q[2];
rz(-0.89322007) q[2];
sx q[2];
rz(0.42579892) q[2];
rz(2.9567772) q[3];
sx q[3];
rz(-2.5033689) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7040779) q[0];
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
rz(-2.9976224) q[3];
sx q[3];
rz(-0.96481779) q[3];
sx q[3];
rz(2.6044566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
