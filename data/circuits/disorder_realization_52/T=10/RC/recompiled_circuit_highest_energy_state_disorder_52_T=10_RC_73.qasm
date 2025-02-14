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
rz(3.8654397) q[0];
sx q[0];
rz(1.2011733) q[0];
sx q[0];
rz(6.7772128) q[0];
rz(-1.5863034) q[1];
sx q[1];
rz(-2.2145693) q[1];
sx q[1];
rz(2.291099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9181831) q[0];
sx q[0];
rz(-0.69108686) q[0];
sx q[0];
rz(1.8331493) q[0];
rz(-pi) q[1];
rz(2.565704) q[2];
sx q[2];
rz(-2.401899) q[2];
sx q[2];
rz(1.0583371) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32341084) q[1];
sx q[1];
rz(-1.8935303) q[1];
sx q[1];
rz(-1.3363911) q[1];
rz(-1.4981573) q[3];
sx q[3];
rz(-1.2056672) q[3];
sx q[3];
rz(2.3453494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4832619) q[2];
sx q[2];
rz(-2.5265145) q[2];
sx q[2];
rz(0.94493803) q[2];
rz(-3.1350709) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(0.25750461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1061123) q[0];
sx q[0];
rz(-0.98863125) q[0];
sx q[0];
rz(-0.60580564) q[0];
rz(2.2094191) q[1];
sx q[1];
rz(-1.7180387) q[1];
sx q[1];
rz(3.0359643) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6208134) q[0];
sx q[0];
rz(-0.70292369) q[0];
sx q[0];
rz(-2.9018794) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95870734) q[2];
sx q[2];
rz(-2.5664133) q[2];
sx q[2];
rz(1.5633595) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2534598) q[1];
sx q[1];
rz(-1.2769298) q[1];
sx q[1];
rz(2.6685358) q[1];
rz(0.83189957) q[3];
sx q[3];
rz(-2.6100592) q[3];
sx q[3];
rz(-3.1011776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85670829) q[2];
sx q[2];
rz(-1.363089) q[2];
sx q[2];
rz(-1.523783) q[2];
rz(-0.81426042) q[3];
sx q[3];
rz(-1.3596478) q[3];
sx q[3];
rz(2.2566569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098323671) q[0];
sx q[0];
rz(-2.3461778) q[0];
sx q[0];
rz(1.5765618) q[0];
rz(2.1463429) q[1];
sx q[1];
rz(-0.97882706) q[1];
sx q[1];
rz(-2.3547122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9691447) q[0];
sx q[0];
rz(-0.24227628) q[0];
sx q[0];
rz(2.1680225) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7943775) q[2];
sx q[2];
rz(-1.6499398) q[2];
sx q[2];
rz(-1.7007507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9287368) q[1];
sx q[1];
rz(-2.6784705) q[1];
sx q[1];
rz(-2.379851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3555879) q[3];
sx q[3];
rz(-0.084298221) q[3];
sx q[3];
rz(0.5072197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83871049) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(0.6380471) q[2];
rz(2.6436515) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(2.0518484) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8968673) q[0];
sx q[0];
rz(-1.7843972) q[0];
sx q[0];
rz(0.063752739) q[0];
rz(-1.6925192) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(1.8316899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44353911) q[0];
sx q[0];
rz(-1.6081282) q[0];
sx q[0];
rz(-3.061767) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3856349) q[2];
sx q[2];
rz(-1.3755535) q[2];
sx q[2];
rz(2.756292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79015649) q[1];
sx q[1];
rz(-2.3871941) q[1];
sx q[1];
rz(2.960207) q[1];
rz(-pi) q[2];
rz(2.2579262) q[3];
sx q[3];
rz(-1.851463) q[3];
sx q[3];
rz(-1.3786045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.070907585) q[2];
sx q[2];
rz(-1.1068152) q[2];
sx q[2];
rz(3.0359388) q[2];
rz(-1.1834772) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(-0.67160523) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79528177) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(1.7972535) q[0];
rz(2.4686939) q[1];
sx q[1];
rz(-1.3105323) q[1];
sx q[1];
rz(0.013462822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97082645) q[0];
sx q[0];
rz(-1.544569) q[0];
sx q[0];
rz(0.88084014) q[0];
rz(-pi) q[1];
rz(1.9849987) q[2];
sx q[2];
rz(-1.557781) q[2];
sx q[2];
rz(-2.3682197) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9078482) q[1];
sx q[1];
rz(-1.9861167) q[1];
sx q[1];
rz(1.5563346) q[1];
rz(-2.0125403) q[3];
sx q[3];
rz(-0.75296445) q[3];
sx q[3];
rz(0.66679614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54231918) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(1.8947961) q[2];
rz(-3.0689734) q[3];
sx q[3];
rz(-2.9473372) q[3];
sx q[3];
rz(-0.3092002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4319864) q[0];
sx q[0];
rz(-2.015634) q[0];
sx q[0];
rz(-0.11269888) q[0];
rz(-0.20507774) q[1];
sx q[1];
rz(-1.8094742) q[1];
sx q[1];
rz(-2.233706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31602177) q[0];
sx q[0];
rz(-2.1793723) q[0];
sx q[0];
rz(-0.84644239) q[0];
rz(2.5656469) q[2];
sx q[2];
rz(-2.4436946) q[2];
sx q[2];
rz(0.18651785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7026313) q[1];
sx q[1];
rz(-1.6621894) q[1];
sx q[1];
rz(-0.63321873) q[1];
rz(2.870918) q[3];
sx q[3];
rz(-2.4587016) q[3];
sx q[3];
rz(-2.0840933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1418566) q[2];
sx q[2];
rz(-2.8047968) q[2];
sx q[2];
rz(-0.067961819) q[2];
rz(-1.147602) q[3];
sx q[3];
rz(-1.8736898) q[3];
sx q[3];
rz(1.2609153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481908) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(2.6780658) q[0];
rz(2.9947128) q[1];
sx q[1];
rz(-1.2356707) q[1];
sx q[1];
rz(0.99259496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79560876) q[0];
sx q[0];
rz(-2.3844686) q[0];
sx q[0];
rz(-0.36619314) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3456773) q[2];
sx q[2];
rz(-0.29307191) q[2];
sx q[2];
rz(-0.64426433) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9001683) q[1];
sx q[1];
rz(-1.1196152) q[1];
sx q[1];
rz(-0.74189775) q[1];
x q[2];
rz(-0.73294183) q[3];
sx q[3];
rz(-1.9187201) q[3];
sx q[3];
rz(2.3331353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85840449) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(-0.44537133) q[2];
rz(1.8177659) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(-2.0557192) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0367947) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(-0.15583663) q[0];
rz(1.8644631) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(3.1256622) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30980047) q[0];
sx q[0];
rz(-0.81912097) q[0];
sx q[0];
rz(2.8099634) q[0];
x q[1];
rz(-2.3492947) q[2];
sx q[2];
rz(-1.8323261) q[2];
sx q[2];
rz(-0.025394414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5585352) q[1];
sx q[1];
rz(-1.0638405) q[1];
sx q[1];
rz(-0.81113775) q[1];
rz(-pi) q[2];
x q[2];
rz(1.995581) q[3];
sx q[3];
rz(-0.97633712) q[3];
sx q[3];
rz(0.068258523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5126123) q[2];
sx q[2];
rz(-0.11233687) q[2];
sx q[2];
rz(-0.12574276) q[2];
rz(-0.9564774) q[3];
sx q[3];
rz(-1.7185017) q[3];
sx q[3];
rz(-0.3392578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823031) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-2.6851728) q[0];
rz(-1.5688815) q[1];
sx q[1];
rz(-2.1379037) q[1];
sx q[1];
rz(-2.454954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4149041) q[0];
sx q[0];
rz(-2.3285638) q[0];
sx q[0];
rz(0.70720478) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3977658) q[2];
sx q[2];
rz(-1.695172) q[2];
sx q[2];
rz(-1.591452) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9514044) q[1];
sx q[1];
rz(-1.4238289) q[1];
sx q[1];
rz(1.6920056) q[1];
x q[2];
rz(-2.7620662) q[3];
sx q[3];
rz(-0.18774334) q[3];
sx q[3];
rz(1.8462586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75866428) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(-2.0056966) q[2];
rz(-0.9797594) q[3];
sx q[3];
rz(-1.8085248) q[3];
sx q[3];
rz(2.4433344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6887688) q[0];
sx q[0];
rz(-1.8166421) q[0];
sx q[0];
rz(2.685637) q[0];
rz(0.042757209) q[1];
sx q[1];
rz(-1.2084992) q[1];
sx q[1];
rz(-2.4483689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15007818) q[0];
sx q[0];
rz(-0.39798361) q[0];
sx q[0];
rz(-1.8282169) q[0];
x q[1];
rz(-1.1666388) q[2];
sx q[2];
rz(-1.1956788) q[2];
sx q[2];
rz(-0.98304316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0574368) q[1];
sx q[1];
rz(-0.098271253) q[1];
sx q[1];
rz(1.6429423) q[1];
rz(-pi) q[2];
rz(-1.2537122) q[3];
sx q[3];
rz(-1.1389995) q[3];
sx q[3];
rz(-0.82796873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91528714) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(2.477296) q[2];
rz(-1.213446) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(2.2693995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.450347) q[0];
sx q[0];
rz(-0.84083122) q[0];
sx q[0];
rz(-1.3837411) q[0];
rz(1.6830403) q[1];
sx q[1];
rz(-2.8589307) q[1];
sx q[1];
rz(0.9160441) q[1];
rz(-3.0677879) q[2];
sx q[2];
rz(-0.37141411) q[2];
sx q[2];
rz(0.85167533) q[2];
rz(-1.0091598) q[3];
sx q[3];
rz(-0.94579332) q[3];
sx q[3];
rz(1.6413064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
