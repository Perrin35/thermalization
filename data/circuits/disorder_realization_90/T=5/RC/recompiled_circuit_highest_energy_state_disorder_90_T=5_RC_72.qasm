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
rz(2.4861205) q[0];
sx q[0];
rz(-0.95215005) q[0];
sx q[0];
rz(-0.093753554) q[0];
rz(1.3682415) q[1];
sx q[1];
rz(-1.5803087) q[1];
sx q[1];
rz(-0.66722792) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.102948) q[0];
sx q[0];
rz(-1.1955452) q[0];
sx q[0];
rz(2.1041591) q[0];
rz(-0.66114963) q[2];
sx q[2];
rz(-2.8183658) q[2];
sx q[2];
rz(-1.8913392) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2161795) q[1];
sx q[1];
rz(-0.88757463) q[1];
sx q[1];
rz(-2.776078) q[1];
rz(2.2556858) q[3];
sx q[3];
rz(-1.0647756) q[3];
sx q[3];
rz(-0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61624709) q[2];
sx q[2];
rz(-1.5364001) q[2];
sx q[2];
rz(3.1180535) q[2];
rz(0.25447887) q[3];
sx q[3];
rz(-1.9813709) q[3];
sx q[3];
rz(-2.3833073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1940521) q[0];
sx q[0];
rz(-1.3964615) q[0];
sx q[0];
rz(-2.5260455) q[0];
rz(2.5415892) q[1];
sx q[1];
rz(-1.871385) q[1];
sx q[1];
rz(0.34872762) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.00086) q[0];
sx q[0];
rz(-1.6475186) q[0];
sx q[0];
rz(0.64312913) q[0];
rz(-pi) q[1];
rz(-2.8519657) q[2];
sx q[2];
rz(-1.2420601) q[2];
sx q[2];
rz(-1.5266974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24838994) q[1];
sx q[1];
rz(-0.55596184) q[1];
sx q[1];
rz(-0.979579) q[1];
x q[2];
rz(-2.0754966) q[3];
sx q[3];
rz(-1.3196006) q[3];
sx q[3];
rz(-0.78758729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6500924) q[2];
sx q[2];
rz(-1.3050175) q[2];
sx q[2];
rz(0.65518641) q[2];
rz(-2.9789467) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(-2.9338525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96838897) q[0];
sx q[0];
rz(-1.116331) q[0];
sx q[0];
rz(-3.047347) q[0];
rz(1.3000129) q[1];
sx q[1];
rz(-0.899122) q[1];
sx q[1];
rz(-1.947044) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0604977) q[0];
sx q[0];
rz(-1.0733685) q[0];
sx q[0];
rz(-1.917385) q[0];
x q[1];
rz(-3.0897806) q[2];
sx q[2];
rz(-1.9345043) q[2];
sx q[2];
rz(-1.6339982) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49385168) q[1];
sx q[1];
rz(-1.6147227) q[1];
sx q[1];
rz(-1.3141339) q[1];
x q[2];
rz(-2.3787375) q[3];
sx q[3];
rz(-1.8486946) q[3];
sx q[3];
rz(-2.7369473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3208348) q[2];
sx q[2];
rz(-1.6944378) q[2];
sx q[2];
rz(-2.7336332) q[2];
rz(-1.3784493) q[3];
sx q[3];
rz(-2.2429376) q[3];
sx q[3];
rz(0.11817008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57678643) q[0];
sx q[0];
rz(-1.3917568) q[0];
sx q[0];
rz(2.3396662) q[0];
rz(2.7469514) q[1];
sx q[1];
rz(-2.092974) q[1];
sx q[1];
rz(-3.1065497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3010522) q[0];
sx q[0];
rz(-0.6834712) q[0];
sx q[0];
rz(2.1816064) q[0];
rz(-0.50582992) q[2];
sx q[2];
rz(-2.7625045) q[2];
sx q[2];
rz(-2.068678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38742018) q[1];
sx q[1];
rz(-1.6129478) q[1];
sx q[1];
rz(1.6142963) q[1];
rz(-2.6529516) q[3];
sx q[3];
rz(-0.39791574) q[3];
sx q[3];
rz(0.15109381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0617712) q[2];
sx q[2];
rz(-1.062919) q[2];
sx q[2];
rz(1.4351832) q[2];
rz(-1.4052514) q[3];
sx q[3];
rz(-2.0515714) q[3];
sx q[3];
rz(1.1100356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6768796) q[0];
sx q[0];
rz(-1.1486624) q[0];
sx q[0];
rz(1.9997464) q[0];
rz(-0.51072085) q[1];
sx q[1];
rz(-1.9074214) q[1];
sx q[1];
rz(1.4265192) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0240117) q[0];
sx q[0];
rz(-0.80285751) q[0];
sx q[0];
rz(-0.76501655) q[0];
x q[1];
rz(-2.9147706) q[2];
sx q[2];
rz(-1.87543) q[2];
sx q[2];
rz(-2.2351928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.63663) q[1];
sx q[1];
rz(-2.1607724) q[1];
sx q[1];
rz(3.0061199) q[1];
rz(-pi) q[2];
rz(-2.9742091) q[3];
sx q[3];
rz(-1.1684844) q[3];
sx q[3];
rz(1.917812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4003754) q[2];
sx q[2];
rz(-1.8343238) q[2];
sx q[2];
rz(-0.2571787) q[2];
rz(0.87279618) q[3];
sx q[3];
rz(-1.8200487) q[3];
sx q[3];
rz(-2.0040373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1207101) q[0];
sx q[0];
rz(-2.6873984) q[0];
sx q[0];
rz(-1.0783476) q[0];
rz(-2.2348166) q[1];
sx q[1];
rz(-0.6858784) q[1];
sx q[1];
rz(-3.0996941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9020136) q[0];
sx q[0];
rz(-1.5348832) q[0];
sx q[0];
rz(1.0980868) q[0];
rz(1.4153775) q[2];
sx q[2];
rz(-0.67314429) q[2];
sx q[2];
rz(-2.239925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4087569) q[1];
sx q[1];
rz(-0.46161595) q[1];
sx q[1];
rz(2.2240586) q[1];
rz(-pi) q[2];
rz(0.83591977) q[3];
sx q[3];
rz(-1.984388) q[3];
sx q[3];
rz(-1.6294162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28356734) q[2];
sx q[2];
rz(-0.70245409) q[2];
sx q[2];
rz(0.89135998) q[2];
rz(0.71409613) q[3];
sx q[3];
rz(-2.4024506) q[3];
sx q[3];
rz(2.0785296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5080268) q[0];
sx q[0];
rz(-1.897568) q[0];
sx q[0];
rz(2.4666393) q[0];
rz(-1.630111) q[1];
sx q[1];
rz(-0.62050301) q[1];
sx q[1];
rz(-2.564548) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16923184) q[0];
sx q[0];
rz(-1.0857538) q[0];
sx q[0];
rz(-2.8718456) q[0];
x q[1];
rz(-0.21930947) q[2];
sx q[2];
rz(-0.44249924) q[2];
sx q[2];
rz(-2.2377917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8220209) q[1];
sx q[1];
rz(-2.4720925) q[1];
sx q[1];
rz(1.7029525) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2636599) q[3];
sx q[3];
rz(-2.0725277) q[3];
sx q[3];
rz(2.2449765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.453489) q[2];
sx q[2];
rz(-2.0192396) q[2];
sx q[2];
rz(0.92448676) q[2];
rz(-1.9214572) q[3];
sx q[3];
rz(-2.852738) q[3];
sx q[3];
rz(0.060062241) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7500551) q[0];
sx q[0];
rz(-2.4704762) q[0];
sx q[0];
rz(1.7975988) q[0];
rz(1.2229819) q[1];
sx q[1];
rz(-1.5593301) q[1];
sx q[1];
rz(-1.5498243) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61448288) q[0];
sx q[0];
rz(-1.7507995) q[0];
sx q[0];
rz(1.8480145) q[0];
rz(-pi) q[1];
rz(2.7666367) q[2];
sx q[2];
rz(-1.4658615) q[2];
sx q[2];
rz(-0.31601672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0613411) q[1];
sx q[1];
rz(-2.3480573) q[1];
sx q[1];
rz(-0.52498753) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0430085) q[3];
sx q[3];
rz(-2.3752593) q[3];
sx q[3];
rz(-2.6232972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.082077114) q[2];
sx q[2];
rz(-2.1166708) q[2];
sx q[2];
rz(0.22641851) q[2];
rz(-0.039479937) q[3];
sx q[3];
rz(-1.4791146) q[3];
sx q[3];
rz(-1.1994908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3969642) q[0];
sx q[0];
rz(-1.9292984) q[0];
sx q[0];
rz(0.82897559) q[0];
rz(0.12116155) q[1];
sx q[1];
rz(-0.96210259) q[1];
sx q[1];
rz(-1.4519579) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054497) q[0];
sx q[0];
rz(-2.2847957) q[0];
sx q[0];
rz(1.3881486) q[0];
rz(2.0214861) q[2];
sx q[2];
rz(-0.7996489) q[2];
sx q[2];
rz(0.19935184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.02579298) q[1];
sx q[1];
rz(-2.3069972) q[1];
sx q[1];
rz(2.2195312) q[1];
rz(-pi) q[2];
rz(2.1690458) q[3];
sx q[3];
rz(-2.498811) q[3];
sx q[3];
rz(-1.169426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58763233) q[2];
sx q[2];
rz(-2.3162737) q[2];
sx q[2];
rz(0.86453214) q[2];
rz(2.4899321) q[3];
sx q[3];
rz(-1.0385907) q[3];
sx q[3];
rz(0.017875044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5038576) q[0];
sx q[0];
rz(-0.88541579) q[0];
sx q[0];
rz(-2.6813685) q[0];
rz(-1.3453311) q[1];
sx q[1];
rz(-0.63667744) q[1];
sx q[1];
rz(1.3705137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2945098) q[0];
sx q[0];
rz(-0.50343236) q[0];
sx q[0];
rz(-2.2398021) q[0];
rz(-pi) q[1];
rz(-2.7659594) q[2];
sx q[2];
rz(-2.4690095) q[2];
sx q[2];
rz(-3.1194558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9784097) q[1];
sx q[1];
rz(-1.8272361) q[1];
sx q[1];
rz(-0.24602166) q[1];
rz(-pi) q[2];
rz(-0.053456177) q[3];
sx q[3];
rz(-1.7776907) q[3];
sx q[3];
rz(-2.2446936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7728077) q[2];
sx q[2];
rz(-1.7926755) q[2];
sx q[2];
rz(-2.9985912) q[2];
rz(0.16608206) q[3];
sx q[3];
rz(-2.4487285) q[3];
sx q[3];
rz(0.79296976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3731257) q[0];
sx q[0];
rz(-1.65092) q[0];
sx q[0];
rz(-1.7478818) q[0];
rz(-1.1557747) q[1];
sx q[1];
rz(-2.0230237) q[1];
sx q[1];
rz(0.3442234) q[1];
rz(-0.73489462) q[2];
sx q[2];
rz(-1.1558888) q[2];
sx q[2];
rz(0.75239858) q[2];
rz(2.803346) q[3];
sx q[3];
rz(-2.1460642) q[3];
sx q[3];
rz(1.2301302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
