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
rz(2.565413) q[0];
sx q[0];
rz(-1.6771069) q[0];
sx q[0];
rz(0.65281868) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(2.519156) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718423) q[0];
sx q[0];
rz(-1.4542219) q[0];
sx q[0];
rz(1.8844834) q[0];
rz(-1.2369701) q[2];
sx q[2];
rz(-0.99572748) q[2];
sx q[2];
rz(2.0890369) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7622432) q[1];
sx q[1];
rz(-2.6062638) q[1];
sx q[1];
rz(-1.204727) q[1];
rz(-pi) q[2];
rz(0.962019) q[3];
sx q[3];
rz(-1.7005657) q[3];
sx q[3];
rz(0.4421086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6301253) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(-2.5959065) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(0.95664501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0551374) q[0];
sx q[0];
rz(-2.4392023) q[0];
sx q[0];
rz(-1.0954683) q[0];
rz(-0.99758863) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(1.9445317) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124312) q[0];
sx q[0];
rz(-2.1502156) q[0];
sx q[0];
rz(-2.8377735) q[0];
rz(-pi) q[1];
rz(-0.62792553) q[2];
sx q[2];
rz(-0.88191477) q[2];
sx q[2];
rz(-0.21180001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45586387) q[1];
sx q[1];
rz(-1.2867754) q[1];
sx q[1];
rz(-1.0953254) q[1];
x q[2];
rz(1.7577111) q[3];
sx q[3];
rz(-1.8946365) q[3];
sx q[3];
rz(-1.327654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41584388) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(-1.3843298) q[2];
rz(1.7978801) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(-1.2725007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62190732) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(-0.4271048) q[0];
rz(-1.9097795) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(-0.61418358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012497525) q[0];
sx q[0];
rz(-1.84082) q[0];
sx q[0];
rz(0.34059033) q[0];
x q[1];
rz(1.7730447) q[2];
sx q[2];
rz(-0.83421153) q[2];
sx q[2];
rz(-1.6088161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78313561) q[1];
sx q[1];
rz(-1.6880717) q[1];
sx q[1];
rz(3.0774119) q[1];
x q[2];
rz(2.7216629) q[3];
sx q[3];
rz(-2.5961317) q[3];
sx q[3];
rz(0.88014102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4517639) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(2.8705987) q[2];
rz(2.6594035) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(-1.8577925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86376205) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(-2.7226287) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.9326262) q[1];
sx q[1];
rz(3.0640501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2231295) q[0];
sx q[0];
rz(-1.0498199) q[0];
sx q[0];
rz(-2.8023802) q[0];
rz(-pi) q[1];
rz(-2.3691142) q[2];
sx q[2];
rz(-1.8260406) q[2];
sx q[2];
rz(-2.8359063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0667013) q[1];
sx q[1];
rz(-0.46240004) q[1];
sx q[1];
rz(2.8179864) q[1];
rz(-pi) q[2];
rz(-0.61309149) q[3];
sx q[3];
rz(-2.7889851) q[3];
sx q[3];
rz(2.4168454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0858687) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(-2.870765) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643352) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(1.3539535) q[0];
rz(1.0352146) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.530102) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852762) q[0];
sx q[0];
rz(-2.0977328) q[0];
sx q[0];
rz(-0.88165347) q[0];
rz(-pi) q[1];
rz(-1.4064838) q[2];
sx q[2];
rz(-1.2466058) q[2];
sx q[2];
rz(2.6099082) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.93188796) q[1];
sx q[1];
rz(-1.6242289) q[1];
sx q[1];
rz(0.40219743) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13134457) q[3];
sx q[3];
rz(-0.41413158) q[3];
sx q[3];
rz(1.9589748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3349541) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(3.0777001) q[2];
rz(-2.4334) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(-1.3341058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719249) q[0];
sx q[0];
rz(-0.023119211) q[0];
sx q[0];
rz(1.0759906) q[0];
rz(0.05323449) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(-0.3784953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9747978) q[0];
sx q[0];
rz(-2.85436) q[0];
sx q[0];
rz(0.47606456) q[0];
rz(-pi) q[1];
rz(-1.6588062) q[2];
sx q[2];
rz(-1.529379) q[2];
sx q[2];
rz(-0.27829042) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6798181) q[1];
sx q[1];
rz(-1.8289806) q[1];
sx q[1];
rz(-1.3854909) q[1];
rz(-pi) q[2];
rz(-2.1355576) q[3];
sx q[3];
rz(-0.88712403) q[3];
sx q[3];
rz(-0.61236741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9516912) q[2];
sx q[2];
rz(-1.6804164) q[2];
sx q[2];
rz(1.0339197) q[2];
rz(2.2377491) q[3];
sx q[3];
rz(-2.240182) q[3];
sx q[3];
rz(-3.097539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251625) q[0];
sx q[0];
rz(-1.4465541) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(3.016839) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(0.4932901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3083444) q[0];
sx q[0];
rz(-1.8878636) q[0];
sx q[0];
rz(-0.29219101) q[0];
rz(0.10693018) q[2];
sx q[2];
rz(-2.9122346) q[2];
sx q[2];
rz(1.3170674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1715162) q[1];
sx q[1];
rz(-2.070721) q[1];
sx q[1];
rz(2.058567) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8658537) q[3];
sx q[3];
rz(-1.6746124) q[3];
sx q[3];
rz(-2.8732515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7571681) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(-1.1473131) q[2];
rz(-2.3186963) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(-2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.5230474) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(-0.4185032) q[0];
rz(-3.0283527) q[1];
sx q[1];
rz(-1.4694045) q[1];
sx q[1];
rz(-1.07771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61138701) q[0];
sx q[0];
rz(-2.8525087) q[0];
sx q[0];
rz(1.5075839) q[0];
x q[1];
rz(-1.244721) q[2];
sx q[2];
rz(-0.80279826) q[2];
sx q[2];
rz(-1.7434415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4227161) q[1];
sx q[1];
rz(-2.7759429) q[1];
sx q[1];
rz(-1.7083733) q[1];
rz(-pi) q[2];
rz(0.65656482) q[3];
sx q[3];
rz(-1.6180919) q[3];
sx q[3];
rz(1.8790203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8354127) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(-1.6861247) q[2];
rz(-2.9742187) q[3];
sx q[3];
rz(-2.6264329) q[3];
sx q[3];
rz(2.0696056) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587332) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(-2.7400548) q[0];
rz(0.43831476) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(-1.6848791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7244723) q[0];
sx q[0];
rz(-0.03532413) q[0];
sx q[0];
rz(-1.2958584) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28172617) q[2];
sx q[2];
rz(-0.67517074) q[2];
sx q[2];
rz(0.074568579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46692013) q[1];
sx q[1];
rz(-2.2385257) q[1];
sx q[1];
rz(0.93514644) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0450122) q[3];
sx q[3];
rz(-2.0353531) q[3];
sx q[3];
rz(-1.7024226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7598286) q[2];
sx q[2];
rz(-2.9324053) q[2];
sx q[2];
rz(-0.71167243) q[2];
rz(3.107374) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3592247) q[0];
sx q[0];
rz(-1.490626) q[0];
sx q[0];
rz(0.79886287) q[0];
rz(2.0955739) q[1];
sx q[1];
rz(-1.9452399) q[1];
sx q[1];
rz(-2.9050713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4170161) q[0];
sx q[0];
rz(-1.5484705) q[0];
sx q[0];
rz(1.5100368) q[0];
x q[1];
rz(0.92774074) q[2];
sx q[2];
rz(-2.1143713) q[2];
sx q[2];
rz(1.7801931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1139086) q[1];
sx q[1];
rz(-0.56237537) q[1];
sx q[1];
rz(2.4706868) q[1];
x q[2];
rz(0.71098401) q[3];
sx q[3];
rz(-1.4812638) q[3];
sx q[3];
rz(2.3481365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.025042621) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(-0.64209783) q[2];
rz(-2.5721278) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(-3.0552982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(-1.9500465) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(-0.12267648) q[2];
sx q[2];
rz(-2.3458866) q[2];
sx q[2];
rz(2.2550368) q[2];
rz(0.91611923) q[3];
sx q[3];
rz(-0.51894938) q[3];
sx q[3];
rz(2.9872158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
