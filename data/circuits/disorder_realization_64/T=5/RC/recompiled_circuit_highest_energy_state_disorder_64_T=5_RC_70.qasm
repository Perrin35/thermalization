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
rz(2.3447073) q[0];
sx q[0];
rz(-1.7303884) q[0];
sx q[0];
rz(2.222173) q[0];
rz(-1.8384276) q[1];
sx q[1];
rz(-3.0923831) q[1];
sx q[1];
rz(-2.4278909) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588736) q[0];
sx q[0];
rz(-1.8451515) q[0];
sx q[0];
rz(0.15299847) q[0];
x q[1];
rz(2.0242516) q[2];
sx q[2];
rz(-2.9988117) q[2];
sx q[2];
rz(0.26856865) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9956869) q[1];
sx q[1];
rz(-0.45098454) q[1];
sx q[1];
rz(-0.0031975071) q[1];
x q[2];
rz(1.0994751) q[3];
sx q[3];
rz(-1.6013535) q[3];
sx q[3];
rz(1.9224642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0539703) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(2.7719882) q[2];
rz(-1.1450279) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0358148) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(2.5685487) q[0];
rz(1.0967968) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(0.47164741) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1435463) q[0];
sx q[0];
rz(-2.5295288) q[0];
sx q[0];
rz(1.7630811) q[0];
rz(-pi) q[1];
rz(0.69219442) q[2];
sx q[2];
rz(-1.0613777) q[2];
sx q[2];
rz(2.8329527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14395606) q[1];
sx q[1];
rz(-2.156971) q[1];
sx q[1];
rz(2.2464941) q[1];
rz(-0.92949683) q[3];
sx q[3];
rz(-2.5166582) q[3];
sx q[3];
rz(0.24111023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(1.088885) q[2];
rz(1.0645083) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0629145) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(-2.4315289) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(-2.5753218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44646548) q[0];
sx q[0];
rz(-0.74120616) q[0];
sx q[0];
rz(-1.5873853) q[0];
rz(-pi) q[1];
x q[1];
rz(2.429661) q[2];
sx q[2];
rz(-1.997479) q[2];
sx q[2];
rz(2.349145) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0349087) q[1];
sx q[1];
rz(-1.7572173) q[1];
sx q[1];
rz(-1.6572464) q[1];
rz(-pi) q[2];
rz(-2.8763883) q[3];
sx q[3];
rz(-1.5062766) q[3];
sx q[3];
rz(-1.2539188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8038586) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(-3.1324978) q[2];
rz(3.005262) q[3];
sx q[3];
rz(-0.74458849) q[3];
sx q[3];
rz(-1.2135308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71753865) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(0.16615443) q[0];
rz(2.0280139) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(0.16214935) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1095162) q[0];
sx q[0];
rz(-1.6639532) q[0];
sx q[0];
rz(1.7514125) q[0];
rz(3.1012898) q[2];
sx q[2];
rz(-0.54931123) q[2];
sx q[2];
rz(1.2072762) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.6238193) q[1];
sx q[1];
rz(-0.86408593) q[1];
sx q[1];
rz(2.049974) q[1];
x q[2];
rz(1.9982073) q[3];
sx q[3];
rz(-0.95164883) q[3];
sx q[3];
rz(2.038508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98264155) q[2];
sx q[2];
rz(-2.0621767) q[2];
sx q[2];
rz(0.52784935) q[2];
rz(-0.85150254) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(-0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.270179) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(2.340509) q[0];
rz(0.78394765) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(0.98091006) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0410515) q[0];
sx q[0];
rz(-1.1642191) q[0];
sx q[0];
rz(-1.340614) q[0];
rz(-pi) q[1];
rz(-2.7655914) q[2];
sx q[2];
rz(-2.0472976) q[2];
sx q[2];
rz(-0.27793542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.09771654) q[1];
sx q[1];
rz(-1.6642769) q[1];
sx q[1];
rz(0.059475617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9526236) q[3];
sx q[3];
rz(-2.5562048) q[3];
sx q[3];
rz(2.8849734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(-0.8832461) q[2];
rz(-2.9412681) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9969295) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(-2.1494179) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(-1.7209524) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29692263) q[0];
sx q[0];
rz(-1.8886744) q[0];
sx q[0];
rz(-2.4852089) q[0];
rz(-pi) q[1];
rz(0.10481609) q[2];
sx q[2];
rz(-1.040852) q[2];
sx q[2];
rz(-1.3111834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1375422) q[1];
sx q[1];
rz(-0.95064771) q[1];
sx q[1];
rz(-1.1637079) q[1];
rz(-pi) q[2];
rz(-1.4735953) q[3];
sx q[3];
rz(-1.5335113) q[3];
sx q[3];
rz(0.71747045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27986032) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(-0.53708616) q[2];
rz(-0.01072695) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456197) q[0];
sx q[0];
rz(-0.27181044) q[0];
sx q[0];
rz(-0.33988345) q[0];
rz(-1.3019568) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(1.1713015) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3070967) q[0];
sx q[0];
rz(-2.8423873) q[0];
sx q[0];
rz(-0.12339573) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7843855) q[2];
sx q[2];
rz(-1.3554975) q[2];
sx q[2];
rz(2.4395296) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9885607) q[1];
sx q[1];
rz(-1.5383729) q[1];
sx q[1];
rz(0.02213636) q[1];
x q[2];
rz(-2.9006935) q[3];
sx q[3];
rz(-0.63848513) q[3];
sx q[3];
rz(-2.1158259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12477144) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-2.4701414) q[2];
rz(-0.5365544) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(-2.0898537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3143828) q[0];
sx q[0];
rz(-2.0874513) q[0];
sx q[0];
rz(2.2453454) q[0];
rz(0.25262901) q[1];
sx q[1];
rz(-2.2815506) q[1];
sx q[1];
rz(-1.0505229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3327961) q[0];
sx q[0];
rz(-2.0967297) q[0];
sx q[0];
rz(-1.6288847) q[0];
rz(-2.4765313) q[2];
sx q[2];
rz(-1.6740693) q[2];
sx q[2];
rz(-2.2879083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0493361) q[1];
sx q[1];
rz(-0.4886371) q[1];
sx q[1];
rz(1.8948003) q[1];
x q[2];
rz(-1.6231617) q[3];
sx q[3];
rz(-2.3279421) q[3];
sx q[3];
rz(-2.6770075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5184021) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(-1.5790348) q[2];
rz(1.3671499) q[3];
sx q[3];
rz(-2.0953777) q[3];
sx q[3];
rz(0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52687454) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(2.7516464) q[0];
rz(-2.3155616) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(1.0521851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1601801) q[0];
sx q[0];
rz(-0.76428586) q[0];
sx q[0];
rz(2.0212964) q[0];
rz(-pi) q[1];
rz(1.7462535) q[2];
sx q[2];
rz(-1.1793841) q[2];
sx q[2];
rz(3.1188426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3996399) q[1];
sx q[1];
rz(-1.1832848) q[1];
sx q[1];
rz(-0.31419973) q[1];
rz(-0.47635079) q[3];
sx q[3];
rz(-0.99925502) q[3];
sx q[3];
rz(2.0921354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(-0.52919394) q[2];
rz(-0.75096327) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(-0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879409) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(-1.9770812) q[0];
rz(1.8544082) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(-0.50416344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34216322) q[0];
sx q[0];
rz(-1.829188) q[0];
sx q[0];
rz(2.9151205) q[0];
rz(-pi) q[1];
rz(1.590056) q[2];
sx q[2];
rz(-1.9474721) q[2];
sx q[2];
rz(0.45805675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.530624) q[1];
sx q[1];
rz(-1.305831) q[1];
sx q[1];
rz(0.19719657) q[1];
rz(-pi) q[2];
rz(-2.9243117) q[3];
sx q[3];
rz(-0.8978399) q[3];
sx q[3];
rz(-0.5289883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9265499) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(3.1070993) q[2];
rz(-1.6132332) q[3];
sx q[3];
rz(-2.4210763) q[3];
sx q[3];
rz(-0.31149402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16019776) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(-1.002671) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(-0.38242292) q[2];
sx q[2];
rz(-2.6527846) q[2];
sx q[2];
rz(-2.6058886) q[2];
rz(-2.814449) q[3];
sx q[3];
rz(-1.5113304) q[3];
sx q[3];
rz(-2.3979901) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
