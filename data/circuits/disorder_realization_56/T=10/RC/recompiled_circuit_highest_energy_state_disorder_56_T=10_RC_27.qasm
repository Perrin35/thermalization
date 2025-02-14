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
rz(2.7409878) q[0];
sx q[0];
rz(-0.40979835) q[0];
sx q[0];
rz(1.0709437) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(3.7096042) q[1];
sx q[1];
rz(8.5811442) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7535665) q[0];
sx q[0];
rz(-0.26615289) q[0];
sx q[0];
rz(1.4583336) q[0];
rz(-pi) q[1];
rz(2.8461692) q[2];
sx q[2];
rz(-1.7082126) q[2];
sx q[2];
rz(0.58770056) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93701332) q[1];
sx q[1];
rz(-0.76256231) q[1];
sx q[1];
rz(-1.0916876) q[1];
rz(2.6671404) q[3];
sx q[3];
rz(-1.1233002) q[3];
sx q[3];
rz(1.8008658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3598651) q[2];
sx q[2];
rz(-2.3024776) q[2];
sx q[2];
rz(-0.016999379) q[2];
rz(1.4012236) q[3];
sx q[3];
rz(-2.5798116) q[3];
sx q[3];
rz(1.9248272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725937) q[0];
sx q[0];
rz(-1.3324791) q[0];
sx q[0];
rz(0.78729415) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(1.9040727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2923861) q[0];
sx q[0];
rz(-0.9885177) q[0];
sx q[0];
rz(1.4968064) q[0];
rz(1.9745047) q[2];
sx q[2];
rz(-2.4879527) q[2];
sx q[2];
rz(2.6622651) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.552265) q[1];
sx q[1];
rz(-0.54937148) q[1];
sx q[1];
rz(2.2087847) q[1];
x q[2];
rz(-0.54333289) q[3];
sx q[3];
rz(-2.6381734) q[3];
sx q[3];
rz(1.0375298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7201207) q[2];
sx q[2];
rz(-1.4482435) q[2];
sx q[2];
rz(0.063701542) q[2];
rz(1.8675768) q[3];
sx q[3];
rz(-0.84309045) q[3];
sx q[3];
rz(-1.5773076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.098123) q[0];
sx q[0];
rz(-1.1935357) q[0];
sx q[0];
rz(2.3980339) q[0];
rz(0.45383635) q[1];
sx q[1];
rz(-1.3498787) q[1];
sx q[1];
rz(2.7226864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009368816) q[0];
sx q[0];
rz(-1.1282578) q[0];
sx q[0];
rz(-0.21712961) q[0];
x q[1];
rz(2.3055196) q[2];
sx q[2];
rz(-1.3385069) q[2];
sx q[2];
rz(0.79233263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6082904) q[1];
sx q[1];
rz(-0.69606298) q[1];
sx q[1];
rz(-1.1185557) q[1];
rz(2.2499372) q[3];
sx q[3];
rz(-1.4640995) q[3];
sx q[3];
rz(1.5859491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0539187) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(1.9107001) q[2];
rz(-0.538921) q[3];
sx q[3];
rz(-1.475622) q[3];
sx q[3];
rz(-1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8777799) q[0];
sx q[0];
rz(-2.3872264) q[0];
sx q[0];
rz(0.60212773) q[0];
rz(-1.4471794) q[1];
sx q[1];
rz(-2.4592631) q[1];
sx q[1];
rz(0.86300659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857638) q[0];
sx q[0];
rz(-1.6842972) q[0];
sx q[0];
rz(-0.72300006) q[0];
rz(-0.30089897) q[2];
sx q[2];
rz(-1.2105807) q[2];
sx q[2];
rz(0.27733251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80490404) q[1];
sx q[1];
rz(-2.4987578) q[1];
sx q[1];
rz(-2.0175319) q[1];
rz(-0.3163655) q[3];
sx q[3];
rz(-2.2834407) q[3];
sx q[3];
rz(1.0567088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.065923125) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(2.2389331) q[2];
rz(-2.3265808) q[3];
sx q[3];
rz(-0.60324001) q[3];
sx q[3];
rz(2.7284315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0937423) q[0];
sx q[0];
rz(-3.0939026) q[0];
sx q[0];
rz(0.93511859) q[0];
rz(-1.6085666) q[1];
sx q[1];
rz(-0.32564274) q[1];
sx q[1];
rz(2.8757222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1810721) q[0];
sx q[0];
rz(-1.4645394) q[0];
sx q[0];
rz(-1.8828859) q[0];
x q[1];
rz(2.8433617) q[2];
sx q[2];
rz(-0.77285337) q[2];
sx q[2];
rz(-1.4082343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8740377) q[1];
sx q[1];
rz(-1.4103762) q[1];
sx q[1];
rz(-2.7880413) q[1];
x q[2];
rz(-1.3060027) q[3];
sx q[3];
rz(-1.4715428) q[3];
sx q[3];
rz(0.81391993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68044668) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(-2.6908596) q[2];
rz(1.1905253) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(0.43959555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615237) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(-0.044483749) q[0];
rz(-2.8728409) q[1];
sx q[1];
rz(-0.93695295) q[1];
sx q[1];
rz(0.43089795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.160153) q[0];
sx q[0];
rz(-2.4593076) q[0];
sx q[0];
rz(-0.42295608) q[0];
rz(-1.0111647) q[2];
sx q[2];
rz(-0.75087386) q[2];
sx q[2];
rz(-2.4023285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84316501) q[1];
sx q[1];
rz(-0.80363217) q[1];
sx q[1];
rz(-2.1433565) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5650125) q[3];
sx q[3];
rz(-0.91051379) q[3];
sx q[3];
rz(-0.74951142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0253133) q[2];
sx q[2];
rz(-2.7677324) q[2];
sx q[2];
rz(0.6244134) q[2];
rz(-1.1011018) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(-0.98744923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91200149) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(-0.68369317) q[0];
rz(1.6992441) q[1];
sx q[1];
rz(-0.78386274) q[1];
sx q[1];
rz(2.9396465) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.604902) q[0];
sx q[0];
rz(-0.32819191) q[0];
sx q[0];
rz(1.5668014) q[0];
rz(-pi) q[1];
rz(-1.6228637) q[2];
sx q[2];
rz(-0.74801842) q[2];
sx q[2];
rz(-1.5128653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4836639) q[1];
sx q[1];
rz(-0.52573181) q[1];
sx q[1];
rz(-3.0083198) q[1];
x q[2];
rz(-0.087183909) q[3];
sx q[3];
rz(-2.6823061) q[3];
sx q[3];
rz(0.82827866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.086229) q[2];
sx q[2];
rz(-1.4198885) q[2];
sx q[2];
rz(1.281338) q[2];
rz(2.4751439) q[3];
sx q[3];
rz(-1.0969578) q[3];
sx q[3];
rz(2.5950529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619038) q[0];
sx q[0];
rz(-2.5894916) q[0];
sx q[0];
rz(0.51396489) q[0];
rz(-0.95298302) q[1];
sx q[1];
rz(-2.516808) q[1];
sx q[1];
rz(2.0228588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3866769) q[0];
sx q[0];
rz(-2.710452) q[0];
sx q[0];
rz(0.19622959) q[0];
rz(-pi) q[1];
rz(1.1360815) q[2];
sx q[2];
rz(-0.82821199) q[2];
sx q[2];
rz(-1.0880053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36597628) q[1];
sx q[1];
rz(-1.9785641) q[1];
sx q[1];
rz(2.9080176) q[1];
rz(-pi) q[2];
rz(0.54688485) q[3];
sx q[3];
rz(-2.5207479) q[3];
sx q[3];
rz(1.8016694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.755456) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(2.238671) q[2];
rz(-1.4593982) q[3];
sx q[3];
rz(-1.7995588) q[3];
sx q[3];
rz(0.82227796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4626386) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(-1.951304) q[0];
rz(-0.32084385) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(-1.2215337) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71305481) q[0];
sx q[0];
rz(-0.24745169) q[0];
sx q[0];
rz(2.5525981) q[0];
rz(-pi) q[1];
rz(-0.83619976) q[2];
sx q[2];
rz(-1.3635907) q[2];
sx q[2];
rz(1.5839603) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22449422) q[1];
sx q[1];
rz(-2.3116744) q[1];
sx q[1];
rz(1.8486449) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.753958) q[3];
sx q[3];
rz(-1.7622928) q[3];
sx q[3];
rz(0.93246704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(0.3332738) q[2];
rz(-1.0171558) q[3];
sx q[3];
rz(-0.95484304) q[3];
sx q[3];
rz(-2.8216383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3453813) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(-2.2228125) q[0];
rz(2.9504919) q[1];
sx q[1];
rz(-0.42593503) q[1];
sx q[1];
rz(2.3474615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.213906) q[0];
sx q[0];
rz(-1.0451259) q[0];
sx q[0];
rz(2.3679683) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3599858) q[2];
sx q[2];
rz(-1.2173869) q[2];
sx q[2];
rz(1.7164053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5743374) q[1];
sx q[1];
rz(-0.87346389) q[1];
sx q[1];
rz(0.35178565) q[1];
x q[2];
rz(2.8201032) q[3];
sx q[3];
rz(-1.3703385) q[3];
sx q[3];
rz(-3.0908302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0309151) q[2];
sx q[2];
rz(-0.6157178) q[2];
sx q[2];
rz(0.75882971) q[2];
rz(2.2714254) q[3];
sx q[3];
rz(-0.78799677) q[3];
sx q[3];
rz(1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0017241521) q[0];
sx q[0];
rz(-1.682946) q[0];
sx q[0];
rz(-1.2276822) q[0];
rz(2.7036746) q[1];
sx q[1];
rz(-1.7057849) q[1];
sx q[1];
rz(1.5954856) q[1];
rz(1.154806) q[2];
sx q[2];
rz(-2.0218973) q[2];
sx q[2];
rz(0.99662957) q[2];
rz(2.6043456) q[3];
sx q[3];
rz(-2.2840847) q[3];
sx q[3];
rz(1.7084697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
