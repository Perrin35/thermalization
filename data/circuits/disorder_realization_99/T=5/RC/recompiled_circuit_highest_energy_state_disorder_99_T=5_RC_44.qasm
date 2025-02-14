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
rz(-1.7412269) q[0];
sx q[0];
rz(-2.9018612) q[0];
sx q[0];
rz(2.8170407) q[0];
rz(-0.95207721) q[1];
sx q[1];
rz(-0.2258741) q[1];
sx q[1];
rz(1.3083375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150779) q[0];
sx q[0];
rz(-2.240772) q[0];
sx q[0];
rz(1.7448533) q[0];
rz(-pi) q[1];
rz(-0.4842437) q[2];
sx q[2];
rz(-1.4560369) q[2];
sx q[2];
rz(1.6387303) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5950299) q[1];
sx q[1];
rz(-1.1076704) q[1];
sx q[1];
rz(2.8519277) q[1];
rz(-2.9746858) q[3];
sx q[3];
rz(-0.70611533) q[3];
sx q[3];
rz(2.0049068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7294881) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(2.7903902) q[2];
rz(1.8515733) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(1.2235519) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(-2.2112041) q[0];
rz(1.8366086) q[1];
sx q[1];
rz(-2.3001859) q[1];
sx q[1];
rz(-2.5321541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6822426) q[0];
sx q[0];
rz(-1.1455904) q[0];
sx q[0];
rz(-0.11328477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8687042) q[2];
sx q[2];
rz(-1.4845843) q[2];
sx q[2];
rz(-0.49315587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48991443) q[1];
sx q[1];
rz(-0.53333542) q[1];
sx q[1];
rz(2.3467605) q[1];
x q[2];
rz(-1.3136765) q[3];
sx q[3];
rz(-2.5559396) q[3];
sx q[3];
rz(2.3975092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0448407) q[2];
sx q[2];
rz(-2.1614306) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(-2.6212202) q[3];
sx q[3];
rz(-0.64295355) q[3];
sx q[3];
rz(-0.61409942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60270131) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-0.30211788) q[0];
rz(0.27613861) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(1.709323) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0841504) q[0];
sx q[0];
rz(-1.0330811) q[0];
sx q[0];
rz(-2.6938963) q[0];
rz(-pi) q[1];
x q[1];
rz(1.309607) q[2];
sx q[2];
rz(-2.1495594) q[2];
sx q[2];
rz(-2.8566993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4319182) q[1];
sx q[1];
rz(-2.1248105) q[1];
sx q[1];
rz(0.18504237) q[1];
rz(-pi) q[2];
rz(-2.0366482) q[3];
sx q[3];
rz(-1.5072952) q[3];
sx q[3];
rz(-2.7257763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3064731) q[2];
sx q[2];
rz(-2.6758631) q[2];
sx q[2];
rz(1.2960557) q[2];
rz(-0.45201388) q[3];
sx q[3];
rz(-1.1072423) q[3];
sx q[3];
rz(-0.38255102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6119659) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(1.3324598) q[0];
rz(-2.0924856) q[1];
sx q[1];
rz(-2.0499947) q[1];
sx q[1];
rz(-1.5528991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72399607) q[0];
sx q[0];
rz(-1.231047) q[0];
sx q[0];
rz(-1.8906361) q[0];
x q[1];
rz(1.5459486) q[2];
sx q[2];
rz(-2.6406807) q[2];
sx q[2];
rz(0.76297578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9178333) q[1];
sx q[1];
rz(-2.4661015) q[1];
sx q[1];
rz(1.89066) q[1];
rz(0.53557204) q[3];
sx q[3];
rz(-2.2102997) q[3];
sx q[3];
rz(-2.4909693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66854746) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(-2.891053) q[2];
rz(2.1446832) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(2.7610049) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8119891) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(1.3478152) q[0];
rz(2.7373121) q[1];
sx q[1];
rz(-0.75899044) q[1];
sx q[1];
rz(1.1291198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4816138) q[0];
sx q[0];
rz(-0.82757512) q[0];
sx q[0];
rz(-1.0218234) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8724832) q[2];
sx q[2];
rz(-2.2186612) q[2];
sx q[2];
rz(1.3307759) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80218643) q[1];
sx q[1];
rz(-1.6755584) q[1];
sx q[1];
rz(2.2441007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4604767) q[3];
sx q[3];
rz(-1.2577783) q[3];
sx q[3];
rz(2.4576994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3471442) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(-1.8974737) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4853972) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(0.42718497) q[0];
rz(-3.1071013) q[1];
sx q[1];
rz(-1.5692312) q[1];
sx q[1];
rz(-0.010206612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6319868) q[0];
sx q[0];
rz(-1.5684396) q[0];
sx q[0];
rz(-1.5725333) q[0];
x q[1];
rz(1.7882657) q[2];
sx q[2];
rz(-0.64369338) q[2];
sx q[2];
rz(0.74483192) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2040703) q[1];
sx q[1];
rz(-2.2559705) q[1];
sx q[1];
rz(-3.0584945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56513803) q[3];
sx q[3];
rz(-1.1332773) q[3];
sx q[3];
rz(2.0034307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61937845) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(0.79356066) q[2];
rz(2.2788952) q[3];
sx q[3];
rz(-1.5746652) q[3];
sx q[3];
rz(-0.70393744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70235395) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(-2.0080361) q[0];
rz(-1.0776862) q[1];
sx q[1];
rz(-0.51529854) q[1];
sx q[1];
rz(-0.30212197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568475) q[0];
sx q[0];
rz(-2.6805566) q[0];
sx q[0];
rz(2.4088944) q[0];
rz(-pi) q[1];
rz(-2.5446645) q[2];
sx q[2];
rz(-1.9608627) q[2];
sx q[2];
rz(-2.3204107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1715917) q[1];
sx q[1];
rz(-1.0495032) q[1];
sx q[1];
rz(1.820351) q[1];
rz(-pi) q[2];
rz(-2.1572635) q[3];
sx q[3];
rz(-2.1696089) q[3];
sx q[3];
rz(0.58517712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9083531) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(-1.7944149) q[2];
rz(-1.2498648) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5477448) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(0.10840848) q[0];
rz(-1.0477061) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(-1.0188867) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86025809) q[0];
sx q[0];
rz(-2.0906013) q[0];
sx q[0];
rz(2.5297574) q[0];
rz(-pi) q[1];
x q[1];
rz(1.564304) q[2];
sx q[2];
rz(-0.80950709) q[2];
sx q[2];
rz(-2.6463375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26886235) q[1];
sx q[1];
rz(-2.514317) q[1];
sx q[1];
rz(-0.71898798) q[1];
rz(1.4298444) q[3];
sx q[3];
rz(-0.8356072) q[3];
sx q[3];
rz(2.6648389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63742796) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(2.7039779) q[2];
rz(-0.49312433) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(1.5776618) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86579943) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(2.8637874) q[0];
rz(-1.5996784) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(0.42617282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8086373) q[0];
sx q[0];
rz(-1.5530968) q[0];
sx q[0];
rz(-1.5874392) q[0];
x q[1];
rz(0.73400273) q[2];
sx q[2];
rz(-2.2224769) q[2];
sx q[2];
rz(2.3190448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81282367) q[1];
sx q[1];
rz(-1.4534833) q[1];
sx q[1];
rz(-1.3552865) q[1];
x q[2];
rz(-0.79094633) q[3];
sx q[3];
rz(-2.2344347) q[3];
sx q[3];
rz(-2.5928465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49491945) q[2];
sx q[2];
rz(-2.9743331) q[2];
sx q[2];
rz(-1.0529998) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.6152629) q[3];
sx q[3];
rz(1.0927965) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36295715) q[0];
sx q[0];
rz(-2.9243922) q[0];
sx q[0];
rz(-0.92078513) q[0];
rz(-2.6329363) q[1];
sx q[1];
rz(-2.3743036) q[1];
sx q[1];
rz(2.7210534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79179278) q[0];
sx q[0];
rz(-2.0689488) q[0];
sx q[0];
rz(2.34116) q[0];
rz(-pi) q[1];
rz(1.8968209) q[2];
sx q[2];
rz(-2.4851126) q[2];
sx q[2];
rz(-0.028353779) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1433761) q[1];
sx q[1];
rz(-0.80165473) q[1];
sx q[1];
rz(2.2364535) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0402571) q[3];
sx q[3];
rz(-2.2592415) q[3];
sx q[3];
rz(-2.7197414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4113808) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(-2.8249557) q[2];
rz(-0.5591048) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(0.81812286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.485514) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(2.6651233) q[1];
sx q[1];
rz(-2.10119) q[1];
sx q[1];
rz(1.4269921) q[1];
rz(-3.0546679) q[2];
sx q[2];
rz(-2.3200421) q[2];
sx q[2];
rz(2.5210597) q[2];
rz(1.8000525) q[3];
sx q[3];
rz(-2.7490902) q[3];
sx q[3];
rz(0.31828087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
