OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(2.3072825) q[0];
sx q[0];
rz(6.351525) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(0.037820427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67536394) q[0];
sx q[0];
rz(-0.20875202) q[0];
sx q[0];
rz(-1.1842313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4191364) q[2];
sx q[2];
rz(-2.4745387) q[2];
sx q[2];
rz(2.5487713) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.5303858) q[1];
sx q[1];
rz(-1.621765) q[1];
sx q[1];
rz(-2.3947122) q[1];
rz(-pi) q[2];
rz(-1.2245523) q[3];
sx q[3];
rz(-1.7131221) q[3];
sx q[3];
rz(-3.0179253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7906856) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(-1.2834056) q[2];
rz(-3.0154199) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44777563) q[0];
sx q[0];
rz(-2.2651146) q[0];
sx q[0];
rz(-0.59666657) q[0];
rz(1.5555351) q[1];
sx q[1];
rz(-1.7998453) q[1];
sx q[1];
rz(-1.3635051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024433) q[0];
sx q[0];
rz(-1.9958695) q[0];
sx q[0];
rz(-2.5545679) q[0];
rz(-1.5200137) q[2];
sx q[2];
rz(-2.1261458) q[2];
sx q[2];
rz(2.3107049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5982358) q[1];
sx q[1];
rz(-2.1992866) q[1];
sx q[1];
rz(-1.223982) q[1];
rz(-pi) q[2];
rz(0.76916839) q[3];
sx q[3];
rz(-2.2393919) q[3];
sx q[3];
rz(-1.9070966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3229225) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(-2.4070516) q[2];
rz(-0.62720403) q[3];
sx q[3];
rz(-1.6278798) q[3];
sx q[3];
rz(2.396615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4194141) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(-0.27221361) q[0];
rz(2.294337) q[1];
sx q[1];
rz(-1.8224742) q[1];
sx q[1];
rz(2.8289657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659098) q[0];
sx q[0];
rz(-1.624375) q[0];
sx q[0];
rz(0.21955325) q[0];
x q[1];
rz(-2.0414511) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(0.91980308) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9885013) q[1];
sx q[1];
rz(-1.830849) q[1];
sx q[1];
rz(-3.0696763) q[1];
rz(0.76287855) q[3];
sx q[3];
rz(-1.2925914) q[3];
sx q[3];
rz(-2.0335846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0212038) q[2];
sx q[2];
rz(-0.52336064) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(0.96873823) q[3];
sx q[3];
rz(-2.0961943) q[3];
sx q[3];
rz(2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75354904) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(0.62336212) q[0];
rz(2.3240044) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-0.049291074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76694861) q[0];
sx q[0];
rz(-0.76061237) q[0];
sx q[0];
rz(1.6409372) q[0];
rz(-2.8231499) q[2];
sx q[2];
rz(-2.2568984) q[2];
sx q[2];
rz(2.7903914) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99187461) q[1];
sx q[1];
rz(-1.8135934) q[1];
sx q[1];
rz(-1.2319831) q[1];
rz(-pi) q[2];
rz(-2.0972581) q[3];
sx q[3];
rz(-1.1757848) q[3];
sx q[3];
rz(-2.928424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7307044) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(0.98497406) q[2];
rz(-2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(0.27339098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34489283) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(0.087619089) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(0.87096754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317469) q[0];
sx q[0];
rz(-1.5994659) q[0];
sx q[0];
rz(-0.28733758) q[0];
rz(1.2147374) q[2];
sx q[2];
rz(-2.7241754) q[2];
sx q[2];
rz(0.88127121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37967967) q[1];
sx q[1];
rz(-1.2907791) q[1];
sx q[1];
rz(0.67358394) q[1];
rz(-pi) q[2];
rz(2.3351339) q[3];
sx q[3];
rz(-1.4440923) q[3];
sx q[3];
rz(0.63092953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0481723) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(-0.86587632) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(-2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46309328) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(-3.034806) q[0];
rz(-1.9550025) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(-2.1170763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2183285) q[0];
sx q[0];
rz(-0.54696333) q[0];
sx q[0];
rz(-2.5640019) q[0];
rz(-pi) q[1];
rz(-1.1149939) q[2];
sx q[2];
rz(-1.015128) q[2];
sx q[2];
rz(1.232604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66631324) q[1];
sx q[1];
rz(-2.3716455) q[1];
sx q[1];
rz(2.2837) q[1];
rz(0.49196524) q[3];
sx q[3];
rz(-1.8987738) q[3];
sx q[3];
rz(-0.63719751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(-0.8708896) q[2];
rz(1.332256) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(-1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-0.55091888) q[0];
rz(2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(2.904772) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91009249) q[0];
sx q[0];
rz(-1.5331368) q[0];
sx q[0];
rz(-1.670027) q[0];
rz(-pi) q[1];
rz(-2.4004164) q[2];
sx q[2];
rz(-0.20222649) q[2];
sx q[2];
rz(-1.9066332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9501818) q[1];
sx q[1];
rz(-1.6457874) q[1];
sx q[1];
rz(1.300315) q[1];
x q[2];
rz(-0.3053815) q[3];
sx q[3];
rz(-0.85994342) q[3];
sx q[3];
rz(1.4399432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6992496) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(2.701475) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(-1.7949036) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(-0.029504689) q[0];
rz(-0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(0.11880076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7293852) q[0];
sx q[0];
rz(-2.8280624) q[0];
sx q[0];
rz(0.14051147) q[0];
rz(-2.115909) q[2];
sx q[2];
rz(-1.1720177) q[2];
sx q[2];
rz(1.1023956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.48092948) q[1];
sx q[1];
rz(-0.56023635) q[1];
sx q[1];
rz(-0.57628298) q[1];
x q[2];
rz(-0.16468594) q[3];
sx q[3];
rz(-1.0822923) q[3];
sx q[3];
rz(0.96795852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7416731) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(-1.5681533) q[2];
rz(-1.9741156) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-0.82505834) q[0];
sx q[0];
rz(0.15326823) q[0];
rz(1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.9877888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.518084) q[0];
sx q[0];
rz(-2.6797446) q[0];
sx q[0];
rz(-2.5005546) q[0];
rz(2.2527163) q[2];
sx q[2];
rz(-1.1317012) q[2];
sx q[2];
rz(-0.097749226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.02281636) q[1];
sx q[1];
rz(-2.9252242) q[1];
sx q[1];
rz(1.3092036) q[1];
rz(-0.12113573) q[3];
sx q[3];
rz(-0.88705685) q[3];
sx q[3];
rz(-0.63431206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29685059) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(1.5273757) q[2];
rz(-2.1447003) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(-2.2629288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24755724) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(0.19009185) q[0];
rz(-0.67063531) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(-1.6533096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0310693) q[0];
sx q[0];
rz(-3.0123683) q[0];
sx q[0];
rz(0.24137012) q[0];
rz(1.8411631) q[2];
sx q[2];
rz(-1.1425945) q[2];
sx q[2];
rz(-0.92744517) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59939811) q[1];
sx q[1];
rz(-0.8510667) q[1];
sx q[1];
rz(0.69516121) q[1];
rz(-0.90115746) q[3];
sx q[3];
rz(-0.81448758) q[3];
sx q[3];
rz(-2.8709473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1378479) q[2];
sx q[2];
rz(-2.2403084) q[2];
sx q[2];
rz(-2.2424662) q[2];
rz(-2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(-0.17764828) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0777733) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(0.28941119) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(-0.32439705) q[2];
sx q[2];
rz(-2.4051721) q[2];
sx q[2];
rz(0.13798513) q[2];
rz(-1.297847) q[3];
sx q[3];
rz(-1.6932586) q[3];
sx q[3];
rz(0.3441588) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
