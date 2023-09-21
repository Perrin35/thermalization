OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(5.4260317) q[0];
sx q[0];
rz(9.5572588) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(5.6981882) q[1];
sx q[1];
rz(11.873801) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1231175) q[0];
sx q[0];
rz(-1.0536061) q[0];
sx q[0];
rz(0.11488199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5508467) q[2];
sx q[2];
rz(-0.735539) q[2];
sx q[2];
rz(0.22434805) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9308656) q[1];
sx q[1];
rz(-2.9949246) q[1];
sx q[1];
rz(2.0861097) q[1];
rz(-pi) q[2];
rz(1.1107221) q[3];
sx q[3];
rz(-0.51957909) q[3];
sx q[3];
rz(-0.67668623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0341558) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(1.5365323) q[2];
rz(1.6202554) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(0.03604123) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(1.2492299) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(2.5610279) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3378355) q[0];
sx q[0];
rz(-2.8135186) q[0];
sx q[0];
rz(2.8003545) q[0];
rz(-pi) q[1];
rz(-1.7891907) q[2];
sx q[2];
rz(-2.1973655) q[2];
sx q[2];
rz(-0.44109694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9241582) q[1];
sx q[1];
rz(-2.1719451) q[1];
sx q[1];
rz(-2.2648328) q[1];
x q[2];
rz(1.1167691) q[3];
sx q[3];
rz(-1.8899711) q[3];
sx q[3];
rz(1.9516731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(2.2632329) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(-0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(-0.032827854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093124495) q[0];
sx q[0];
rz(-1.3402407) q[0];
sx q[0];
rz(-2.8218517) q[0];
x q[1];
rz(-2.122934) q[2];
sx q[2];
rz(-1.820192) q[2];
sx q[2];
rz(0.92555911) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58363885) q[1];
sx q[1];
rz(-1.4065521) q[1];
sx q[1];
rz(-0.79973952) q[1];
rz(-2.114931) q[3];
sx q[3];
rz(-1.2715724) q[3];
sx q[3];
rz(2.2194089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(-0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70401496) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-0.26279703) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(2.3707726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17731006) q[0];
sx q[0];
rz(-1.5543803) q[0];
sx q[0];
rz(1.55127) q[0];
rz(-pi) q[1];
rz(-2.0609444) q[2];
sx q[2];
rz(-1.4785826) q[2];
sx q[2];
rz(-0.25445081) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.091688823) q[1];
sx q[1];
rz(-0.65444512) q[1];
sx q[1];
rz(0.62033886) q[1];
rz(0.869107) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(1.3645736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(-0.7129933) q[2];
rz(1.0130079) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(-0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825496) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(-1.4404526) q[0];
rz(-3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48430303) q[0];
sx q[0];
rz(-0.7938677) q[0];
sx q[0];
rz(-2.8110709) q[0];
rz(1.1414358) q[2];
sx q[2];
rz(-2.6300207) q[2];
sx q[2];
rz(-0.19829743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.070378455) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(-0.94707625) q[1];
x q[2];
rz(2.3159536) q[3];
sx q[3];
rz(-2.2908483) q[3];
sx q[3];
rz(2.839098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99469441) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(-1.7374932) q[2];
rz(-1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-0.45853841) q[0];
rz(2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(-2.4564254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4165669) q[0];
sx q[0];
rz(-1.5013114) q[0];
sx q[0];
rz(2.2216703) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7148758) q[2];
sx q[2];
rz(-0.79877582) q[2];
sx q[2];
rz(-3.0720403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9381147) q[1];
sx q[1];
rz(-2.1739829) q[1];
sx q[1];
rz(-2.5342062) q[1];
x q[2];
rz(-3.1232883) q[3];
sx q[3];
rz(-0.98494512) q[3];
sx q[3];
rz(-0.5154807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(2.4122453) q[0];
rz(0.29306456) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(-1.1475295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3353951) q[0];
sx q[0];
rz(-1.1578387) q[0];
sx q[0];
rz(-1.6523244) q[0];
rz(-pi) q[1];
rz(-0.36979923) q[2];
sx q[2];
rz(-2.1714006) q[2];
sx q[2];
rz(-0.72593867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7123375) q[1];
sx q[1];
rz(-2.0550248) q[1];
sx q[1];
rz(1.5244563) q[1];
x q[2];
rz(2.6461584) q[3];
sx q[3];
rz(-1.2261651) q[3];
sx q[3];
rz(-2.7652094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78836936) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(0.64368147) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(2.3102982) q[0];
rz(-1.3759026) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9996223) q[0];
sx q[0];
rz(-0.49641434) q[0];
sx q[0];
rz(-1.2312141) q[0];
rz(-pi) q[1];
rz(2.5404846) q[2];
sx q[2];
rz(-1.483344) q[2];
sx q[2];
rz(0.89278883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41259137) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(2.1680135) q[1];
x q[2];
rz(0.83491915) q[3];
sx q[3];
rz(-0.69838006) q[3];
sx q[3];
rz(-1.2554807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9514256) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(-1.1249582) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(-1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(1.0983889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5441355) q[0];
sx q[0];
rz(-1.8034593) q[0];
sx q[0];
rz(-2.068589) q[0];
rz(-pi) q[1];
rz(0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(-2.2965477) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92022773) q[1];
sx q[1];
rz(-0.5792633) q[1];
sx q[1];
rz(-0.00097640493) q[1];
x q[2];
rz(-0.58795712) q[3];
sx q[3];
rz(-2.8391317) q[3];
sx q[3];
rz(-2.1289701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(1.8064921) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(-1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(0.46646068) q[0];
rz(0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-0.62896532) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53101978) q[0];
sx q[0];
rz(-1.5218698) q[0];
sx q[0];
rz(-1.945709) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92832698) q[2];
sx q[2];
rz(-1.3767585) q[2];
sx q[2];
rz(0.59552586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65044636) q[1];
sx q[1];
rz(-1.5728587) q[1];
sx q[1];
rz(-0.44209977) q[1];
rz(-pi) q[2];
rz(-0.38871308) q[3];
sx q[3];
rz(-0.75838381) q[3];
sx q[3];
rz(3.0384516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8978867) q[2];
sx q[2];
rz(-1.0906929) q[2];
sx q[2];
rz(-2.0824599) q[2];
rz(-0.6774261) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(-0.89390755) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28329904) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(-0.25390608) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(-1.614744) q[2];
sx q[2];
rz(-0.84779253) q[2];
sx q[2];
rz(0.070889125) q[2];
rz(-0.068594882) q[3];
sx q[3];
rz(-2.1566236) q[3];
sx q[3];
rz(-1.9038283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
