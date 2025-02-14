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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(-2.6068249) q[0];
rz(0.98400247) q[1];
sx q[1];
rz(0.50552955) q[1];
sx q[1];
rz(11.304392) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0875435) q[0];
sx q[0];
rz(-2.4719878) q[0];
sx q[0];
rz(2.7101507) q[0];
rz(2.6017351) q[2];
sx q[2];
rz(-1.015402) q[2];
sx q[2];
rz(0.29796165) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83916503) q[1];
sx q[1];
rz(-1.622597) q[1];
sx q[1];
rz(2.5940345) q[1];
rz(-pi) q[2];
rz(-2.6747236) q[3];
sx q[3];
rz(-0.43064603) q[3];
sx q[3];
rz(-2.1809585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24903211) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(-1.6585635) q[2];
rz(-2.856971) q[3];
sx q[3];
rz(-0.62323815) q[3];
sx q[3];
rz(1.0641789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677218) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(0.04100767) q[0];
rz(-1.2076591) q[1];
sx q[1];
rz(-0.227808) q[1];
sx q[1];
rz(1.6809195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631993) q[0];
sx q[0];
rz(-1.5413056) q[0];
sx q[0];
rz(-0.14568744) q[0];
rz(0.76717214) q[2];
sx q[2];
rz(-1.9843272) q[2];
sx q[2];
rz(0.77595475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8833958) q[1];
sx q[1];
rz(-0.25942311) q[1];
sx q[1];
rz(1.7218462) q[1];
x q[2];
rz(2.9462141) q[3];
sx q[3];
rz(-2.2535705) q[3];
sx q[3];
rz(-1.9548137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4760806) q[2];
sx q[2];
rz(-1.4613232) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(0.0811854) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7471033) q[0];
sx q[0];
rz(-0.012270027) q[0];
sx q[0];
rz(-2.5315206) q[0];
rz(1.3129781) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(-2.5909766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9351106) q[0];
sx q[0];
rz(-1.2121823) q[0];
sx q[0];
rz(0.7048084) q[0];
x q[1];
rz(2.8189893) q[2];
sx q[2];
rz(-1.3122953) q[2];
sx q[2];
rz(-2.5532364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0455835) q[1];
sx q[1];
rz(-2.7852045) q[1];
sx q[1];
rz(2.136904) q[1];
rz(-pi) q[2];
rz(-1.7254225) q[3];
sx q[3];
rz(-2.3904413) q[3];
sx q[3];
rz(0.42645833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6637471) q[2];
sx q[2];
rz(-2.9691594) q[2];
sx q[2];
rz(1.0663859) q[2];
rz(-1.1586698) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(0.042044736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8313507) q[0];
sx q[0];
rz(-2.5284335) q[0];
sx q[0];
rz(2.2245275) q[0];
rz(0.65545583) q[1];
sx q[1];
rz(-0.90924811) q[1];
sx q[1];
rz(0.38958946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87990078) q[0];
sx q[0];
rz(-1.3845056) q[0];
sx q[0];
rz(-3.0584719) q[0];
x q[1];
rz(0.88574804) q[2];
sx q[2];
rz(-2.0221524) q[2];
sx q[2];
rz(1.7471755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0181737) q[1];
sx q[1];
rz(-1.5999854) q[1];
sx q[1];
rz(1.1694714) q[1];
rz(2.2993777) q[3];
sx q[3];
rz(-0.73916736) q[3];
sx q[3];
rz(0.30876866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3098844) q[2];
sx q[2];
rz(-0.94253057) q[2];
sx q[2];
rz(1.7892828) q[2];
rz(-0.74812198) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(1.5266533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-2.5717773) q[0];
sx q[0];
rz(-1.8001528) q[0];
rz(-2.1603284) q[1];
sx q[1];
rz(-1.860268) q[1];
sx q[1];
rz(0.70995465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34559743) q[0];
sx q[0];
rz(-1.5899204) q[0];
sx q[0];
rz(-2.8157793) q[0];
rz(-2.1469028) q[2];
sx q[2];
rz(-1.369595) q[2];
sx q[2];
rz(0.76801571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1011184) q[1];
sx q[1];
rz(-1.3583364) q[1];
sx q[1];
rz(0.17417769) q[1];
rz(1.665713) q[3];
sx q[3];
rz(-2.1664985) q[3];
sx q[3];
rz(-1.9323521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95825163) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(-2.23488) q[2];
rz(-0.18236154) q[3];
sx q[3];
rz(-1.7014818) q[3];
sx q[3];
rz(-0.85842925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(1.3458503) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(2.9398651) q[0];
rz(-1.7851104) q[1];
sx q[1];
rz(-0.67142612) q[1];
sx q[1];
rz(1.9904402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7089823) q[0];
sx q[0];
rz(-1.994619) q[0];
sx q[0];
rz(0.2291542) q[0];
rz(-pi) q[1];
rz(0.62080748) q[2];
sx q[2];
rz(-2.4352031) q[2];
sx q[2];
rz(-2.8548129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5804739) q[1];
sx q[1];
rz(-2.4137133) q[1];
sx q[1];
rz(0.62904686) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39252354) q[3];
sx q[3];
rz(-2.5965434) q[3];
sx q[3];
rz(2.1529993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1543033) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(-0.34745535) q[2];
rz(-0.82184982) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(-1.6177026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3407985) q[0];
sx q[0];
rz(-0.64685416) q[0];
sx q[0];
rz(1.1232173) q[0];
rz(-0.42876354) q[1];
sx q[1];
rz(-0.88974297) q[1];
sx q[1];
rz(-2.9926328) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4377047) q[0];
sx q[0];
rz(-2.0876679) q[0];
sx q[0];
rz(0.025512841) q[0];
x q[1];
rz(-2.5507418) q[2];
sx q[2];
rz(-2.789641) q[2];
sx q[2];
rz(0.9443493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80854446) q[1];
sx q[1];
rz(-0.77189779) q[1];
sx q[1];
rz(0.21265642) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9765719) q[3];
sx q[3];
rz(-1.5386484) q[3];
sx q[3];
rz(-3.0253077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1296156) q[2];
sx q[2];
rz(-2.7836383) q[2];
sx q[2];
rz(2.820106) q[2];
rz(1.1164411) q[3];
sx q[3];
rz(-0.8588841) q[3];
sx q[3];
rz(1.727625) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9880992) q[0];
sx q[0];
rz(-1.3445925) q[0];
sx q[0];
rz(-3.1186812) q[0];
rz(-1.5921536) q[1];
sx q[1];
rz(-2.0982845) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1003636) q[0];
sx q[0];
rz(-2.1080906) q[0];
sx q[0];
rz(2.8365718) q[0];
x q[1];
rz(-1.171807) q[2];
sx q[2];
rz(-0.771847) q[2];
sx q[2];
rz(0.10153025) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.072619) q[1];
sx q[1];
rz(-2.6311462) q[1];
sx q[1];
rz(2.4859395) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.096209196) q[3];
sx q[3];
rz(-2.2819977) q[3];
sx q[3];
rz(2.8195153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93078485) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(-2.1442294) q[2];
rz(-1.8831683) q[3];
sx q[3];
rz(-1.9505898) q[3];
sx q[3];
rz(-1.2303801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1579943) q[0];
sx q[0];
rz(-2.4857434) q[0];
sx q[0];
rz(-0.47437814) q[0];
rz(2.5899218) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(-0.65753585) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0484587) q[0];
sx q[0];
rz(-3.0952929) q[0];
sx q[0];
rz(0.53673021) q[0];
rz(-pi) q[1];
rz(1.9507061) q[2];
sx q[2];
rz(-0.43336855) q[2];
sx q[2];
rz(-0.72335183) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8847981) q[1];
sx q[1];
rz(-1.9391283) q[1];
sx q[1];
rz(0.017466768) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10513427) q[3];
sx q[3];
rz(-1.4270367) q[3];
sx q[3];
rz(1.2105699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9115596) q[2];
sx q[2];
rz(-2.7342789) q[2];
sx q[2];
rz(1.3113021) q[2];
rz(2.4274872) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1353726) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(-2.5332992) q[0];
rz(0.81781203) q[1];
sx q[1];
rz(-1.9655971) q[1];
sx q[1];
rz(-0.83293319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8084106) q[0];
sx q[0];
rz(-1.5191374) q[0];
sx q[0];
rz(-0.42148659) q[0];
rz(-pi) q[1];
rz(-1.8274177) q[2];
sx q[2];
rz(-1.875281) q[2];
sx q[2];
rz(1.2198795) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2257833) q[1];
sx q[1];
rz(-1.438101) q[1];
sx q[1];
rz(-2.5578294) q[1];
x q[2];
rz(2.5485431) q[3];
sx q[3];
rz(-1.5447953) q[3];
sx q[3];
rz(3.0700695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1828764) q[2];
sx q[2];
rz(-0.35280886) q[2];
sx q[2];
rz(-0.58615169) q[2];
rz(2.2385249) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027973) q[0];
sx q[0];
rz(-2.0113404) q[0];
sx q[0];
rz(2.197862) q[0];
rz(-1.0974274) q[1];
sx q[1];
rz(-2.8891017) q[1];
sx q[1];
rz(-0.15695922) q[1];
rz(0.13931724) q[2];
sx q[2];
rz(-1.9092154) q[2];
sx q[2];
rz(-0.42001324) q[2];
rz(2.1395352) q[3];
sx q[3];
rz(-0.78227038) q[3];
sx q[3];
rz(-3.0113358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
