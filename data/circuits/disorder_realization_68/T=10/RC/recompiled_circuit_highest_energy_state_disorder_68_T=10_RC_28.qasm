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
rz(0.23213586) q[0];
sx q[0];
rz(-0.27096662) q[0];
sx q[0];
rz(1.969307) q[0];
rz(1.5341893) q[1];
sx q[1];
rz(-1.6521896) q[1];
sx q[1];
rz(0.54816562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0294103) q[0];
sx q[0];
rz(-1.8444841) q[0];
sx q[0];
rz(0.11464439) q[0];
rz(-0.51082533) q[2];
sx q[2];
rz(-1.0721954) q[2];
sx q[2];
rz(-3.0650578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89828789) q[1];
sx q[1];
rz(-2.5410748) q[1];
sx q[1];
rz(-0.55796786) q[1];
x q[2];
rz(-0.062730363) q[3];
sx q[3];
rz(-2.0196805) q[3];
sx q[3];
rz(0.99454885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24070172) q[2];
sx q[2];
rz(-1.2071995) q[2];
sx q[2];
rz(1.594225) q[2];
rz(-2.2089925) q[3];
sx q[3];
rz(-0.91351944) q[3];
sx q[3];
rz(2.3368321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43657434) q[0];
sx q[0];
rz(-0.60282928) q[0];
sx q[0];
rz(-2.9233209) q[0];
rz(-1.5087347) q[1];
sx q[1];
rz(-2.2684596) q[1];
sx q[1];
rz(1.1223209) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22161417) q[0];
sx q[0];
rz(-1.99619) q[0];
sx q[0];
rz(-0.90817389) q[0];
rz(-2.400945) q[2];
sx q[2];
rz(-1.7193931) q[2];
sx q[2];
rz(-2.1987178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6186445) q[1];
sx q[1];
rz(-2.1174333) q[1];
sx q[1];
rz(1.0376105) q[1];
x q[2];
rz(0.25457766) q[3];
sx q[3];
rz(-1.8768594) q[3];
sx q[3];
rz(1.7382415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0108769) q[2];
sx q[2];
rz(-1.9515832) q[2];
sx q[2];
rz(0.38401628) q[2];
rz(0.555641) q[3];
sx q[3];
rz(-0.58755392) q[3];
sx q[3];
rz(-2.2491992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71690774) q[0];
sx q[0];
rz(-2.6838344) q[0];
sx q[0];
rz(-2.6440788) q[0];
rz(2.6566907) q[1];
sx q[1];
rz(-1.5919911) q[1];
sx q[1];
rz(-0.44152322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45831523) q[0];
sx q[0];
rz(-0.9181058) q[0];
sx q[0];
rz(-2.9291332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1782568) q[2];
sx q[2];
rz(-2.6148893) q[2];
sx q[2];
rz(-1.8369499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5680629) q[1];
sx q[1];
rz(-0.90319777) q[1];
sx q[1];
rz(2.1993162) q[1];
x q[2];
rz(2.2763292) q[3];
sx q[3];
rz(-2.7332716) q[3];
sx q[3];
rz(2.7029519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32450822) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(0.95477742) q[2];
rz(3.1214516) q[3];
sx q[3];
rz(-2.3432799) q[3];
sx q[3];
rz(-2.8153343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7666053) q[0];
sx q[0];
rz(-0.52244455) q[0];
sx q[0];
rz(-3.1368384) q[0];
rz(0.44113723) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(1.7316679) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49214464) q[0];
sx q[0];
rz(-0.36055598) q[0];
sx q[0];
rz(1.1900405) q[0];
rz(-pi) q[1];
rz(-0.27556117) q[2];
sx q[2];
rz(-2.2803734) q[2];
sx q[2];
rz(1.0312652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9121352) q[1];
sx q[1];
rz(-1.7963103) q[1];
sx q[1];
rz(1.0646125) q[1];
rz(-pi) q[2];
rz(0.027214931) q[3];
sx q[3];
rz(-1.3821954) q[3];
sx q[3];
rz(2.1299429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8498174) q[2];
sx q[2];
rz(-1.3128277) q[2];
sx q[2];
rz(2.4449091) q[2];
rz(-1.5289395) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(-2.4618885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1729537) q[0];
sx q[0];
rz(-0.30690673) q[0];
sx q[0];
rz(2.2343743) q[0];
rz(3.0061159) q[1];
sx q[1];
rz(-1.2716581) q[1];
sx q[1];
rz(2.4438593) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013779) q[0];
sx q[0];
rz(-1.4411949) q[0];
sx q[0];
rz(-1.1762397) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1393797) q[2];
sx q[2];
rz(-1.7656823) q[2];
sx q[2];
rz(0.28838487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5144044) q[1];
sx q[1];
rz(-1.3980734) q[1];
sx q[1];
rz(0.33337969) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1118513) q[3];
sx q[3];
rz(-1.5248767) q[3];
sx q[3];
rz(2.8950402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1873143) q[2];
sx q[2];
rz(-2.4772187) q[2];
sx q[2];
rz(0.3453671) q[2];
rz(-2.6464388) q[3];
sx q[3];
rz(-0.59585714) q[3];
sx q[3];
rz(-3.0143747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0879764) q[0];
sx q[0];
rz(-1.477366) q[0];
sx q[0];
rz(1.9867058) q[0];
rz(-2.3467973) q[1];
sx q[1];
rz(-1.2966917) q[1];
sx q[1];
rz(2.6577139) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.799494) q[0];
sx q[0];
rz(-2.3364107) q[0];
sx q[0];
rz(2.1495887) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15911289) q[2];
sx q[2];
rz(-2.0622353) q[2];
sx q[2];
rz(3.0895777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3920073) q[1];
sx q[1];
rz(-1.3806496) q[1];
sx q[1];
rz(-2.0431678) q[1];
rz(2.5577066) q[3];
sx q[3];
rz(-0.66388452) q[3];
sx q[3];
rz(1.3308457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38857073) q[2];
sx q[2];
rz(-1.9066255) q[2];
sx q[2];
rz(-1.1080326) q[2];
rz(-1.612251) q[3];
sx q[3];
rz(-0.11950167) q[3];
sx q[3];
rz(2.748238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18688467) q[0];
sx q[0];
rz(-2.055838) q[0];
sx q[0];
rz(0.67436522) q[0];
rz(0.87633324) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(0.032329917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040213765) q[0];
sx q[0];
rz(-1.4353936) q[0];
sx q[0];
rz(3.0187294) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59835513) q[2];
sx q[2];
rz(-2.8718635) q[2];
sx q[2];
rz(0.21141569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8994042) q[1];
sx q[1];
rz(-2.5061878) q[1];
sx q[1];
rz(-2.5492378) q[1];
rz(0.19200872) q[3];
sx q[3];
rz(-2.1550094) q[3];
sx q[3];
rz(-1.2835128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9548367) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(2.1046861) q[2];
rz(-1.2284651) q[3];
sx q[3];
rz(-0.90932536) q[3];
sx q[3];
rz(2.9508446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4036338) q[0];
sx q[0];
rz(-0.027733138) q[0];
sx q[0];
rz(-2.8926358) q[0];
rz(-0.20949334) q[1];
sx q[1];
rz(-1.341235) q[1];
sx q[1];
rz(2.478821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4512206) q[0];
sx q[0];
rz(-1.5589542) q[0];
sx q[0];
rz(1.580834) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0515819) q[2];
sx q[2];
rz(-0.86746806) q[2];
sx q[2];
rz(2.472773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1472305) q[1];
sx q[1];
rz(-0.77795815) q[1];
sx q[1];
rz(1.2176355) q[1];
rz(0.39318496) q[3];
sx q[3];
rz(-0.69614702) q[3];
sx q[3];
rz(-2.7659211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2304307) q[2];
sx q[2];
rz(-1.570236) q[2];
sx q[2];
rz(-0.39361185) q[2];
rz(-2.7502934) q[3];
sx q[3];
rz(-0.53520447) q[3];
sx q[3];
rz(-0.75214255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-1.1289718) q[0];
sx q[0];
rz(-0.18652815) q[0];
sx q[0];
rz(0.57998002) q[0];
rz(-2.773556) q[1];
sx q[1];
rz(-2.3263704) q[1];
sx q[1];
rz(-0.11485242) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2942685) q[0];
sx q[0];
rz(-2.2912187) q[0];
sx q[0];
rz(-0.17981932) q[0];
rz(-pi) q[1];
rz(0.1505974) q[2];
sx q[2];
rz(-2.0080749) q[2];
sx q[2];
rz(0.33971805) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99913524) q[1];
sx q[1];
rz(-1.8977381) q[1];
sx q[1];
rz(-1.8349951) q[1];
x q[2];
rz(0.87646342) q[3];
sx q[3];
rz(-1.8871347) q[3];
sx q[3];
rz(-2.9522459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7241235) q[2];
sx q[2];
rz(-2.8026411) q[2];
sx q[2];
rz(2.4453956) q[2];
rz(-1.1160858) q[3];
sx q[3];
rz(-2.3510272) q[3];
sx q[3];
rz(2.4000786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29257947) q[0];
sx q[0];
rz(-2.1722023) q[0];
sx q[0];
rz(0.27266362) q[0];
rz(-0.69372454) q[1];
sx q[1];
rz(-1.1169746) q[1];
sx q[1];
rz(-2.5781217) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54105896) q[0];
sx q[0];
rz(-0.1002914) q[0];
sx q[0];
rz(1.1931399) q[0];
rz(-pi) q[1];
rz(-1.4239975) q[2];
sx q[2];
rz(-1.8052535) q[2];
sx q[2];
rz(-0.80658972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.882) q[1];
sx q[1];
rz(-1.2597076) q[1];
sx q[1];
rz(0.085373665) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2056687) q[3];
sx q[3];
rz(-2.1800123) q[3];
sx q[3];
rz(0.09906957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1571265) q[2];
sx q[2];
rz(-1.0447341) q[2];
sx q[2];
rz(-2.5926479) q[2];
rz(2.1460311) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(0.63001776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41848771) q[0];
sx q[0];
rz(-2.1514819) q[0];
sx q[0];
rz(2.6958382) q[0];
rz(2.2743629) q[1];
sx q[1];
rz(-1.1031311) q[1];
sx q[1];
rz(1.1581609) q[1];
rz(-2.6643999) q[2];
sx q[2];
rz(-1.7765316) q[2];
sx q[2];
rz(2.3892567) q[2];
rz(-2.5936962) q[3];
sx q[3];
rz(-1.6954281) q[3];
sx q[3];
rz(2.2466698) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
