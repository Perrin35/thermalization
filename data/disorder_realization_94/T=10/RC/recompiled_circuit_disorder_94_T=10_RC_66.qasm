OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(1.6341524) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(2.5453321) q[1];
sx q[1];
rz(8.8095713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1546254) q[0];
sx q[0];
rz(-1.1095424) q[0];
sx q[0];
rz(2.2485562) q[0];
rz(-0.46703672) q[2];
sx q[2];
rz(-2.7454498) q[2];
sx q[2];
rz(3.0774088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8021009) q[1];
sx q[1];
rz(-1.1420982) q[1];
sx q[1];
rz(0.93356737) q[1];
rz(-pi) q[2];
rz(-0.30652133) q[3];
sx q[3];
rz(-1.74311) q[3];
sx q[3];
rz(-0.70787187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(-2.8033076) q[2];
rz(-1.4398549) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.671316) q[0];
sx q[0];
rz(-0.19965262) q[0];
sx q[0];
rz(1.5792363) q[0];
rz(-pi) q[1];
rz(-1.6127869) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(-3.0595879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46088947) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(-0.99980385) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1094692) q[3];
sx q[3];
rz(-0.20878775) q[3];
sx q[3];
rz(1.8148282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-2.8895203) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7092428) q[0];
sx q[0];
rz(-0.81626695) q[0];
sx q[0];
rz(2.4791251) q[0];
rz(2.4263072) q[2];
sx q[2];
rz(-2.2429357) q[2];
sx q[2];
rz(2.8298024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0059800681) q[1];
sx q[1];
rz(-2.1826535) q[1];
sx q[1];
rz(-0.23000418) q[1];
rz(-pi) q[2];
rz(0.39450816) q[3];
sx q[3];
rz(-1.2393701) q[3];
sx q[3];
rz(3.1098207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1217653) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(1.6148286) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-2.4345051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83311659) q[0];
sx q[0];
rz(-1.1261228) q[0];
sx q[0];
rz(2.0067257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6962887) q[2];
sx q[2];
rz(-1.1271994) q[2];
sx q[2];
rz(-2.0656297) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45338079) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(1.3761671) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9261758) q[3];
sx q[3];
rz(-1.8724752) q[3];
sx q[3];
rz(2.8043583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(-0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(1.7472349) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8481962) q[0];
sx q[0];
rz(-3.0612429) q[0];
sx q[0];
rz(-1.3991762) q[0];
rz(-2.1483634) q[2];
sx q[2];
rz(-1.5830056) q[2];
sx q[2];
rz(0.73405594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9533206) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(2.7642247) q[1];
x q[2];
rz(-0.8615287) q[3];
sx q[3];
rz(-1.9083175) q[3];
sx q[3];
rz(2.5634114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(1.3767892) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(0.29944637) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(2.9350231) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5305938) q[0];
sx q[0];
rz(-1.3337413) q[0];
sx q[0];
rz(0.80942746) q[0];
rz(-pi) q[1];
rz(-0.95025392) q[2];
sx q[2];
rz(-0.62180078) q[2];
sx q[2];
rz(-1.3602464) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5094604) q[1];
sx q[1];
rz(-1.3941947) q[1];
sx q[1];
rz(2.280974) q[1];
rz(-pi) q[2];
rz(-2.8387186) q[3];
sx q[3];
rz(-2.2040963) q[3];
sx q[3];
rz(-0.57297046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151826) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(-1.3279703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2845924) q[0];
sx q[0];
rz(-0.45495957) q[0];
sx q[0];
rz(-1.3083463) q[0];
rz(-pi) q[1];
rz(0.94857256) q[2];
sx q[2];
rz(-1.9554536) q[2];
sx q[2];
rz(2.4228061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2139637) q[1];
sx q[1];
rz(-0.8968401) q[1];
sx q[1];
rz(-2.9272635) q[1];
rz(0.51289576) q[3];
sx q[3];
rz(-1.9478056) q[3];
sx q[3];
rz(-0.60037724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.7810812) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-2.9597136) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(2.1906733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(-2.8911203) q[0];
x q[1];
rz(1.7958926) q[2];
sx q[2];
rz(-2.8689119) q[2];
sx q[2];
rz(-0.30579145) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2086522) q[1];
sx q[1];
rz(-1.3622074) q[1];
sx q[1];
rz(-1.637146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1675646) q[3];
sx q[3];
rz(-1.9208761) q[3];
sx q[3];
rz(1.10266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2237079) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(0.075909464) q[2];
rz(-2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.52371812) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10044554) q[0];
sx q[0];
rz(-0.99248306) q[0];
sx q[0];
rz(-2.5006177) q[0];
x q[1];
rz(-1.3525891) q[2];
sx q[2];
rz(-2.3404684) q[2];
sx q[2];
rz(1.6831236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0162504) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(2.7793105) q[1];
rz(-1.6660059) q[3];
sx q[3];
rz(-0.54064893) q[3];
sx q[3];
rz(1.8298061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24348564) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(1.8359258) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0730597) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(0.43138327) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87558526) q[2];
sx q[2];
rz(-2.0250118) q[2];
sx q[2];
rz(2.9241965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73951057) q[1];
sx q[1];
rz(-1.2716736) q[1];
sx q[1];
rz(-1.8226536) q[1];
rz(0.28966784) q[3];
sx q[3];
rz(-2.3710459) q[3];
sx q[3];
rz(-2.88248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-0.097461854) q[2];
sx q[2];
rz(-0.41257358) q[2];
sx q[2];
rz(-1.67795) q[2];
rz(-0.61836615) q[3];
sx q[3];
rz(-1.5636087) q[3];
sx q[3];
rz(-2.9611361) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];