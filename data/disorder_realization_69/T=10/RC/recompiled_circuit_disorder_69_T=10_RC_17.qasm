OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(5.2678582) q[0];
sx q[0];
rz(9.8922748) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8898534) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(2.2962909) q[0];
x q[1];
rz(0.99256398) q[2];
sx q[2];
rz(-0.34878584) q[2];
sx q[2];
rz(2.0860096) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80402741) q[1];
sx q[1];
rz(-2.0731888) q[1];
sx q[1];
rz(2.875209) q[1];
x q[2];
rz(2.8075571) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9900069) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(-2.4123689) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4366348) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(-1.0013642) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(2.2569879) q[0];
x q[1];
rz(0.87191138) q[2];
sx q[2];
rz(-1.9053835) q[2];
sx q[2];
rz(-1.3427693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3327649) q[1];
sx q[1];
rz(-1.9462703) q[1];
sx q[1];
rz(1.3198225) q[1];
x q[2];
rz(-2.8855521) q[3];
sx q[3];
rz(-2.5898993) q[3];
sx q[3];
rz(1.8191169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(1.8015507) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(-0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(-0.62477338) q[0];
rz(-0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(-2.952081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288085) q[0];
sx q[0];
rz(-1.3002987) q[0];
sx q[0];
rz(-0.74032797) q[0];
rz(-pi) q[1];
rz(1.465221) q[2];
sx q[2];
rz(-1.1974679) q[2];
sx q[2];
rz(2.2764652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1358007) q[1];
sx q[1];
rz(-1.6264919) q[1];
sx q[1];
rz(2.9018351) q[1];
rz(-pi) q[2];
rz(2.965851) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(2.9807846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(1.754388) q[2];
rz(-0.43131367) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4819734) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(1.1384456) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(1.0901573) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91256234) q[0];
sx q[0];
rz(-1.9711442) q[0];
sx q[0];
rz(-1.1926665) q[0];
rz(-pi) q[1];
rz(1.0341187) q[2];
sx q[2];
rz(-2.3148429) q[2];
sx q[2];
rz(0.85988753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5174597) q[1];
sx q[1];
rz(-2.1412666) q[1];
sx q[1];
rz(-2.0681417) q[1];
rz(-0.43153901) q[3];
sx q[3];
rz(-2.6772237) q[3];
sx q[3];
rz(2.4114256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1426992) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(0.85582716) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4677508) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(2.2348316) q[0];
rz(-pi) q[1];
rz(-1.1779551) q[2];
sx q[2];
rz(-2.422214) q[2];
sx q[2];
rz(2.9658085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7829195) q[1];
sx q[1];
rz(-1.1766608) q[1];
sx q[1];
rz(-0.47788099) q[1];
rz(-pi) q[2];
rz(1.7647469) q[3];
sx q[3];
rz(-1.0411106) q[3];
sx q[3];
rz(0.84996163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.871792) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(-3.0043547) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653172) q[0];
sx q[0];
rz(-0.99940171) q[0];
sx q[0];
rz(0.12302834) q[0];
x q[1];
rz(0.85645533) q[2];
sx q[2];
rz(-1.5894784) q[2];
sx q[2];
rz(-0.088353889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7676047) q[1];
sx q[1];
rz(-1.3629706) q[1];
sx q[1];
rz(0.61088224) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7298922) q[3];
sx q[3];
rz(-1.6414335) q[3];
sx q[3];
rz(0.27765805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(-2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(2.8939261) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(1.0621747) q[0];
rz(1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(1.4020845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37004334) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(2.7605961) q[0];
rz(-2.8962171) q[2];
sx q[2];
rz(-0.78005314) q[2];
sx q[2];
rz(1.9666372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2981616) q[1];
sx q[1];
rz(-1.2097675) q[1];
sx q[1];
rz(1.7346738) q[1];
x q[2];
rz(-1.7756895) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(0.3860592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.359882) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(2.6002398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4810825) q[0];
sx q[0];
rz(-2.5773002) q[0];
sx q[0];
rz(3.1206467) q[0];
x q[1];
rz(-1.2301684) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(-2.9609749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.036877) q[1];
sx q[1];
rz(-0.86537213) q[1];
sx q[1];
rz(-2.5887262) q[1];
x q[2];
rz(-1.5726611) q[3];
sx q[3];
rz(-1.1904753) q[3];
sx q[3];
rz(-2.5170381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(0.77735916) q[2];
rz(-0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.6962956) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9706668) q[0];
sx q[0];
rz(-2.2609841) q[0];
sx q[0];
rz(-0.70113457) q[0];
rz(-pi) q[1];
rz(-0.89500918) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(0.5627788) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6362793) q[1];
sx q[1];
rz(-1.5863824) q[1];
sx q[1];
rz(2.511335) q[1];
rz(-0.89740885) q[3];
sx q[3];
rz(-1.0703147) q[3];
sx q[3];
rz(2.9874779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-2.5342069) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(-0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(2.0954258) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(0.39224958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0961571) q[0];
sx q[0];
rz(-2.3491612) q[0];
sx q[0];
rz(-2.5655454) q[0];
x q[1];
rz(2.1456111) q[2];
sx q[2];
rz(-1.0258342) q[2];
sx q[2];
rz(-1.121322) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1019198) q[1];
sx q[1];
rz(-2.3452248) q[1];
sx q[1];
rz(2.2069195) q[1];
rz(1.5997821) q[3];
sx q[3];
rz(-1.3658804) q[3];
sx q[3];
rz(-2.8857185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(0.60047853) q[2];
rz(-2.0742119) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(-2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(2.375013) q[2];
sx q[2];
rz(-1.8324413) q[2];
sx q[2];
rz(1.2698297) q[2];
rz(2.2157833) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
