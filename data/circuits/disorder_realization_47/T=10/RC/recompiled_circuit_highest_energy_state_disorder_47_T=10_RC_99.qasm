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
rz(-2.0464309) q[0];
sx q[0];
rz(-0.14552966) q[0];
sx q[0];
rz(-2.6382883) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(1.4454747) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5531122) q[0];
sx q[0];
rz(-2.7403826) q[0];
sx q[0];
rz(-3.0305639) q[0];
rz(-2.2022871) q[2];
sx q[2];
rz(-2.21647) q[2];
sx q[2];
rz(1.948057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9564285) q[1];
sx q[1];
rz(-1.899029) q[1];
sx q[1];
rz(0.64394711) q[1];
x q[2];
rz(2.8347391) q[3];
sx q[3];
rz(-0.05159353) q[3];
sx q[3];
rz(1.2449698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4989) q[2];
sx q[2];
rz(-1.2314726) q[2];
sx q[2];
rz(1.5042245) q[2];
rz(-2.4046992) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(-0.98114291) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40642834) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(-0.51097393) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-0.32646349) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8998979) q[0];
sx q[0];
rz(-2.7621671) q[0];
sx q[0];
rz(-2.3510758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9316177) q[2];
sx q[2];
rz(-2.2433503) q[2];
sx q[2];
rz(0.58174101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.105868) q[1];
sx q[1];
rz(-2.1505023) q[1];
sx q[1];
rz(-0.3649664) q[1];
rz(-pi) q[2];
rz(-2.191698) q[3];
sx q[3];
rz(-0.46484676) q[3];
sx q[3];
rz(-0.45528938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9713126) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(1.8453321) q[2];
rz(2.7235459) q[3];
sx q[3];
rz(-0.97801912) q[3];
sx q[3];
rz(-0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068785) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(-2.5991154) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(-0.03104041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2834969) q[0];
sx q[0];
rz(-1.265321) q[0];
sx q[0];
rz(-1.7683755) q[0];
x q[1];
rz(-2.228187) q[2];
sx q[2];
rz(-1.7511466) q[2];
sx q[2];
rz(-0.099152891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6185304) q[1];
sx q[1];
rz(-1.3253322) q[1];
sx q[1];
rz(-2.1851468) q[1];
x q[2];
rz(-2.5865745) q[3];
sx q[3];
rz(-1.5848098) q[3];
sx q[3];
rz(2.009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62640181) q[2];
sx q[2];
rz(-0.73651892) q[2];
sx q[2];
rz(-2.314563) q[2];
rz(-0.51860297) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14722918) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(-1.2116785) q[0];
rz(2.0934824) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.3406219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3163534) q[0];
sx q[0];
rz(-1.2636501) q[0];
sx q[0];
rz(0.65748837) q[0];
rz(-2.5848992) q[2];
sx q[2];
rz(-0.27678267) q[2];
sx q[2];
rz(-2.3900696) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1010434) q[1];
sx q[1];
rz(-1.2770997) q[1];
sx q[1];
rz(2.5459176) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9817438) q[3];
sx q[3];
rz(-1.2839497) q[3];
sx q[3];
rz(-2.8887859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2010605) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(2.000467) q[2];
rz(-1.5308258) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56213266) q[0];
sx q[0];
rz(-1.9714404) q[0];
sx q[0];
rz(-0.60254565) q[0];
rz(-2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(-0.7811195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3042219) q[0];
sx q[0];
rz(-1.2835555) q[0];
sx q[0];
rz(-0.17673136) q[0];
rz(0.76483043) q[2];
sx q[2];
rz(-2.0553608) q[2];
sx q[2];
rz(-2.7861905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8210817) q[1];
sx q[1];
rz(-2.2901504) q[1];
sx q[1];
rz(-3.0266552) q[1];
rz(-pi) q[2];
rz(-2.4158998) q[3];
sx q[3];
rz(-1.4376493) q[3];
sx q[3];
rz(-2.468688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96131229) q[2];
sx q[2];
rz(-2.0013516) q[2];
sx q[2];
rz(-2.9126634) q[2];
rz(-1.3198352) q[3];
sx q[3];
rz(-1.4819375) q[3];
sx q[3];
rz(2.2975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22885403) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(-1.9629021) q[0];
rz(1.9121869) q[1];
sx q[1];
rz(-2.7472159) q[1];
sx q[1];
rz(1.2961402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79149198) q[0];
sx q[0];
rz(-1.7069567) q[0];
sx q[0];
rz(-1.8541992) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77028124) q[2];
sx q[2];
rz(-2.2438223) q[2];
sx q[2];
rz(1.9865303) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4001105) q[1];
sx q[1];
rz(-2.0026836) q[1];
sx q[1];
rz(-2.6071965) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48719897) q[3];
sx q[3];
rz(-1.6604843) q[3];
sx q[3];
rz(-1.439217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65427762) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(-0.26228341) q[2];
rz(-0.30019635) q[3];
sx q[3];
rz(-1.9191977) q[3];
sx q[3];
rz(-0.20849553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4406141) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(-1.3025008) q[0];
rz(0.66849661) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(1.0350593) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9205017) q[0];
sx q[0];
rz(-1.971444) q[0];
sx q[0];
rz(-0.23200881) q[0];
rz(-1.9403753) q[2];
sx q[2];
rz(-1.3962922) q[2];
sx q[2];
rz(-2.7480887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8866402) q[1];
sx q[1];
rz(-1.6463841) q[1];
sx q[1];
rz(-1.2264538) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74372282) q[3];
sx q[3];
rz(-1.111314) q[3];
sx q[3];
rz(2.2032025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11789007) q[2];
sx q[2];
rz(-1.002545) q[2];
sx q[2];
rz(-1.4855851) q[2];
rz(0.28247908) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(-2.7070467) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78541237) q[0];
sx q[0];
rz(-2.0996576) q[0];
sx q[0];
rz(2.8670512) q[0];
rz(1.9369594) q[1];
sx q[1];
rz(-2.5614673) q[1];
sx q[1];
rz(2.5801632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1875136) q[0];
sx q[0];
rz(-2.2642379) q[0];
sx q[0];
rz(2.2827143) q[0];
rz(-2.4653685) q[2];
sx q[2];
rz(-1.3315787) q[2];
sx q[2];
rz(-0.70331159) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9778487) q[1];
sx q[1];
rz(-1.9789961) q[1];
sx q[1];
rz(0.331649) q[1];
x q[2];
rz(2.2116978) q[3];
sx q[3];
rz(-0.42503438) q[3];
sx q[3];
rz(2.2512521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5795035) q[2];
sx q[2];
rz(-2.0359437) q[2];
sx q[2];
rz(-1.4037464) q[2];
rz(1.3459407) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(-2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3812934) q[0];
sx q[0];
rz(-0.60307044) q[0];
sx q[0];
rz(2.2316933) q[0];
rz(1.9445885) q[1];
sx q[1];
rz(-1.0531813) q[1];
sx q[1];
rz(-2.2534456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.600888) q[0];
sx q[0];
rz(-1.6645169) q[0];
sx q[0];
rz(-1.4613495) q[0];
rz(-pi) q[1];
rz(-2.4817564) q[2];
sx q[2];
rz(-1.6183231) q[2];
sx q[2];
rz(-2.8328676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.55938086) q[1];
sx q[1];
rz(-2.191698) q[1];
sx q[1];
rz(2.0526753) q[1];
rz(-0.03753438) q[3];
sx q[3];
rz(-1.2753295) q[3];
sx q[3];
rz(0.42643828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1278648) q[2];
sx q[2];
rz(-1.4415386) q[2];
sx q[2];
rz(-2.0590651) q[2];
rz(2.9108289) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(2.6575991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.3103264) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(-2.5600774) q[0];
rz(-0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(1.8448578) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0727507) q[0];
sx q[0];
rz(-1.5525991) q[0];
sx q[0];
rz(-0.010404603) q[0];
x q[1];
rz(0.026785568) q[2];
sx q[2];
rz(-2.7973632) q[2];
sx q[2];
rz(-2.8452498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2263068) q[1];
sx q[1];
rz(-1.1476333) q[1];
sx q[1];
rz(2.8598815) q[1];
rz(-pi) q[2];
rz(-2.4695685) q[3];
sx q[3];
rz(-1.250953) q[3];
sx q[3];
rz(-0.32957382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7221308) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(2.6922743) q[2];
rz(0.32014534) q[3];
sx q[3];
rz(-1.3195427) q[3];
sx q[3];
rz(1.8951353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768059) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(1.7569348) q[1];
sx q[1];
rz(-2.90381) q[1];
sx q[1];
rz(-0.79868383) q[1];
rz(1.2133219) q[2];
sx q[2];
rz(-2.4926179) q[2];
sx q[2];
rz(0.31821584) q[2];
rz(-2.7560227) q[3];
sx q[3];
rz(-2.0115934) q[3];
sx q[3];
rz(2.6075076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
