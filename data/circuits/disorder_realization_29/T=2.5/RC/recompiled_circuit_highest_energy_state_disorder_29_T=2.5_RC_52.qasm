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
rz(1.2443378) q[0];
sx q[0];
rz(2.3490348) q[0];
sx q[0];
rz(9.6777182) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(0.3869431) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3186424) q[0];
sx q[0];
rz(-1.8272663) q[0];
sx q[0];
rz(-1.2366813) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57780452) q[2];
sx q[2];
rz(-1.4981604) q[2];
sx q[2];
rz(-0.60608038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9305715) q[1];
sx q[1];
rz(-1.4912349) q[1];
sx q[1];
rz(-0.08427517) q[1];
rz(0.67220848) q[3];
sx q[3];
rz(-2.432193) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0953377) q[2];
sx q[2];
rz(-2.3825808) q[2];
sx q[2];
rz(-1.4026027) q[2];
rz(-1.4143573) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(-2.6426017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614748) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(0.71037355) q[0];
rz(-2.12517) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(0.76510915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0046526) q[0];
sx q[0];
rz(-2.1015133) q[0];
sx q[0];
rz(1.1999446) q[0];
rz(0.69074735) q[2];
sx q[2];
rz(-2.1155069) q[2];
sx q[2];
rz(-0.70554698) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73280947) q[1];
sx q[1];
rz(-1.8955232) q[1];
sx q[1];
rz(-0.88692437) q[1];
x q[2];
rz(1.160687) q[3];
sx q[3];
rz(-0.36204673) q[3];
sx q[3];
rz(2.1459879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0244828) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(-2.4369241) q[2];
rz(-0.17413983) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(2.7599938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.0891377) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(0.88743368) q[0];
rz(-1.2499836) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(-1.7807622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.783238) q[0];
sx q[0];
rz(-1.947442) q[0];
sx q[0];
rz(0.54623099) q[0];
rz(-pi) q[1];
rz(-1.7086231) q[2];
sx q[2];
rz(-1.4915474) q[2];
sx q[2];
rz(-0.85286959) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.683663) q[1];
sx q[1];
rz(-0.81240801) q[1];
sx q[1];
rz(3.0034566) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0988337) q[3];
sx q[3];
rz(-1.1819289) q[3];
sx q[3];
rz(-1.0145018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7736194) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(-1.6486637) q[2];
rz(0.44935539) q[3];
sx q[3];
rz(-1.8466693) q[3];
sx q[3];
rz(1.9213283) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0538977) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.4087403) q[0];
rz(-0.98211163) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(-1.7353479) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3262451) q[0];
sx q[0];
rz(-0.24213386) q[0];
sx q[0];
rz(1.0905488) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1020145) q[2];
sx q[2];
rz(-1.6235895) q[2];
sx q[2];
rz(0.18597183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2212053) q[1];
sx q[1];
rz(-1.3246598) q[1];
sx q[1];
rz(1.0369151) q[1];
x q[2];
rz(1.9191756) q[3];
sx q[3];
rz(-0.67871415) q[3];
sx q[3];
rz(1.410759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1534319) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(-0.15596685) q[2];
rz(2.4912452) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-0.11421886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0735737) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(2.9408348) q[0];
rz(-1.1772032) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(-1.1600561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61134385) q[0];
sx q[0];
rz(-2.8599116) q[0];
sx q[0];
rz(0.92667093) q[0];
x q[1];
rz(0.54520901) q[2];
sx q[2];
rz(-0.13449796) q[2];
sx q[2];
rz(0.70869499) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1694035) q[1];
sx q[1];
rz(-2.6693925) q[1];
sx q[1];
rz(-2.4257823) q[1];
rz(-2.1006159) q[3];
sx q[3];
rz(-1.3966832) q[3];
sx q[3];
rz(-1.8903738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7572299) q[2];
sx q[2];
rz(-0.94567662) q[2];
sx q[2];
rz(-0.45822701) q[2];
rz(-2.5107757) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(-2.6423776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.683627) q[0];
sx q[0];
rz(-0.85245913) q[0];
sx q[0];
rz(3.1347347) q[0];
rz(-2.9099756) q[1];
sx q[1];
rz(-1.7624785) q[1];
sx q[1];
rz(-1.0622271) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214913) q[0];
sx q[0];
rz(-1.694223) q[0];
sx q[0];
rz(2.8272259) q[0];
x q[1];
rz(0.49906667) q[2];
sx q[2];
rz(-1.7573414) q[2];
sx q[2];
rz(-1.3908433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1072709) q[1];
sx q[1];
rz(-0.91718972) q[1];
sx q[1];
rz(1.2493709) q[1];
rz(-pi) q[2];
rz(2.2799479) q[3];
sx q[3];
rz(-1.6966736) q[3];
sx q[3];
rz(-1.5365841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14083938) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(-1.2678649) q[2];
rz(0.86152348) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(-1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93290257) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(-0.21555899) q[0];
rz(0.49628273) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(-2.9821679) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8474583) q[0];
sx q[0];
rz(-0.93592866) q[0];
sx q[0];
rz(3.1383187) q[0];
rz(1.764545) q[2];
sx q[2];
rz(-1.1074578) q[2];
sx q[2];
rz(-0.090557052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9558952) q[1];
sx q[1];
rz(-1.61191) q[1];
sx q[1];
rz(-1.3354882) q[1];
rz(-pi) q[2];
rz(-1.0219021) q[3];
sx q[3];
rz(-1.9733784) q[3];
sx q[3];
rz(-2.2233913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40519199) q[2];
sx q[2];
rz(-1.1242194) q[2];
sx q[2];
rz(-0.51652017) q[2];
rz(1.2107595) q[3];
sx q[3];
rz(-1.6885933) q[3];
sx q[3];
rz(-1.3794544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.4844168) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(0.54022378) q[0];
rz(1.6172488) q[1];
sx q[1];
rz(-1.0018307) q[1];
sx q[1];
rz(-1.2670955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5634573) q[0];
sx q[0];
rz(-1.7892013) q[0];
sx q[0];
rz(-2.6294623) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77719633) q[2];
sx q[2];
rz(-1.0269477) q[2];
sx q[2];
rz(-1.0842619) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46173501) q[1];
sx q[1];
rz(-1.9377794) q[1];
sx q[1];
rz(1.7253897) q[1];
x q[2];
rz(-2.2121939) q[3];
sx q[3];
rz(-1.2817973) q[3];
sx q[3];
rz(0.81833392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8871062) q[2];
sx q[2];
rz(-1.1595414) q[2];
sx q[2];
rz(-0.17733388) q[2];
rz(-1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(-2.9593318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808481) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(-3.0408707) q[0];
rz(1.1489457) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(-1.9675868) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35648221) q[0];
sx q[0];
rz(-0.63725162) q[0];
sx q[0];
rz(-2.8204945) q[0];
rz(-0.76624932) q[2];
sx q[2];
rz(-1.6403997) q[2];
sx q[2];
rz(-1.0468259) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8619949) q[1];
sx q[1];
rz(-1.5902385) q[1];
sx q[1];
rz(-1.6085546) q[1];
x q[2];
rz(1.4467008) q[3];
sx q[3];
rz(-0.69487232) q[3];
sx q[3];
rz(1.2563039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.164087) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(0.40317765) q[2];
rz(-0.66344231) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(-3.1373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4106301) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(3.0626815) q[0];
rz(-2.2321189) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(-0.39631072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27061227) q[0];
sx q[0];
rz(-1.5617234) q[0];
sx q[0];
rz(-1.5551644) q[0];
rz(2.3478823) q[2];
sx q[2];
rz(-1.5657445) q[2];
sx q[2];
rz(1.2271301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8974298) q[1];
sx q[1];
rz(-1.3363991) q[1];
sx q[1];
rz(0.48903709) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5509866) q[3];
sx q[3];
rz(-2.1175623) q[3];
sx q[3];
rz(-1.8282229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4497455) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(2.3425102) q[2];
rz(-3.0633022) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(-0.6453132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.66452022) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(1.9480582) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(0.1706201) q[2];
sx q[2];
rz(-2.0224051) q[2];
sx q[2];
rz(-0.97013459) q[2];
rz(0.54936784) q[3];
sx q[3];
rz(-0.80437575) q[3];
sx q[3];
rz(1.5301216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
