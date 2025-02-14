OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.032441) q[0];
sx q[0];
rz(-1.1530131) q[0];
sx q[0];
rz(3.0670526) q[0];
rz(-1.5300765) q[1];
sx q[1];
rz(-0.059377436) q[1];
sx q[1];
rz(-2.6174954) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3976879) q[0];
sx q[0];
rz(-0.70262229) q[0];
sx q[0];
rz(2.1545059) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8139548) q[2];
sx q[2];
rz(-1.8983525) q[2];
sx q[2];
rz(-1.0741247) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14088881) q[1];
sx q[1];
rz(-2.2743164) q[1];
sx q[1];
rz(0.8013273) q[1];
rz(-pi) q[2];
rz(-3.1007441) q[3];
sx q[3];
rz(-1.4788879) q[3];
sx q[3];
rz(2.330454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3268299) q[2];
sx q[2];
rz(-2.6430898) q[2];
sx q[2];
rz(2.2461183) q[2];
rz(0.26252663) q[3];
sx q[3];
rz(-0.34590507) q[3];
sx q[3];
rz(-3.053022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.217591) q[0];
sx q[0];
rz(-1.1955248) q[0];
sx q[0];
rz(-3.0254645) q[0];
rz(2.5092292) q[1];
sx q[1];
rz(-2.875681) q[1];
sx q[1];
rz(1.6557453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53609849) q[0];
sx q[0];
rz(-1.5935531) q[0];
sx q[0];
rz(-3.0748585) q[0];
rz(2.112442) q[2];
sx q[2];
rz(-1.7746762) q[2];
sx q[2];
rz(-1.1612877) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0318438) q[1];
sx q[1];
rz(-2.9753775) q[1];
sx q[1];
rz(2.4480661) q[1];
x q[2];
rz(0.13074517) q[3];
sx q[3];
rz(-2.2669889) q[3];
sx q[3];
rz(-1.8896904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31132013) q[2];
sx q[2];
rz(-2.1612284) q[2];
sx q[2];
rz(2.2118528) q[2];
rz(-0.43944198) q[3];
sx q[3];
rz(-1.6012871) q[3];
sx q[3];
rz(0.75696993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0759401) q[0];
sx q[0];
rz(-2.3891698) q[0];
sx q[0];
rz(-2.2546076) q[0];
rz(1.8807962) q[1];
sx q[1];
rz(-1.1075243) q[1];
sx q[1];
rz(-1.4874123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5607779) q[0];
sx q[0];
rz(-0.11753035) q[0];
sx q[0];
rz(-1.5752026) q[0];
x q[1];
rz(1.4753129) q[2];
sx q[2];
rz(-1.1794568) q[2];
sx q[2];
rz(3.0961159) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9076242) q[1];
sx q[1];
rz(-0.32999906) q[1];
sx q[1];
rz(3.0872295) q[1];
rz(-pi) q[2];
rz(0.58165197) q[3];
sx q[3];
rz(-1.8466788) q[3];
sx q[3];
rz(1.1946071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4999353) q[2];
sx q[2];
rz(-0.96330088) q[2];
sx q[2];
rz(2.7063043) q[2];
rz(0.38517243) q[3];
sx q[3];
rz(-1.2110854) q[3];
sx q[3];
rz(-2.1729573) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4061072) q[0];
sx q[0];
rz(-0.084097363) q[0];
sx q[0];
rz(-3.0870068) q[0];
rz(-1.5054043) q[1];
sx q[1];
rz(-1.2137698) q[1];
sx q[1];
rz(-0.41753599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62359257) q[0];
sx q[0];
rz(-1.6290278) q[0];
sx q[0];
rz(-2.7043155) q[0];
rz(2.5092441) q[2];
sx q[2];
rz(-1.8397619) q[2];
sx q[2];
rz(2.7346345) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8167953) q[1];
sx q[1];
rz(-2.9831605) q[1];
sx q[1];
rz(-3.0976803) q[1];
rz(-pi) q[2];
rz(-1.8154816) q[3];
sx q[3];
rz(-0.26991329) q[3];
sx q[3];
rz(-0.0072561023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85228449) q[2];
sx q[2];
rz(-2.4384629) q[2];
sx q[2];
rz(-0.0011778041) q[2];
rz(-1.7581455) q[3];
sx q[3];
rz(-2.8774084) q[3];
sx q[3];
rz(1.0435102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99172878) q[0];
sx q[0];
rz(-2.9809451) q[0];
sx q[0];
rz(0.37586656) q[0];
rz(-2.1916892) q[1];
sx q[1];
rz(-1.4774731) q[1];
sx q[1];
rz(2.4163767) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8729018) q[0];
sx q[0];
rz(-0.16624545) q[0];
sx q[0];
rz(2.3670271) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64230772) q[2];
sx q[2];
rz(-2.1644582) q[2];
sx q[2];
rz(0.35333179) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9053697) q[1];
sx q[1];
rz(-3.0439483) q[1];
sx q[1];
rz(0.059448551) q[1];
x q[2];
rz(-1.3106724) q[3];
sx q[3];
rz(-2.2571466) q[3];
sx q[3];
rz(2.7174866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3788562) q[2];
sx q[2];
rz(-1.4521658) q[2];
sx q[2];
rz(2.6689996) q[2];
rz(3.0224814) q[3];
sx q[3];
rz(-2.5178858) q[3];
sx q[3];
rz(-0.81010336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8238207) q[0];
sx q[0];
rz(-2.9404984) q[0];
sx q[0];
rz(-0.066545181) q[0];
rz(-2.9464974) q[1];
sx q[1];
rz(-1.0761484) q[1];
sx q[1];
rz(-1.1544352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5148564) q[0];
sx q[0];
rz(-1.3382698) q[0];
sx q[0];
rz(2.0621081) q[0];
x q[1];
rz(-0.17627861) q[2];
sx q[2];
rz(-0.46277324) q[2];
sx q[2];
rz(2.0001992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6007541) q[1];
sx q[1];
rz(-2.7620865) q[1];
sx q[1];
rz(1.2787766) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1514417) q[3];
sx q[3];
rz(-1.5816474) q[3];
sx q[3];
rz(-0.059338245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2164312) q[2];
sx q[2];
rz(-1.602828) q[2];
sx q[2];
rz(0.72539854) q[2];
rz(1.409449) q[3];
sx q[3];
rz(-1.1812482) q[3];
sx q[3];
rz(0.60666549) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0856536) q[0];
sx q[0];
rz(-0.5760718) q[0];
sx q[0];
rz(2.2357909) q[0];
rz(-0.55465758) q[1];
sx q[1];
rz(-0.83811086) q[1];
sx q[1];
rz(-0.14403266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0222826) q[0];
sx q[0];
rz(-0.57252175) q[0];
sx q[0];
rz(1.0767816) q[0];
x q[1];
rz(1.0582744) q[2];
sx q[2];
rz(-1.6450168) q[2];
sx q[2];
rz(-2.3942301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1293948) q[1];
sx q[1];
rz(-1.8499814) q[1];
sx q[1];
rz(-1.1711981) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1062076) q[3];
sx q[3];
rz(-2.852172) q[3];
sx q[3];
rz(1.1076224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4732699) q[2];
sx q[2];
rz(-0.39445764) q[2];
sx q[2];
rz(0.99297601) q[2];
rz(0.18443491) q[3];
sx q[3];
rz(-0.95576972) q[3];
sx q[3];
rz(-0.9930281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1297146) q[0];
sx q[0];
rz(-2.5373902) q[0];
sx q[0];
rz(0.39464828) q[0];
rz(2.6036085) q[1];
sx q[1];
rz(-0.6901651) q[1];
sx q[1];
rz(1.1916377) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25090835) q[0];
sx q[0];
rz(-0.60870752) q[0];
sx q[0];
rz(-2.6842791) q[0];
rz(1.3778119) q[2];
sx q[2];
rz(-1.9375083) q[2];
sx q[2];
rz(0.16736469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.887693) q[1];
sx q[1];
rz(-2.116345) q[1];
sx q[1];
rz(-2.6274526) q[1];
rz(-2.0268782) q[3];
sx q[3];
rz(-2.1563765) q[3];
sx q[3];
rz(0.28156137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.789232) q[2];
sx q[2];
rz(-1.5287986) q[2];
sx q[2];
rz(2.1996876) q[2];
rz(2.3801129) q[3];
sx q[3];
rz(-1.1250291) q[3];
sx q[3];
rz(-3.1010845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2725459) q[0];
sx q[0];
rz(-1.7750374) q[0];
sx q[0];
rz(-1.0719365) q[0];
rz(-0.6140703) q[1];
sx q[1];
rz(-0.89659381) q[1];
sx q[1];
rz(2.695172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9337685) q[0];
sx q[0];
rz(-2.38777) q[0];
sx q[0];
rz(-1.1541744) q[0];
rz(-pi) q[1];
rz(-0.62320407) q[2];
sx q[2];
rz(-0.7936306) q[2];
sx q[2];
rz(0.49937427) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.804396) q[1];
sx q[1];
rz(-2.2257518) q[1];
sx q[1];
rz(1.7095997) q[1];
x q[2];
rz(-1.4037973) q[3];
sx q[3];
rz(-2.0701968) q[3];
sx q[3];
rz(1.4016408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3093695) q[2];
sx q[2];
rz(-0.95755219) q[2];
sx q[2];
rz(1.2584125) q[2];
rz(1.2331412) q[3];
sx q[3];
rz(-0.62658739) q[3];
sx q[3];
rz(1.3174177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719139) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(1.7183787) q[0];
rz(2.8990959) q[1];
sx q[1];
rz(-2.6624694) q[1];
sx q[1];
rz(1.2538145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72287382) q[0];
sx q[0];
rz(-2.6986045) q[0];
sx q[0];
rz(2.0110058) q[0];
rz(-pi) q[1];
rz(3.0331633) q[2];
sx q[2];
rz(-3.0745818) q[2];
sx q[2];
rz(2.5052414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.032728) q[1];
sx q[1];
rz(-2.0833587) q[1];
sx q[1];
rz(-2.7973919) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1997211) q[3];
sx q[3];
rz(-1.365445) q[3];
sx q[3];
rz(-1.273231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8157876) q[2];
sx q[2];
rz(-2.3141404) q[2];
sx q[2];
rz(0.023836689) q[2];
rz(1.2787974) q[3];
sx q[3];
rz(-2.8102504) q[3];
sx q[3];
rz(-2.3648025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1635638) q[0];
sx q[0];
rz(-1.4373056) q[0];
sx q[0];
rz(1.9594255) q[0];
rz(-2.4143746) q[1];
sx q[1];
rz(-1.4677508) q[1];
sx q[1];
rz(-0.8263091) q[1];
rz(-1.4197311) q[2];
sx q[2];
rz(-1.551566) q[2];
sx q[2];
rz(-2.8788064) q[2];
rz(1.942523) q[3];
sx q[3];
rz(-3.0109497) q[3];
sx q[3];
rz(2.0554832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
