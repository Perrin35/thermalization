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
rz(-1.545067) q[1];
sx q[1];
rz(-2.5453321) q[1];
sx q[1];
rz(2.526386) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98696729) q[0];
sx q[0];
rz(-2.0320503) q[0];
sx q[0];
rz(0.89303645) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7569321) q[2];
sx q[2];
rz(-1.2190483) q[2];
sx q[2];
rz(-0.56439161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7429744) q[1];
sx q[1];
rz(-0.75098872) q[1];
sx q[1];
rz(-0.91575925) q[1];
rz(-pi) q[2];
rz(-0.30652133) q[3];
sx q[3];
rz(-1.3984826) q[3];
sx q[3];
rz(0.70787187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(2.8033076) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(3.1112444) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.5240086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.679927) q[0];
sx q[0];
rz(-1.3711509) q[0];
sx q[0];
rz(3.1398849) q[0];
rz(-1.6127869) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(0.08200478) q[2];
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
rz(-1.6398805) q[1];
sx q[1];
rz(0.99980385) q[1];
rz(-pi) q[2];
rz(1.7582943) q[3];
sx q[3];
rz(-1.6631931) q[3];
sx q[3];
rz(0.20860162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(-1.1478109) q[2];
rz(-1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(-2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(-0.31578627) q[0];
rz(0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-2.8895203) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35226563) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(0.69791039) q[0];
rz(-2.4263072) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(-0.31179024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7604916) q[1];
sx q[1];
rz(-0.64844202) q[1];
sx q[1];
rz(1.2566503) q[1];
rz(-0.39450816) q[3];
sx q[3];
rz(-1.2393701) q[3];
sx q[3];
rz(0.031771914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(-0.034051731) q[2];
rz(-3.1241336) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(2.0461369) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83311659) q[0];
sx q[0];
rz(-2.0154698) q[0];
sx q[0];
rz(-2.0067257) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(0.69869631) q[1];
sx q[1];
rz(-2.2203608) q[1];
sx q[1];
rz(2.991308) q[1];
rz(-pi) q[2];
rz(0.96890038) q[3];
sx q[3];
rz(-0.36877353) q[3];
sx q[3];
rz(-0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-2.6814931) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(-0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8481962) q[0];
sx q[0];
rz(-0.080349803) q[0];
sx q[0];
rz(-1.3991762) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1483634) q[2];
sx q[2];
rz(-1.5585871) q[2];
sx q[2];
rz(-2.4075367) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6756145) q[1];
sx q[1];
rz(-1.8120159) q[1];
sx q[1];
rz(2.4720008) q[1];
x q[2];
rz(1.0765431) q[3];
sx q[3];
rz(-0.77270618) q[3];
sx q[3];
rz(0.62437526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(1.6453751) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(0.20656955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9402007) q[0];
sx q[0];
rz(-0.79027806) q[0];
sx q[0];
rz(-1.9076365) q[0];
x q[1];
rz(2.1913387) q[2];
sx q[2];
rz(-0.62180078) q[2];
sx q[2];
rz(-1.3602464) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0530015) q[1];
sx q[1];
rz(-0.87391657) q[1];
sx q[1];
rz(-0.23115302) q[1];
x q[2];
rz(0.30287403) q[3];
sx q[3];
rz(-2.2040963) q[3];
sx q[3];
rz(-0.57297046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.841659) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(-2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2151826) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(1.8136224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0921811) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(-1.1294424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46220025) q[2];
sx q[2];
rz(-2.1415347) q[2];
sx q[2];
rz(0.58909033) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3634062) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(-0.88552514) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6789356) q[3];
sx q[3];
rz(-0.62641615) q[3];
sx q[3];
rz(2.7501447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(-0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(2.1906733) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2598341) q[0];
sx q[0];
rz(-1.8108597) q[0];
sx q[0];
rz(-1.2756707) q[0];
x q[1];
rz(1.7958926) q[2];
sx q[2];
rz(-2.8689119) q[2];
sx q[2];
rz(-0.30579145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2086522) q[1];
sx q[1];
rz(-1.7793852) q[1];
sx q[1];
rz(-1.5044466) q[1];
rz(-pi) q[2];
rz(-2.7637134) q[3];
sx q[3];
rz(-1.9482908) q[3];
sx q[3];
rz(-2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-0.075909464) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(2.7245522) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-2.4818647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83803672) q[0];
sx q[0];
rz(-2.3065901) q[0];
sx q[0];
rz(-2.3124218) q[0];
rz(-pi) q[1];
rz(-0.21978901) q[2];
sx q[2];
rz(-0.79384365) q[2];
sx q[2];
rz(1.1500051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8068741) q[1];
sx q[1];
rz(-1.8611307) q[1];
sx q[1];
rz(-0.90805407) q[1];
rz(-pi) q[2];
rz(-0.057007313) q[3];
sx q[3];
rz(-2.1087286) q[3];
sx q[3];
rz(-1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(2.1155604) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-1.0326577) q[3];
sx q[3];
rz(0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.3056668) q[0];
rz(-1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0337692) q[0];
sx q[0];
rz(-0.53350893) q[0];
sx q[0];
rz(-0.67722042) q[0];
rz(-pi) q[1];
rz(-0.91949384) q[2];
sx q[2];
rz(-2.3323625) q[2];
sx q[2];
rz(0.86916718) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1199011) q[1];
sx q[1];
rz(-0.38858116) q[1];
sx q[1];
rz(2.4619224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3922937) q[3];
sx q[3];
rz(-1.7710847) q[3];
sx q[3];
rz(2.0405958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(0.36469665) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37968996) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(1.5079386) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(0.097461854) q[2];
sx q[2];
rz(-2.7290191) q[2];
sx q[2];
rz(1.4636427) q[2];
rz(-0.012398331) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];