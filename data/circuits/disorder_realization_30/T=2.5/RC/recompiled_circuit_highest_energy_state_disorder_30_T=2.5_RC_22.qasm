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
rz(0.27424115) q[0];
sx q[0];
rz(-2.7354136) q[0];
sx q[0];
rz(0.22845536) q[0];
rz(-2.7574381) q[1];
sx q[1];
rz(-1.6698807) q[1];
sx q[1];
rz(2.6928112) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160546) q[0];
sx q[0];
rz(-1.9384697) q[0];
sx q[0];
rz(1.7203164) q[0];
rz(-1.2972791) q[2];
sx q[2];
rz(-0.43028545) q[2];
sx q[2];
rz(-1.6543433) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27707057) q[1];
sx q[1];
rz(-1.6344353) q[1];
sx q[1];
rz(-3.1176104) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9308546) q[3];
sx q[3];
rz(-0.02861261) q[3];
sx q[3];
rz(0.60625854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3936798) q[2];
sx q[2];
rz(-2.9069275) q[2];
sx q[2];
rz(2.5133384) q[2];
rz(1.0660394) q[3];
sx q[3];
rz(-2.1431036) q[3];
sx q[3];
rz(-1.1294686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5265441) q[0];
sx q[0];
rz(-2.5874309) q[0];
sx q[0];
rz(-0.91914415) q[0];
rz(2.9650086) q[1];
sx q[1];
rz(-1.44839) q[1];
sx q[1];
rz(2.6939556) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2528322) q[0];
sx q[0];
rz(-1.5764142) q[0];
sx q[0];
rz(3.1340213) q[0];
rz(-pi) q[1];
rz(-1.8488919) q[2];
sx q[2];
rz(-1.1881042) q[2];
sx q[2];
rz(-1.1130321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7940951) q[1];
sx q[1];
rz(-1.7298051) q[1];
sx q[1];
rz(-0.68557941) q[1];
rz(-0.80609808) q[3];
sx q[3];
rz(-2.4662123) q[3];
sx q[3];
rz(2.4247629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24719079) q[2];
sx q[2];
rz(-1.9220507) q[2];
sx q[2];
rz(0.88376898) q[2];
rz(-0.73924685) q[3];
sx q[3];
rz(-1.913162) q[3];
sx q[3];
rz(1.8993186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3833375) q[0];
sx q[0];
rz(-2.1012335) q[0];
sx q[0];
rz(0.0058767949) q[0];
rz(-2.8351496) q[1];
sx q[1];
rz(-2.0649316) q[1];
sx q[1];
rz(1.2642911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976083) q[0];
sx q[0];
rz(-2.6559444) q[0];
sx q[0];
rz(2.4729054) q[0];
rz(-pi) q[1];
rz(-1.4219079) q[2];
sx q[2];
rz(-1.604594) q[2];
sx q[2];
rz(1.6848909) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2325461) q[1];
sx q[1];
rz(-2.1276444) q[1];
sx q[1];
rz(1.805748) q[1];
rz(-pi) q[2];
rz(0.98250947) q[3];
sx q[3];
rz(-2.4949565) q[3];
sx q[3];
rz(-1.8573295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82187247) q[2];
sx q[2];
rz(-1.7629273) q[2];
sx q[2];
rz(0.6591447) q[2];
rz(1.871073) q[3];
sx q[3];
rz(-2.8489385) q[3];
sx q[3];
rz(1.6856153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6419411) q[0];
sx q[0];
rz(-0.15436509) q[0];
sx q[0];
rz(-2.2041712) q[0];
rz(-1.5501267) q[1];
sx q[1];
rz(-2.4731686) q[1];
sx q[1];
rz(-0.74916565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031627071) q[0];
sx q[0];
rz(-2.3195889) q[0];
sx q[0];
rz(-1.1618932) q[0];
rz(-pi) q[1];
rz(-0.62994434) q[2];
sx q[2];
rz(-1.0039181) q[2];
sx q[2];
rz(-0.34467372) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3377825) q[1];
sx q[1];
rz(-0.69548464) q[1];
sx q[1];
rz(-1.2616322) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9177298) q[3];
sx q[3];
rz(-2.3626257) q[3];
sx q[3];
rz(1.1457486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59430993) q[2];
sx q[2];
rz(-2.3197033) q[2];
sx q[2];
rz(-1.8853356) q[2];
rz(-2.3516583) q[3];
sx q[3];
rz(-0.94090763) q[3];
sx q[3];
rz(-1.1676211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.054852) q[0];
sx q[0];
rz(-2.5331443) q[0];
sx q[0];
rz(-0.21859455) q[0];
rz(0.78544468) q[1];
sx q[1];
rz(-2.0365069) q[1];
sx q[1];
rz(1.8484263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620813) q[0];
sx q[0];
rz(-2.0326699) q[0];
sx q[0];
rz(-2.1681468) q[0];
x q[1];
rz(-1.2720455) q[2];
sx q[2];
rz(-2.8674539) q[2];
sx q[2];
rz(-1.1098286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.433328) q[1];
sx q[1];
rz(-0.86769968) q[1];
sx q[1];
rz(-0.82872699) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1035091) q[3];
sx q[3];
rz(-1.786219) q[3];
sx q[3];
rz(-1.076339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.834632) q[2];
sx q[2];
rz(-1.0647651) q[2];
sx q[2];
rz(1.628423) q[2];
rz(-2.4026134) q[3];
sx q[3];
rz(-1.8578953) q[3];
sx q[3];
rz(-1.6966604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4479248) q[0];
sx q[0];
rz(-1.4883) q[0];
sx q[0];
rz(-0.36201763) q[0];
rz(-1.1297049) q[1];
sx q[1];
rz(-1.874322) q[1];
sx q[1];
rz(0.76090181) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99859649) q[0];
sx q[0];
rz(-0.79299295) q[0];
sx q[0];
rz(-2.8727813) q[0];
x q[1];
rz(0.82693066) q[2];
sx q[2];
rz(-1.2614577) q[2];
sx q[2];
rz(0.29238809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43266549) q[1];
sx q[1];
rz(-0.58307465) q[1];
sx q[1];
rz(-1.1770171) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6057556) q[3];
sx q[3];
rz(-0.69147516) q[3];
sx q[3];
rz(2.7590318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1165983) q[2];
sx q[2];
rz(-2.8572539) q[2];
sx q[2];
rz(-1.1762168) q[2];
rz(-2.1182012) q[3];
sx q[3];
rz(-2.0659645) q[3];
sx q[3];
rz(0.37337506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1171595) q[0];
sx q[0];
rz(-2.7130782) q[0];
sx q[0];
rz(1.3065216) q[0];
rz(-0.10082968) q[1];
sx q[1];
rz(-1.3918624) q[1];
sx q[1];
rz(-0.94591013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3358766) q[0];
sx q[0];
rz(-1.44084) q[0];
sx q[0];
rz(0.30673271) q[0];
rz(2.628075) q[2];
sx q[2];
rz(-2.0198722) q[2];
sx q[2];
rz(0.30563909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57179831) q[1];
sx q[1];
rz(-1.6299953) q[1];
sx q[1];
rz(-3.0059391) q[1];
x q[2];
rz(-1.9705371) q[3];
sx q[3];
rz(-2.3869166) q[3];
sx q[3];
rz(-0.012059742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0239608) q[2];
sx q[2];
rz(-1.5317711) q[2];
sx q[2];
rz(-2.7839933) q[2];
rz(-0.050497342) q[3];
sx q[3];
rz(-2.7183618) q[3];
sx q[3];
rz(-0.38664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015775) q[0];
sx q[0];
rz(-1.3192663) q[0];
sx q[0];
rz(2.5460119) q[0];
rz(-1.2254084) q[1];
sx q[1];
rz(-1.775368) q[1];
sx q[1];
rz(2.7092686) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5670938) q[0];
sx q[0];
rz(-1.4050975) q[0];
sx q[0];
rz(-2.6461547) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7244699) q[2];
sx q[2];
rz(-1.6343079) q[2];
sx q[2];
rz(2.2907298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22031517) q[1];
sx q[1];
rz(-2.5743458) q[1];
sx q[1];
rz(-3.1392643) q[1];
x q[2];
rz(-0.8882167) q[3];
sx q[3];
rz(-1.2176664) q[3];
sx q[3];
rz(-2.5251561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7633742) q[2];
sx q[2];
rz(-1.8209527) q[2];
sx q[2];
rz(2.1136368) q[2];
rz(-1.6592735) q[3];
sx q[3];
rz(-1.6927398) q[3];
sx q[3];
rz(-2.3894374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1886275) q[0];
sx q[0];
rz(-1.231622) q[0];
sx q[0];
rz(0.46284437) q[0];
rz(-1.2340744) q[1];
sx q[1];
rz(-1.9201098) q[1];
sx q[1];
rz(5*pi/8) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1723768) q[0];
sx q[0];
rz(-1.5961884) q[0];
sx q[0];
rz(1.7459016) q[0];
rz(-0.47204702) q[2];
sx q[2];
rz(-0.85523048) q[2];
sx q[2];
rz(-0.89821076) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2888711) q[1];
sx q[1];
rz(-0.72119547) q[1];
sx q[1];
rz(1.3836963) q[1];
x q[2];
rz(-2.7998802) q[3];
sx q[3];
rz(-1.4572506) q[3];
sx q[3];
rz(-2.7734203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50596333) q[2];
sx q[2];
rz(-1.7545981) q[2];
sx q[2];
rz(0.76816922) q[2];
rz(0.65028894) q[3];
sx q[3];
rz(-2.8322329) q[3];
sx q[3];
rz(-0.55997854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3637417) q[0];
sx q[0];
rz(-0.44054458) q[0];
sx q[0];
rz(-1.995723) q[0];
rz(2.5762985) q[1];
sx q[1];
rz(-0.71563131) q[1];
sx q[1];
rz(0.40564767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0019819) q[0];
sx q[0];
rz(-2.2898556) q[0];
sx q[0];
rz(-2.0400892) q[0];
x q[1];
rz(2.8641607) q[2];
sx q[2];
rz(-1.9346969) q[2];
sx q[2];
rz(0.24574797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5989117) q[1];
sx q[1];
rz(-0.57334954) q[1];
sx q[1];
rz(2.663055) q[1];
rz(-pi) q[2];
rz(1.6784366) q[3];
sx q[3];
rz(-1.1336796) q[3];
sx q[3];
rz(0.48278615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.046935) q[2];
sx q[2];
rz(-1.2311225) q[2];
sx q[2];
rz(2.6222353) q[2];
rz(-2.6804067) q[3];
sx q[3];
rz(-2.0824194) q[3];
sx q[3];
rz(-0.69380277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8984579) q[0];
sx q[0];
rz(-1.5512137) q[0];
sx q[0];
rz(2.2358538) q[0];
rz(-1.657919) q[1];
sx q[1];
rz(-1.1653733) q[1];
sx q[1];
rz(1.9001874) q[1];
rz(-2.9979669) q[2];
sx q[2];
rz(-0.96114071) q[2];
sx q[2];
rz(-3.0428934) q[2];
rz(-2.0679071) q[3];
sx q[3];
rz(-1.1563672) q[3];
sx q[3];
rz(1.4856114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
