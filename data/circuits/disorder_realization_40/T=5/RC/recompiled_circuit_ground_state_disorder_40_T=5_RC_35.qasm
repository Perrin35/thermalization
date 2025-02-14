OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52387828) q[0];
sx q[0];
rz(3.7359306) q[0];
sx q[0];
rz(9.115968) q[0];
rz(0.29769695) q[1];
sx q[1];
rz(-1.8841789) q[1];
sx q[1];
rz(0.6119734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966105) q[0];
sx q[0];
rz(-2.5685446) q[0];
sx q[0];
rz(-1.7368421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8582552) q[2];
sx q[2];
rz(-1.4170215) q[2];
sx q[2];
rz(2.2361627) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0667324) q[1];
sx q[1];
rz(-1.43044) q[1];
sx q[1];
rz(1.9060777) q[1];
x q[2];
rz(-0.060617491) q[3];
sx q[3];
rz(-2.4481886) q[3];
sx q[3];
rz(-1.7934679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1674126) q[2];
sx q[2];
rz(-1.9638502) q[2];
sx q[2];
rz(-0.56498945) q[2];
rz(-2.6307093) q[3];
sx q[3];
rz(-0.23879819) q[3];
sx q[3];
rz(1.6142982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329472) q[0];
sx q[0];
rz(-2.8605509) q[0];
sx q[0];
rz(-0.99579048) q[0];
rz(0.58798724) q[1];
sx q[1];
rz(-2.7016787) q[1];
sx q[1];
rz(1.3194552) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4026227) q[0];
sx q[0];
rz(-1.3330632) q[0];
sx q[0];
rz(2.0031702) q[0];
rz(0.075602268) q[2];
sx q[2];
rz(-1.2443064) q[2];
sx q[2];
rz(1.9734427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.022835613) q[1];
sx q[1];
rz(-2.3110227) q[1];
sx q[1];
rz(0.28121314) q[1];
rz(-pi) q[2];
rz(-0.48254099) q[3];
sx q[3];
rz(-1.4630894) q[3];
sx q[3];
rz(1.5272702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0413282) q[2];
sx q[2];
rz(-1.0039971) q[2];
sx q[2];
rz(-3.1190994) q[2];
rz(3.0371173) q[3];
sx q[3];
rz(-1.5375117) q[3];
sx q[3];
rz(2.2711066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81258881) q[0];
sx q[0];
rz(-1.0306232) q[0];
sx q[0];
rz(-1.6382244) q[0];
rz(-0.080987856) q[1];
sx q[1];
rz(-2.4826725) q[1];
sx q[1];
rz(1.9714877) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2572944) q[0];
sx q[0];
rz(-2.5405875) q[0];
sx q[0];
rz(-2.8544507) q[0];
rz(-2.8798772) q[2];
sx q[2];
rz(-2.4804575) q[2];
sx q[2];
rz(2.9950855) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2418211) q[1];
sx q[1];
rz(-1.5994497) q[1];
sx q[1];
rz(-2.0709527) q[1];
x q[2];
rz(1.0634138) q[3];
sx q[3];
rz(-1.3510002) q[3];
sx q[3];
rz(-0.45377094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4795503) q[2];
sx q[2];
rz(-0.68632555) q[2];
sx q[2];
rz(-0.043206841) q[2];
rz(0.19733812) q[3];
sx q[3];
rz(-1.03136) q[3];
sx q[3];
rz(-2.7092547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2546805) q[0];
sx q[0];
rz(-1.0858902) q[0];
sx q[0];
rz(0.32549724) q[0];
rz(-2.2858641) q[1];
sx q[1];
rz(-2.7266462) q[1];
sx q[1];
rz(-2.0554481) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153692) q[0];
sx q[0];
rz(-1.6994711) q[0];
sx q[0];
rz(2.6475304) q[0];
rz(2.3065673) q[2];
sx q[2];
rz(-0.97423282) q[2];
sx q[2];
rz(1.2460097) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1126744) q[1];
sx q[1];
rz(-2.5956545) q[1];
sx q[1];
rz(1.8914936) q[1];
x q[2];
rz(0.086768199) q[3];
sx q[3];
rz(-1.7132856) q[3];
sx q[3];
rz(0.87228197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93422741) q[2];
sx q[2];
rz(-2.1659329) q[2];
sx q[2];
rz(-2.9577241) q[2];
rz(1.81987) q[3];
sx q[3];
rz(-0.70122856) q[3];
sx q[3];
rz(3.0085861) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4399399) q[0];
sx q[0];
rz(-2.8337605) q[0];
sx q[0];
rz(-2.2046748) q[0];
rz(0.91494715) q[1];
sx q[1];
rz(-2.2888384) q[1];
sx q[1];
rz(1.459704) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383612) q[0];
sx q[0];
rz(-2.438926) q[0];
sx q[0];
rz(-1.42704) q[0];
rz(-pi) q[1];
rz(-2.1906846) q[2];
sx q[2];
rz(-1.0703329) q[2];
sx q[2];
rz(-1.2318947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.416977) q[1];
sx q[1];
rz(-1.4733581) q[1];
sx q[1];
rz(-2.3259054) q[1];
x q[2];
rz(2.567148) q[3];
sx q[3];
rz(-0.46009053) q[3];
sx q[3];
rz(2.121832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0437643) q[2];
sx q[2];
rz(-1.8057258) q[2];
sx q[2];
rz(-2.9975927) q[2];
rz(2.0740267) q[3];
sx q[3];
rz(-0.34393603) q[3];
sx q[3];
rz(0.6914773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095793) q[0];
sx q[0];
rz(-2.3342275) q[0];
sx q[0];
rz(0.76835865) q[0];
rz(-2.2840624) q[1];
sx q[1];
rz(-1.0464959) q[1];
sx q[1];
rz(-0.55364496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2087649) q[0];
sx q[0];
rz(-2.0220482) q[0];
sx q[0];
rz(0.28214595) q[0];
rz(-1.3411936) q[2];
sx q[2];
rz(-0.78289778) q[2];
sx q[2];
rz(1.7097434) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2170797) q[1];
sx q[1];
rz(-1.9929152) q[1];
sx q[1];
rz(-0.98084992) q[1];
rz(-pi) q[2];
rz(0.72154273) q[3];
sx q[3];
rz(-2.0399722) q[3];
sx q[3];
rz(2.814722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.386261) q[2];
sx q[2];
rz(-2.7109881) q[2];
sx q[2];
rz(-2.6957896) q[2];
rz(-1.2465994) q[3];
sx q[3];
rz(-1.3695025) q[3];
sx q[3];
rz(-2.2915452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069139473) q[0];
sx q[0];
rz(-0.9599762) q[0];
sx q[0];
rz(0.14007105) q[0];
rz(-2.2135997) q[1];
sx q[1];
rz(-1.3368139) q[1];
sx q[1];
rz(1.450052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04270794) q[0];
sx q[0];
rz(-1.6804114) q[0];
sx q[0];
rz(1.4102458) q[0];
rz(-1.9533402) q[2];
sx q[2];
rz(-1.3556644) q[2];
sx q[2];
rz(2.6184788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19181025) q[1];
sx q[1];
rz(-0.78652421) q[1];
sx q[1];
rz(-2.5108348) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4387162) q[3];
sx q[3];
rz(-0.23628391) q[3];
sx q[3];
rz(0.98881665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.037584) q[2];
sx q[2];
rz(-0.17825492) q[2];
sx q[2];
rz(0.083871052) q[2];
rz(0.9907848) q[3];
sx q[3];
rz(-1.1889941) q[3];
sx q[3];
rz(-1.8535463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8952119) q[0];
sx q[0];
rz(-1.3061433) q[0];
sx q[0];
rz(0.35686785) q[0];
rz(1.4419979) q[1];
sx q[1];
rz(-2.5394963) q[1];
sx q[1];
rz(-0.009036202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9401902) q[0];
sx q[0];
rz(-1.3497258) q[0];
sx q[0];
rz(-1.3371435) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6918654) q[2];
sx q[2];
rz(-2.3894261) q[2];
sx q[2];
rz(-2.8945685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7459834) q[1];
sx q[1];
rz(-1.8425178) q[1];
sx q[1];
rz(-3.0467669) q[1];
rz(-pi) q[2];
rz(-1.610454) q[3];
sx q[3];
rz(-1.1278099) q[3];
sx q[3];
rz(-2.7970527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.25702) q[2];
sx q[2];
rz(-1.7947861) q[2];
sx q[2];
rz(-0.70551562) q[2];
rz(0.26933119) q[3];
sx q[3];
rz(-1.0764542) q[3];
sx q[3];
rz(-2.1721325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89328289) q[0];
sx q[0];
rz(-0.80535424) q[0];
sx q[0];
rz(-2.4321108) q[0];
rz(-2.4995038) q[1];
sx q[1];
rz(-2.700192) q[1];
sx q[1];
rz(-0.4253687) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.294153) q[0];
sx q[0];
rz(-2.3340539) q[0];
sx q[0];
rz(-2.5320413) q[0];
x q[1];
rz(0.9832731) q[2];
sx q[2];
rz(-0.66907489) q[2];
sx q[2];
rz(-2.4226505) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3166362) q[1];
sx q[1];
rz(-1.6610089) q[1];
sx q[1];
rz(-2.6949203) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6900979) q[3];
sx q[3];
rz(-0.62314829) q[3];
sx q[3];
rz(-0.54766207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2702668) q[2];
sx q[2];
rz(-2.9623803) q[2];
sx q[2];
rz(0.47879177) q[2];
rz(0.35017961) q[3];
sx q[3];
rz(-1.9637354) q[3];
sx q[3];
rz(0.42696264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-0.66788524) q[0];
sx q[0];
rz(-2.844664) q[0];
sx q[0];
rz(-0.83734751) q[0];
rz(0.26652023) q[1];
sx q[1];
rz(-1.3198677) q[1];
sx q[1];
rz(1.7841608) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11649179) q[0];
sx q[0];
rz(-3.0550346) q[0];
sx q[0];
rz(-1.5309912) q[0];
rz(-1.6896115) q[2];
sx q[2];
rz(-1.6336361) q[2];
sx q[2];
rz(-2.6667537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19465654) q[1];
sx q[1];
rz(-1.5583073) q[1];
sx q[1];
rz(-0.65232531) q[1];
x q[2];
rz(2.8164163) q[3];
sx q[3];
rz(-2.0132228) q[3];
sx q[3];
rz(0.57462245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9799161) q[2];
sx q[2];
rz(-2.6075173) q[2];
sx q[2];
rz(-2.7034289) q[2];
rz(-2.2257889) q[3];
sx q[3];
rz(-1.5235498) q[3];
sx q[3];
rz(0.80317909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5220779) q[0];
sx q[0];
rz(-1.2418455) q[0];
sx q[0];
rz(-1.0157304) q[0];
rz(0.10454128) q[1];
sx q[1];
rz(-1.8395945) q[1];
sx q[1];
rz(1.3855388) q[1];
rz(-1.2404493) q[2];
sx q[2];
rz(-1.5120244) q[2];
sx q[2];
rz(2.8373847) q[2];
rz(-0.85827479) q[3];
sx q[3];
rz(-1.0714053) q[3];
sx q[3];
rz(-1.3965931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
