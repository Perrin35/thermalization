OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57269078) q[0];
sx q[0];
rz(-2.1532018) q[0];
sx q[0];
rz(0.44901499) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(0.13154498) q[1];
sx q[1];
rz(8.2934525) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97206632) q[0];
sx q[0];
rz(-2.457805) q[0];
sx q[0];
rz(1.0975361) q[0];
x q[1];
rz(-2.4822818) q[2];
sx q[2];
rz(-0.74987312) q[2];
sx q[2];
rz(1.7323959) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1234963) q[1];
sx q[1];
rz(-2.5275748) q[1];
sx q[1];
rz(-1.8753002) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9692124) q[3];
sx q[3];
rz(-2.1076492) q[3];
sx q[3];
rz(2.7011724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.955287) q[2];
sx q[2];
rz(-0.39262843) q[2];
sx q[2];
rz(3.0379831) q[2];
rz(0.3668395) q[3];
sx q[3];
rz(-1.6678526) q[3];
sx q[3];
rz(1.3175255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1845301) q[0];
sx q[0];
rz(-2.6657031) q[0];
sx q[0];
rz(2.6153508) q[0];
rz(-1.2435675) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(-0.19613656) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945902) q[0];
sx q[0];
rz(-2.6795912) q[0];
sx q[0];
rz(-1.4205971) q[0];
rz(2.4320658) q[2];
sx q[2];
rz(-1.4440184) q[2];
sx q[2];
rz(-1.4122054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96207075) q[1];
sx q[1];
rz(-1.4286433) q[1];
sx q[1];
rz(-2.0217871) q[1];
x q[2];
rz(0.72153458) q[3];
sx q[3];
rz(-0.01577687) q[3];
sx q[3];
rz(1.1392913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5940932) q[2];
sx q[2];
rz(-2.8671691) q[2];
sx q[2];
rz(-2.7653232) q[2];
rz(-0.88092342) q[3];
sx q[3];
rz(-1.8062402) q[3];
sx q[3];
rz(2.3295565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6905717) q[0];
sx q[0];
rz(-0.32830992) q[0];
sx q[0];
rz(2.1300533) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-1.9790383) q[1];
sx q[1];
rz(0.54723251) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1062463) q[0];
sx q[0];
rz(-1.5900471) q[0];
sx q[0];
rz(-0.20142844) q[0];
rz(-pi) q[1];
rz(0.34788068) q[2];
sx q[2];
rz(-0.82150412) q[2];
sx q[2];
rz(-0.88319983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32729748) q[1];
sx q[1];
rz(-1.3298099) q[1];
sx q[1];
rz(-2.2404033) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7748479) q[3];
sx q[3];
rz(-2.3746852) q[3];
sx q[3];
rz(2.6605061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4554567) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(-1.6395462) q[2];
rz(-0.57602588) q[3];
sx q[3];
rz(-1.6582158) q[3];
sx q[3];
rz(0.36147931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107002) q[0];
sx q[0];
rz(-0.80131131) q[0];
sx q[0];
rz(-2.8019688) q[0];
rz(-2.6009808) q[1];
sx q[1];
rz(-0.70053354) q[1];
sx q[1];
rz(2.2115754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41893533) q[0];
sx q[0];
rz(-1.6412853) q[0];
sx q[0];
rz(1.427729) q[0];
rz(2.5581237) q[2];
sx q[2];
rz(-1.8975048) q[2];
sx q[2];
rz(-2.938156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2423289) q[1];
sx q[1];
rz(-1.5643969) q[1];
sx q[1];
rz(2.0501627) q[1];
rz(-2.2131683) q[3];
sx q[3];
rz(-1.1517715) q[3];
sx q[3];
rz(-2.2862072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12863079) q[2];
sx q[2];
rz(-1.7065115) q[2];
sx q[2];
rz(1.5677412) q[2];
rz(-0.74357998) q[3];
sx q[3];
rz(-2.3502374) q[3];
sx q[3];
rz(-0.96405205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4670694) q[0];
sx q[0];
rz(-2.2860797) q[0];
sx q[0];
rz(-0.056644406) q[0];
rz(-1.4757587) q[1];
sx q[1];
rz(-1.1328127) q[1];
sx q[1];
rz(-1.3585565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751793) q[0];
sx q[0];
rz(-0.73690276) q[0];
sx q[0];
rz(-0.76216682) q[0];
rz(-0.26014056) q[2];
sx q[2];
rz(-0.67408757) q[2];
sx q[2];
rz(-2.0345794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0322235) q[1];
sx q[1];
rz(-2.0521161) q[1];
sx q[1];
rz(-0.66415031) q[1];
x q[2];
rz(0.56157063) q[3];
sx q[3];
rz(-1.6706781) q[3];
sx q[3];
rz(0.20250721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8654827) q[2];
sx q[2];
rz(-1.7296187) q[2];
sx q[2];
rz(-2.0152246) q[2];
rz(1.3364835) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(-2.0520463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4514076) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(0.13352808) q[0];
rz(0.97767699) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(2.8547063) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3505976) q[0];
sx q[0];
rz(-2.4387601) q[0];
sx q[0];
rz(-1.6446471) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8989766) q[2];
sx q[2];
rz(-1.5300473) q[2];
sx q[2];
rz(0.894067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0485018) q[1];
sx q[1];
rz(-1.347692) q[1];
sx q[1];
rz(-2.1483509) q[1];
x q[2];
rz(-2.7703076) q[3];
sx q[3];
rz(-0.6273191) q[3];
sx q[3];
rz(3.1255699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7202683) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.486091) q[2];
rz(0.36763516) q[3];
sx q[3];
rz(-2.1143819) q[3];
sx q[3];
rz(1.0953085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4246282) q[0];
sx q[0];
rz(-0.76911887) q[0];
sx q[0];
rz(-2.6677483) q[0];
rz(2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(-1.7582105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6598005) q[0];
sx q[0];
rz(-2.0542025) q[0];
sx q[0];
rz(0.76235911) q[0];
rz(-pi) q[1];
rz(-1.9888617) q[2];
sx q[2];
rz(-1.5961313) q[2];
sx q[2];
rz(2.096614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7935087) q[1];
sx q[1];
rz(-1.3474819) q[1];
sx q[1];
rz(0.23621724) q[1];
rz(-pi) q[2];
rz(0.9035191) q[3];
sx q[3];
rz(-0.28378962) q[3];
sx q[3];
rz(-1.2354148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9845) q[2];
sx q[2];
rz(-1.7326771) q[2];
sx q[2];
rz(-2.816693) q[2];
rz(0.38069185) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(1.0977753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710881) q[0];
sx q[0];
rz(-0.66362137) q[0];
sx q[0];
rz(-0.93884236) q[0];
rz(-0.70010575) q[1];
sx q[1];
rz(-0.63799262) q[1];
sx q[1];
rz(-0.28900388) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2719921) q[0];
sx q[0];
rz(-1.648252) q[0];
sx q[0];
rz(-1.3415706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2174066) q[2];
sx q[2];
rz(-1.1614387) q[2];
sx q[2];
rz(-1.6381303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71017337) q[1];
sx q[1];
rz(-1.050921) q[1];
sx q[1];
rz(-0.10837491) q[1];
x q[2];
rz(-0.17811505) q[3];
sx q[3];
rz(-2.7550972) q[3];
sx q[3];
rz(2.9535171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8812022) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(-0.41268665) q[2];
rz(-0.9203426) q[3];
sx q[3];
rz(-0.50655443) q[3];
sx q[3];
rz(-1.9119561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.6270139) q[0];
sx q[0];
rz(-0.98015061) q[0];
sx q[0];
rz(-2.8421616) q[0];
rz(-2.0191655) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(-1.0312414) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5297197) q[0];
sx q[0];
rz(-2.7752989) q[0];
sx q[0];
rz(-1.5407015) q[0];
x q[1];
rz(-0.38487969) q[2];
sx q[2];
rz(-0.49321929) q[2];
sx q[2];
rz(-0.2156336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.448062) q[1];
sx q[1];
rz(-1.281257) q[1];
sx q[1];
rz(2.2953643) q[1];
x q[2];
rz(-2.5599285) q[3];
sx q[3];
rz(-2.1278893) q[3];
sx q[3];
rz(2.795199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1624182) q[2];
sx q[2];
rz(-1.3434429) q[2];
sx q[2];
rz(0.24270414) q[2];
rz(-1.8414712) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(0.21568957) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31371394) q[0];
sx q[0];
rz(-2.9254881) q[0];
sx q[0];
rz(1.4208273) q[0];
rz(-0.31006649) q[1];
sx q[1];
rz(-0.87691751) q[1];
sx q[1];
rz(1.5440595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89429796) q[0];
sx q[0];
rz(-1.9820807) q[0];
sx q[0];
rz(-0.572834) q[0];
rz(-pi) q[1];
rz(1.3256959) q[2];
sx q[2];
rz(-2.7900006) q[2];
sx q[2];
rz(2.0577672) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4793933) q[1];
sx q[1];
rz(-1.7879722) q[1];
sx q[1];
rz(0.26611664) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15940729) q[3];
sx q[3];
rz(-1.7091572) q[3];
sx q[3];
rz(2.0970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0073504) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(1.3664112) q[2];
rz(0.8775231) q[3];
sx q[3];
rz(-1.6759422) q[3];
sx q[3];
rz(-0.88959488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14961814) q[0];
sx q[0];
rz(-1.5562973) q[0];
sx q[0];
rz(1.656116) q[0];
rz(-2.2772475) q[1];
sx q[1];
rz(-1.0135916) q[1];
sx q[1];
rz(2.7085173) q[1];
rz(-2.7769185) q[2];
sx q[2];
rz(-2.6119947) q[2];
sx q[2];
rz(0.15577015) q[2];
rz(-1.0917615) q[3];
sx q[3];
rz(-1.619966) q[3];
sx q[3];
rz(1.0109284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
