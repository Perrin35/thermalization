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
rz(-2.4401378) q[0];
sx q[0];
rz(-0.48818809) q[0];
sx q[0];
rz(0.39946431) q[0];
rz(2.0692628) q[1];
sx q[1];
rz(-0.89858276) q[1];
sx q[1];
rz(-2.8314765) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7785491) q[0];
sx q[0];
rz(-1.1828711) q[0];
sx q[0];
rz(-1.2004195) q[0];
rz(-pi) q[1];
rz(-2.1793876) q[2];
sx q[2];
rz(-0.96675379) q[2];
sx q[2];
rz(-0.60985987) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18067828) q[1];
sx q[1];
rz(-0.49094683) q[1];
sx q[1];
rz(0.21638201) q[1];
x q[2];
rz(3.0324494) q[3];
sx q[3];
rz(-1.1880842) q[3];
sx q[3];
rz(-2.7963137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8862137) q[2];
sx q[2];
rz(-1.1905406) q[2];
sx q[2];
rz(1.0211241) q[2];
rz(-1.6518263) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(0.79045734) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0806231) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(-0.27394295) q[0];
rz(-0.25320369) q[1];
sx q[1];
rz(-0.72414032) q[1];
sx q[1];
rz(-0.23987548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016099361) q[0];
sx q[0];
rz(-0.81917995) q[0];
sx q[0];
rz(-2.3045834) q[0];
rz(-0.48440964) q[2];
sx q[2];
rz(-2.5229215) q[2];
sx q[2];
rz(-1.8687488) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66544724) q[1];
sx q[1];
rz(-2.7283333) q[1];
sx q[1];
rz(-0.59935092) q[1];
x q[2];
rz(-1.3680458) q[3];
sx q[3];
rz(-0.79907954) q[3];
sx q[3];
rz(0.54982027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46951744) q[2];
sx q[2];
rz(-1.3152452) q[2];
sx q[2];
rz(-1.7496656) q[2];
rz(-0.64432708) q[3];
sx q[3];
rz(-1.143012) q[3];
sx q[3];
rz(-2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7740087) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(-0.17502633) q[0];
rz(-1.4352098) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(-2.3451436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38742313) q[0];
sx q[0];
rz(-1.8909374) q[0];
sx q[0];
rz(2.8693136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.286722) q[2];
sx q[2];
rz(-0.98978087) q[2];
sx q[2];
rz(0.43339455) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0092344) q[1];
sx q[1];
rz(-1.9446484) q[1];
sx q[1];
rz(-2.641885) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17603342) q[3];
sx q[3];
rz(-2.6280132) q[3];
sx q[3];
rz(-1.5777509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.986787) q[2];
sx q[2];
rz(-3.0974168) q[2];
sx q[2];
rz(-0.046048306) q[2];
rz(-1.1151399) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(-2.0980289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030877) q[0];
sx q[0];
rz(-2.5512888) q[0];
sx q[0];
rz(0.16511551) q[0];
rz(-1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(1.4547691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6957729) q[0];
sx q[0];
rz(-2.8759888) q[0];
sx q[0];
rz(0.37047343) q[0];
rz(-pi) q[1];
rz(-2.7023849) q[2];
sx q[2];
rz(-1.3813218) q[2];
sx q[2];
rz(2.2486692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36697179) q[1];
sx q[1];
rz(-2.4200388) q[1];
sx q[1];
rz(-1.5040892) q[1];
rz(0.28819542) q[3];
sx q[3];
rz(-1.1600947) q[3];
sx q[3];
rz(-2.2626143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7016865) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(-0.17882027) q[2];
rz(1.6245101) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4894067) q[0];
sx q[0];
rz(-1.711373) q[0];
sx q[0];
rz(1.408761) q[0];
rz(0.57213712) q[1];
sx q[1];
rz(-1.9913543) q[1];
sx q[1];
rz(0.22452721) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2437203) q[0];
sx q[0];
rz(-2.2099054) q[0];
sx q[0];
rz(-2.0452818) q[0];
rz(-pi) q[1];
rz(-0.83518402) q[2];
sx q[2];
rz(-1.3408644) q[2];
sx q[2];
rz(2.0809157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71345879) q[1];
sx q[1];
rz(-1.1468219) q[1];
sx q[1];
rz(0.78935539) q[1];
x q[2];
rz(-1.8711617) q[3];
sx q[3];
rz(-1.9711509) q[3];
sx q[3];
rz(-2.3107236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81395927) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(2.853788) q[2];
rz(0.52590251) q[3];
sx q[3];
rz(-1.9764427) q[3];
sx q[3];
rz(-2.5743918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7228058) q[0];
sx q[0];
rz(-2.4961508) q[0];
sx q[0];
rz(3.1268815) q[0];
rz(0.22444935) q[1];
sx q[1];
rz(-1.4597471) q[1];
sx q[1];
rz(-0.86123484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6447198) q[0];
sx q[0];
rz(-1.5278491) q[0];
sx q[0];
rz(-1.5687902) q[0];
rz(-pi) q[1];
rz(1.9056896) q[2];
sx q[2];
rz(-1.6097203) q[2];
sx q[2];
rz(0.30920497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6622195) q[1];
sx q[1];
rz(-1.8623317) q[1];
sx q[1];
rz(0.71034224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9771044) q[3];
sx q[3];
rz(-0.93014088) q[3];
sx q[3];
rz(1.759884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8031249) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(-3.0306446) q[2];
rz(-2.0153996) q[3];
sx q[3];
rz(-1.6323099) q[3];
sx q[3];
rz(0.66669983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.2647389) q[0];
sx q[0];
rz(-0.56346098) q[0];
sx q[0];
rz(-1.2891084) q[0];
rz(1.3339169) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(1.8905554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97020075) q[0];
sx q[0];
rz(-1.3193865) q[0];
sx q[0];
rz(2.2455477) q[0];
rz(1.2463322) q[2];
sx q[2];
rz(-2.1212541) q[2];
sx q[2];
rz(2.7389527) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.928628) q[1];
sx q[1];
rz(-1.8142833) q[1];
sx q[1];
rz(-2.2446471) q[1];
x q[2];
rz(-2.2381225) q[3];
sx q[3];
rz(-2.7674737) q[3];
sx q[3];
rz(1.7382517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0879849) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(1.4245707) q[2];
rz(0.21814957) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(2.654259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3394534) q[0];
sx q[0];
rz(-2.5886783) q[0];
sx q[0];
rz(0.46808991) q[0];
rz(-0.83174902) q[1];
sx q[1];
rz(-0.91107285) q[1];
sx q[1];
rz(-2.1591878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8881852) q[0];
sx q[0];
rz(-1.6352075) q[0];
sx q[0];
rz(3.0712434) q[0];
rz(2.1929412) q[2];
sx q[2];
rz(-1.7402044) q[2];
sx q[2];
rz(-2.624315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47226631) q[1];
sx q[1];
rz(-1.2756747) q[1];
sx q[1];
rz(-0.2448611) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3221075) q[3];
sx q[3];
rz(-2.9277955) q[3];
sx q[3];
rz(-2.4219861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9437647) q[2];
sx q[2];
rz(-2.6524537) q[2];
sx q[2];
rz(2.2042507) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.849544) q[3];
sx q[3];
rz(2.9924801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8775403) q[0];
sx q[0];
rz(-0.5101246) q[0];
sx q[0];
rz(-0.75793761) q[0];
rz(2.0856048) q[1];
sx q[1];
rz(-2.2426558) q[1];
sx q[1];
rz(-2.9008289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42007459) q[0];
sx q[0];
rz(-1.2596247) q[0];
sx q[0];
rz(-0.56367409) q[0];
rz(-pi) q[1];
rz(-0.40177138) q[2];
sx q[2];
rz(-0.86454287) q[2];
sx q[2];
rz(1.0313977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9537462) q[1];
sx q[1];
rz(-2.3930379) q[1];
sx q[1];
rz(-1.1331609) q[1];
x q[2];
rz(2.3189697) q[3];
sx q[3];
rz(-2.6634376) q[3];
sx q[3];
rz(-1.9229696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9869438) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(1.1692952) q[2];
rz(-2.5985006) q[3];
sx q[3];
rz(-1.8648059) q[3];
sx q[3];
rz(2.5663466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52694046) q[0];
sx q[0];
rz(-1.696538) q[0];
sx q[0];
rz(-2.8618983) q[0];
rz(-1.1138629) q[1];
sx q[1];
rz(-2.0852456) q[1];
sx q[1];
rz(-2.9934771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0954656) q[0];
sx q[0];
rz(-1.7610504) q[0];
sx q[0];
rz(-2.4664761) q[0];
rz(0.19377562) q[2];
sx q[2];
rz(-2.4420718) q[2];
sx q[2];
rz(1.4137088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.095299) q[1];
sx q[1];
rz(-2.110184) q[1];
sx q[1];
rz(2.1124243) q[1];
rz(2.8612265) q[3];
sx q[3];
rz(-0.70322817) q[3];
sx q[3];
rz(-0.91152465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0327882) q[2];
sx q[2];
rz(-2.2970436) q[2];
sx q[2];
rz(0.27842251) q[2];
rz(-1.4480715) q[3];
sx q[3];
rz(-1.4689987) q[3];
sx q[3];
rz(1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096810452) q[0];
sx q[0];
rz(-0.7881931) q[0];
sx q[0];
rz(-1.0544554) q[0];
rz(2.0432368) q[1];
sx q[1];
rz(-1.8041246) q[1];
sx q[1];
rz(0.95380797) q[1];
rz(-2.5658812) q[2];
sx q[2];
rz(-1.4288708) q[2];
sx q[2];
rz(2.5484011) q[2];
rz(1.2145417) q[3];
sx q[3];
rz(-0.14394017) q[3];
sx q[3];
rz(-2.3059358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
