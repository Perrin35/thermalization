OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5476721) q[0];
sx q[0];
rz(-0.5991109) q[0];
sx q[0];
rz(2.9648018) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(-2.63201) q[1];
sx q[1];
rz(-0.086960763) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3172042) q[0];
sx q[0];
rz(-1.1939216) q[0];
sx q[0];
rz(2.5254842) q[0];
x q[1];
rz(2.5702417) q[2];
sx q[2];
rz(-1.4077912) q[2];
sx q[2];
rz(0.75955078) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0217758) q[1];
sx q[1];
rz(-1.6476742) q[1];
sx q[1];
rz(-1.6989811) q[1];
rz(-pi) q[2];
rz(1.8636892) q[3];
sx q[3];
rz(-0.97145069) q[3];
sx q[3];
rz(1.0224316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1842492) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(0.43821487) q[3];
sx q[3];
rz(-0.44132909) q[3];
sx q[3];
rz(-0.87619585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1205207) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(1.0907809) q[0];
rz(-2.4202994) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(0.24478197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2782841) q[0];
sx q[0];
rz(-1.1296891) q[0];
sx q[0];
rz(0.37190227) q[0];
rz(1.1058979) q[2];
sx q[2];
rz(-0.80543488) q[2];
sx q[2];
rz(-1.9400846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3192057) q[1];
sx q[1];
rz(-0.97619769) q[1];
sx q[1];
rz(-0.90984224) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29676389) q[3];
sx q[3];
rz(-1.7113842) q[3];
sx q[3];
rz(2.1641017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(-2.2954693) q[2];
rz(2.7766679) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038079809) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(1.6733276) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(0.066468261) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398427) q[0];
sx q[0];
rz(-2.2494613) q[0];
sx q[0];
rz(1.4561653) q[0];
rz(-1.1622381) q[2];
sx q[2];
rz(-0.58149946) q[2];
sx q[2];
rz(1.3571908) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.089515226) q[1];
sx q[1];
rz(-1.6680535) q[1];
sx q[1];
rz(-0.11701028) q[1];
rz(2.9598438) q[3];
sx q[3];
rz(-0.74859339) q[3];
sx q[3];
rz(-2.4570217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0299783) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(0.31271333) q[2];
rz(-2.8800268) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(-0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12544352) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(-0.087015986) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5785297) q[1];
sx q[1];
rz(0.048390128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780086) q[0];
sx q[0];
rz(-2.3391294) q[0];
sx q[0];
rz(1.3924696) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15026413) q[2];
sx q[2];
rz(-1.3925526) q[2];
sx q[2];
rz(1.5709637) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9632513) q[1];
sx q[1];
rz(-1.5164485) q[1];
sx q[1];
rz(-0.69596325) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88416962) q[3];
sx q[3];
rz(-2.3034952) q[3];
sx q[3];
rz(0.17894408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4309569) q[2];
sx q[2];
rz(-0.29971665) q[2];
sx q[2];
rz(-2.7002913) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(1.8080447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271673) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(-1.4432888) q[1];
sx q[1];
rz(-1.8836421) q[1];
sx q[1];
rz(-1.4220994) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073339065) q[0];
sx q[0];
rz(-1.4311106) q[0];
sx q[0];
rz(-0.54648593) q[0];
rz(-2.431972) q[2];
sx q[2];
rz(-2.8292252) q[2];
sx q[2];
rz(1.720495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.369097) q[1];
sx q[1];
rz(-0.14892347) q[1];
sx q[1];
rz(2.5293674) q[1];
rz(0.57797076) q[3];
sx q[3];
rz(-2.1506566) q[3];
sx q[3];
rz(-0.65366369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(-1.5809853) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(0.54792255) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1230028) q[0];
sx q[0];
rz(-0.81273166) q[0];
sx q[0];
rz(-2.4216968) q[0];
rz(0.24049354) q[1];
sx q[1];
rz(-1.980282) q[1];
sx q[1];
rz(-0.24615157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1937597) q[0];
sx q[0];
rz(-1.5475377) q[0];
sx q[0];
rz(-2.860753) q[0];
rz(-pi) q[1];
rz(-1.9754378) q[2];
sx q[2];
rz(-2.2863467) q[2];
sx q[2];
rz(-0.095794769) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0689905) q[1];
sx q[1];
rz(-1.700987) q[1];
sx q[1];
rz(-0.73579244) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9897377) q[3];
sx q[3];
rz(-2.2139858) q[3];
sx q[3];
rz(-2.1758428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64421946) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(2.6175446) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(-0.016949765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84435695) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(-2.2204087) q[0];
rz(1.4211897) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(0.99501077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667082) q[0];
sx q[0];
rz(-2.6972983) q[0];
sx q[0];
rz(2.1106476) q[0];
x q[1];
rz(-3.0757881) q[2];
sx q[2];
rz(-1.1514246) q[2];
sx q[2];
rz(2.9247627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.095671818) q[1];
sx q[1];
rz(-1.0532044) q[1];
sx q[1];
rz(-2.059883) q[1];
x q[2];
rz(2.7217676) q[3];
sx q[3];
rz(-1.1683162) q[3];
sx q[3];
rz(-2.4703006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2034188) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(2.7040238) q[2];
rz(-0.82018745) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(-2.1831128) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(1.0108277) q[0];
rz(-3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(0.54862499) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4806145) q[0];
sx q[0];
rz(-1.3950431) q[0];
sx q[0];
rz(-3.1094167) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8993511) q[2];
sx q[2];
rz(-1.9827843) q[2];
sx q[2];
rz(3.0871473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49017845) q[1];
sx q[1];
rz(-1.7140577) q[1];
sx q[1];
rz(-0.68639042) q[1];
rz(-pi) q[2];
rz(-0.083666936) q[3];
sx q[3];
rz(-1.4691969) q[3];
sx q[3];
rz(-1.6179832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7944472) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(-2.6084206) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(-0.78224283) q[0];
rz(-0.070146322) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(2.9152962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057581456) q[0];
sx q[0];
rz(-1.5001443) q[0];
sx q[0];
rz(1.6397301) q[0];
x q[1];
rz(-2.2489297) q[2];
sx q[2];
rz(-1.7177267) q[2];
sx q[2];
rz(0.28444296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.037776064) q[1];
sx q[1];
rz(-2.1063707) q[1];
sx q[1];
rz(-1.6984254) q[1];
x q[2];
rz(-3.1191344) q[3];
sx q[3];
rz(-1.8216672) q[3];
sx q[3];
rz(2.89507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0349064) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.7314343) q[2];
rz(-2.8790224) q[3];
sx q[3];
rz(-1.5841443) q[3];
sx q[3];
rz(-0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866289) q[0];
sx q[0];
rz(-0.2121191) q[0];
sx q[0];
rz(0.18375272) q[0];
rz(-2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(-0.62350887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6717017) q[0];
sx q[0];
rz(-2.5755223) q[0];
sx q[0];
rz(0.59622391) q[0];
rz(2.170606) q[2];
sx q[2];
rz(-2.1040593) q[2];
sx q[2];
rz(2.2704934) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19621721) q[1];
sx q[1];
rz(-1.5480642) q[1];
sx q[1];
rz(0.10597056) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4922114) q[3];
sx q[3];
rz(-2.7659263) q[3];
sx q[3];
rz(-2.5933655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3925675) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(-0.17879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9380209) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(0.45166311) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(2.6091433) q[2];
sx q[2];
rz(-1.3146853) q[2];
sx q[2];
rz(2.5934861) q[2];
rz(0.27576294) q[3];
sx q[3];
rz(-1.388474) q[3];
sx q[3];
rz(0.17920517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
