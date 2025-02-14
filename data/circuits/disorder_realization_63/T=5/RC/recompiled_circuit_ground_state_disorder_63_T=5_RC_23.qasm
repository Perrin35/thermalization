OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7665793) q[0];
sx q[0];
rz(-1.2393247) q[0];
sx q[0];
rz(1.0898606) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(2.2383632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3147557) q[0];
sx q[0];
rz(-0.96393782) q[0];
sx q[0];
rz(-2.3640102) q[0];
rz(-pi) q[1];
rz(1.5888693) q[2];
sx q[2];
rz(-0.68328349) q[2];
sx q[2];
rz(2.9073213) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5168239) q[1];
sx q[1];
rz(-0.98691578) q[1];
sx q[1];
rz(1.1560239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7890268) q[3];
sx q[3];
rz(-0.50980836) q[3];
sx q[3];
rz(0.17579432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.696306) q[2];
sx q[2];
rz(-0.13662766) q[2];
sx q[2];
rz(-2.7083) q[2];
rz(0.28996921) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(2.8210631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9669773) q[0];
sx q[0];
rz(-0.86242914) q[0];
sx q[0];
rz(-1.8509266) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(0.89250934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9226869) q[0];
sx q[0];
rz(-1.7480047) q[0];
sx q[0];
rz(-3.1161581) q[0];
rz(-1.9694445) q[2];
sx q[2];
rz(-2.1855693) q[2];
sx q[2];
rz(-0.022567011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.635517) q[1];
sx q[1];
rz(-2.2137224) q[1];
sx q[1];
rz(1.9159813) q[1];
rz(1.1970911) q[3];
sx q[3];
rz(-1.804816) q[3];
sx q[3];
rz(1.4732052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14900011) q[2];
sx q[2];
rz(-1.0789824) q[2];
sx q[2];
rz(2.0089669) q[2];
rz(-2.4752786) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(-1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742974) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(-0.98180109) q[0];
rz(2.1994195) q[1];
sx q[1];
rz(-1.8323545) q[1];
sx q[1];
rz(-0.91032496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2766586) q[0];
sx q[0];
rz(-0.66682839) q[0];
sx q[0];
rz(2.4163321) q[0];
rz(1.6907755) q[2];
sx q[2];
rz(-1.2960805) q[2];
sx q[2];
rz(-0.071823013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.040346216) q[1];
sx q[1];
rz(-1.2381499) q[1];
sx q[1];
rz(0.31942792) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0542473) q[3];
sx q[3];
rz(-2.6539565) q[3];
sx q[3];
rz(2.3211562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94074774) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(-0.18690898) q[2];
rz(-2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(2.4439268) q[0];
rz(-0.54667073) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(-1.9465416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9126667) q[0];
sx q[0];
rz(-0.9879092) q[0];
sx q[0];
rz(-2.7253662) q[0];
x q[1];
rz(0.19104345) q[2];
sx q[2];
rz(-1.7037539) q[2];
sx q[2];
rz(-1.3933448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4801712) q[1];
sx q[1];
rz(-1.464559) q[1];
sx q[1];
rz(1.7376971) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0680858) q[3];
sx q[3];
rz(-1.2871082) q[3];
sx q[3];
rz(-2.7105041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(-2.9441693) q[2];
rz(-0.10489634) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(-1.0092258) q[0];
rz(0.028845305) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(-2.4829594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4761676) q[0];
sx q[0];
rz(-1.9840249) q[0];
sx q[0];
rz(-0.40979235) q[0];
x q[1];
rz(0.20435996) q[2];
sx q[2];
rz(-1.0931226) q[2];
sx q[2];
rz(-0.43276873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5168136) q[1];
sx q[1];
rz(-0.1923696) q[1];
sx q[1];
rz(-0.50290795) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4464408) q[3];
sx q[3];
rz(-1.3926695) q[3];
sx q[3];
rz(-0.069308829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7625526) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(-2.9075882) q[2];
rz(2.0512569) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(-1.0696629) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2284018) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(-1.2983904) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(1.7826805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8131065) q[0];
sx q[0];
rz(-0.66603249) q[0];
sx q[0];
rz(-0.12909992) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4491389) q[2];
sx q[2];
rz(-0.34231774) q[2];
sx q[2];
rz(2.3794425) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1383776) q[1];
sx q[1];
rz(-1.5724465) q[1];
sx q[1];
rz(1.5695851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6796713) q[3];
sx q[3];
rz(-1.2826903) q[3];
sx q[3];
rz(1.896331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4141554) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(0.1725014) q[2];
rz(2.6532459) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(-0.33997047) q[0];
rz(-2.8151457) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(1.9065769) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8064855) q[0];
sx q[0];
rz(-0.77101427) q[0];
sx q[0];
rz(-0.54698995) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65151524) q[2];
sx q[2];
rz(-0.55547248) q[2];
sx q[2];
rz(3.0998203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.75720471) q[1];
sx q[1];
rz(-1.5761307) q[1];
sx q[1];
rz(0.47949201) q[1];
rz(-pi) q[2];
rz(2.7478479) q[3];
sx q[3];
rz(-1.3662158) q[3];
sx q[3];
rz(2.1285299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0299224) q[2];
sx q[2];
rz(-2.3387574) q[2];
sx q[2];
rz(-2.5862582) q[2];
rz(-2.9739144) q[3];
sx q[3];
rz(-1.6043112) q[3];
sx q[3];
rz(0.47784561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70917201) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(-2.0322556) q[0];
rz(2.0093911) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(-2.1999377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50186673) q[0];
sx q[0];
rz(-1.7302156) q[0];
sx q[0];
rz(-1.0109883) q[0];
rz(-1.2852816) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(-0.51644737) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1600765) q[1];
sx q[1];
rz(-1.6972516) q[1];
sx q[1];
rz(2.2156348) q[1];
x q[2];
rz(0.38576084) q[3];
sx q[3];
rz(-1.3732135) q[3];
sx q[3];
rz(-2.945874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38810101) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(1.6027742) q[2];
rz(1.7704891) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(-0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082315363) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(2.1737461) q[0];
rz(-1.8702501) q[1];
sx q[1];
rz(-1.2920734) q[1];
sx q[1];
rz(3.079792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7907352) q[0];
sx q[0];
rz(-2.4567607) q[0];
sx q[0];
rz(-2.5032296) q[0];
rz(0.25904663) q[2];
sx q[2];
rz(-1.7648089) q[2];
sx q[2];
rz(-1.8679179) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7605412) q[1];
sx q[1];
rz(-2.2534721) q[1];
sx q[1];
rz(-1.6939916) q[1];
rz(-pi) q[2];
rz(0.7155719) q[3];
sx q[3];
rz(-0.4288097) q[3];
sx q[3];
rz(0.174086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1900078) q[2];
sx q[2];
rz(-2.9374359) q[2];
sx q[2];
rz(0.07587138) q[2];
rz(-1.9742981) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(-2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9645914) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(-0.28468537) q[0];
rz(-3.0668861) q[1];
sx q[1];
rz(-1.2295405) q[1];
sx q[1];
rz(1.4178735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41499785) q[0];
sx q[0];
rz(-2.3624067) q[0];
sx q[0];
rz(-0.27611538) q[0];
rz(-pi) q[1];
rz(-2.9921164) q[2];
sx q[2];
rz(-2.1227395) q[2];
sx q[2];
rz(1.9665444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26894179) q[1];
sx q[1];
rz(-1.9851908) q[1];
sx q[1];
rz(-0.68434836) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6154434) q[3];
sx q[3];
rz(-2.4952125) q[3];
sx q[3];
rz(0.97585362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8031926) q[2];
sx q[2];
rz(-2.000838) q[2];
sx q[2];
rz(-1.2791862) q[2];
rz(0.98215669) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(-0.88472432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2411156) q[0];
sx q[0];
rz(-1.3237088) q[0];
sx q[0];
rz(1.5644912) q[0];
rz(2.1263532) q[1];
sx q[1];
rz(-1.1493586) q[1];
sx q[1];
rz(-1.1372067) q[1];
rz(1.954409) q[2];
sx q[2];
rz(-0.28228471) q[2];
sx q[2];
rz(-1.6105035) q[2];
rz(-0.99526631) q[3];
sx q[3];
rz(-0.91777663) q[3];
sx q[3];
rz(-1.4025626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
