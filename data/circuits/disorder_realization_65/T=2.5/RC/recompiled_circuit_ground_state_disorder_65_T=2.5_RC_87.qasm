OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.559691) q[0];
sx q[0];
rz(3.0540967) q[0];
sx q[0];
rz(9.3257186) q[0];
rz(2.8506408) q[1];
sx q[1];
rz(-1.7506316) q[1];
sx q[1];
rz(-1.5150392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5387548) q[0];
sx q[0];
rz(-1.3405717) q[0];
sx q[0];
rz(1.099596) q[0];
rz(2.6886772) q[2];
sx q[2];
rz(-1.1057111) q[2];
sx q[2];
rz(1.0207435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8352055) q[1];
sx q[1];
rz(-0.82535078) q[1];
sx q[1];
rz(-0.21791509) q[1];
rz(2.4060449) q[3];
sx q[3];
rz(-0.41829073) q[3];
sx q[3];
rz(2.4052704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0809975) q[2];
sx q[2];
rz(-0.011757714) q[2];
sx q[2];
rz(0.11059977) q[2];
rz(-0.009875385) q[3];
sx q[3];
rz(-0.014009352) q[3];
sx q[3];
rz(-0.88147718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1572384) q[0];
sx q[0];
rz(-3.053559) q[0];
sx q[0];
rz(1.7617759) q[0];
rz(-0.012160483) q[1];
sx q[1];
rz(-0.98406839) q[1];
sx q[1];
rz(1.5252569) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1943239) q[0];
sx q[0];
rz(-2.8229694) q[0];
sx q[0];
rz(-0.18100321) q[0];
rz(-pi) q[1];
rz(1.3697769) q[2];
sx q[2];
rz(-1.9819385) q[2];
sx q[2];
rz(2.0697921) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5926275) q[1];
sx q[1];
rz(-1.4969259) q[1];
sx q[1];
rz(0.17358257) q[1];
x q[2];
rz(0.44604519) q[3];
sx q[3];
rz(-1.4592663) q[3];
sx q[3];
rz(-2.1457637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3297367) q[2];
sx q[2];
rz(-0.024710329) q[2];
sx q[2];
rz(3.109566) q[2];
rz(0.61102593) q[3];
sx q[3];
rz(-1.150584) q[3];
sx q[3];
rz(1.7648296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585612) q[0];
sx q[0];
rz(-2.2295781) q[0];
sx q[0];
rz(1.4645905) q[0];
rz(-1.5088082) q[1];
sx q[1];
rz(-1.8784411) q[1];
sx q[1];
rz(3.0806697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2476544) q[0];
sx q[0];
rz(-1.1919855) q[0];
sx q[0];
rz(1.3734067) q[0];
rz(-pi) q[1];
rz(1.6834358) q[2];
sx q[2];
rz(-1.4667744) q[2];
sx q[2];
rz(-2.011781) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4767392) q[1];
sx q[1];
rz(-1.7283943) q[1];
sx q[1];
rz(-1.3955302) q[1];
x q[2];
rz(-2.7512128) q[3];
sx q[3];
rz(-1.2900616) q[3];
sx q[3];
rz(-2.8685307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9691201) q[2];
sx q[2];
rz(-0.0027593297) q[2];
sx q[2];
rz(0.74392444) q[2];
rz(2.7562691) q[3];
sx q[3];
rz(-3.139956) q[3];
sx q[3];
rz(-0.90414944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0124403) q[0];
sx q[0];
rz(-0.04239447) q[0];
sx q[0];
rz(-0.84268919) q[0];
rz(2.1878302) q[1];
sx q[1];
rz(-1.5019491) q[1];
sx q[1];
rz(-1.5488254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208389) q[0];
sx q[0];
rz(-0.11891236) q[0];
sx q[0];
rz(-0.17158385) q[0];
rz(0.15797961) q[2];
sx q[2];
rz(-1.176774) q[2];
sx q[2];
rz(-1.3713149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5332414) q[1];
sx q[1];
rz(-2.5035739) q[1];
sx q[1];
rz(-1.4916496) q[1];
rz(0.78068818) q[3];
sx q[3];
rz(-1.4997291) q[3];
sx q[3];
rz(2.4960757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10161764) q[2];
sx q[2];
rz(-3.1246694) q[2];
sx q[2];
rz(0.42538154) q[2];
rz(-1.3500805) q[3];
sx q[3];
rz(-1.5964419) q[3];
sx q[3];
rz(2.8369956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.72163433) q[0];
sx q[0];
rz(-0.0036792734) q[0];
sx q[0];
rz(0.71596181) q[0];
rz(-1.6250027) q[1];
sx q[1];
rz(-2.6874028) q[1];
sx q[1];
rz(-2.9862278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778188) q[0];
sx q[0];
rz(-1.4646104) q[0];
sx q[0];
rz(-3.1255577) q[0];
rz(-1.4808319) q[2];
sx q[2];
rz(-1.1954006) q[2];
sx q[2];
rz(-1.1590955) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1137312) q[1];
sx q[1];
rz(-1.3044954) q[1];
sx q[1];
rz(1.6990425) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.243216) q[3];
sx q[3];
rz(-1.9168609) q[3];
sx q[3];
rz(-2.4806674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6905489) q[2];
sx q[2];
rz(-2.0205708) q[2];
sx q[2];
rz(1.6070018) q[2];
rz(-1.0891886) q[3];
sx q[3];
rz(-3.1179699) q[3];
sx q[3];
rz(1.7781809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71517336) q[0];
sx q[0];
rz(-0.024113163) q[0];
sx q[0];
rz(0.56060785) q[0];
rz(2.2757065) q[1];
sx q[1];
rz(-0.0090323369) q[1];
sx q[1];
rz(-0.83585709) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2888907) q[0];
sx q[0];
rz(-1.2749854) q[0];
sx q[0];
rz(-1.5519141) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1634665) q[2];
sx q[2];
rz(-1.5390804) q[2];
sx q[2];
rz(-0.18975604) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.071358368) q[1];
sx q[1];
rz(-2.8615983) q[1];
sx q[1];
rz(-2.7482949) q[1];
x q[2];
rz(-2.4847651) q[3];
sx q[3];
rz(-0.39077911) q[3];
sx q[3];
rz(0.093599044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5777638) q[2];
sx q[2];
rz(-1.3059629) q[2];
sx q[2];
rz(1.5888265) q[2];
rz(2.0896437) q[3];
sx q[3];
rz(-2.6237539) q[3];
sx q[3];
rz(0.53798211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.8727259) q[0];
sx q[0];
rz(-0.59393847) q[0];
sx q[0];
rz(-1.8178222) q[0];
rz(-0.84953228) q[1];
sx q[1];
rz(-3.140026) q[1];
sx q[1];
rz(0.82226396) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182151) q[0];
sx q[0];
rz(-1.4750622) q[0];
sx q[0];
rz(1.7507919) q[0];
rz(-pi) q[1];
rz(-1.6258114) q[2];
sx q[2];
rz(-1.9085247) q[2];
sx q[2];
rz(1.6134451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9692053) q[1];
sx q[1];
rz(-1.6403461) q[1];
sx q[1];
rz(0.0065253536) q[1];
rz(-1.7674462) q[3];
sx q[3];
rz(-1.8314284) q[3];
sx q[3];
rz(-2.9109241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5620586) q[2];
sx q[2];
rz(-1.9066533) q[2];
sx q[2];
rz(0.84260064) q[2];
rz(-0.54251999) q[3];
sx q[3];
rz(-3.123057) q[3];
sx q[3];
rz(-0.59244573) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056996718) q[0];
sx q[0];
rz(-1.5942986) q[0];
sx q[0];
rz(-1.0513167) q[0];
rz(1.7683138) q[1];
sx q[1];
rz(-3.1128502) q[1];
sx q[1];
rz(-1.9017259) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9997885) q[0];
sx q[0];
rz(-1.3281137) q[0];
sx q[0];
rz(1.5304327) q[0];
rz(-pi) q[1];
rz(-1.6592024) q[2];
sx q[2];
rz(-2.7475221) q[2];
sx q[2];
rz(-1.5752156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4989325) q[1];
sx q[1];
rz(-2.1524384) q[1];
sx q[1];
rz(-10/(13*pi)) q[1];
rz(-pi) q[2];
rz(2.9536581) q[3];
sx q[3];
rz(-1.3085164) q[3];
sx q[3];
rz(-0.89406195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54590589) q[2];
sx q[2];
rz(-3.1133856) q[2];
sx q[2];
rz(2.9941881) q[2];
rz(1.5386511) q[3];
sx q[3];
rz(-1.7037062) q[3];
sx q[3];
rz(-0.182972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1786757) q[0];
sx q[0];
rz(-2.8067532) q[0];
sx q[0];
rz(1.8033002) q[0];
rz(-1.6637404) q[1];
sx q[1];
rz(-1.1344974) q[1];
sx q[1];
rz(-0.17325625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4497183) q[0];
sx q[0];
rz(-1.0074179) q[0];
sx q[0];
rz(-2.0823576) q[0];
x q[1];
rz(0.40399) q[2];
sx q[2];
rz(-2.1232018) q[2];
sx q[2];
rz(-0.45726038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.063364) q[1];
sx q[1];
rz(-1.8174108) q[1];
sx q[1];
rz(1.7298841) q[1];
rz(-pi) q[2];
rz(1.162649) q[3];
sx q[3];
rz(-0.014685304) q[3];
sx q[3];
rz(2.5774308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0080537) q[2];
sx q[2];
rz(-3.1296788) q[2];
sx q[2];
rz(-2.1214205) q[2];
rz(-0.080848761) q[3];
sx q[3];
rz(-0.0011708502) q[3];
sx q[3];
rz(-2.0118058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6218277) q[0];
sx q[0];
rz(-2.8563359) q[0];
sx q[0];
rz(-1.6602302) q[0];
rz(-0.18352428) q[1];
sx q[1];
rz(-2.8158975) q[1];
sx q[1];
rz(-3.0537925) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0073951) q[0];
sx q[0];
rz(-2.0124276) q[0];
sx q[0];
rz(1.6820119) q[0];
rz(-2.8351458) q[2];
sx q[2];
rz(-1.4922835) q[2];
sx q[2];
rz(-1.3456089) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1488249) q[1];
sx q[1];
rz(-1.449905) q[1];
sx q[1];
rz(-2.8962027) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2689104) q[3];
sx q[3];
rz(-2.2175237) q[3];
sx q[3];
rz(-0.069096126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8537019) q[2];
sx q[2];
rz(-3.1292858) q[2];
sx q[2];
rz(0.88407174) q[2];
rz(0.70841241) q[3];
sx q[3];
rz(-3.1413779) q[3];
sx q[3];
rz(2.8475672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6207599) q[0];
sx q[0];
rz(-1.5694869) q[0];
sx q[0];
rz(-1.6318305) q[0];
rz(3.1132501) q[1];
sx q[1];
rz(-2.5181073) q[1];
sx q[1];
rz(0.070652031) q[1];
rz(1.2030884) q[2];
sx q[2];
rz(-1.9919507) q[2];
sx q[2];
rz(0.19797355) q[2];
rz(0.72851758) q[3];
sx q[3];
rz(-1.3195724) q[3];
sx q[3];
rz(2.6825678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
