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
rz(-0.17679086) q[0];
rz(-1.085936) q[1];
sx q[1];
rz(-0.50958264) q[1];
sx q[1];
rz(-3.0546319) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3172042) q[0];
sx q[0];
rz(-1.1939216) q[0];
sx q[0];
rz(-2.5254842) q[0];
x q[1];
rz(-1.3777138) q[2];
sx q[2];
rz(-1.0079441) q[2];
sx q[2];
rz(2.4342997) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0551974) q[1];
sx q[1];
rz(-2.9922303) q[1];
sx q[1];
rz(-1.0284852) q[1];
rz(-2.7417408) q[3];
sx q[3];
rz(-2.4824871) q[3];
sx q[3];
rz(0.53149283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95734346) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(0.43821487) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(-2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1205207) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(-0.72129321) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(-0.24478197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0185623) q[0];
sx q[0];
rz(-2.5726312) q[0];
sx q[0];
rz(2.2267692) q[0];
rz(2.7049644) q[2];
sx q[2];
rz(-0.87030137) q[2];
sx q[2];
rz(-1.8281405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3192057) q[1];
sx q[1];
rz(-2.165395) q[1];
sx q[1];
rz(0.90984224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8448288) q[3];
sx q[3];
rz(-1.4302084) q[3];
sx q[3];
rz(0.97749099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(-2.2954693) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(2.164771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1035128) q[0];
sx q[0];
rz(-0.80779034) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(-1.6733276) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-0.066468261) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9212355) q[0];
sx q[0];
rz(-0.68676239) q[0];
sx q[0];
rz(-3.0007017) q[0];
rz(-1.0280175) q[2];
sx q[2];
rz(-1.3508056) q[2];
sx q[2];
rz(2.5808711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3507562) q[1];
sx q[1];
rz(-0.15201026) q[1];
sx q[1];
rz(2.445447) q[1];
rz(-1.4044365) q[3];
sx q[3];
rz(-2.3041953) q[3];
sx q[3];
rz(2.211253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0299783) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(2.8288793) q[2];
rz(-0.2615658) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0161491) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(0.087015986) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5785297) q[1];
sx q[1];
rz(-3.0932025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2242369) q[0];
sx q[0];
rz(-2.356987) q[0];
sx q[0];
rz(-0.18152256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.264123) q[2];
sx q[2];
rz(-0.23261586) q[2];
sx q[2];
rz(2.2777429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9632513) q[1];
sx q[1];
rz(-1.6251441) q[1];
sx q[1];
rz(0.69596325) q[1];
x q[2];
rz(2.2807924) q[3];
sx q[3];
rz(-1.080092) q[3];
sx q[3];
rz(-1.8933637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4309569) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(-0.44130138) q[2];
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
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271673) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(1.4432888) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(1.7194933) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.559645) q[0];
sx q[0];
rz(-2.1113681) q[0];
sx q[0];
rz(1.4076884) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70962064) q[2];
sx q[2];
rz(-2.8292252) q[2];
sx q[2];
rz(-1.4210977) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7362914) q[1];
sx q[1];
rz(-1.4854238) q[1];
sx q[1];
rz(-3.0194204) q[1];
x q[2];
rz(2.2660013) q[3];
sx q[3];
rz(-0.79447047) q[3];
sx q[3];
rz(-2.9231185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(1.5606073) q[2];
rz(-1.3382781) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(-2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(3.1230028) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(-2.4216968) q[0];
rz(-0.24049354) q[1];
sx q[1];
rz(-1.980282) q[1];
sx q[1];
rz(0.24615157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7712646) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(-1.5465897) q[0];
rz(-pi) q[1];
rz(-0.75743875) q[2];
sx q[2];
rz(-1.2691109) q[2];
sx q[2];
rz(1.2011004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38118874) q[1];
sx q[1];
rz(-0.84263984) q[1];
sx q[1];
rz(-1.3959753) q[1];
x q[2];
rz(1.9897377) q[3];
sx q[3];
rz(-0.92760689) q[3];
sx q[3];
rz(0.96574984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4973732) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(-0.52404809) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2972357) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(0.921184) q[0];
rz(-1.720403) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(0.99501077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5748844) q[0];
sx q[0];
rz(-0.44429438) q[0];
sx q[0];
rz(-2.1106476) q[0];
rz(-3.0757881) q[2];
sx q[2];
rz(-1.9901681) q[2];
sx q[2];
rz(-2.9247627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2240963) q[1];
sx q[1];
rz(-0.69643785) q[1];
sx q[1];
rz(-2.4516979) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80713804) q[3];
sx q[3];
rz(-0.57315956) q[3];
sx q[3];
rz(-1.5218211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.93817389) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(-0.43756884) q[2];
rz(-2.3214052) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95847982) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(1.0108277) q[0];
rz(-0.084835947) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(-2.5929677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91544596) q[0];
sx q[0];
rz(-1.5391162) q[0];
sx q[0];
rz(-1.7466387) q[0];
rz(1.8993511) q[2];
sx q[2];
rz(-1.1588084) q[2];
sx q[2];
rz(-3.0871473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6514142) q[1];
sx q[1];
rz(-1.427535) q[1];
sx q[1];
rz(0.68639042) q[1];
x q[2];
rz(2.2574378) q[3];
sx q[3];
rz(-0.13152371) q[3];
sx q[3];
rz(2.2145074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34714547) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(0.53317201) q[2];
rz(-1.0779856) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(0.50895154) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0549523) q[0];
sx q[0];
rz(-2.7821879) q[0];
sx q[0];
rz(-2.3593498) q[0];
rz(3.0714463) q[1];
sx q[1];
rz(-0.47824305) q[1];
sx q[1];
rz(-2.9152962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5180888) q[0];
sx q[0];
rz(-1.6395578) q[0];
sx q[0];
rz(-3.070773) q[0];
x q[1];
rz(-1.802472) q[2];
sx q[2];
rz(-0.69139987) q[2];
sx q[2];
rz(-2.0350128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1038166) q[1];
sx q[1];
rz(-2.1063707) q[1];
sx q[1];
rz(1.4431672) q[1];
x q[2];
rz(1.6582011) q[3];
sx q[3];
rz(-2.8897396) q[3];
sx q[3];
rz(0.15628584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0349064) q[2];
sx q[2];
rz(-0.28768134) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(2.8790224) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(-0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(-2.9578399) q[0];
rz(-0.19206583) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(-0.62350887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.469891) q[0];
sx q[0];
rz(-0.56607038) q[0];
sx q[0];
rz(-0.59622391) q[0];
rz(-2.170606) q[2];
sx q[2];
rz(-1.0375334) q[2];
sx q[2];
rz(2.2704934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3721613) q[1];
sx q[1];
rz(-1.4648533) q[1];
sx q[1];
rz(1.547936) q[1];
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
rz(1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(-2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9380209) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(-2.6899295) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(2.6091433) q[2];
sx q[2];
rz(-1.3146853) q[2];
sx q[2];
rz(2.5934861) q[2];
rz(-0.59521159) q[3];
sx q[3];
rz(-2.8122936) q[3];
sx q[3];
rz(-0.82174792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
