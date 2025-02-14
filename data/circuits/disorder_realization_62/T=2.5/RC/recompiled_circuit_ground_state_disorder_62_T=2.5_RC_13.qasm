OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.59392053) q[0];
sx q[0];
rz(3.7407036) q[0];
sx q[0];
rz(9.6015688) q[0];
rz(-1.085936) q[1];
sx q[1];
rz(-0.50958264) q[1];
sx q[1];
rz(-3.0546319) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22573839) q[0];
sx q[0];
rz(-2.432352) q[0];
sx q[0];
rz(0.60053696) q[0];
rz(-pi) q[1];
rz(-1.7638788) q[2];
sx q[2];
rz(-1.0079441) q[2];
sx q[2];
rz(-2.4342997) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1198169) q[1];
sx q[1];
rz(-1.6476742) q[1];
sx q[1];
rz(-1.4426115) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7417408) q[3];
sx q[3];
rz(-2.4824871) q[3];
sx q[3];
rz(2.6100998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95734346) q[2];
sx q[2];
rz(-0.49746305) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(2.7033778) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-1.483016) q[1];
sx q[1];
rz(0.24478197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1230303) q[0];
sx q[0];
rz(-2.5726312) q[0];
sx q[0];
rz(-2.2267692) q[0];
rz(-pi) q[1];
rz(-0.4366283) q[2];
sx q[2];
rz(-2.2712913) q[2];
sx q[2];
rz(-1.3134522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.822387) q[1];
sx q[1];
rz(-2.165395) q[1];
sx q[1];
rz(0.90984224) q[1];
rz(-0.29676389) q[3];
sx q[3];
rz(-1.4302084) q[3];
sx q[3];
rz(-0.97749099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(2.2954693) q[2];
rz(2.7766679) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(2.164771) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038079809) q[0];
sx q[0];
rz(-0.80779034) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(1.6733276) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-3.0751244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9212355) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(0.14089091) q[0];
rz(-pi) q[1];
rz(2.8861553) q[2];
sx q[2];
rz(-1.0424926) q[2];
sx q[2];
rz(-2.262399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.089515226) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(0.11701028) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9598438) q[3];
sx q[3];
rz(-0.74859339) q[3];
sx q[3];
rz(2.4570217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0299783) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(-2.8288793) q[2];
rz(-0.2615658) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(-0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0161491) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(-3.0545767) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(-0.048390128) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2242369) q[0];
sx q[0];
rz(-2.356987) q[0];
sx q[0];
rz(-2.9600701) q[0];
x q[1];
rz(-0.15026413) q[2];
sx q[2];
rz(-1.74904) q[2];
sx q[2];
rz(1.5706289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4573867) q[1];
sx q[1];
rz(-2.4438639) q[1];
sx q[1];
rz(0.084650234) q[1];
x q[2];
rz(-0.88416962) q[3];
sx q[3];
rz(-2.3034952) q[3];
sx q[3];
rz(2.9626486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71063572) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(-0.44130138) q[2];
rz(0.86853164) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5271673) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(1.6983039) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(1.4220994) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5819477) q[0];
sx q[0];
rz(-1.0302246) q[0];
sx q[0];
rz(-1.7339043) q[0];
x q[1];
rz(-1.7781813) q[2];
sx q[2];
rz(-1.8060914) q[2];
sx q[2];
rz(2.4547142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.369097) q[1];
sx q[1];
rz(-2.9926692) q[1];
sx q[1];
rz(0.61222525) q[1];
rz(-pi) q[2];
rz(-0.57797076) q[3];
sx q[3];
rz(-2.1506566) q[3];
sx q[3];
rz(0.65366369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(-1.5606073) q[2];
rz(-1.8033146) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018589858) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(2.4216968) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(0.24615157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3703281) q[0];
sx q[0];
rz(-1.851558) q[0];
sx q[0];
rz(1.595003) q[0];
rz(-pi) q[1];
rz(-0.75743875) q[2];
sx q[2];
rz(-1.2691109) q[2];
sx q[2];
rz(1.2011004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64061058) q[1];
sx q[1];
rz(-2.3965008) q[1];
sx q[1];
rz(-2.9489345) q[1];
rz(1.9897377) q[3];
sx q[3];
rz(-2.2139858) q[3];
sx q[3];
rz(2.1758428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4973732) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(-2.6175446) q[3];
sx q[3];
rz(-2.7553813) q[3];
sx q[3];
rz(3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84435695) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(2.2204087) q[0];
rz(-1.4211897) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(0.99501077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5000347) q[0];
sx q[0];
rz(-1.7935658) q[0];
sx q[0];
rz(1.1831229) q[0];
rz(-pi) q[1];
rz(-1.424355) q[2];
sx q[2];
rz(-2.7173923) q[2];
sx q[2];
rz(2.7643124) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2240963) q[1];
sx q[1];
rz(-0.69643785) q[1];
sx q[1];
rz(-0.68989474) q[1];
rz(-pi) q[2];
rz(1.1345484) q[3];
sx q[3];
rz(-1.9552257) q[3];
sx q[3];
rz(-0.72641295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93817389) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(-2.7040238) q[2];
rz(2.3214052) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95847982) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(-1.0108277) q[0];
rz(3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(-0.54862499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47893229) q[0];
sx q[0];
rz(-2.9629483) q[0];
sx q[0];
rz(1.7500072) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8993511) q[2];
sx q[2];
rz(-1.9827843) q[2];
sx q[2];
rz(-0.054445353) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8884436) q[1];
sx q[1];
rz(-2.4427919) q[1];
sx q[1];
rz(-2.9177925) q[1];
x q[2];
rz(1.6727499) q[3];
sx q[3];
rz(-1.6540308) q[3];
sx q[3];
rz(0.038681313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34714547) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(2.6084206) q[2];
rz(-2.063607) q[3];
sx q[3];
rz(-2.2205133) q[3];
sx q[3];
rz(0.50895154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0549523) q[0];
sx q[0];
rz(-2.7821879) q[0];
sx q[0];
rz(2.3593498) q[0];
rz(-3.0714463) q[1];
sx q[1];
rz(-0.47824305) q[1];
sx q[1];
rz(-0.22629647) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4248764) q[0];
sx q[0];
rz(-0.098669395) q[0];
sx q[0];
rz(-0.77186056) q[0];
x q[1];
rz(-2.2489297) q[2];
sx q[2];
rz(-1.423866) q[2];
sx q[2];
rz(-0.28444296) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6739686) q[1];
sx q[1];
rz(-1.4611164) q[1];
sx q[1];
rz(-2.6024271) q[1];
rz(-1.6582011) q[3];
sx q[3];
rz(-2.8897396) q[3];
sx q[3];
rz(2.9853068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0349064) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.7314343) q[2];
rz(-2.8790224) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(3.0866289) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(0.18375272) q[0];
rz(2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(-2.5180838) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348891) q[0];
sx q[0];
rz(-1.1109795) q[0];
sx q[0];
rz(-1.9135273) q[0];
rz(0.97098668) q[2];
sx q[2];
rz(-2.1040593) q[2];
sx q[2];
rz(0.87109921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7694313) q[1];
sx q[1];
rz(-1.6767394) q[1];
sx q[1];
rz(-1.547936) q[1];
rz(-pi) q[2];
rz(1.9454089) q[3];
sx q[3];
rz(-1.5419898) q[3];
sx q[3];
rz(2.0458986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(3.1039216) q[2];
rz(-2.7231349) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(0.17879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(1.2757318) q[2];
sx q[2];
rz(-2.0841334) q[2];
sx q[2];
rz(1.1708553) q[2];
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
