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
rz(2.0556567) q[1];
sx q[1];
rz(-2.63201) q[1];
sx q[1];
rz(-0.086960763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6330948) q[0];
sx q[0];
rz(-1.0035536) q[0];
sx q[0];
rz(2.022341) q[0];
x q[1];
rz(0.571351) q[2];
sx q[2];
rz(-1.4077912) q[2];
sx q[2];
rz(2.3820419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6807144) q[1];
sx q[1];
rz(-1.6986004) q[1];
sx q[1];
rz(3.0640814) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5218202) q[3];
sx q[3];
rz(-1.3300782) q[3];
sx q[3];
rz(-2.4247269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1842492) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(2.5498665) q[2];
rz(2.7033778) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0210719) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(-2.8968107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1230303) q[0];
sx q[0];
rz(-0.5689615) q[0];
sx q[0];
rz(-0.91482343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4366283) q[2];
sx q[2];
rz(-0.87030137) q[2];
sx q[2];
rz(1.8281405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3192057) q[1];
sx q[1];
rz(-2.165395) q[1];
sx q[1];
rz(-0.90984224) q[1];
rz(-pi) q[2];
rz(-2.8448288) q[3];
sx q[3];
rz(-1.4302084) q[3];
sx q[3];
rz(-2.1641017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(0.84612334) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.038079809) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(-1.6733276) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-0.066468261) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9212355) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(-0.14089091) q[0];
rz(-pi) q[1];
rz(1.9793545) q[2];
sx q[2];
rz(-0.58149946) q[2];
sx q[2];
rz(1.3571908) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.089515226) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(3.0245824) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9598438) q[3];
sx q[3];
rz(-2.3929993) q[3];
sx q[3];
rz(0.68457097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1116144) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(-2.8288793) q[2];
rz(0.2615658) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(-0.79791445) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0161491) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(-0.087015986) q[0];
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
rz(0.78260471) q[0];
sx q[0];
rz(-1.442897) q[0];
sx q[0];
rz(-0.77632287) q[0];
x q[1];
rz(-0.15026413) q[2];
sx q[2];
rz(-1.3925526) q[2];
sx q[2];
rz(1.5709637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.794487) q[1];
sx q[1];
rz(-0.87606591) q[1];
sx q[1];
rz(1.6415651) q[1];
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
rz(-0.71063572) q[2];
sx q[2];
rz(-0.29971665) q[2];
sx q[2];
rz(0.44130138) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.3683616) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61442536) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(1.4432888) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(-1.4220994) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5819477) q[0];
sx q[0];
rz(-1.0302246) q[0];
sx q[0];
rz(-1.4076884) q[0];
rz(-pi) q[1];
rz(-2.9013394) q[2];
sx q[2];
rz(-1.3692055) q[2];
sx q[2];
rz(-0.83490419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9865668) q[1];
sx q[1];
rz(-1.4490713) q[1];
sx q[1];
rz(-1.6568068) q[1];
rz(-pi) q[2];
x q[2];
rz(2.234455) q[3];
sx q[3];
rz(-1.0961514) q[3];
sx q[3];
rz(-1.8812219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9299341) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(1.5809853) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(-2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-3.1230028) q[0];
sx q[0];
rz(-0.81273166) q[0];
sx q[0];
rz(0.71989584) q[0];
rz(0.24049354) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(0.24615157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6841078) q[0];
sx q[0];
rz(-2.8598169) q[0];
sx q[0];
rz(-3.0578567) q[0];
rz(-0.75743875) q[2];
sx q[2];
rz(-1.8724818) q[2];
sx q[2];
rz(1.9404922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5009821) q[1];
sx q[1];
rz(-2.3965008) q[1];
sx q[1];
rz(2.9489345) q[1];
rz(0.4972552) q[3];
sx q[3];
rz(-0.75102931) q[3];
sx q[3];
rz(2.8145144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64421946) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(-2.3419044) q[2];
rz(-2.6175446) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(0.016949765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84435695) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(-2.2204087) q[0];
rz(1.720403) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667082) q[0];
sx q[0];
rz(-2.6972983) q[0];
sx q[0];
rz(1.0309451) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.424355) q[2];
sx q[2];
rz(-0.42420039) q[2];
sx q[2];
rz(0.37728024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0459208) q[1];
sx q[1];
rz(-1.0532044) q[1];
sx q[1];
rz(1.0817097) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1345484) q[3];
sx q[3];
rz(-1.1863669) q[3];
sx q[3];
rz(-2.4151797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2034188) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(-2.7040238) q[2];
rz(0.82018745) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1831128) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(1.0108277) q[0];
rz(-3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(-2.5929677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91544596) q[0];
sx q[0];
rz(-1.5391162) q[0];
sx q[0];
rz(1.394954) q[0];
rz(-pi) q[1];
rz(0.43253501) q[2];
sx q[2];
rz(-1.2706332) q[2];
sx q[2];
rz(1.6520239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2531491) q[1];
sx q[1];
rz(-0.69880077) q[1];
sx q[1];
rz(-0.22380016) q[1];
rz(-pi) q[2];
rz(0.083666936) q[3];
sx q[3];
rz(-1.6723958) q[3];
sx q[3];
rz(1.5236095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34714547) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(2.6084206) q[2];
rz(-1.0779856) q[3];
sx q[3];
rz(-2.2205133) q[3];
sx q[3];
rz(2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6235038) q[0];
sx q[0];
rz(-1.5020348) q[0];
sx q[0];
rz(-0.070819617) q[0];
x q[1];
rz(1.3391206) q[2];
sx q[2];
rz(-0.69139987) q[2];
sx q[2];
rz(1.1065799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4676241) q[1];
sx q[1];
rz(-1.6804763) q[1];
sx q[1];
rz(0.53916559) q[1];
x q[2];
rz(-1.4833916) q[3];
sx q[3];
rz(-2.8897396) q[3];
sx q[3];
rz(-2.9853068) q[3];
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
rz(-0.26257026) q[3];
sx q[3];
rz(-1.5841443) q[3];
sx q[3];
rz(-3.0098651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(-0.18375272) q[0];
rz(2.9495268) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(2.5180838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348891) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(-1.2280653) q[0];
rz(-pi) q[1];
rz(0.97098668) q[2];
sx q[2];
rz(-1.0375334) q[2];
sx q[2];
rz(2.2704934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9453754) q[1];
sx q[1];
rz(-1.5480642) q[1];
sx q[1];
rz(-3.0356221) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1961837) q[3];
sx q[3];
rz(-1.5996029) q[3];
sx q[3];
rz(-1.0956941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7490251) q[2];
sx q[2];
rz(-0.6124658) q[2];
sx q[2];
rz(0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(0.17879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(-1.2757318) q[2];
sx q[2];
rz(-1.0574592) q[2];
sx q[2];
rz(-1.9707373) q[2];
rz(-1.3814817) q[3];
sx q[3];
rz(-1.8418722) q[3];
sx q[3];
rz(1.6987396) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
