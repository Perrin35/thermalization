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
rz(-1.085936) q[1];
sx q[1];
rz(-0.50958264) q[1];
sx q[1];
rz(-3.0546319) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50849781) q[0];
sx q[0];
rz(-2.1380391) q[0];
sx q[0];
rz(2.022341) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8463507) q[2];
sx q[2];
rz(-0.59165162) q[2];
sx q[2];
rz(-1.0585143) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1198169) q[1];
sx q[1];
rz(-1.6476742) q[1];
sx q[1];
rz(1.6989811) q[1];
rz(2.5218202) q[3];
sx q[3];
rz(-1.8115145) q[3];
sx q[3];
rz(0.71686577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95734346) q[2];
sx q[2];
rz(-0.49746305) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(2.7033778) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(-0.87619585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0210719) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(0.72129321) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(2.8968107) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0185623) q[0];
sx q[0];
rz(-2.5726312) q[0];
sx q[0];
rz(2.2267692) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7049644) q[2];
sx q[2];
rz(-2.2712913) q[2];
sx q[2];
rz(1.3134522) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2657425) q[1];
sx q[1];
rz(-0.85803723) q[1];
sx q[1];
rz(0.73709388) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4238722) q[3];
sx q[3];
rz(-1.8645446) q[3];
sx q[3];
rz(2.5054641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.81405866) q[2];
sx q[2];
rz(-0.56698292) q[2];
sx q[2];
rz(-0.84612334) q[2];
rz(0.36492473) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(2.164771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1035128) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(2.9981151) q[0];
rz(-1.4682651) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(3.0751244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398427) q[0];
sx q[0];
rz(-2.2494613) q[0];
sx q[0];
rz(1.6854273) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9793545) q[2];
sx q[2];
rz(-0.58149946) q[2];
sx q[2];
rz(-1.7844019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0520774) q[1];
sx q[1];
rz(-1.6680535) q[1];
sx q[1];
rz(3.0245824) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74031728) q[3];
sx q[3];
rz(-1.6941287) q[3];
sx q[3];
rz(-2.3892059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0299783) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(2.8288793) q[2];
rz(2.8800268) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(-0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0161491) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(0.087015986) q[0];
rz(1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(-3.0932025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91735578) q[0];
sx q[0];
rz(-2.356987) q[0];
sx q[0];
rz(-0.18152256) q[0];
rz(-1.751028) q[2];
sx q[2];
rz(-1.7186621) q[2];
sx q[2];
rz(-3.1145873) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1783414) q[1];
sx q[1];
rz(-1.5164485) q[1];
sx q[1];
rz(0.69596325) q[1];
rz(-pi) q[2];
rz(-2.2807924) q[3];
sx q[3];
rz(-1.080092) q[3];
sx q[3];
rz(-1.2482289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4309569) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5271673) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(2.4145678) q[0];
rz(-1.4432888) q[1];
sx q[1];
rz(-1.8836421) q[1];
sx q[1];
rz(-1.4220994) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5819477) q[0];
sx q[0];
rz(-1.0302246) q[0];
sx q[0];
rz(1.4076884) q[0];
rz(1.3634113) q[2];
sx q[2];
rz(-1.8060914) q[2];
sx q[2];
rz(-0.6868785) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15502587) q[1];
sx q[1];
rz(-1.6925214) q[1];
sx q[1];
rz(-1.6568068) q[1];
x q[2];
rz(0.90713769) q[3];
sx q[3];
rz(-1.0961514) q[3];
sx q[3];
rz(-1.2603708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21165851) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(-1.5606073) q[2];
rz(1.3382781) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(-2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1230028) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(-2.4216968) q[0];
rz(0.24049354) q[1];
sx q[1];
rz(-1.980282) q[1];
sx q[1];
rz(-0.24615157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3703281) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(1.5465897) q[0];
rz(1.1661548) q[2];
sx q[2];
rz(-2.2863467) q[2];
sx q[2];
rz(-0.095794769) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5009821) q[1];
sx q[1];
rz(-2.3965008) q[1];
sx q[1];
rz(-0.19265811) q[1];
x q[2];
rz(2.6443374) q[3];
sx q[3];
rz(-0.75102931) q[3];
sx q[3];
rz(-2.8145144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.64421946) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(2.3419044) q[2];
rz(-0.52404809) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(-0.016949765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84435695) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(-2.2204087) q[0];
rz(-1.720403) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5000347) q[0];
sx q[0];
rz(-1.7935658) q[0];
sx q[0];
rz(1.1831229) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.065804577) q[2];
sx q[2];
rz(-1.9901681) q[2];
sx q[2];
rz(2.9247627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9239569) q[1];
sx q[1];
rz(-1.991365) q[1];
sx q[1];
rz(-2.5687508) q[1];
x q[2];
rz(1.1345484) q[3];
sx q[3];
rz(-1.1863669) q[3];
sx q[3];
rz(-2.4151797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2034188) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(0.43756884) q[2];
rz(0.82018745) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(3.0062413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1831128) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(1.0108277) q[0];
rz(0.084835947) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(2.5929677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6626604) q[0];
sx q[0];
rz(-2.9629483) q[0];
sx q[0];
rz(1.3915855) q[0];
rz(1.2422416) q[2];
sx q[2];
rz(-1.1588084) q[2];
sx q[2];
rz(3.0871473) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1774166) q[1];
sx q[1];
rz(-0.89277041) q[1];
sx q[1];
rz(1.3864338) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0579257) q[3];
sx q[3];
rz(-1.6723958) q[3];
sx q[3];
rz(1.5236095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7944472) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(0.53317201) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-2.7821879) q[0];
sx q[0];
rz(0.78224283) q[0];
rz(0.070146322) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(-2.9152962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71671623) q[0];
sx q[0];
rz(-3.0429233) q[0];
sx q[0];
rz(-2.3697321) q[0];
rz(-pi) q[1];
rz(0.18780577) q[2];
sx q[2];
rz(-2.2402798) q[2];
sx q[2];
rz(-1.4037496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8574684) q[1];
sx q[1];
rz(-2.5924666) q[1];
sx q[1];
rz(-2.930307) q[1];
rz(-1.3198648) q[3];
sx q[3];
rz(-1.5490412) q[3];
sx q[3];
rz(1.3298498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(0.26257026) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866289) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(-2.9578399) q[0];
rz(-2.9495268) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(-2.5180838) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348891) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(1.9135273) q[0];
x q[1];
rz(-0.76304014) q[2];
sx q[2];
rz(-2.3614778) q[2];
sx q[2];
rz(-3.0811276) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.19621721) q[1];
sx q[1];
rz(-1.5935285) q[1];
sx q[1];
rz(-0.10597056) q[1];
rz(-1.4922114) q[3];
sx q[3];
rz(-2.7659263) q[3];
sx q[3];
rz(-0.54822719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3925675) q[2];
sx q[2];
rz(-0.6124658) q[2];
sx q[2];
rz(-3.1039216) q[2];
rz(-0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(-2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2035718) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(0.45166311) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(0.47623604) q[2];
sx q[2];
rz(-2.5561437) q[2];
sx q[2];
rz(-2.5249425) q[2];
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
