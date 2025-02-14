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
rz(2.0556567) q[1];
sx q[1];
rz(3.6511753) q[1];
sx q[1];
rz(9.3378172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6330948) q[0];
sx q[0];
rz(-1.0035536) q[0];
sx q[0];
rz(1.1192516) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5702417) q[2];
sx q[2];
rz(-1.4077912) q[2];
sx q[2];
rz(-2.3820419) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.086395212) q[1];
sx q[1];
rz(-2.9922303) q[1];
sx q[1];
rz(1.0284852) q[1];
rz(1.8636892) q[3];
sx q[3];
rz(-2.170142) q[3];
sx q[3];
rz(-1.0224316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95734346) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(2.5498665) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(-0.87619585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0210719) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(1.0907809) q[0];
rz(-2.4202994) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(-2.8968107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0185623) q[0];
sx q[0];
rz(-2.5726312) q[0];
sx q[0];
rz(2.2267692) q[0];
rz(-pi) q[1];
rz(-2.7049644) q[2];
sx q[2];
rz(-0.87030137) q[2];
sx q[2];
rz(1.8281405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9823581) q[1];
sx q[1];
rz(-1.0373384) q[1];
sx q[1];
rz(-0.70833556) q[1];
x q[2];
rz(2.6908633) q[3];
sx q[3];
rz(-2.8141032) q[3];
sx q[3];
rz(2.9779676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(-0.84612334) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(-0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1035128) q[0];
sx q[0];
rz(-0.80779034) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(1.4682651) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(-3.0751244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398427) q[0];
sx q[0];
rz(-0.89213138) q[0];
sx q[0];
rz(-1.4561653) q[0];
x q[1];
rz(-0.25543737) q[2];
sx q[2];
rz(-2.0991) q[2];
sx q[2];
rz(-0.87919368) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.089515226) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(0.11701028) q[1];
rz(-pi) q[2];
rz(-1.4044365) q[3];
sx q[3];
rz(-2.3041953) q[3];
sx q[3];
rz(2.211253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0299783) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(-0.31271333) q[2];
rz(2.8800268) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(2.3436782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(3.0161491) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(-3.0545767) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(-0.048390128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66358405) q[0];
sx q[0];
rz(-0.80246325) q[0];
sx q[0];
rz(-1.3924696) q[0];
rz(0.15026413) q[2];
sx q[2];
rz(-1.3925526) q[2];
sx q[2];
rz(-1.5709637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1783414) q[1];
sx q[1];
rz(-1.6251441) q[1];
sx q[1];
rz(2.4456294) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6137621) q[3];
sx q[3];
rz(-0.9582954) q[3];
sx q[3];
rz(-0.70741725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.71063572) q[2];
sx q[2];
rz(-0.29971665) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61442536) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(0.72702485) q[0];
rz(-1.6983039) q[1];
sx q[1];
rz(-1.8836421) q[1];
sx q[1];
rz(1.4220994) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0682536) q[0];
sx q[0];
rz(-1.4311106) q[0];
sx q[0];
rz(-2.5951067) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.431972) q[2];
sx q[2];
rz(-0.31236744) q[2];
sx q[2];
rz(1.4210977) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9865668) q[1];
sx q[1];
rz(-1.4490713) q[1];
sx q[1];
rz(1.4847859) q[1];
rz(2.234455) q[3];
sx q[3];
rz(-2.0454413) q[3];
sx q[3];
rz(1.8812219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9299341) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(1.5606073) q[2];
rz(1.3382781) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1230028) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(-0.71989584) q[0];
rz(-0.24049354) q[1];
sx q[1];
rz(-1.980282) q[1];
sx q[1];
rz(-2.8954411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.947833) q[0];
sx q[0];
rz(-1.5475377) q[0];
sx q[0];
rz(0.2808397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9754378) q[2];
sx q[2];
rz(-0.85524594) q[2];
sx q[2];
rz(0.095794769) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7604039) q[1];
sx q[1];
rz(-0.84263984) q[1];
sx q[1];
rz(-1.3959753) q[1];
x q[2];
rz(-1.9897377) q[3];
sx q[3];
rz(-0.92760689) q[3];
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
rz(-2.5770598) q[2];
sx q[2];
rz(-0.79968828) q[2];
rz(-2.6175446) q[3];
sx q[3];
rz(-2.7553813) q[3];
sx q[3];
rz(3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(-2.4589296) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.641558) q[0];
sx q[0];
rz(-1.7935658) q[0];
sx q[0];
rz(1.9584697) q[0];
rz(-pi) q[1];
rz(1.424355) q[2];
sx q[2];
rz(-0.42420039) q[2];
sx q[2];
rz(-0.37728024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9239569) q[1];
sx q[1];
rz(-1.991365) q[1];
sx q[1];
rz(2.5687508) q[1];
rz(-pi) q[2];
rz(0.80713804) q[3];
sx q[3];
rz(-2.5684331) q[3];
sx q[3];
rz(1.5218211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93817389) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(-2.7040238) q[2];
rz(-0.82018745) q[3];
sx q[3];
rz(-1.5929675) q[3];
sx q[3];
rz(3.0062413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1831128) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(-1.0108277) q[0];
rz(3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(-0.54862499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2261467) q[0];
sx q[0];
rz(-1.5391162) q[0];
sx q[0];
rz(1.394954) q[0];
rz(-pi) q[1];
rz(-0.63603129) q[2];
sx q[2];
rz(-0.52101982) q[2];
sx q[2];
rz(-0.65083671) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2531491) q[1];
sx q[1];
rz(-2.4427919) q[1];
sx q[1];
rz(-2.9177925) q[1];
rz(0.083666936) q[3];
sx q[3];
rz(-1.6723958) q[3];
sx q[3];
rz(-1.6179832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34714547) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(2.6084206) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-2.2205133) q[3];
sx q[3];
rz(-0.50895154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.47824305) q[1];
sx q[1];
rz(0.22629647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6235038) q[0];
sx q[0];
rz(-1.5020348) q[0];
sx q[0];
rz(-3.070773) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.802472) q[2];
sx q[2];
rz(-0.69139987) q[2];
sx q[2];
rz(1.1065799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28412425) q[1];
sx q[1];
rz(-0.54912607) q[1];
sx q[1];
rz(2.930307) q[1];
rz(1.4833916) q[3];
sx q[3];
rz(-0.25185302) q[3];
sx q[3];
rz(0.15628584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1066863) q[2];
sx q[2];
rz(-0.28768134) q[2];
sx q[2];
rz(-1.4101583) q[2];
rz(-0.26257026) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(3.0098651) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(2.9578399) q[0];
rz(-0.19206583) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(-0.62350887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62105084) q[0];
sx q[0];
rz(-1.2648925) q[0];
sx q[0];
rz(0.48407475) q[0];
rz(-pi) q[1];
rz(-2.5207768) q[2];
sx q[2];
rz(-2.0784272) q[2];
sx q[2];
rz(2.1073282) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3721613) q[1];
sx q[1];
rz(-1.6767394) q[1];
sx q[1];
rz(1.547936) q[1];
rz(-0.030951854) q[3];
sx q[3];
rz(-1.1963468) q[3];
sx q[3];
rz(0.46377814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2035718) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(-0.45166311) q[1];
sx q[1];
rz(-1.3126806) q[1];
sx q[1];
rz(-1.5246593) q[1];
rz(-2.6653566) q[2];
sx q[2];
rz(-2.5561437) q[2];
sx q[2];
rz(-2.5249425) q[2];
rz(1.3814817) q[3];
sx q[3];
rz(-1.2997205) q[3];
sx q[3];
rz(-1.442853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
