OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(4.0840277) q[1];
sx q[1];
rz(10.066636) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9489884) q[0];
sx q[0];
rz(-1.3697213) q[0];
sx q[0];
rz(1.6866392) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4840706) q[2];
sx q[2];
rz(-3.1397592) q[2];
sx q[2];
rz(2.4065903) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.042965021) q[1];
sx q[1];
rz(-2.7471099) q[1];
sx q[1];
rz(0.41542713) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1987207) q[3];
sx q[3];
rz(-2.298405) q[3];
sx q[3];
rz(0.6641297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1610819) q[2];
sx q[2];
rz(-1.0970205) q[2];
sx q[2];
rz(-2.3316627) q[2];
rz(0.03446456) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(-0.1035498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.506839) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(0.65938812) q[0];
rz(-1.4393282) q[1];
sx q[1];
rz(-1.6292452) q[1];
sx q[1];
rz(-0.48092458) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5372807) q[0];
sx q[0];
rz(-1.7780281) q[0];
sx q[0];
rz(1.1234493) q[0];
x q[1];
rz(-2.2873053) q[2];
sx q[2];
rz(-2.169974) q[2];
sx q[2];
rz(2.6024659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2175423) q[1];
sx q[1];
rz(-2.9311507) q[1];
sx q[1];
rz(-1.7173052) q[1];
rz(-pi) q[2];
x q[2];
rz(1.102785) q[3];
sx q[3];
rz(-1.0429405) q[3];
sx q[3];
rz(0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-2.3048293) q[2];
sx q[2];
rz(0.62977201) q[2];
rz(1.1791641) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(-1.111697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513007) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(2.1731398) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.7882971) q[1];
sx q[1];
rz(2.5586939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4752858) q[0];
sx q[0];
rz(-2.9687442) q[0];
sx q[0];
rz(-2.8768507) q[0];
x q[1];
rz(1.8331362) q[2];
sx q[2];
rz(-2.3916187) q[2];
sx q[2];
rz(1.394608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.041545871) q[1];
sx q[1];
rz(-2.1963781) q[1];
sx q[1];
rz(0.85432521) q[1];
x q[2];
rz(-0.85286136) q[3];
sx q[3];
rz(-2.7282469) q[3];
sx q[3];
rz(2.1242382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53050238) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(0.21406847) q[2];
rz(-0.30238447) q[3];
sx q[3];
rz(-2.3979135) q[3];
sx q[3];
rz(-0.42588699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18836235) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(-1.0444214) q[0];
rz(2.2638679) q[1];
sx q[1];
rz(-1.5526155) q[1];
sx q[1];
rz(-2.9728319) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8498189) q[0];
sx q[0];
rz(-2.5892604) q[0];
sx q[0];
rz(-2.0122347) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9364768) q[2];
sx q[2];
rz(-1.2502708) q[2];
sx q[2];
rz(-1.4687302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9738439) q[1];
sx q[1];
rz(-2.1653321) q[1];
sx q[1];
rz(-1.3092625) q[1];
x q[2];
rz(-3.0996347) q[3];
sx q[3];
rz(-2.2760512) q[3];
sx q[3];
rz(-2.369997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7793444) q[2];
sx q[2];
rz(-1.8579973) q[2];
sx q[2];
rz(2.0812422) q[2];
rz(2.849071) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(0.084107548) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5523858) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(-0.86828434) q[0];
rz(2.5153416) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(1.3583604) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6407961) q[0];
sx q[0];
rz(-1.2828886) q[0];
sx q[0];
rz(-1.663371) q[0];
rz(2.6399459) q[2];
sx q[2];
rz(-2.7124462) q[2];
sx q[2];
rz(2.8535064) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7029801) q[1];
sx q[1];
rz(-1.7963406) q[1];
sx q[1];
rz(1.3425407) q[1];
rz(-pi) q[2];
x q[2];
rz(1.129398) q[3];
sx q[3];
rz(-2.1611161) q[3];
sx q[3];
rz(0.92720997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9354349) q[2];
sx q[2];
rz(-0.81659603) q[2];
sx q[2];
rz(2.6182776) q[2];
rz(-2.7715136) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10941457) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(0.35414645) q[0];
rz(2.1996563) q[1];
sx q[1];
rz(-1.6862005) q[1];
sx q[1];
rz(-1.3166434) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958651) q[0];
sx q[0];
rz(-2.6168129) q[0];
sx q[0];
rz(-2.9324233) q[0];
x q[1];
rz(-3.1194341) q[2];
sx q[2];
rz(-2.5922734) q[2];
sx q[2];
rz(-1.4469128) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.020654708) q[1];
sx q[1];
rz(-1.6319425) q[1];
sx q[1];
rz(2.967359) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8527669) q[3];
sx q[3];
rz(-0.94579711) q[3];
sx q[3];
rz(0.19240141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3256623) q[2];
sx q[2];
rz(-2.7553835) q[2];
sx q[2];
rz(-0.33829921) q[2];
rz(2.6650688) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4395831) q[0];
sx q[0];
rz(-1.5832573) q[0];
sx q[0];
rz(0.53771341) q[0];
rz(-1.733755) q[1];
sx q[1];
rz(-0.45999637) q[1];
sx q[1];
rz(-0.62526155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940373) q[0];
sx q[0];
rz(-0.68935822) q[0];
sx q[0];
rz(-2.5563142) q[0];
rz(-3.0187716) q[2];
sx q[2];
rz(-2.5993957) q[2];
sx q[2];
rz(1.922883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-9*pi/10) q[1];
sx q[1];
rz(-2.4971967) q[1];
sx q[1];
rz(0.36233904) q[1];
rz(1.9373364) q[3];
sx q[3];
rz(-0.82895422) q[3];
sx q[3];
rz(2.3693565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.16767804) q[2];
sx q[2];
rz(-0.46678552) q[2];
sx q[2];
rz(1.3803049) q[2];
rz(-0.4365094) q[3];
sx q[3];
rz(-1.0218388) q[3];
sx q[3];
rz(-2.4280587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.9647144) q[0];
sx q[0];
rz(-2.5202993) q[0];
sx q[0];
rz(-2.6253413) q[0];
rz(-2.3618354) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(-2.0947184) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488432) q[0];
sx q[0];
rz(-1.7147281) q[0];
sx q[0];
rz(0.30152614) q[0];
rz(2.5604232) q[2];
sx q[2];
rz(-2.0072674) q[2];
sx q[2];
rz(-1.798686) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88491625) q[1];
sx q[1];
rz(-0.11029989) q[1];
sx q[1];
rz(1.0156471) q[1];
x q[2];
rz(0.47886301) q[3];
sx q[3];
rz(-1.2148464) q[3];
sx q[3];
rz(1.009481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20389916) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(2.1843809) q[2];
rz(-2.8474478) q[3];
sx q[3];
rz(-2.011994) q[3];
sx q[3];
rz(2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1421563) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(-0.18705046) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-0.43771935) q[1];
sx q[1];
rz(-1.4923219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0566643) q[0];
sx q[0];
rz(-1.6414375) q[0];
sx q[0];
rz(1.4991374) q[0];
x q[1];
rz(1.9052986) q[2];
sx q[2];
rz(-1.750947) q[2];
sx q[2];
rz(-2.1366773) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88384151) q[1];
sx q[1];
rz(-1.7353936) q[1];
sx q[1];
rz(1.1562892) q[1];
x q[2];
rz(1.4142597) q[3];
sx q[3];
rz(-1.6633342) q[3];
sx q[3];
rz(0.7660367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46818587) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(-0.51472384) q[3];
sx q[3];
rz(-2.616021) q[3];
sx q[3];
rz(1.2334067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24432261) q[0];
sx q[0];
rz(-1.6615302) q[0];
sx q[0];
rz(2.3829714) q[0];
rz(-1.9482535) q[1];
sx q[1];
rz(-1.9494282) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0139931) q[0];
sx q[0];
rz(-1.7993449) q[0];
sx q[0];
rz(2.9319113) q[0];
x q[1];
rz(-3.0994268) q[2];
sx q[2];
rz(-1.773196) q[2];
sx q[2];
rz(0.037994904) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7566061) q[1];
sx q[1];
rz(-0.95545095) q[1];
sx q[1];
rz(0.032757515) q[1];
rz(2.2281353) q[3];
sx q[3];
rz(-2.6413915) q[3];
sx q[3];
rz(-0.72266912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54185581) q[2];
sx q[2];
rz(-2.2208322) q[2];
sx q[2];
rz(0.80121458) q[2];
rz(0.93252212) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(-2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5781317) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(-0.32147944) q[1];
sx q[1];
rz(-0.97578661) q[1];
sx q[1];
rz(-1.6937561) q[1];
rz(-1.6952487) q[2];
sx q[2];
rz(-1.4067408) q[2];
sx q[2];
rz(-0.066700145) q[2];
rz(2.911261) q[3];
sx q[3];
rz(-1.8445704) q[3];
sx q[3];
rz(-0.45382378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
