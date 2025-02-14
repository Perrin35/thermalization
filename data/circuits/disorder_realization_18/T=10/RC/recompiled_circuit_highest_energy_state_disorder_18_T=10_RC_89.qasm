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
rz(1.7914766) q[0];
sx q[0];
rz(-0.67576367) q[0];
sx q[0];
rz(3.1337373) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(1.6176728) q[1];
sx q[1];
rz(9.6674506) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033256) q[0];
sx q[0];
rz(-1.4446065) q[0];
sx q[0];
rz(0.85889205) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2062293) q[2];
sx q[2];
rz(-1.1375712) q[2];
sx q[2];
rz(1.3490167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2897799) q[1];
sx q[1];
rz(-0.35408005) q[1];
sx q[1];
rz(3.0733527) q[1];
x q[2];
rz(-2.8318066) q[3];
sx q[3];
rz(-2.4521146) q[3];
sx q[3];
rz(-0.30287095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1050038) q[2];
sx q[2];
rz(-1.9052817) q[2];
sx q[2];
rz(-2.429602) q[2];
rz(-2.8365734) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(-3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104093) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(1.7741868) q[0];
rz(-0.039904682) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(2.5423999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0586226) q[0];
sx q[0];
rz(-0.48958594) q[0];
sx q[0];
rz(-0.24607308) q[0];
rz(0.97936432) q[2];
sx q[2];
rz(-2.2970846) q[2];
sx q[2];
rz(-1.756211) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89394648) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(1.9878073) q[1];
x q[2];
rz(-1.2892339) q[3];
sx q[3];
rz(-0.90471876) q[3];
sx q[3];
rz(0.56872845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.082531884) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(1.3085636) q[2];
rz(3.1176873) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(-1.5628975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6858653) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(-0.84079963) q[0];
rz(-2.1127545) q[1];
sx q[1];
rz(-0.40887555) q[1];
sx q[1];
rz(-1.6811446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90636364) q[0];
sx q[0];
rz(-1.1113941) q[0];
sx q[0];
rz(1.4109341) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1991791) q[2];
sx q[2];
rz(-2.016168) q[2];
sx q[2];
rz(-0.98486131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19059316) q[1];
sx q[1];
rz(-1.237932) q[1];
sx q[1];
rz(0.69575633) q[1];
rz(2.5342503) q[3];
sx q[3];
rz(-1.6173956) q[3];
sx q[3];
rz(0.49639116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2568405) q[2];
sx q[2];
rz(-1.2105056) q[2];
sx q[2];
rz(2.9849226) q[2];
rz(-1.5001851) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0750065) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(2.0196041) q[1];
sx q[1];
rz(-2.5009218) q[1];
sx q[1];
rz(1.5544308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6268353) q[0];
sx q[0];
rz(-2.3957402) q[0];
sx q[0];
rz(1.4738073) q[0];
rz(-pi) q[1];
rz(0.3360668) q[2];
sx q[2];
rz(-3.0185351) q[2];
sx q[2];
rz(-0.50257909) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69103855) q[1];
sx q[1];
rz(-0.86133707) q[1];
sx q[1];
rz(2.4616509) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5347363) q[3];
sx q[3];
rz(-1.292406) q[3];
sx q[3];
rz(-0.75137072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88177219) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(-2.1176977) q[2];
rz(-0.77670589) q[3];
sx q[3];
rz(-1.9909765) q[3];
sx q[3];
rz(-2.1385433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.40989947) q[0];
sx q[0];
rz(-2.2875146) q[0];
sx q[0];
rz(2.5140629) q[0];
rz(0.69951406) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(-0.90739179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2058682) q[0];
sx q[0];
rz(-0.75399929) q[0];
sx q[0];
rz(0.058490965) q[0];
rz(-pi) q[1];
rz(0.19030119) q[2];
sx q[2];
rz(-1.7277328) q[2];
sx q[2];
rz(-2.4630594) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4470565) q[1];
sx q[1];
rz(-1.2061822) q[1];
sx q[1];
rz(2.1237808) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058919546) q[3];
sx q[3];
rz(-2.1061961) q[3];
sx q[3];
rz(-0.31188341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0044535) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(0.84490204) q[2];
rz(-2.4431303) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(0.79160488) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38248211) q[0];
sx q[0];
rz(-1.4610721) q[0];
sx q[0];
rz(1.2680898) q[0];
rz(-1.5269439) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(0.39438927) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1295107) q[0];
sx q[0];
rz(-2.2144268) q[0];
sx q[0];
rz(-0.056931007) q[0];
rz(3.0494681) q[2];
sx q[2];
rz(-0.73706223) q[2];
sx q[2];
rz(1.2110405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2339237) q[1];
sx q[1];
rz(-0.65429293) q[1];
sx q[1];
rz(-0.073542417) q[1];
rz(1.4930356) q[3];
sx q[3];
rz(-1.2111411) q[3];
sx q[3];
rz(-0.80870562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4322728) q[2];
sx q[2];
rz(-0.061881438) q[2];
sx q[2];
rz(2.5800932) q[2];
rz(-2.0105441) q[3];
sx q[3];
rz(-2.3384422) q[3];
sx q[3];
rz(-2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5973709) q[0];
sx q[0];
rz(-0.18246305) q[0];
sx q[0];
rz(-0.23319787) q[0];
rz(0.11416642) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(0.17359576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9812935) q[0];
sx q[0];
rz(-0.7795142) q[0];
sx q[0];
rz(1.0407838) q[0];
x q[1];
rz(2.0162705) q[2];
sx q[2];
rz(-1.2230754) q[2];
sx q[2];
rz(-1.7537774) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0376589) q[1];
sx q[1];
rz(-1.3500431) q[1];
sx q[1];
rz(-2.9751865) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0775142) q[3];
sx q[3];
rz(-0.18905583) q[3];
sx q[3];
rz(2.1259675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0357828) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(0.85476533) q[2];
rz(3.0586976) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(-0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6524803) q[0];
sx q[0];
rz(-1.2615477) q[0];
sx q[0];
rz(1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.2636355) q[1];
sx q[1];
rz(2.9875535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73418854) q[0];
sx q[0];
rz(-1.2449236) q[0];
sx q[0];
rz(-2.4209321) q[0];
rz(-2.4392022) q[2];
sx q[2];
rz(-2.1419814) q[2];
sx q[2];
rz(2.1399501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8785868) q[1];
sx q[1];
rz(-1.1060426) q[1];
sx q[1];
rz(-0.33585637) q[1];
rz(2.4064111) q[3];
sx q[3];
rz(-2.7594341) q[3];
sx q[3];
rz(2.2459386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68086326) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(0.66780773) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(2.5051129) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77968303) q[0];
sx q[0];
rz(-1.4754262) q[0];
sx q[0];
rz(2.5478126) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.8437754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0827328) q[0];
sx q[0];
rz(-1.4987336) q[0];
sx q[0];
rz(-0.58499344) q[0];
x q[1];
rz(0.62784451) q[2];
sx q[2];
rz(-1.825001) q[2];
sx q[2];
rz(1.1140299) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1480963) q[1];
sx q[1];
rz(-1.2132753) q[1];
sx q[1];
rz(-1.4361678) q[1];
rz(-1.7421977) q[3];
sx q[3];
rz(-0.77091445) q[3];
sx q[3];
rz(-0.41219974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7456776) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(-1.6321261) q[2];
rz(3.0294026) q[3];
sx q[3];
rz(-2.1144919) q[3];
sx q[3];
rz(1.3794911) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141948) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(-2.8759586) q[0];
rz(2.1521125) q[1];
sx q[1];
rz(-1.8786636) q[1];
sx q[1];
rz(2.8094453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6393747) q[0];
sx q[0];
rz(-1.9042339) q[0];
sx q[0];
rz(0.0040472814) q[0];
rz(-pi) q[1];
rz(0.29075622) q[2];
sx q[2];
rz(-1.0117048) q[2];
sx q[2];
rz(1.5955995) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0451584) q[1];
sx q[1];
rz(-1.2306884) q[1];
sx q[1];
rz(-0.051326871) q[1];
rz(-3.0933558) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(-1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0997448) q[2];
sx q[2];
rz(-1.0724649) q[2];
sx q[2];
rz(2.9912046) q[2];
rz(-2.2930875) q[3];
sx q[3];
rz(-2.7917807) q[3];
sx q[3];
rz(1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0583508) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(1.3311483) q[1];
sx q[1];
rz(-0.7970627) q[1];
sx q[1];
rz(2.0373559) q[1];
rz(1.642557) q[2];
sx q[2];
rz(-2.4007779) q[2];
sx q[2];
rz(3.130198) q[2];
rz(-1.5450618) q[3];
sx q[3];
rz(-2.2774057) q[3];
sx q[3];
rz(-2.0792815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
