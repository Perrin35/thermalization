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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(-2.057743) q[0];
rz(2.2639182) q[1];
sx q[1];
rz(4.35507) q[1];
sx q[1];
rz(9.9014643) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2798729) q[0];
sx q[0];
rz(-1.3127767) q[0];
sx q[0];
rz(-1.1610384) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79390766) q[2];
sx q[2];
rz(-2.53144) q[2];
sx q[2];
rz(0.33294233) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.013416524) q[1];
sx q[1];
rz(-0.25826301) q[1];
sx q[1];
rz(1.5034666) q[1];
rz(1.7909369) q[3];
sx q[3];
rz(-1.1486037) q[3];
sx q[3];
rz(3.0289502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.013926355) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(2.9909383) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-0.25130513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0852614) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(-0.37013176) q[0];
rz(-2.6689957) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(-2.1841689) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3963378) q[0];
sx q[0];
rz(-2.6364713) q[0];
sx q[0];
rz(-1.4046679) q[0];
rz(-pi) q[1];
rz(1.2217058) q[2];
sx q[2];
rz(-2.6654976) q[2];
sx q[2];
rz(-2.8698336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53471662) q[1];
sx q[1];
rz(-0.83321111) q[1];
sx q[1];
rz(1.6859158) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4647275) q[3];
sx q[3];
rz(-1.3086623) q[3];
sx q[3];
rz(-1.8755975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1364253) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(1.0235419) q[2];
rz(-1.7510022) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(-1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904974) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(2.5787831) q[0];
rz(1.5142745) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-0.079158457) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1437627) q[0];
sx q[0];
rz(-1.0799066) q[0];
sx q[0];
rz(1.2116371) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5299893) q[2];
sx q[2];
rz(-1.0804324) q[2];
sx q[2];
rz(2.9029569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.99874931) q[1];
sx q[1];
rz(-0.85341893) q[1];
sx q[1];
rz(-0.99695506) q[1];
rz(-pi) q[2];
rz(0.8757365) q[3];
sx q[3];
rz(-2.1453224) q[3];
sx q[3];
rz(0.72573001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2341653) q[2];
sx q[2];
rz(-1.558446) q[2];
sx q[2];
rz(-1.8094212) q[2];
rz(-0.35587707) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-1.1112777) q[0];
rz(-0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(2.8880602) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64719114) q[0];
sx q[0];
rz(-1.8888942) q[0];
sx q[0];
rz(2.6107016) q[0];
rz(-pi) q[1];
rz(0.16040921) q[2];
sx q[2];
rz(-1.6840076) q[2];
sx q[2];
rz(-0.28217523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82099709) q[1];
sx q[1];
rz(-1.7554469) q[1];
sx q[1];
rz(-2.3750633) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0757789) q[3];
sx q[3];
rz(-1.2974713) q[3];
sx q[3];
rz(-1.2697288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(1.4595855) q[2];
rz(0.5736351) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(-3.0972163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(-1.4264872) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41353696) q[0];
sx q[0];
rz(-1.1794834) q[0];
sx q[0];
rz(1.1831468) q[0];
rz(2.5002062) q[2];
sx q[2];
rz(-1.8872627) q[2];
sx q[2];
rz(-3.0015025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0795143) q[1];
sx q[1];
rz(-2.8833564) q[1];
sx q[1];
rz(-2.9596427) q[1];
x q[2];
rz(-2.9699202) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(-1.6115007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9023989) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(1.1725461) q[2];
rz(-3.0555365) q[3];
sx q[3];
rz(-2.1890169) q[3];
sx q[3];
rz(2.9854767) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657144) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(2.8942096) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(0.47201306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9111689) q[0];
sx q[0];
rz(-2.2448475) q[0];
sx q[0];
rz(-2.0087025) q[0];
rz(1.3026779) q[2];
sx q[2];
rz(-0.52137085) q[2];
sx q[2];
rz(0.15586317) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5322383) q[1];
sx q[1];
rz(-1.8731469) q[1];
sx q[1];
rz(-1.3855893) q[1];
rz(0.21228822) q[3];
sx q[3];
rz(-2.1822896) q[3];
sx q[3];
rz(-1.4820198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(1.3989353) q[2];
rz(-1.4553962) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(-2.5808891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6550734) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(2.6229677) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(-1.8124883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5547253) q[0];
sx q[0];
rz(-1.9828078) q[0];
sx q[0];
rz(2.7464965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05032866) q[2];
sx q[2];
rz(-1.2939046) q[2];
sx q[2];
rz(-2.8406539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3336528) q[1];
sx q[1];
rz(-0.83631714) q[1];
sx q[1];
rz(-3.0972788) q[1];
rz(1.5409327) q[3];
sx q[3];
rz(-2.5549508) q[3];
sx q[3];
rz(-1.7778974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3507877) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(-3.1034071) q[2];
rz(2.3112678) q[3];
sx q[3];
rz(-1.7540951) q[3];
sx q[3];
rz(3.0939046) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4627948) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(-1.3702673) q[0];
rz(-0.38617745) q[1];
sx q[1];
rz(-1.736707) q[1];
sx q[1];
rz(0.6699627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46028194) q[0];
sx q[0];
rz(-0.71397266) q[0];
sx q[0];
rz(-2.332469) q[0];
rz(-0.092548056) q[2];
sx q[2];
rz(-0.71882788) q[2];
sx q[2];
rz(-2.3859442) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.252287) q[1];
sx q[1];
rz(-2.5355967) q[1];
sx q[1];
rz(1.7701946) q[1];
rz(-pi) q[2];
rz(1.7827598) q[3];
sx q[3];
rz(-1.6073213) q[3];
sx q[3];
rz(-1.9741457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92114821) q[2];
sx q[2];
rz(-2.8577652) q[2];
sx q[2];
rz(1.05668) q[2];
rz(1.4165261) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(-1.8511124) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0678299) q[0];
sx q[0];
rz(-1.0404328) q[0];
sx q[0];
rz(1.006806) q[0];
rz(2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(-0.83008343) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070979764) q[0];
sx q[0];
rz(-1.5375397) q[0];
sx q[0];
rz(-0.279056) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30980457) q[2];
sx q[2];
rz(-2.5199157) q[2];
sx q[2];
rz(-2.4182188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3398847) q[1];
sx q[1];
rz(-0.23708585) q[1];
sx q[1];
rz(0.70357844) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1277783) q[3];
sx q[3];
rz(-1.5823871) q[3];
sx q[3];
rz(3.1168337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81862187) q[2];
sx q[2];
rz(-1.8399723) q[2];
sx q[2];
rz(-2.0446365) q[2];
rz(1.8638301) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.291236) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(0.37503234) q[0];
rz(3.0229783) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(-0.92432252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40709201) q[0];
sx q[0];
rz(-1.0043) q[0];
sx q[0];
rz(1.2248216) q[0];
rz(-pi) q[1];
rz(-2.0807939) q[2];
sx q[2];
rz(-0.96335232) q[2];
sx q[2];
rz(0.80112544) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3872747) q[1];
sx q[1];
rz(-1.5994834) q[1];
sx q[1];
rz(2.7806675) q[1];
rz(1.9162991) q[3];
sx q[3];
rz(-1.9959108) q[3];
sx q[3];
rz(2.2540852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88860005) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(0.89317733) q[2];
rz(-2.0231953) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(0.18192667) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(-0.091808783) q[2];
sx q[2];
rz(-1.257007) q[2];
sx q[2];
rz(-1.9329482) q[2];
rz(-2.6245194) q[3];
sx q[3];
rz(-1.9319677) q[3];
sx q[3];
rz(-2.464184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
