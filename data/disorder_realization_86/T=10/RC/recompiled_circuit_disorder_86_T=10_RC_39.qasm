OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(-1.7295184) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1887814) q[0];
sx q[0];
rz(-0.45438284) q[0];
sx q[0];
rz(-2.7812468) q[0];
x q[1];
rz(1.4380768) q[2];
sx q[2];
rz(-1.3157529) q[2];
sx q[2];
rz(1.0058574) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9285779) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(1.0282474) q[1];
rz(-2.1336018) q[3];
sx q[3];
rz(-1.6562914) q[3];
sx q[3];
rz(0.29719719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98510629) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.2877864) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433379) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(1.6014618) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-2.1781133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5945157) q[0];
sx q[0];
rz(-1.5684959) q[0];
sx q[0];
rz(0.95675795) q[0];
x q[1];
rz(-2.0630402) q[2];
sx q[2];
rz(-0.62765593) q[2];
sx q[2];
rz(0.17222675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0429093) q[1];
sx q[1];
rz(-0.60815647) q[1];
sx q[1];
rz(0.98867464) q[1];
rz(2.1727932) q[3];
sx q[3];
rz(-0.41799823) q[3];
sx q[3];
rz(-2.9434162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(3.0070686) q[2];
rz(2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(2.581596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13941923) q[0];
sx q[0];
rz(-1.9529441) q[0];
sx q[0];
rz(2.9234773) q[0];
x q[1];
rz(-1.5921002) q[2];
sx q[2];
rz(-1.8739803) q[2];
sx q[2];
rz(1.4959178) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8289889) q[1];
sx q[1];
rz(-1.7692411) q[1];
sx q[1];
rz(-0.36638422) q[1];
x q[2];
rz(-2.6651354) q[3];
sx q[3];
rz(-1.1343079) q[3];
sx q[3];
rz(-0.95526327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-2.0836232) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.2987312) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643698) q[0];
sx q[0];
rz(-1.3457314) q[0];
sx q[0];
rz(-1.0610915) q[0];
rz(-pi) q[1];
rz(-1.5694593) q[2];
sx q[2];
rz(-2.3751811) q[2];
sx q[2];
rz(3.1198451) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.270077) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(-2.1427076) q[1];
rz(-pi) q[2];
rz(0.11233791) q[3];
sx q[3];
rz(-1.897052) q[3];
sx q[3];
rz(-0.60409594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(2.3875333) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6032747) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-0.2968266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083858) q[0];
sx q[0];
rz(-2.1413681) q[0];
sx q[0];
rz(1.617336) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74283959) q[2];
sx q[2];
rz(-1.683871) q[2];
sx q[2];
rz(0.46087056) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8991124) q[1];
sx q[1];
rz(-2.50933) q[1];
sx q[1];
rz(2.2805023) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0389901) q[3];
sx q[3];
rz(-0.71628621) q[3];
sx q[3];
rz(2.7308381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(-2.3902067) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(2.5352535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4088926) q[0];
sx q[0];
rz(-0.048763976) q[0];
sx q[0];
rz(-0.34838895) q[0];
x q[1];
rz(1.0688071) q[2];
sx q[2];
rz(-0.74337465) q[2];
sx q[2];
rz(2.5943287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6341083) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(-2.9299111) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0665928) q[3];
sx q[3];
rz(-1.3538176) q[3];
sx q[3];
rz(-2.0892339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(1.6566488) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(-3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(-0.79024822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3797487) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(-0.6167114) q[0];
rz(-pi) q[1];
rz(-1.5590645) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(-2.8649883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33752791) q[1];
sx q[1];
rz(-1.3009326) q[1];
sx q[1];
rz(-0.40562628) q[1];
x q[2];
rz(2.4942057) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(-2.2341773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(-2.0943663) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(1.75288) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7097276) q[0];
sx q[0];
rz(-1.0114397) q[0];
sx q[0];
rz(0.035168408) q[0];
rz(-pi) q[1];
rz(1.1461166) q[2];
sx q[2];
rz(-2.0600024) q[2];
sx q[2];
rz(0.81791544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3988406) q[1];
sx q[1];
rz(-1.5594149) q[1];
sx q[1];
rz(-0.091817261) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9854529) q[3];
sx q[3];
rz(-0.64838833) q[3];
sx q[3];
rz(0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1887112) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(-2.2593373) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.4260938) q[0];
rz(-0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-0.55823278) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0575858) q[0];
sx q[0];
rz(-1.9965729) q[0];
sx q[0];
rz(2.7140679) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0874861) q[2];
sx q[2];
rz(-1.6781085) q[2];
sx q[2];
rz(3.0215614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3840752) q[1];
sx q[1];
rz(-1.6346524) q[1];
sx q[1];
rz(2.3149895) q[1];
rz(2.6551412) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(2.6182981) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035318035) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(-1.6465228) q[0];
rz(-pi) q[1];
rz(1.6388355) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(0.21493658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8955252) q[1];
sx q[1];
rz(-1.2125891) q[1];
sx q[1];
rz(-1.9073061) q[1];
rz(-pi) q[2];
rz(-2.8966122) q[3];
sx q[3];
rz(-1.7953201) q[3];
sx q[3];
rz(0.45104879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(-2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8158648) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.4670463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(1.5291443) q[2];
sx q[2];
rz(-2.8786082) q[2];
sx q[2];
rz(-1.0515121) q[2];
rz(-1.8105103) q[3];
sx q[3];
rz(-0.78835434) q[3];
sx q[3];
rz(2.2314856) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
