OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1971561) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(0.42874254) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6719195) q[2];
sx q[2];
rz(-0.28684068) q[2];
sx q[2];
rz(-1.4925721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9285779) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(-1.0282474) q[1];
rz(-pi) q[2];
rz(3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(-1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-0.86581725) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(3.1153733) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(2.1781133) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5945157) q[0];
sx q[0];
rz(-1.5730968) q[0];
sx q[0];
rz(0.95675795) q[0];
x q[1];
rz(-0.33032592) q[2];
sx q[2];
rz(-1.0268372) q[2];
sx q[2];
rz(-0.75737539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0429093) q[1];
sx q[1];
rz(-0.60815647) q[1];
sx q[1];
rz(0.98867464) q[1];
rz(2.8951737) q[3];
sx q[3];
rz(-1.911947) q[3];
sx q[3];
rz(0.44647549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2117675) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0021734) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(2.9234773) q[0];
x q[1];
rz(-2.838344) q[2];
sx q[2];
rz(-1.5504642) q[2];
sx q[2];
rz(-3.0603527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8078976) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(1.3586504) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79389823) q[3];
sx q[3];
rz(-2.5069935) q[3];
sx q[3];
rz(-1.3018228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(-0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.8428615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239149) q[0];
sx q[0];
rz(-2.0664584) q[0];
sx q[0];
rz(2.8850874) q[0];
rz(-pi) q[1];
rz(-1.5694593) q[2];
sx q[2];
rz(-2.3751811) q[2];
sx q[2];
rz(-0.02174755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.270077) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(2.1427076) q[1];
x q[2];
rz(-3.0292547) q[3];
sx q[3];
rz(-1.897052) q[3];
sx q[3];
rz(2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5359042) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(-0.38468012) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6032747) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(0.2968266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332069) q[0];
sx q[0];
rz(-1.0002245) q[0];
sx q[0];
rz(-1.617336) q[0];
rz(0.74283959) q[2];
sx q[2];
rz(-1.683871) q[2];
sx q[2];
rz(-2.6807221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0592812) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(2.6962198) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4818707) q[3];
sx q[3];
rz(-2.2825135) q[3];
sx q[3];
rz(-2.5951648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(-0.39247593) q[2];
rz(-1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(-0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(2.3902067) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-0.60633916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0814708) q[0];
sx q[0];
rz(-1.6166286) q[0];
sx q[0];
rz(-1.5541374) q[0];
rz(-2.2491127) q[2];
sx q[2];
rz(-1.9024897) q[2];
sx q[2];
rz(-1.4075556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6341083) q[1];
sx q[1];
rz(-1.023523) q[1];
sx q[1];
rz(2.9299111) q[1];
rz(-1.7883676) q[3];
sx q[3];
rz(-1.4975582) q[3];
sx q[3];
rz(-0.53461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63885826) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-0.10425723) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095603) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-0.50810057) q[0];
rz(1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.761844) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(2.5248812) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5825282) q[2];
sx q[2];
rz(-1.3571697) q[2];
sx q[2];
rz(-0.27660433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0223169) q[1];
sx q[1];
rz(-1.9609309) q[1];
sx q[1];
rz(1.8632061) q[1];
rz(-pi) q[2];
rz(-0.64738691) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(-2.2341773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4247894) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(1.75288) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213288) q[0];
sx q[0];
rz(-1.5409894) q[0];
sx q[0];
rz(2.1304312) q[0];
x q[1];
rz(0.52877229) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(0.96226529) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4365303) q[1];
sx q[1];
rz(-3.0490746) q[1];
sx q[1];
rz(0.1235048) q[1];
rz(-pi) q[2];
rz(-1.1561398) q[3];
sx q[3];
rz(-2.4932043) q[3];
sx q[3];
rz(0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1887112) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700478) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(0.55823278) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687913) q[0];
sx q[0];
rz(-1.9580012) q[0];
sx q[0];
rz(2.0331435) q[0];
rz(-3.0874861) q[2];
sx q[2];
rz(-1.6781085) q[2];
sx q[2];
rz(-3.0215614) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0135632) q[1];
sx q[1];
rz(-0.82847825) q[1];
sx q[1];
rz(3.0548884) q[1];
rz(-2.6551412) q[3];
sx q[3];
rz(-1.4064186) q[3];
sx q[3];
rz(2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-2.9294087) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(-2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(0.5232946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4921724) q[0];
sx q[0];
rz(-1.5086552) q[0];
sx q[0];
rz(0.6092351) q[0];
rz(-pi) q[1];
rz(-1.5027572) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(0.21493658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.20269468) q[1];
sx q[1];
rz(-1.2564066) q[1];
sx q[1];
rz(2.7640192) q[1];
rz(-pi) q[2];
rz(-1.3396026) q[3];
sx q[3];
rz(-1.3320859) q[3];
sx q[3];
rz(2.0774487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3045197) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.5291443) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(0.79694637) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];