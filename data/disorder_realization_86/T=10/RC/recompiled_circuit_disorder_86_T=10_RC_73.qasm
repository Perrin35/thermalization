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
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(-1.4821948) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1971561) q[0];
sx q[0];
rz(-1.4154139) q[0];
sx q[0];
rz(2.7128501) q[0];
x q[1];
rz(-1.7035159) q[2];
sx q[2];
rz(-1.8258397) q[2];
sx q[2];
rz(2.1357352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60018051) q[1];
sx q[1];
rz(-1.0804847) q[1];
sx q[1];
rz(2.6580826) q[1];
rz(-pi) q[2];
rz(-3.0406038) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(-1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433379) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(3.1153733) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-0.96347934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026982633) q[0];
sx q[0];
rz(-0.61404213) q[0];
sx q[0];
rz(1.5747889) q[0];
rz(-1.0019148) q[2];
sx q[2];
rz(-1.8520253) q[2];
sx q[2];
rz(2.1525454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0986833) q[1];
sx q[1];
rz(-2.5334362) q[1];
sx q[1];
rz(2.152918) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6271237) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0021734) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(2.9234773) q[0];
rz(-1.5494924) q[2];
sx q[2];
rz(-1.2676123) q[2];
sx q[2];
rz(-1.6456749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8078976) q[1];
sx q[1];
rz(-1.2119319) q[1];
sx q[1];
rz(1.7829423) q[1];
x q[2];
rz(1.0873763) q[3];
sx q[3];
rz(-1.1421575) q[3];
sx q[3];
rz(-2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-2.0836232) q[3];
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
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.8428615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772229) q[0];
sx q[0];
rz(-1.3457314) q[0];
sx q[0];
rz(1.0610915) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0012871731) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(3.121701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58060927) q[1];
sx q[1];
rz(-2.1316075) q[1];
sx q[1];
rz(-0.22052712) q[1];
rz(1.8989765) q[3];
sx q[3];
rz(-1.6771852) q[3];
sx q[3];
rz(2.1387517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(-0.38468012) q[2];
rz(-2.3875333) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.6872905) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-2.8447661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332069) q[0];
sx q[0];
rz(-2.1413681) q[0];
sx q[0];
rz(1.5242566) q[0];
x q[1];
rz(1.7237687) q[2];
sx q[2];
rz(-0.83380552) q[2];
sx q[2];
rz(-1.0066777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0592812) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(0.44537284) q[1];
x q[2];
rz(1.6597219) q[3];
sx q[3];
rz(-0.85907912) q[3];
sx q[3];
rz(-2.5951648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(0.39247593) q[2];
rz(-1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(0.75138599) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-0.60633916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060121814) q[0];
sx q[0];
rz(-1.524964) q[0];
sx q[0];
rz(-1.5541374) q[0];
x q[1];
rz(0.89247993) q[2];
sx q[2];
rz(-1.2391029) q[2];
sx q[2];
rz(-1.734037) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11583466) q[1];
sx q[1];
rz(-0.58287139) q[1];
sx q[1];
rz(-1.238766) q[1];
x q[2];
rz(1.243152) q[3];
sx q[3];
rz(-0.22938211) q[3];
sx q[3];
rz(-2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(1.6566488) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5320324) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(2.3513444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3797487) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(-0.6167114) q[0];
rz(-pi) q[1];
rz(2.9279518) q[2];
sx q[2];
rz(-1.5593312) q[2];
sx q[2];
rz(1.2917047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8040647) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(2.7359664) q[1];
x q[2];
rz(0.4831794) q[3];
sx q[3];
rz(-2.4348767) q[3];
sx q[3];
rz(0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0885075) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(1.3887127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.9430338) q[2];
sx q[2];
rz(2.1793274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3125004) q[1];
sx q[1];
rz(-1.478985) q[1];
sx q[1];
rz(1.5822259) q[1];
rz(1.9854529) q[3];
sx q[3];
rz(-2.4932043) q[3];
sx q[3];
rz(0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9528815) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(-3.0601314) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(-0.55823278) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0575858) q[0];
sx q[0];
rz(-1.9965729) q[0];
sx q[0];
rz(-2.7140679) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0360557) q[2];
sx q[2];
rz(-0.12013398) q[2];
sx q[2];
rz(2.5533954) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12802943) q[1];
sx q[1];
rz(-0.82847825) q[1];
sx q[1];
rz(-3.0548884) q[1];
rz(-1.3853119) q[3];
sx q[3];
rz(-2.0501325) q[3];
sx q[3];
rz(-2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50487173) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(-1.6037534) q[0];
rz(-0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(0.5232946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494203) q[0];
sx q[0];
rz(-1.6329375) q[0];
sx q[0];
rz(-0.6092351) q[0];
rz(2.142749) q[2];
sx q[2];
rz(-1.6076644) q[2];
sx q[2];
rz(1.8429304) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8955252) q[1];
sx q[1];
rz(-1.2125891) q[1];
sx q[1];
rz(1.9073061) q[1];
rz(-0.2449805) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(-2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8158648) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-3.1303828) q[2];
sx q[2];
rz(-1.8335473) q[2];
sx q[2];
rz(2.0469472) q[2];
rz(-1.3310824) q[3];
sx q[3];
rz(-2.3532383) q[3];
sx q[3];
rz(-0.91010703) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
