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
rz(-1.00296) q[0];
sx q[0];
rz(3.6149663) q[0];
sx q[0];
rz(10.511721) q[0];
rz(2.3776157) q[1];
sx q[1];
rz(-2.6978701) q[1];
sx q[1];
rz(0.8313764) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1241123) q[0];
sx q[0];
rz(-2.3907733) q[0];
sx q[0];
rz(-0.73877891) q[0];
rz(-pi) q[1];
rz(2.0152807) q[2];
sx q[2];
rz(-0.51156564) q[2];
sx q[2];
rz(-1.4064521) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30637851) q[1];
sx q[1];
rz(-1.6740546) q[1];
sx q[1];
rz(-1.6914961) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3545019) q[3];
sx q[3];
rz(-3.0648764) q[3];
sx q[3];
rz(1.8767954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6143371) q[2];
sx q[2];
rz(-1.273512) q[2];
sx q[2];
rz(-2.7296076) q[2];
rz(-0.24474239) q[3];
sx q[3];
rz(-2.5089846) q[3];
sx q[3];
rz(2.8542724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43993846) q[0];
sx q[0];
rz(-2.169862) q[0];
sx q[0];
rz(-2.7213851) q[0];
rz(-2.6523051) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(-1.18321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0121883) q[0];
sx q[0];
rz(-0.79372915) q[0];
sx q[0];
rz(2.1994526) q[0];
rz(-pi) q[1];
rz(-1.8941882) q[2];
sx q[2];
rz(-1.2719435) q[2];
sx q[2];
rz(-3.1355203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.816531) q[1];
sx q[1];
rz(-1.1899714) q[1];
sx q[1];
rz(-2.517609) q[1];
rz(-pi) q[2];
rz(-1.7443827) q[3];
sx q[3];
rz(-1.4549482) q[3];
sx q[3];
rz(2.4811452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6861787) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(3.0012644) q[2];
rz(-0.27648196) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(0.28237835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33271933) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(1.3746388) q[0];
rz(2.2287492) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(1.0309781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46488097) q[0];
sx q[0];
rz(-1.8959989) q[0];
sx q[0];
rz(-0.19285481) q[0];
x q[1];
rz(2.1616251) q[2];
sx q[2];
rz(-0.64729948) q[2];
sx q[2];
rz(1.574263) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46315868) q[1];
sx q[1];
rz(-2.7292876) q[1];
sx q[1];
rz(1.5957455) q[1];
rz(-pi) q[2];
rz(0.52521962) q[3];
sx q[3];
rz(-0.37074836) q[3];
sx q[3];
rz(-1.4300089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.012152) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(0.63808179) q[2];
rz(-0.10073999) q[3];
sx q[3];
rz(-2.1295857) q[3];
sx q[3];
rz(0.49873275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8345555) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(1.3822973) q[0];
rz(0.31940976) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(-2.3841948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6721181) q[0];
sx q[0];
rz(-1.4961157) q[0];
sx q[0];
rz(-1.5793708) q[0];
rz(-0.73586936) q[2];
sx q[2];
rz(-0.19093787) q[2];
sx q[2];
rz(-1.5519993) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0590225) q[1];
sx q[1];
rz(-2.3772514) q[1];
sx q[1];
rz(-0.90103006) q[1];
rz(0.69587995) q[3];
sx q[3];
rz(-2.4081488) q[3];
sx q[3];
rz(-2.2242198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3697529) q[2];
sx q[2];
rz(-2.2660793) q[2];
sx q[2];
rz(2.3337951) q[2];
rz(1.4098903) q[3];
sx q[3];
rz(-1.8818703) q[3];
sx q[3];
rz(0.65214777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659502) q[0];
sx q[0];
rz(-2.763651) q[0];
sx q[0];
rz(-0.98528969) q[0];
rz(1.4957042) q[1];
sx q[1];
rz(-1.138405) q[1];
sx q[1];
rz(2.2122502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083223) q[0];
sx q[0];
rz(-0.65289557) q[0];
sx q[0];
rz(-2.6920431) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55589755) q[2];
sx q[2];
rz(-0.68624845) q[2];
sx q[2];
rz(2.2874119) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4080491) q[1];
sx q[1];
rz(-1.0449303) q[1];
sx q[1];
rz(-0.2307363) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1272917) q[3];
sx q[3];
rz(-1.0096483) q[3];
sx q[3];
rz(2.1545269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.052496584) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(1.0864786) q[2];
rz(-0.9440445) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(-1.121678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.1679967) q[0];
sx q[0];
rz(-2.7041628) q[0];
sx q[0];
rz(-0.22015372) q[0];
rz(1.6379697) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(-2.5880609) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93388825) q[0];
sx q[0];
rz(-1.7676864) q[0];
sx q[0];
rz(-1.6879199) q[0];
x q[1];
rz(-1.3936744) q[2];
sx q[2];
rz(-1.9290886) q[2];
sx q[2];
rz(2.9492593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11343918) q[1];
sx q[1];
rz(-0.69481297) q[1];
sx q[1];
rz(-3.0656612) q[1];
x q[2];
rz(2.4559654) q[3];
sx q[3];
rz(-1.3925288) q[3];
sx q[3];
rz(-0.36742154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0963875) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(-2.9333147) q[2];
rz(-2.0221209) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(0.82818762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646362) q[0];
sx q[0];
rz(-1.3883075) q[0];
sx q[0];
rz(-2.0679423) q[0];
rz(3.011009) q[1];
sx q[1];
rz(-1.5675631) q[1];
sx q[1];
rz(-1.0098339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70390095) q[0];
sx q[0];
rz(-2.7006375) q[0];
sx q[0];
rz(-2.6954014) q[0];
rz(-pi) q[1];
rz(-2.9931941) q[2];
sx q[2];
rz(-1.4981393) q[2];
sx q[2];
rz(1.0277445) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1853789) q[1];
sx q[1];
rz(-1.662743) q[1];
sx q[1];
rz(2.5673709) q[1];
rz(-pi) q[2];
rz(-2.6920986) q[3];
sx q[3];
rz(-2.426034) q[3];
sx q[3];
rz(0.64330949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69663557) q[2];
sx q[2];
rz(-2.770165) q[2];
sx q[2];
rz(2.8182287) q[2];
rz(-0.67445406) q[3];
sx q[3];
rz(-2.2489397) q[3];
sx q[3];
rz(-3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0167639) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(-1.0892518) q[0];
rz(-2.543154) q[1];
sx q[1];
rz(-1.8611703) q[1];
sx q[1];
rz(0.79900297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75507823) q[0];
sx q[0];
rz(-1.2909954) q[0];
sx q[0];
rz(-2.8248019) q[0];
rz(2.2369821) q[2];
sx q[2];
rz(-1.6450201) q[2];
sx q[2];
rz(2.7192309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1027229) q[1];
sx q[1];
rz(-0.33867044) q[1];
sx q[1];
rz(-2.6790957) q[1];
x q[2];
rz(2.899762) q[3];
sx q[3];
rz(-1.5655091) q[3];
sx q[3];
rz(1.9092803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4716855) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(1.8890107) q[2];
rz(-1.0505098) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85225409) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(0.5156714) q[0];
rz(-1.6853261) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(-2.8175443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1103678) q[0];
sx q[0];
rz(-0.47981167) q[0];
sx q[0];
rz(-2.5672002) q[0];
rz(-pi) q[1];
rz(0.91941763) q[2];
sx q[2];
rz(-1.2044665) q[2];
sx q[2];
rz(1.9969783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4676592) q[1];
sx q[1];
rz(-2.0083462) q[1];
sx q[1];
rz(-2.9872341) q[1];
rz(-pi) q[2];
rz(2.1284655) q[3];
sx q[3];
rz(-1.5032167) q[3];
sx q[3];
rz(-0.24012071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1590165) q[2];
sx q[2];
rz(-2.0197208) q[2];
sx q[2];
rz(-2.4328361) q[2];
rz(2.3580264) q[3];
sx q[3];
rz(-1.2015899) q[3];
sx q[3];
rz(-1.5172575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6056972) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(-2.9851483) q[0];
rz(0.63977301) q[1];
sx q[1];
rz(-0.6548869) q[1];
sx q[1];
rz(2.1515062) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4626083) q[0];
sx q[0];
rz(-2.1741674) q[0];
sx q[0];
rz(2.4067448) q[0];
rz(-pi) q[1];
rz(1.3285816) q[2];
sx q[2];
rz(-2.7822084) q[2];
sx q[2];
rz(0.49312544) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1875547) q[1];
sx q[1];
rz(-1.0196263) q[1];
sx q[1];
rz(1.7269082) q[1];
rz(0.57164945) q[3];
sx q[3];
rz(-0.23861966) q[3];
sx q[3];
rz(-2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5618374) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(0.61335316) q[2];
rz(2.8776045) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(-2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092125208) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(2.3650305) q[2];
sx q[2];
rz(-0.53755098) q[2];
sx q[2];
rz(2.9992486) q[2];
rz(2.522884) q[3];
sx q[3];
rz(-1.6029463) q[3];
sx q[3];
rz(-1.1986986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
