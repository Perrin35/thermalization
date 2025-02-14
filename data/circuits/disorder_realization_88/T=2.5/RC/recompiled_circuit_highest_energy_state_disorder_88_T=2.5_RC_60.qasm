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
rz(-3.1336194) q[0];
sx q[0];
rz(-2.70533) q[0];
sx q[0];
rz(-2.1313957) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(-2.3354524) q[1];
sx q[1];
rz(-1.0239209) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9925594) q[0];
sx q[0];
rz(-1.5776538) q[0];
sx q[0];
rz(2.0586781) q[0];
rz(-pi) q[1];
rz(0.3918565) q[2];
sx q[2];
rz(-1.7855682) q[2];
sx q[2];
rz(1.2607167) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1271474) q[1];
sx q[1];
rz(-1.5902441) q[1];
sx q[1];
rz(-1.4102742) q[1];
rz(-pi) q[2];
rz(2.8918582) q[3];
sx q[3];
rz(-1.6563376) q[3];
sx q[3];
rz(1.6289323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94850928) q[2];
sx q[2];
rz(-0.32863363) q[2];
sx q[2];
rz(2.9060717) q[2];
rz(-2.113302) q[3];
sx q[3];
rz(-1.2582658) q[3];
sx q[3];
rz(-1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.704772) q[0];
sx q[0];
rz(-1.3587767) q[0];
sx q[0];
rz(0.92192465) q[0];
rz(3.0251265) q[1];
sx q[1];
rz(-1.7200229) q[1];
sx q[1];
rz(0.88752735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5796611) q[0];
sx q[0];
rz(-1.4401739) q[0];
sx q[0];
rz(-1.3829253) q[0];
rz(-pi) q[1];
rz(-2.791128) q[2];
sx q[2];
rz(-1.5157358) q[2];
sx q[2];
rz(2.6907235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4749467) q[1];
sx q[1];
rz(-1.4578867) q[1];
sx q[1];
rz(0.571588) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13963066) q[3];
sx q[3];
rz(-0.88362304) q[3];
sx q[3];
rz(0.72134127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5827215) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(0.23898807) q[2];
rz(0.72238266) q[3];
sx q[3];
rz(-2.3014849) q[3];
sx q[3];
rz(1.9766138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31579414) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(2.353299) q[0];
rz(-1.3901419) q[1];
sx q[1];
rz(-2.7056521) q[1];
sx q[1];
rz(1.699126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4460488) q[0];
sx q[0];
rz(-2.9835577) q[0];
sx q[0];
rz(-1.1084854) q[0];
rz(-pi) q[1];
rz(0.29533889) q[2];
sx q[2];
rz(-1.9304233) q[2];
sx q[2];
rz(-2.4725898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6953814) q[1];
sx q[1];
rz(-1.3331116) q[1];
sx q[1];
rz(-2.2878245) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0736819) q[3];
sx q[3];
rz(-1.4137795) q[3];
sx q[3];
rz(-2.5105421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73551377) q[2];
sx q[2];
rz(-1.1512076) q[2];
sx q[2];
rz(-2.8847412) q[2];
rz(-0.86366051) q[3];
sx q[3];
rz(-1.0830027) q[3];
sx q[3];
rz(-0.5736205) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094512016) q[0];
sx q[0];
rz(-3.1086476) q[0];
sx q[0];
rz(-1.5091913) q[0];
rz(-0.31653658) q[1];
sx q[1];
rz(-1.5455952) q[1];
sx q[1];
rz(2.1260156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8636531) q[0];
sx q[0];
rz(-0.67623338) q[0];
sx q[0];
rz(-0.77267026) q[0];
x q[1];
rz(-1.847083) q[2];
sx q[2];
rz(-0.7965379) q[2];
sx q[2];
rz(-0.73819064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16962651) q[1];
sx q[1];
rz(-1.3463273) q[1];
sx q[1];
rz(1.1826993) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82364453) q[3];
sx q[3];
rz(-1.6948786) q[3];
sx q[3];
rz(-1.9923576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.29869002) q[2];
sx q[2];
rz(-1.0516673) q[2];
sx q[2];
rz(1.451937) q[2];
rz(1.2339833) q[3];
sx q[3];
rz(-2.0472287) q[3];
sx q[3];
rz(-0.88172495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3877617) q[0];
sx q[0];
rz(-1.3085288) q[0];
sx q[0];
rz(-2.0552788) q[0];
rz(1.7408675) q[1];
sx q[1];
rz(-0.81175214) q[1];
sx q[1];
rz(-0.0033671826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30154341) q[0];
sx q[0];
rz(-0.97629428) q[0];
sx q[0];
rz(1.5778753) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0515767) q[2];
sx q[2];
rz(-0.64489105) q[2];
sx q[2];
rz(-2.1457248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0778179) q[1];
sx q[1];
rz(-1.5146202) q[1];
sx q[1];
rz(-0.50775524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.037854511) q[3];
sx q[3];
rz(-2.541171) q[3];
sx q[3];
rz(2.5021533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1668261) q[2];
sx q[2];
rz(-2.4415015) q[2];
sx q[2];
rz(0.33494803) q[2];
rz(2.8454928) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(3.0193442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5122546) q[0];
sx q[0];
rz(-2.6423995) q[0];
sx q[0];
rz(0.41803023) q[0];
rz(0.66608518) q[1];
sx q[1];
rz(-1.2434375) q[1];
sx q[1];
rz(0.24806771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1676511) q[0];
sx q[0];
rz(-0.60567666) q[0];
sx q[0];
rz(2.1781237) q[0];
x q[1];
rz(-2.0204629) q[2];
sx q[2];
rz(-1.5213335) q[2];
sx q[2];
rz(-0.63724697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0679761) q[1];
sx q[1];
rz(-0.66367894) q[1];
sx q[1];
rz(-0.90866088) q[1];
x q[2];
rz(0.94132386) q[3];
sx q[3];
rz(-1.7793312) q[3];
sx q[3];
rz(-1.0104313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69998133) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(-1.9786394) q[2];
rz(-2.5476088) q[3];
sx q[3];
rz(-1.5409527) q[3];
sx q[3];
rz(-2.1514413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.537259) q[0];
sx q[0];
rz(-0.27892497) q[0];
sx q[0];
rz(0.064924031) q[0];
rz(0.36554947) q[1];
sx q[1];
rz(-1.431501) q[1];
sx q[1];
rz(2.4117267) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.176031) q[0];
sx q[0];
rz(-2.4906127) q[0];
sx q[0];
rz(1.7528379) q[0];
rz(-pi) q[1];
rz(-2.9403749) q[2];
sx q[2];
rz(-0.7732174) q[2];
sx q[2];
rz(-1.7380444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.262558) q[1];
sx q[1];
rz(-2.1260736) q[1];
sx q[1];
rz(-3.0295275) q[1];
rz(0.49298005) q[3];
sx q[3];
rz(-2.3928309) q[3];
sx q[3];
rz(-1.4593339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87021747) q[2];
sx q[2];
rz(-1.368618) q[2];
sx q[2];
rz(-2.0591056) q[2];
rz(1.4879976) q[3];
sx q[3];
rz(-0.36181417) q[3];
sx q[3];
rz(0.35058072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38700405) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(-0.3311232) q[0];
rz(1.6638727) q[1];
sx q[1];
rz(-0.96550566) q[1];
sx q[1];
rz(0.32041916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5950111) q[0];
sx q[0];
rz(-0.10215952) q[0];
sx q[0];
rz(-0.25545303) q[0];
rz(-1.4244979) q[2];
sx q[2];
rz(-1.4182036) q[2];
sx q[2];
rz(-0.998978) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72057331) q[1];
sx q[1];
rz(-1.58632) q[1];
sx q[1];
rz(2.7153354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.057366296) q[3];
sx q[3];
rz(-0.88729492) q[3];
sx q[3];
rz(-0.43903109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42935818) q[2];
sx q[2];
rz(-1.4171436) q[2];
sx q[2];
rz(-2.7799535) q[2];
rz(2.2278191) q[3];
sx q[3];
rz(-0.29699609) q[3];
sx q[3];
rz(2.698212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8906032) q[0];
sx q[0];
rz(-0.78762233) q[0];
sx q[0];
rz(3.0255764) q[0];
rz(-1.8138255) q[1];
sx q[1];
rz(-1.9783744) q[1];
sx q[1];
rz(-1.9414925) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2105149) q[0];
sx q[0];
rz(-1.3793945) q[0];
sx q[0];
rz(-1.5899508) q[0];
rz(-pi) q[1];
rz(-1.7723254) q[2];
sx q[2];
rz(-2.4459165) q[2];
sx q[2];
rz(1.3112549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7197693) q[1];
sx q[1];
rz(-2.1323816) q[1];
sx q[1];
rz(1.5535068) q[1];
x q[2];
rz(2.9755244) q[3];
sx q[3];
rz(-1.8072053) q[3];
sx q[3];
rz(2.0061559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8414595) q[2];
sx q[2];
rz(-1.4347142) q[2];
sx q[2];
rz(-0.33162281) q[2];
rz(0.2837818) q[3];
sx q[3];
rz(-1.1016568) q[3];
sx q[3];
rz(-0.42154977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3109741) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(2.9750138) q[0];
rz(2.2919948) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(3.0679682) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66921418) q[0];
sx q[0];
rz(-1.3166691) q[0];
sx q[0];
rz(-0.73331244) q[0];
x q[1];
rz(-1.0208836) q[2];
sx q[2];
rz(-1.9859909) q[2];
sx q[2];
rz(-2.6299494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8371997) q[1];
sx q[1];
rz(-1.0719258) q[1];
sx q[1];
rz(1.1863329) q[1];
x q[2];
rz(-1.0898548) q[3];
sx q[3];
rz(-1.0006128) q[3];
sx q[3];
rz(0.4057623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9061884) q[2];
sx q[2];
rz(-2.3222458) q[2];
sx q[2];
rz(-2.4833615) q[2];
rz(-1.8002347) q[3];
sx q[3];
rz(-2.1737289) q[3];
sx q[3];
rz(-2.2906176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3157208) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(0.52147621) q[1];
sx q[1];
rz(-1.5745402) q[1];
sx q[1];
rz(1.5706617) q[1];
rz(2.1076219) q[2];
sx q[2];
rz(-0.85672094) q[2];
sx q[2];
rz(-2.232205) q[2];
rz(2.9110649) q[3];
sx q[3];
rz(-0.63481019) q[3];
sx q[3];
rz(-1.732375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
