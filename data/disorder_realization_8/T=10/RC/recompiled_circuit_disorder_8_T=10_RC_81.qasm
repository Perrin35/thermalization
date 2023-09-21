OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(4.0868563) q[0];
sx q[0];
rz(9.950369) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63431595) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(1.0213724) q[0];
rz(-pi) q[1];
rz(-0.22613871) q[2];
sx q[2];
rz(-1.3790501) q[2];
sx q[2];
rz(0.5118256) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.875784) q[1];
sx q[1];
rz(-1.8059397) q[1];
sx q[1];
rz(-1.1633412) q[1];
rz(-pi) q[2];
rz(2.9526887) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(2.4019935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(-0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7006943) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(-2.3764215) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(2.4786425) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239319) q[0];
sx q[0];
rz(-1.6370602) q[0];
sx q[0];
rz(-0.06225417) q[0];
rz(-pi) q[1];
rz(0.40132482) q[2];
sx q[2];
rz(-2.1839645) q[2];
sx q[2];
rz(0.80712986) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.502683) q[1];
sx q[1];
rz(-1.1306292) q[1];
sx q[1];
rz(-1.6619976) q[1];
rz(-pi) q[2];
rz(2.8421721) q[3];
sx q[3];
rz(-0.50588183) q[3];
sx q[3];
rz(-0.29380709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(-2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(0.51868784) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9496574) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(2.3200672) q[0];
rz(-pi) q[1];
rz(2.8027595) q[2];
sx q[2];
rz(-2.860184) q[2];
sx q[2];
rz(2.0237405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0842523) q[1];
sx q[1];
rz(-0.80576128) q[1];
sx q[1];
rz(1.3303824) q[1];
rz(-pi) q[2];
rz(-3.1268901) q[3];
sx q[3];
rz(-0.054617453) q[3];
sx q[3];
rz(2.2811449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(2.5391501) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(0.13901916) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.65072) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78631567) q[0];
sx q[0];
rz(-0.62505165) q[0];
sx q[0];
rz(3.0062208) q[0];
rz(-pi) q[1];
rz(-0.94167534) q[2];
sx q[2];
rz(-2.2721014) q[2];
sx q[2];
rz(2.6710682) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7375813) q[1];
sx q[1];
rz(-1.1176795) q[1];
sx q[1];
rz(1.416942) q[1];
rz(-pi) q[2];
rz(-0.030839132) q[3];
sx q[3];
rz(-1.7770924) q[3];
sx q[3];
rz(-1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47485581) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693562) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-2.4647734) q[0];
rz(-0.49304402) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(2.5255323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0115833) q[0];
sx q[0];
rz(-1.6066215) q[0];
sx q[0];
rz(-1.623276) q[0];
rz(-pi) q[1];
rz(1.6904171) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(-0.56064831) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39561135) q[1];
sx q[1];
rz(-0.99318722) q[1];
sx q[1];
rz(-1.5773768) q[1];
rz(-pi) q[2];
rz(-2.8354007) q[3];
sx q[3];
rz(-0.91727835) q[3];
sx q[3];
rz(-2.9432076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(-2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42246321) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(-2.2568259) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17381829) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.7675179) q[0];
rz(-pi) q[1];
rz(0.066263513) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(-1.2103684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6217475) q[1];
sx q[1];
rz(-1.2832844) q[1];
sx q[1];
rz(-1.9985755) q[1];
rz(-pi) q[2];
rz(2.8956036) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(-0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(0.3113783) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(3.080522) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-0.73289245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.830025) q[0];
sx q[0];
rz(-2.4550779) q[0];
sx q[0];
rz(-0.65450432) q[0];
rz(-pi) q[1];
x q[1];
rz(2.272846) q[2];
sx q[2];
rz(-1.947543) q[2];
sx q[2];
rz(-1.6422611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7072308) q[1];
sx q[1];
rz(-1.5642628) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(3.1002058) q[3];
sx q[3];
rz(-1.4635411) q[3];
sx q[3];
rz(-1.0556575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0682893) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(2.627009) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(-1.6363232) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(-0.27871305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.404971) q[0];
sx q[0];
rz(-1.4355037) q[0];
sx q[0];
rz(2.9297329) q[0];
rz(-1.8144572) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(-1.7903763) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7292494) q[1];
sx q[1];
rz(-1.3530429) q[1];
sx q[1];
rz(-1.8885683) q[1];
x q[2];
rz(1.2392427) q[3];
sx q[3];
rz(-0.23570508) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(0.056079496) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99496019) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(-2.6682207) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(-1.999058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24507095) q[0];
sx q[0];
rz(-0.30884305) q[0];
sx q[0];
rz(1.8737428) q[0];
rz(-pi) q[1];
rz(-0.0028902729) q[2];
sx q[2];
rz(-2.4382466) q[2];
sx q[2];
rz(0.10290111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6856319) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(0.23994069) q[1];
rz(-pi) q[2];
rz(-0.20299083) q[3];
sx q[3];
rz(-2.2309125) q[3];
sx q[3];
rz(-2.4329894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.9030301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9311213) q[0];
sx q[0];
rz(-0.56715542) q[0];
sx q[0];
rz(-1.9803067) q[0];
rz(-pi) q[1];
rz(1.5194703) q[2];
sx q[2];
rz(-0.81846279) q[2];
sx q[2];
rz(-0.78879702) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4148256) q[1];
sx q[1];
rz(-1.6884202) q[1];
sx q[1];
rz(-0.91314258) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0276887) q[3];
sx q[3];
rz(-2.1602727) q[3];
sx q[3];
rz(-1.1978428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(3.0977541) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(0.025370601) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-2.8364137) q[2];
sx q[2];
rz(-1.9532433) q[2];
sx q[2];
rz(0.80079186) q[2];
rz(2.3425441) q[3];
sx q[3];
rz(-1.4481164) q[3];
sx q[3];
rz(-1.8275402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];