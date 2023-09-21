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
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(-1.9089729) q[1];
sx q[1];
rz(0.90484172) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68569505) q[0];
sx q[0];
rz(-1.0766317) q[0];
sx q[0];
rz(-2.6463638) q[0];
x q[1];
rz(-2.4279847) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(-1.7507391) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20476725) q[1];
sx q[1];
rz(-1.966404) q[1];
sx q[1];
rz(-0.25524615) q[1];
x q[2];
rz(0.18890394) q[3];
sx q[3];
rz(-1.2652664) q[3];
sx q[3];
rz(-0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(0.096244372) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(-0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(-2.4786425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5040341) q[0];
sx q[0];
rz(-3.0507037) q[0];
sx q[0];
rz(-2.3239517) q[0];
x q[1];
rz(2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4274802) q[1];
sx q[1];
rz(-0.44891) q[1];
sx q[1];
rz(2.950579) q[1];
x q[2];
rz(2.6547673) q[3];
sx q[3];
rz(-1.7142222) q[3];
sx q[3];
rz(-1.5407345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(-2.6229048) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1919353) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(2.3200672) q[0];
x q[1];
rz(-2.8027595) q[2];
sx q[2];
rz(-2.860184) q[2];
sx q[2];
rz(-2.0237405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8232302) q[1];
sx q[1];
rz(-1.743411) q[1];
sx q[1];
rz(2.361972) q[1];
rz(-pi) q[2];
rz(3.0869811) q[3];
sx q[3];
rz(-1.5699937) q[3];
sx q[3];
rz(2.4459248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(-0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(2.6285016) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78631567) q[0];
sx q[0];
rz(-0.62505165) q[0];
sx q[0];
rz(-0.13537188) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60855234) q[2];
sx q[2];
rz(-2.2366479) q[2];
sx q[2];
rz(0.37492875) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7375813) q[1];
sx q[1];
rz(-2.0239132) q[1];
sx q[1];
rz(1.416942) q[1];
rz(-3.1107535) q[3];
sx q[3];
rz(-1.7770924) q[3];
sx q[3];
rz(-1.4269958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-2.9117057) q[2];
rz(2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-2.5255323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426684) q[0];
sx q[0];
rz(-1.6232423) q[0];
sx q[0];
rz(0.035874493) q[0];
rz(-2.8370503) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(-0.234137) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38356009) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(-3.1314965) q[1];
rz(-2.8354007) q[3];
sx q[3];
rz(-2.2243143) q[3];
sx q[3];
rz(-0.19838504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(-2.8862254) q[2];
rz(-1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(-2.6690924) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(0.88476673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17381829) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.7675179) q[0];
rz(-2.9315345) q[2];
sx q[2];
rz(-1.5846328) q[2];
sx q[2];
rz(-2.7163598) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2192167) q[1];
sx q[1];
rz(-1.1616542) q[1];
sx q[1];
rz(0.31422305) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24598908) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(-0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(-2.8302144) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(-2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(-2.4087002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3467348) q[0];
sx q[0];
rz(-1.96694) q[0];
sx q[0];
rz(-0.57647716) q[0];
rz(-pi) q[1];
rz(-2.6635025) q[2];
sx q[2];
rz(-0.92667898) q[2];
sx q[2];
rz(-2.9113876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1454754) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(-1.5819509) q[1];
x q[2];
rz(-1.2039127) q[3];
sx q[3];
rz(-3.0266579) q[3];
sx q[3];
rz(-2.4550408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0682893) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(-2.627009) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(-1.6363232) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(2.8628796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.404971) q[0];
sx q[0];
rz(-1.706089) q[0];
sx q[0];
rz(-0.21185974) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3271354) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(1.7903763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0540788) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(-0.22884303) q[1];
rz(-1.2392427) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.99496019) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(-2.1726998) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(-1.1425346) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.105285) q[0];
sx q[0];
rz(-1.6616016) q[0];
sx q[0];
rz(-1.8663976) q[0];
rz(-pi) q[1];
rz(0.0028902729) q[2];
sx q[2];
rz(-0.70334607) q[2];
sx q[2];
rz(0.10290111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6856319) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(2.901652) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2409866) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(0.98774324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(-0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(2.2258863) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(-1.9030301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21047132) q[0];
sx q[0];
rz(-0.56715542) q[0];
sx q[0];
rz(-1.161286) q[0];
rz(-0.054758666) q[2];
sx q[2];
rz(-0.75373947) q[2];
sx q[2];
rz(0.71376212) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8348332) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(-1.761761) q[1];
x q[2];
rz(-1.7391316) q[3];
sx q[3];
rz(-2.5424962) q[3];
sx q[3];
rz(-1.7408016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(0.043838538) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(-2.2123443) q[2];
sx q[2];
rz(-2.6570448) q[2];
sx q[2];
rz(-1.6397283) q[2];
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