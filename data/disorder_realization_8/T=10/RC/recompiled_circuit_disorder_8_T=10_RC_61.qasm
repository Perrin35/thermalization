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
rz(2.2367509) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68569505) q[0];
sx q[0];
rz(-1.0766317) q[0];
sx q[0];
rz(-0.49522884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4279847) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(1.7507391) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20476725) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(-2.8863465) q[1];
x q[2];
rz(-1.2600793) q[3];
sx q[3];
rz(-1.3907392) q[3];
sx q[3];
rz(2.3678399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(0.15371418) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7006943) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.457263) q[0];
sx q[0];
rz(-1.6329137) q[0];
sx q[0];
rz(1.5044042) q[0];
rz(0.91815572) q[2];
sx q[2];
rz(-1.2456206) q[2];
sx q[2];
rz(0.52415372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6389097) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(1.6619976) q[1];
rz(1.7327659) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(-3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(2.8125787) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5488141) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(0.51868784) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0654304) q[0];
sx q[0];
rz(-0.87834529) q[0];
sx q[0];
rz(-2.6718219) q[0];
rz(1.475004) q[2];
sx q[2];
rz(-1.8358069) q[2];
sx q[2];
rz(-1.67213) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0842523) q[1];
sx q[1];
rz(-0.80576128) q[1];
sx q[1];
rz(1.8112103) q[1];
rz(1.5716001) q[3];
sx q[3];
rz(-1.5161848) q[3];
sx q[3];
rz(-2.2664203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(2.6285016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1888694) q[0];
sx q[0];
rz(-0.9523305) q[0];
sx q[0];
rz(1.6678715) q[0];
rz(0.80715837) q[2];
sx q[2];
rz(-1.1045189) q[2];
sx q[2];
rz(1.5392787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40401134) q[1];
sx q[1];
rz(-1.1176795) q[1];
sx q[1];
rz(-1.7246507) q[1];
rz(-pi) q[2];
rz(1.4245093) q[3];
sx q[3];
rz(-0.20855599) q[3];
sx q[3];
rz(1.2775161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(-2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(0.67681926) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(-2.5255323) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6989243) q[0];
sx q[0];
rz(-1.5183503) q[0];
sx q[0];
rz(-3.1057182) q[0];
rz(-pi) q[1];
rz(-2.8370503) q[2];
sx q[2];
rz(-2.7603622) q[2];
sx q[2];
rz(-2.9074557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39561135) q[1];
sx q[1];
rz(-2.1484054) q[1];
sx q[1];
rz(-1.5642158) q[1];
rz(-pi) q[2];
rz(-0.89415278) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(-1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(-1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(2.2568259) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17381829) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(-1.3740747) q[0];
x q[1];
rz(0.21005819) q[2];
sx q[2];
rz(-1.5846328) q[2];
sx q[2];
rz(-2.7163598) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92237597) q[1];
sx q[1];
rz(-1.9799385) q[1];
sx q[1];
rz(2.8273696) q[1];
x q[2];
rz(-2.1045477) q[3];
sx q[3];
rz(-1.3579206) q[3];
sx q[3];
rz(0.45149976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(0.73289245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31156763) q[0];
sx q[0];
rz(-0.68651474) q[0];
sx q[0];
rz(-0.65450432) q[0];
rz(-1.0211208) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(-2.6598425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.43436189) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(-2.1965501) q[1];
rz(-pi) q[2];
rz(1.6781428) q[3];
sx q[3];
rz(-1.6119453) q[3];
sx q[3];
rz(-0.51957182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(-0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(2.432166) q[0];
rz(-1.5052694) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-0.27871305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3942791) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(-0.5745116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3271354) q[2];
sx q[2];
rz(-1.0255587) q[2];
sx q[2];
rz(1.7903763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0875138) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(-0.22884303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3474943) q[3];
sx q[3];
rz(-1.4947065) q[3];
sx q[3];
rz(-0.014811024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-0.056079496) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.1425346) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965217) q[0];
sx q[0];
rz(-0.30884305) q[0];
sx q[0];
rz(-1.2678498) q[0];
rz(-pi) q[1];
rz(-2.4382486) q[2];
sx q[2];
rz(-1.5726657) q[2];
sx q[2];
rz(1.4700996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.090230378) q[1];
sx q[1];
rz(-1.6962595) q[1];
sx q[1];
rz(0.54162234) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90060602) q[3];
sx q[3];
rz(-1.4108676) q[3];
sx q[3];
rz(-2.1538494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1322051) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(2.0995887) q[0];
x q[1];
rz(0.054758666) q[2];
sx q[2];
rz(-2.3878532) q[2];
sx q[2];
rz(-2.4278305) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30675948) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(1.761761) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11390399) q[3];
sx q[3];
rz(-0.98131991) q[3];
sx q[3];
rz(1.1978428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4548268) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.025370601) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(2.2123443) q[2];
sx q[2];
rz(-0.48454787) q[2];
sx q[2];
rz(1.5018644) q[2];
rz(-1.3958037) q[3];
sx q[3];
rz(-0.77944118) q[3];
sx q[3];
rz(-0.13164095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
