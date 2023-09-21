OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5358937) q[0];
sx q[0];
rz(-0.68471691) q[0];
sx q[0];
rz(0.84795714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7674255) q[2];
sx q[2];
rz(-1.7927205) q[2];
sx q[2];
rz(-2.1264399) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.875784) q[1];
sx q[1];
rz(-1.8059397) q[1];
sx q[1];
rz(-1.1633412) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0337898) q[3];
sx q[3];
rz(-0.35764965) q[3];
sx q[3];
rz(-1.8358177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(3.0453483) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(-0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-2.4786425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6843296) q[0];
sx q[0];
rz(-1.6329137) q[0];
sx q[0];
rz(1.5044042) q[0];
rz(-pi) q[1];
rz(2.0775954) q[2];
sx q[2];
rz(-0.71841824) q[2];
sx q[2];
rz(1.4425299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.502683) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(1.4795951) q[1];
rz(-pi) q[2];
rz(-2.6547673) q[3];
sx q[3];
rz(-1.4273705) q[3];
sx q[3];
rz(1.6008582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(2.8125787) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5488141) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(-0.51868784) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0761622) q[0];
sx q[0];
rz(-0.87834529) q[0];
sx q[0];
rz(-2.6718219) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8754183) q[2];
sx q[2];
rz(-1.4783579) q[2];
sx q[2];
rz(-3.0150974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0573404) q[1];
sx q[1];
rz(-0.80576128) q[1];
sx q[1];
rz(1.3303824) q[1];
rz(-pi) q[2];
rz(-0.014702602) q[3];
sx q[3];
rz(-0.054617453) q[3];
sx q[3];
rz(0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(3.0730491) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(0.13901916) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(2.6285016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1888694) q[0];
sx q[0];
rz(-2.1892622) q[0];
sx q[0];
rz(1.4737211) q[0];
rz(-pi) q[1];
rz(2.3344343) q[2];
sx q[2];
rz(-1.1045189) q[2];
sx q[2];
rz(1.6023139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.063555524) q[1];
sx q[1];
rz(-2.6647898) q[1];
sx q[1];
rz(-2.8366691) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1107535) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(2.9117057) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(0.61606032) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426684) q[0];
sx q[0];
rz(-1.5183503) q[0];
sx q[0];
rz(3.1057182) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8370503) q[2];
sx q[2];
rz(-2.7603622) q[2];
sx q[2];
rz(-2.9074557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38356009) q[1];
sx q[1];
rz(-0.57764232) q[1];
sx q[1];
rz(3.1314965) q[1];
x q[2];
rz(-2.8354007) q[3];
sx q[3];
rz(-0.91727835) q[3];
sx q[3];
rz(0.19838504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(-0.25536728) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7191294) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(-0.47250026) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(2.2568259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17381829) q[0];
sx q[0];
rz(-2.7533555) q[0];
sx q[0];
rz(1.3740747) q[0];
x q[1];
rz(-0.21005819) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(0.42523281) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2192167) q[1];
sx q[1];
rz(-1.9799385) q[1];
sx q[1];
rz(-0.31422305) q[1];
rz(1.1690508) q[3];
sx q[3];
rz(-2.5707977) q[3];
sx q[3];
rz(1.4626383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(0.3113783) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-0.73289245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31156763) q[0];
sx q[0];
rz(-2.4550779) q[0];
sx q[0];
rz(0.65450432) q[0];
x q[1];
rz(-2.272846) q[2];
sx q[2];
rz(-1.947543) q[2];
sx q[2];
rz(1.6422611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9961173) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(-1.5596418) q[1];
x q[2];
rz(1.4634499) q[3];
sx q[3];
rz(-1.6119453) q[3];
sx q[3];
rz(-2.6220208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(-2.627009) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.085389) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(-1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(0.27871305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80517171) q[0];
sx q[0];
rz(-1.780691) q[0];
sx q[0];
rz(-1.4324485) q[0];
x q[1];
rz(1.3271354) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(-1.7903763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0875138) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(-0.22884303) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2392427) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30148208) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(0.056079496) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(-0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.1425346) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56209598) q[0];
sx q[0];
rz(-1.8651433) q[0];
sx q[0];
rz(0.094898183) q[0];
x q[1];
rz(0.70334401) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(-1.4700996) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6856319) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(-0.23994069) q[1];
rz(-1.3167131) q[3];
sx q[3];
rz(-0.68613201) q[3];
sx q[3];
rz(-0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.2385626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4558444) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(-0.2483764) q[0];
x q[1];
rz(2.3886016) q[2];
sx q[2];
rz(-1.6082616) q[2];
sx q[2];
rz(-2.3245036) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8348332) q[1];
sx q[1];
rz(-2.4750437) q[1];
sx q[1];
rz(1.761761) q[1];
rz(-1.7391316) q[3];
sx q[3];
rz(-0.59909648) q[3];
sx q[3];
rz(1.7408016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4548268) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.745789) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
