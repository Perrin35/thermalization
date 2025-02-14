OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.051209) q[0];
sx q[0];
rz(-0.78855711) q[0];
sx q[0];
rz(2.8330044) q[0];
rz(-2.7170972) q[1];
sx q[1];
rz(-1.0841882) q[1];
sx q[1];
rz(0.37556136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080506026) q[0];
sx q[0];
rz(-1.2895786) q[0];
sx q[0];
rz(3.0628231) q[0];
rz(-1.5386861) q[2];
sx q[2];
rz(-1.8726204) q[2];
sx q[2];
rz(-0.15649199) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5243317) q[1];
sx q[1];
rz(-1.4325805) q[1];
sx q[1];
rz(2.6928965) q[1];
x q[2];
rz(1.5442763) q[3];
sx q[3];
rz(-2.9980368) q[3];
sx q[3];
rz(1.7100818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8033119) q[2];
sx q[2];
rz(-2.8901926) q[2];
sx q[2];
rz(0.54331642) q[2];
rz(1.4079037) q[3];
sx q[3];
rz(-1.7270154) q[3];
sx q[3];
rz(0.89731115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59178418) q[0];
sx q[0];
rz(-1.9928638) q[0];
sx q[0];
rz(-0.59709221) q[0];
rz(1.0626556) q[1];
sx q[1];
rz(-2.5025044) q[1];
sx q[1];
rz(0.39892453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9358109) q[0];
sx q[0];
rz(-1.9535741) q[0];
sx q[0];
rz(-2.1924928) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024256134) q[2];
sx q[2];
rz(-1.0139216) q[2];
sx q[2];
rz(0.66489894) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.058417) q[1];
sx q[1];
rz(-0.99715573) q[1];
sx q[1];
rz(-0.87223335) q[1];
x q[2];
rz(2.2116304) q[3];
sx q[3];
rz(-2.0130664) q[3];
sx q[3];
rz(-1.7785566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4661633) q[2];
sx q[2];
rz(-1.9025981) q[2];
sx q[2];
rz(-2.0002401) q[2];
rz(2.2232248) q[3];
sx q[3];
rz(-0.71056241) q[3];
sx q[3];
rz(-0.58951283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68701768) q[0];
sx q[0];
rz(-2.7641251) q[0];
sx q[0];
rz(2.9662509) q[0];
rz(-1.31458) q[1];
sx q[1];
rz(-1.8914696) q[1];
sx q[1];
rz(-0.66741991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771362) q[0];
sx q[0];
rz(-2.6303997) q[0];
sx q[0];
rz(-1.3335467) q[0];
x q[1];
rz(0.66593093) q[2];
sx q[2];
rz(-0.70901044) q[2];
sx q[2];
rz(-0.60707742) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9248501) q[1];
sx q[1];
rz(-1.5713673) q[1];
sx q[1];
rz(-0.3158098) q[1];
x q[2];
rz(-0.61428689) q[3];
sx q[3];
rz(-1.3516055) q[3];
sx q[3];
rz(3.1295071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1158925) q[2];
sx q[2];
rz(-2.1860055) q[2];
sx q[2];
rz(2.058775) q[2];
rz(-1.3057905) q[3];
sx q[3];
rz(-0.81597733) q[3];
sx q[3];
rz(0.79735565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1378655) q[0];
sx q[0];
rz(-0.19345134) q[0];
sx q[0];
rz(-0.61099148) q[0];
rz(-2.5773279) q[1];
sx q[1];
rz(-1.0287501) q[1];
sx q[1];
rz(0.70704031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5093197) q[0];
sx q[0];
rz(-1.5406666) q[0];
sx q[0];
rz(2.3448461) q[0];
rz(0.66989278) q[2];
sx q[2];
rz(-0.50648738) q[2];
sx q[2];
rz(-0.8669002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8312861) q[1];
sx q[1];
rz(-1.093068) q[1];
sx q[1];
rz(-1.076144) q[1];
rz(-pi) q[2];
rz(-2.6534326) q[3];
sx q[3];
rz(-0.55503856) q[3];
sx q[3];
rz(-2.6223573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45346144) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(1.0180265) q[2];
rz(-2.7536143) q[3];
sx q[3];
rz(-1.4464902) q[3];
sx q[3];
rz(0.33897266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46868789) q[0];
sx q[0];
rz(-1.1274811) q[0];
sx q[0];
rz(0.52184033) q[0];
rz(-2.8537967) q[1];
sx q[1];
rz(-0.58240533) q[1];
sx q[1];
rz(-0.60307455) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1408242) q[0];
sx q[0];
rz(-0.10061564) q[0];
sx q[0];
rz(-0.510143) q[0];
rz(0.65400161) q[2];
sx q[2];
rz(-0.79833657) q[2];
sx q[2];
rz(-2.4545074) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9085616) q[1];
sx q[1];
rz(-2.7082293) q[1];
sx q[1];
rz(0.20058452) q[1];
x q[2];
rz(0.92080812) q[3];
sx q[3];
rz(-0.42144708) q[3];
sx q[3];
rz(-1.8325266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7120984) q[2];
sx q[2];
rz(-0.79271972) q[2];
sx q[2];
rz(-0.73954868) q[2];
rz(-2.126501) q[3];
sx q[3];
rz(-2.6138217) q[3];
sx q[3];
rz(-0.019006193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5571112) q[0];
sx q[0];
rz(-2.8040573) q[0];
sx q[0];
rz(-2.2767516) q[0];
rz(0.64019126) q[1];
sx q[1];
rz(-0.65018153) q[1];
sx q[1];
rz(-0.4943628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579635) q[0];
sx q[0];
rz(-1.3682162) q[0];
sx q[0];
rz(2.0264739) q[0];
x q[1];
rz(1.6310235) q[2];
sx q[2];
rz(-2.5382747) q[2];
sx q[2];
rz(2.0996527) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8714218) q[1];
sx q[1];
rz(-3.0746859) q[1];
sx q[1];
rz(-0.54065458) q[1];
x q[2];
rz(-0.81624372) q[3];
sx q[3];
rz(-2.4640969) q[3];
sx q[3];
rz(1.578581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.87259) q[2];
sx q[2];
rz(-0.63429093) q[2];
sx q[2];
rz(2.1799901) q[2];
rz(1.5087992) q[3];
sx q[3];
rz(-1.6274933) q[3];
sx q[3];
rz(-0.076920286) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85840571) q[0];
sx q[0];
rz(-0.34342331) q[0];
sx q[0];
rz(-1.0431694) q[0];
rz(0.6768325) q[1];
sx q[1];
rz(-2.4969641) q[1];
sx q[1];
rz(2.4136995) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5344369) q[0];
sx q[0];
rz(-0.42947665) q[0];
sx q[0];
rz(-1.8398102) q[0];
rz(-0.0038996242) q[2];
sx q[2];
rz(-0.74206381) q[2];
sx q[2];
rz(-2.8067547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54503578) q[1];
sx q[1];
rz(-1.903773) q[1];
sx q[1];
rz(2.2117028) q[1];
x q[2];
rz(2.5520127) q[3];
sx q[3];
rz(-1.2834719) q[3];
sx q[3];
rz(-0.048487566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4595043) q[2];
sx q[2];
rz(-0.75152087) q[2];
sx q[2];
rz(2.6085243) q[2];
rz(1.8367977) q[3];
sx q[3];
rz(-2.6663836) q[3];
sx q[3];
rz(2.3615725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.099139) q[0];
sx q[0];
rz(-0.071439698) q[0];
sx q[0];
rz(-3.1169917) q[0];
rz(-3.091231) q[1];
sx q[1];
rz(-2.5854526) q[1];
sx q[1];
rz(0.76739001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3022223) q[0];
sx q[0];
rz(-1.5451533) q[0];
sx q[0];
rz(-2.6013732) q[0];
rz(-pi) q[1];
rz(-2.8117958) q[2];
sx q[2];
rz(-1.7539795) q[2];
sx q[2];
rz(-0.51815301) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7473306) q[1];
sx q[1];
rz(-0.74968265) q[1];
sx q[1];
rz(2.5452627) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9629836) q[3];
sx q[3];
rz(-1.0309891) q[3];
sx q[3];
rz(-1.7843877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47160992) q[2];
sx q[2];
rz(-0.47405425) q[2];
sx q[2];
rz(1.6101884) q[2];
rz(2.1286185) q[3];
sx q[3];
rz(-1.6771202) q[3];
sx q[3];
rz(-2.5392635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36039627) q[0];
sx q[0];
rz(-0.66806) q[0];
sx q[0];
rz(0.49041954) q[0];
rz(-0.40052739) q[1];
sx q[1];
rz(-2.782395) q[1];
sx q[1];
rz(-1.5706971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6436764) q[0];
sx q[0];
rz(-2.1093858) q[0];
sx q[0];
rz(-1.9890189) q[0];
x q[1];
rz(0.3272662) q[2];
sx q[2];
rz(-1.9169352) q[2];
sx q[2];
rz(2.1838783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8388673) q[1];
sx q[1];
rz(-1.4701533) q[1];
sx q[1];
rz(1.9221646) q[1];
x q[2];
rz(2.3113046) q[3];
sx q[3];
rz(-2.2133996) q[3];
sx q[3];
rz(2.682529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5639497) q[2];
sx q[2];
rz(-2.9190013) q[2];
sx q[2];
rz(2.8655748) q[2];
rz(-0.66463071) q[3];
sx q[3];
rz(-0.97739995) q[3];
sx q[3];
rz(-2.9450534) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4589602) q[0];
sx q[0];
rz(-0.17550547) q[0];
sx q[0];
rz(-1.5631787) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.1140607) q[1];
sx q[1];
rz(-0.064090699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9784713) q[0];
sx q[0];
rz(-2.6829217) q[0];
sx q[0];
rz(1.5171097) q[0];
x q[1];
rz(0.59696609) q[2];
sx q[2];
rz(-2.4325437) q[2];
sx q[2];
rz(0.89151357) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0371263) q[1];
sx q[1];
rz(-2.0580225) q[1];
sx q[1];
rz(-3.0477316) q[1];
x q[2];
rz(-0.96336295) q[3];
sx q[3];
rz(-2.2775536) q[3];
sx q[3];
rz(2.6092495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2854707) q[2];
sx q[2];
rz(-2.8557114) q[2];
sx q[2];
rz(2.2355283) q[2];
rz(-1.9237349) q[3];
sx q[3];
rz(-1.8702312) q[3];
sx q[3];
rz(-2.0876032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8080407) q[0];
sx q[0];
rz(-1.6186436) q[0];
sx q[0];
rz(2.9343395) q[0];
rz(-2.8605657) q[1];
sx q[1];
rz(-1.749975) q[1];
sx q[1];
rz(-0.71798807) q[1];
rz(2.7396474) q[2];
sx q[2];
rz(-1.4951757) q[2];
sx q[2];
rz(-2.7972372) q[2];
rz(-3.1140399) q[3];
sx q[3];
rz(-0.47712986) q[3];
sx q[3];
rz(-2.5148077) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
