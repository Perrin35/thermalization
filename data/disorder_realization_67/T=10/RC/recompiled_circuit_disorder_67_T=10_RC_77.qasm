OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(2.0818721) q[0];
sx q[0];
rz(11.835397) q[0];
rz(1.641474) q[1];
sx q[1];
rz(-1.0348231) q[1];
sx q[1];
rz(2.1980481) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9991) q[0];
sx q[0];
rz(-1.6377791) q[0];
sx q[0];
rz(-1.6058558) q[0];
x q[1];
rz(-1.4859096) q[2];
sx q[2];
rz(-0.92157084) q[2];
sx q[2];
rz(-2.6390586) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.01602068) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(-1.9858951) q[1];
x q[2];
rz(-0.30480095) q[3];
sx q[3];
rz(-1.3918575) q[3];
sx q[3];
rz(1.5233056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8865108) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(-1.250766) q[2];
rz(1.7154153) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(2.1616518) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0086867) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(1.8006181) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(0.96639955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37106284) q[0];
sx q[0];
rz(-0.75936985) q[0];
sx q[0];
rz(0.66803996) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069300058) q[2];
sx q[2];
rz(-2.5182704) q[2];
sx q[2];
rz(2.9181883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4233154) q[1];
sx q[1];
rz(-1.4916972) q[1];
sx q[1];
rz(0.11421108) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45774777) q[3];
sx q[3];
rz(-0.80349892) q[3];
sx q[3];
rz(-0.97704923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0559343) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(1.9484693) q[3];
sx q[3];
rz(-1.5914702) q[3];
sx q[3];
rz(1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0369204) q[0];
sx q[0];
rz(-1.5043229) q[0];
sx q[0];
rz(-1.1874636) q[0];
rz(1.9056412) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(-1.8240066) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353969) q[0];
sx q[0];
rz(-0.72512308) q[0];
sx q[0];
rz(2.8185185) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74430978) q[2];
sx q[2];
rz(-0.28713206) q[2];
sx q[2];
rz(-2.9922275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31955645) q[1];
sx q[1];
rz(-1.469538) q[1];
sx q[1];
rz(-0.42145573) q[1];
rz(-pi) q[2];
rz(1.7749952) q[3];
sx q[3];
rz(-0.52650982) q[3];
sx q[3];
rz(-0.029475676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(-2.9349566) q[2];
rz(0.7080428) q[3];
sx q[3];
rz(-0.20773023) q[3];
sx q[3];
rz(0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(2.2553717) q[0];
rz(1.0097424) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.2264235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7851631) q[0];
sx q[0];
rz(-2.4117081) q[0];
sx q[0];
rz(-1.6701783) q[0];
x q[1];
rz(2.2573651) q[2];
sx q[2];
rz(-1.9366169) q[2];
sx q[2];
rz(-0.58931749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36818477) q[1];
sx q[1];
rz(-0.28973026) q[1];
sx q[1];
rz(-0.75399953) q[1];
rz(-pi) q[2];
rz(-2.6552116) q[3];
sx q[3];
rz(-1.9035305) q[3];
sx q[3];
rz(-0.90852028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(2.148596) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-2.1195181) q[3];
sx q[3];
rz(-2.657857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(1.1859878) q[0];
rz(-1.4233308) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(-2.5591154) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5929778) q[0];
sx q[0];
rz(-2.5354404) q[0];
sx q[0];
rz(1.6968326) q[0];
rz(-pi) q[1];
rz(-0.8930348) q[2];
sx q[2];
rz(-1.6871916) q[2];
sx q[2];
rz(1.900577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1867265) q[1];
sx q[1];
rz(-0.71471067) q[1];
sx q[1];
rz(-3.1123118) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6378239) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(3.0598818) q[2];
rz(-2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24494568) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(-3.1337877) q[0];
rz(1.7410949) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-2.0369464) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68708006) q[0];
sx q[0];
rz(-2.2283163) q[0];
sx q[0];
rz(-1.7772654) q[0];
x q[1];
rz(0.32128895) q[2];
sx q[2];
rz(-2.4761204) q[2];
sx q[2];
rz(1.2207536) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1659516) q[1];
sx q[1];
rz(-0.8756606) q[1];
sx q[1];
rz(-2.402311) q[1];
x q[2];
rz(1.8940582) q[3];
sx q[3];
rz(-2.6544016) q[3];
sx q[3];
rz(0.67684735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8905939) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884488) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(3.0798262) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2568946) q[0];
sx q[0];
rz(-1.1567133) q[0];
sx q[0];
rz(0.82722442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6484117) q[2];
sx q[2];
rz(-1.5845808) q[2];
sx q[2];
rz(-2.8790561) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55918499) q[1];
sx q[1];
rz(-1.294194) q[1];
sx q[1];
rz(-0.72699593) q[1];
rz(-pi) q[2];
rz(-1.816733) q[3];
sx q[3];
rz(-1.8898367) q[3];
sx q[3];
rz(-0.46003534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3322488) q[2];
sx q[2];
rz(-2.377254) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(1.404473) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5161045) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(1.6199934) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(0.63751784) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.92256) q[0];
sx q[0];
rz(-1.0244644) q[0];
sx q[0];
rz(0.68740293) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72458467) q[2];
sx q[2];
rz(-1.2983592) q[2];
sx q[2];
rz(-2.9582634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70367614) q[1];
sx q[1];
rz(-0.99214593) q[1];
sx q[1];
rz(2.0520567) q[1];
rz(-2.5708837) q[3];
sx q[3];
rz(-1.686704) q[3];
sx q[3];
rz(-1.4012208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1311538) q[2];
sx q[2];
rz(-2.4413979) q[2];
sx q[2];
rz(-2.0054224) q[2];
rz(1.6561967) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5159601) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(0.67614722) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(-2.1264145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43468201) q[0];
sx q[0];
rz(-2.7043531) q[0];
sx q[0];
rz(-0.31243639) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4887772) q[2];
sx q[2];
rz(-1.6932994) q[2];
sx q[2];
rz(-2.0147689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0820513) q[1];
sx q[1];
rz(-0.50652981) q[1];
sx q[1];
rz(-0.7854714) q[1];
rz(2.7916662) q[3];
sx q[3];
rz(-2.3438128) q[3];
sx q[3];
rz(-1.6328904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-2.1742415) q[2];
rz(-1.597065) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5627562) q[0];
sx q[0];
rz(-3.0850947) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(1.4046232) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79849762) q[0];
sx q[0];
rz(-1.4416749) q[0];
sx q[0];
rz(-1.207418) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37168665) q[2];
sx q[2];
rz(-0.32495299) q[2];
sx q[2];
rz(1.5124958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23552588) q[1];
sx q[1];
rz(-2.0350254) q[1];
sx q[1];
rz(-2.3266351) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9910562) q[3];
sx q[3];
rz(-1.5621462) q[3];
sx q[3];
rz(-2.0899525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(-2.8354697) q[2];
rz(-2.7434769) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(-1.9819992) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(-2.878306) q[2];
sx q[2];
rz(-1.2821715) q[2];
sx q[2];
rz(0.53496219) q[2];
rz(-1.0032734) q[3];
sx q[3];
rz(-0.56837396) q[3];
sx q[3];
rz(3.0336998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];