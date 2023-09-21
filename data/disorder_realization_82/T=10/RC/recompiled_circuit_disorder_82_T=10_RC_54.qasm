OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049392603) q[0];
sx q[0];
rz(-0.28486262) q[0];
sx q[0];
rz(-1.0546791) q[0];
x q[1];
rz(0.42135294) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(2.8533964) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0496088) q[1];
sx q[1];
rz(-0.99124747) q[1];
sx q[1];
rz(1.0773354) q[1];
rz(-2.4217371) q[3];
sx q[3];
rz(-2.8350283) q[3];
sx q[3];
rz(-0.43678624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1352284) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(0.41729331) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-0.59536368) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(0.48746902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66747626) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(1.4897896) q[0];
x q[1];
rz(2.0080645) q[2];
sx q[2];
rz(-1.9568866) q[2];
sx q[2];
rz(-1.776945) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65709719) q[1];
sx q[1];
rz(-0.24689455) q[1];
sx q[1];
rz(2.2730278) q[1];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(0.3953735) q[2];
rz(1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(2.9188459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091464) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(1.976165) q[0];
rz(-pi) q[1];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(0.97937102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0157156) q[1];
sx q[1];
rz(-1.5029969) q[1];
sx q[1];
rz(0.61342872) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24554159) q[3];
sx q[3];
rz(-1.89851) q[3];
sx q[3];
rz(-1.6944193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(1.9474691) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5266787) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(2.148669) q[0];
rz(-2.9158981) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(2.8452578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4184119) q[1];
sx q[1];
rz(-0.81058093) q[1];
sx q[1];
rz(2.4852738) q[1];
x q[2];
rz(-2.6357189) q[3];
sx q[3];
rz(-0.40588356) q[3];
sx q[3];
rz(-2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.7368332) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720854) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(2.0399658) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7136953) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(0.75682109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0616152) q[1];
sx q[1];
rz(-1.4887759) q[1];
sx q[1];
rz(-1.3498989) q[1];
x q[2];
rz(2.1702607) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(-0.72367469) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.628081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0244004) q[0];
sx q[0];
rz(-1.8937366) q[0];
sx q[0];
rz(-2.6431503) q[0];
x q[1];
rz(0.061821826) q[2];
sx q[2];
rz(-1.1218058) q[2];
sx q[2];
rz(-2.2668554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7102393) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(2.5297574) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7558277) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(-1.827084) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(0.32062809) q[0];
rz(-pi) q[1];
rz(0.059139472) q[2];
sx q[2];
rz(-1.5416317) q[2];
sx q[2];
rz(1.9524706) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4759051) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(-1.7687294) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39077057) q[3];
sx q[3];
rz(-0.6676538) q[3];
sx q[3];
rz(2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.7857893) q[2];
rz(-1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4521342) q[0];
sx q[0];
rz(-1.9118714) q[0];
sx q[0];
rz(2.336691) q[0];
rz(-pi) q[1];
rz(-1.4872929) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(0.63937843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.028138782) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(-2.0466652) q[1];
rz(-pi) q[2];
rz(-2.1645032) q[3];
sx q[3];
rz(-1.1247375) q[3];
sx q[3];
rz(-2.2462728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(-3.0714152) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(1.1351599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0907744) q[0];
sx q[0];
rz(-0.43681112) q[0];
sx q[0];
rz(0.45849053) q[0];
rz(-pi) q[1];
rz(1.4831545) q[2];
sx q[2];
rz(-2.8200375) q[2];
sx q[2];
rz(-1.9784387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3865349) q[1];
sx q[1];
rz(-2.6549576) q[1];
sx q[1];
rz(-1.6965894) q[1];
rz(1.5376066) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88403945) q[0];
sx q[0];
rz(-2.1243874) q[0];
sx q[0];
rz(-2.2474225) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4354544) q[2];
sx q[2];
rz(-1.6207098) q[2];
sx q[2];
rz(-2.4382255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-2.4073283) q[1];
sx q[1];
rz(-1.5937362) q[1];
rz(-pi) q[2];
rz(-0.65532834) q[3];
sx q[3];
rz(-0.95643759) q[3];
sx q[3];
rz(0.67656803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(1.3148057) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(0.23399227) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(-1.3410939) q[3];
sx q[3];
rz(-1.5256186) q[3];
sx q[3];
rz(-1.516173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];