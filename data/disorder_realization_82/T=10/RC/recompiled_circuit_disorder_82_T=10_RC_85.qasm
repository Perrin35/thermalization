OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(4.1788221) q[0];
sx q[0];
rz(9.0691789) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0227538) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.820178) q[0];
rz(-pi) q[1];
rz(1.3470634) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.3734693) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76525926) q[1];
sx q[1];
rz(-1.9783124) q[1];
sx q[1];
rz(-0.63912649) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9079307) q[3];
sx q[3];
rz(-1.3705001) q[3];
sx q[3];
rz(-1.3113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-2.6541236) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9092642) q[0];
sx q[0];
rz(-3.0380913) q[0];
sx q[0];
rz(-2.2444025) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7198118) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(0.031907206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.060354787) q[1];
sx q[1];
rz(-1.3831257) q[1];
sx q[1];
rz(-2.9802122) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2536212) q[3];
sx q[3];
rz(-1.653169) q[3];
sx q[3];
rz(-1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(0.3953735) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19668002) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-2.9188459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5644154) q[0];
sx q[0];
rz(-0.555493) q[0];
sx q[0];
rz(2.3781858) q[0];
x q[1];
rz(1.0388971) q[2];
sx q[2];
rz(-2.7173018) q[2];
sx q[2];
rz(-1.588856) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5388515) q[1];
sx q[1];
rz(-2.1826084) q[1];
sx q[1];
rz(1.6536504) q[1];
rz(-pi) q[2];
rz(1.907903) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(-0.043134886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738149) q[0];
sx q[0];
rz(-0.91705634) q[0];
sx q[0];
rz(-2.6184404) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6572861) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(0.53777018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4184119) q[1];
sx q[1];
rz(-0.81058093) q[1];
sx q[1];
rz(-2.4852738) q[1];
x q[2];
rz(-0.50587378) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(0.94569262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9437287) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.7368332) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(0.4531025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3081007) q[0];
sx q[0];
rz(-1.5335324) q[0];
sx q[0];
rz(1.4971855) q[0];
x q[1];
rz(-0.034899072) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(0.81894433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0616152) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(1.3498989) q[1];
rz(-pi) q[2];
rz(-0.94282486) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(-3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0742764) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(2.5308454) q[0];
x q[1];
rz(-1.6983301) q[2];
sx q[2];
rz(-2.6886534) q[2];
sx q[2];
rz(2.4085101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7102393) q[1];
sx q[1];
rz(-2.2696886) q[1];
sx q[1];
rz(2.5297574) q[1];
rz(-pi) q[2];
rz(2.7558277) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.05904077) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(0.39757279) q[0];
rz(-1.827084) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.9546753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8531187) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(0.32062809) q[0];
rz(-1.5415807) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(-0.37994775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3992607) q[1];
sx q[1];
rz(-0.84335828) q[1];
sx q[1];
rz(0.17952339) q[1];
rz(-pi) q[2];
rz(-0.62992923) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5671974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521342) q[0];
sx q[0];
rz(-1.9118714) q[0];
sx q[0];
rz(0.80490168) q[0];
rz(0.048642283) q[2];
sx q[2];
rz(-1.0445134) q[2];
sx q[2];
rz(-2.4056048) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(-1.2703018) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1645032) q[3];
sx q[3];
rz(-2.0168552) q[3];
sx q[3];
rz(0.89531985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29256233) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410256) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(2.7450949) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4831545) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(-1.1631539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4371722) q[1];
sx q[1];
rz(-1.51209) q[1];
sx q[1];
rz(-1.0874332) q[1];
rz(-1.5376066) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(-0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-2.9454254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88403945) q[0];
sx q[0];
rz(-2.1243874) q[0];
sx q[0];
rz(2.2474225) q[0];
rz(1.7061383) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(-0.70336715) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0386104) q[1];
sx q[1];
rz(-0.83676941) q[1];
sx q[1];
rz(3.1208913) q[1];
rz(-2.2833061) q[3];
sx q[3];
rz(-0.86601102) q[3];
sx q[3];
rz(-1.603905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(-1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-2.0666163) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186196) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(-2.9076004) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(-0.046394596) q[3];
sx q[3];
rz(-1.3413324) q[3];
sx q[3];
rz(-3.0975292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
