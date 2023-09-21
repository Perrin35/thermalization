OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96072223) q[0];
sx q[0];
rz(-1.7579494) q[0];
sx q[0];
rz(0.10589177) q[0];
x q[1];
rz(2.7772929) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(2.7147164) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6601662) q[1];
sx q[1];
rz(-0.62454849) q[1];
sx q[1];
rz(0.98510965) q[1];
rz(1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(-0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30277006) q[0];
sx q[0];
rz(-1.6478224) q[0];
sx q[0];
rz(-2.1588615) q[0];
x q[1];
rz(-0.73866663) q[2];
sx q[2];
rz(-1.7324442) q[2];
sx q[2];
rz(2.6612298) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6290366) q[1];
sx q[1];
rz(-1.7450383) q[1];
sx q[1];
rz(-2.0577355) q[1];
x q[2];
rz(2.7178571) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(0.37950294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-0.55535299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(-0.81911266) q[0];
rz(-0.31366445) q[2];
sx q[2];
rz(-1.0278388) q[2];
sx q[2];
rz(1.107159) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33298102) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.3036742) q[1];
rz(-2.1159806) q[3];
sx q[3];
rz(-1.7435929) q[3];
sx q[3];
rz(-1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(2.8864158) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2164453) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(-1.2544592) q[0];
x q[1];
rz(1.9011263) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(1.6883048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7123588) q[1];
sx q[1];
rz(-2.7731967) q[1];
sx q[1];
rz(-0.952094) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4820443) q[3];
sx q[3];
rz(-2.6188861) q[3];
sx q[3];
rz(0.50597092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-3.0392652) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027095196) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(-0.61371213) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1115083) q[2];
sx q[2];
rz(-1.8845203) q[2];
sx q[2];
rz(-0.26278652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1069378) q[1];
sx q[1];
rz(-1.5562623) q[1];
sx q[1];
rz(-2.8388139) q[1];
rz(1.0104996) q[3];
sx q[3];
rz(-1.6464525) q[3];
sx q[3];
rz(-1.2922985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-2.1870959) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6765103) q[0];
sx q[0];
rz(-2.0970779) q[0];
sx q[0];
rz(-0.17980534) q[0];
x q[1];
rz(-3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(-0.45331732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94297385) q[1];
sx q[1];
rz(-2.0795515) q[1];
sx q[1];
rz(-2.1224623) q[1];
x q[2];
rz(-0.16320634) q[3];
sx q[3];
rz(-1.8582134) q[3];
sx q[3];
rz(0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(-1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(-2.5315703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689432) q[0];
sx q[0];
rz(-2.5805051) q[0];
sx q[0];
rz(-0.48604301) q[0];
rz(-pi) q[1];
rz(-0.1158175) q[2];
sx q[2];
rz(-2.2216185) q[2];
sx q[2];
rz(-1.7388369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66079084) q[1];
sx q[1];
rz(-1.5958438) q[1];
sx q[1];
rz(0.90344306) q[1];
x q[2];
rz(2.9355572) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(2.2231893) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.6280744) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1153699) q[0];
sx q[0];
rz(-1.6008458) q[0];
sx q[0];
rz(-0.010859246) q[0];
rz(0.10891624) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(-3.0963754) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31863403) q[1];
sx q[1];
rz(-1.4689323) q[1];
sx q[1];
rz(-0.46551367) q[1];
x q[2];
rz(-0.32324507) q[3];
sx q[3];
rz(-1.1439699) q[3];
sx q[3];
rz(-2.1904898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1071876) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(2.396778) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9841524) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-2.0122583) q[0];
x q[1];
rz(-0.79297519) q[2];
sx q[2];
rz(-1.8599531) q[2];
sx q[2];
rz(-0.58794978) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77764952) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(-1.2121483) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39448491) q[3];
sx q[3];
rz(-0.82291616) q[3];
sx q[3];
rz(-2.8701973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.1709447) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5491345) q[0];
sx q[0];
rz(0.42692703) q[0];
x q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(-2.4547581) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78466258) q[1];
sx q[1];
rz(-1.5481879) q[1];
sx q[1];
rz(3.1329586) q[1];
rz(-2.5578299) q[3];
sx q[3];
rz(-0.49920344) q[3];
sx q[3];
rz(-1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.6177572) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82407172) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-2.2291017) q[2];
sx q[2];
rz(-0.9463263) q[2];
sx q[2];
rz(1.8778388) q[2];
rz(0.016146544) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];