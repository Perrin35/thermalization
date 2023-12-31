OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(-2.5807227) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5117447) q[0];
sx q[0];
rz(-1.4667604) q[0];
sx q[0];
rz(1.3826136) q[0];
rz(0.36429976) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(0.42687624) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79599586) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(-0.3791581) q[1];
x q[2];
rz(1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(-0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(1.6024626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3829271) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(-1.4325607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23761959) q[2];
sx q[2];
rz(-0.75287205) q[2];
sx q[2];
rz(1.876229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3746678) q[1];
sx q[1];
rz(-0.51480773) q[1];
sx q[1];
rz(-1.2109846) q[1];
rz(-pi) q[2];
rz(-2.6529796) q[3];
sx q[3];
rz(-2.6693137) q[3];
sx q[3];
rz(-2.3924535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-0.55535299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12268513) q[0];
sx q[0];
rz(-2.2203831) q[0];
sx q[0];
rz(-2.5193549) q[0];
rz(1.0981512) q[2];
sx q[2];
rz(-2.522509) q[2];
sx q[2];
rz(1.667779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.024228) q[1];
sx q[1];
rz(-1.7324565) q[1];
sx q[1];
rz(-2.2092186) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8954574) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(-0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811532) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(2.8864158) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2164453) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(1.2544592) q[0];
rz(2.7093676) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(1.150711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5863122) q[1];
sx q[1];
rz(-1.781207) q[1];
sx q[1];
rz(1.2661238) q[1];
rz(2.0917986) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(1.9998159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5909114) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.3758855) q[0];
rz(1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2984021) q[0];
sx q[0];
rz(-2.139233) q[0];
sx q[0];
rz(-1.1355023) q[0];
rz(2.202583) q[2];
sx q[2];
rz(-2.5917705) q[2];
sx q[2];
rz(0.75013559) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6009112) q[1];
sx q[1];
rz(-1.2680506) q[1];
sx q[1];
rz(1.5555698) q[1];
rz(1.4291184) q[3];
sx q[3];
rz(-2.5767527) q[3];
sx q[3];
rz(-2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(0.9544968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4650824) q[0];
sx q[0];
rz(-2.0970779) q[0];
sx q[0];
rz(-0.17980534) q[0];
rz(-pi) q[1];
rz(0.94021057) q[2];
sx q[2];
rz(-1.5988837) q[2];
sx q[2];
rz(-2.0036151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2973605) q[1];
sx q[1];
rz(-0.73207049) q[1];
sx q[1];
rz(0.75433235) q[1];
rz(-pi) q[2];
rz(2.0733842) q[3];
sx q[3];
rz(-0.32940255) q[3];
sx q[3];
rz(-0.1474895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(1.0423638) q[2];
rz(-2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2110721) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.2852438) q[0];
x q[1];
rz(1.4201944) q[2];
sx q[2];
rz(-0.65956958) q[2];
sx q[2];
rz(-1.9285551) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89027379) q[1];
sx q[1];
rz(-2.2379025) q[1];
sx q[1];
rz(0.031884738) q[1];
rz(2.9355572) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8053749) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.6280744) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-0.73658529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1153699) q[0];
sx q[0];
rz(-1.5407469) q[0];
sx q[0];
rz(3.1307334) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2543711) q[2];
sx q[2];
rz(-0.3379312) q[2];
sx q[2];
rz(-0.28883176) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(-2.9176941) q[1];
rz(-1.1235808) q[3];
sx q[3];
rz(-1.2774602) q[3];
sx q[3];
rz(-0.48188996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.1972315) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-0.7448147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5972827) q[0];
sx q[0];
rz(-2.3047857) q[0];
sx q[0];
rz(-0.44041667) q[0];
rz(-pi) q[1];
rz(-0.3955598) q[2];
sx q[2];
rz(-0.8330847) q[2];
sx q[2];
rz(-1.2566483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77764952) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(1.9294444) q[1];
rz(-pi) q[2];
rz(2.7471077) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(-0.2713954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1054489) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.1709447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02689657) q[0];
sx q[0];
rz(-2.7141502) q[0];
sx q[0];
rz(-0.052274152) q[0];
rz(-pi) q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(0.68683456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3552637) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(-1.5481871) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2788494) q[3];
sx q[3];
rz(-1.1598831) q[3];
sx q[3];
rz(-1.3565855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-2.4026985) q[2];
sx q[2];
rz(-1.0514435) q[2];
sx q[2];
rz(-2.4098868) q[2];
rz(1.5242143) q[3];
sx q[3];
rz(-2.8077879) q[3];
sx q[3];
rz(2.2624364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
