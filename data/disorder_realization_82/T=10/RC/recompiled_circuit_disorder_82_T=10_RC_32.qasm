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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1188388) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.3214146) q[0];
x q[1];
rz(-2.6718164) q[2];
sx q[2];
rz(-0.46576408) q[2];
sx q[2];
rz(2.2848406) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(2.5024662) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9079307) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(-1.83028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(2.485086) q[2];
rz(2.4025829) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2330403) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(-0.064706133) q[0];
rz(-pi) q[1];
rz(2.3357453) q[2];
sx q[2];
rz(-2.5666551) q[2];
sx q[2];
rz(-0.88428674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.060354787) q[1];
sx q[1];
rz(-1.3831257) q[1];
sx q[1];
rz(-0.16138046) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08667605) q[3];
sx q[3];
rz(-1.254734) q[3];
sx q[3];
rz(-2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-2.9188459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6762786) q[0];
sx q[0];
rz(-1.9440117) q[0];
sx q[0];
rz(0.42155427) q[0];
rz(-pi) q[1];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-0.97937102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45914868) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(0.1174121) q[1];
rz(-2.1915216) q[3];
sx q[3];
rz(-2.7347703) q[3];
sx q[3];
rz(2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-2.2213675) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.614914) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(-2.148669) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22569457) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(-0.29633488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4184119) q[1];
sx q[1];
rz(-2.3310117) q[1];
sx q[1];
rz(2.4852738) q[1];
rz(0.50587378) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-2.6884902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720854) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(-1.1016269) q[0];
rz(-1.4278973) q[2];
sx q[2];
rz(-1.5362527) q[2];
sx q[2];
rz(-2.3847716) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.3506538) q[1];
sx q[1];
rz(-0.084053587) q[1];
rz(-0.971332) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(-0.23813716) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3753525) q[0];
sx q[0];
rz(-2.0413114) q[0];
sx q[0];
rz(-1.9348295) q[0];
rz(2.020535) q[2];
sx q[2];
rz(-1.626484) q[2];
sx q[2];
rz(2.4186717) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7102393) q[1];
sx q[1];
rz(-0.87190404) q[1];
sx q[1];
rz(-0.61183521) q[1];
rz(-pi) q[2];
rz(-2.7558277) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(0.39757279) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(1.1869173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531187) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(-0.32062809) q[0];
x q[1];
rz(0.45849623) q[2];
sx q[2];
rz(-0.065932238) q[2];
sx q[2];
rz(0.83930783) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4759051) q[1];
sx q[1];
rz(-2.396282) q[1];
sx q[1];
rz(-1.7687294) q[1];
rz(0.62992923) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(2.2858455) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521342) q[0];
sx q[0];
rz(-1.2297213) q[0];
sx q[0];
rz(2.336691) q[0];
x q[1];
rz(-1.0439992) q[2];
sx q[2];
rz(-1.5287405) q[2];
sx q[2];
rz(0.85925697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6701263) q[1];
sx q[1];
rz(-1.2738436) q[1];
sx q[1];
rz(-2.9832277) q[1];
rz(-pi) q[2];
rz(-2.6183073) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(0.31627396) q[0];
rz(1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(2.0064328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59199698) q[0];
sx q[0];
rz(-1.9599008) q[0];
sx q[0];
rz(-1.367021) q[0];
x q[1];
rz(-0.029149292) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(-2.0707891) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89716298) q[1];
sx q[1];
rz(-1.0883372) q[1];
sx q[1];
rz(3.0753115) q[1];
rz(3.1177403) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(-0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(-0.19616729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88403945) q[0];
sx q[0];
rz(-2.1243874) q[0];
sx q[0];
rz(-0.89417017) q[0];
x q[1];
rz(1.7061383) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(-0.70336715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-0.73426437) q[1];
sx q[1];
rz(-1.5478565) q[1];
rz(-pi) q[2];
rz(-0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(-1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(-0.23399227) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(-1.7668016) q[3];
sx q[3];
rz(-2.9075665) q[3];
sx q[3];
rz(-2.8961765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
