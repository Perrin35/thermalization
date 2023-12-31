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
rz(-0.58334938) q[0];
sx q[0];
rz(-1.8177176) q[0];
sx q[0];
rz(0.14351828) q[0];
x q[1];
rz(-2.6718164) q[2];
sx q[2];
rz(-0.46576408) q[2];
sx q[2];
rz(2.2848406) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0919839) q[1];
sx q[1];
rz(-0.99124747) q[1];
sx q[1];
rz(2.0642573) q[1];
rz(-pi) q[2];
rz(-1.7765331) q[3];
sx q[3];
rz(-1.3418901) q[3];
sx q[3];
rz(-0.3068026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(-0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(2.546229) q[0];
rz(-0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(0.48746902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2330403) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(-3.0768865) q[0];
x q[1];
rz(0.80584731) q[2];
sx q[2];
rz(-0.57493756) q[2];
sx q[2];
rz(-0.88428674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6007874) q[1];
sx q[1];
rz(-1.729319) q[1];
sx q[1];
rz(-1.7608789) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8295733) q[3];
sx q[3];
rz(-0.32734713) q[3];
sx q[3];
rz(-2.5395288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-2.7462192) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-2.9188459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091464) q[0];
sx q[0];
rz(-1.961686) q[0];
sx q[0];
rz(-1.976165) q[0];
rz(2.9163755) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(-2.1622216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45914868) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(3.0241806) q[1];
rz(-1.907903) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(-3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(2.0866701) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(2.2213675) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677778) q[0];
sx q[0];
rz(-0.91705634) q[0];
sx q[0];
rz(-2.6184404) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9158981) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(0.29633488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4184119) q[1];
sx q[1];
rz(-2.3310117) q[1];
sx q[1];
rz(-0.65631887) q[1];
rz(1.7761049) q[3];
sx q[3];
rz(-1.218154) q[3];
sx q[3];
rz(-0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(-0.81930339) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950726) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(-2.0399658) q[0];
rz(0.034899072) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(-0.81894433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9822426) q[1];
sx q[1];
rz(-0.23540007) q[1];
sx q[1];
rz(1.2118641) q[1];
rz(-0.971332) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(-2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-0.23813716) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.5135117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3753525) q[0];
sx q[0];
rz(-1.1002812) q[0];
sx q[0];
rz(1.9348295) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6983301) q[2];
sx q[2];
rz(-0.45293929) q[2];
sx q[2];
rz(-2.4085101) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5781128) q[1];
sx q[1];
rz(-1.1155177) q[1];
sx q[1];
rz(0.77225765) q[1];
rz(-2.2699039) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(0.20760078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(-0.17091621) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(0.39757279) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(1.9546753) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966973) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(-1.0446192) q[0];
rz(0.059139472) q[2];
sx q[2];
rz(-1.5999609) q[2];
sx q[2];
rz(-1.9524706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.742332) q[1];
sx q[1];
rz(-2.2982344) q[1];
sx q[1];
rz(2.9620693) q[1];
rz(-pi) q[2];
rz(2.5116634) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5671974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(0.9986977) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(0.85574714) q[0];
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
rz(-1.94901) q[0];
sx q[0];
rz(-0.858925) q[0];
sx q[0];
rz(0.45759767) q[0];
x q[1];
rz(2.0975935) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(2.2823357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1134539) q[1];
sx q[1];
rz(-2.8061562) q[1];
sx q[1];
rz(1.0949275) q[1];
rz(-2.6183073) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(-0.39213359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29256233) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(0.58644811) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-2.4119222) q[1];
sx q[1];
rz(-1.1351599) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410256) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(-0.39649773) q[0];
rz(-pi) q[1];
rz(-0.029149292) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(2.0707891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.70442048) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(1.0874332) q[1];
x q[2];
rz(-1.6039861) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(-2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(1.5967782) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-2.9454254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575532) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(-0.89417017) q[0];
x q[1];
rz(-1.2162131) q[2];
sx q[2];
rz(-0.14419975) q[2];
sx q[2];
rz(1.2186288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0386104) q[1];
sx q[1];
rz(-0.83676941) q[1];
sx q[1];
rz(-3.1208913) q[1];
rz(-pi) q[2];
rz(0.85828652) q[3];
sx q[3];
rz(-0.86601102) q[3];
sx q[3];
rz(1.5376877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(-2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186196) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(0.20028533) q[2];
sx q[2];
rz(-2.9030048) q[2];
sx q[2];
rz(1.5422274) q[2];
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
