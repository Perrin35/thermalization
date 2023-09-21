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
rz(1.8619327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049392603) q[0];
sx q[0];
rz(-0.28486262) q[0];
sx q[0];
rz(2.0869135) q[0];
rz(-pi) q[1];
rz(1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(1.3734693) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3763334) q[1];
sx q[1];
rz(-1.9783124) q[1];
sx q[1];
rz(0.63912649) q[1];
rz(-2.9079307) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(1.83028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(-0.59536368) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(0.48746902) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90855234) q[0];
sx q[0];
rz(-1.6516343) q[0];
sx q[0];
rz(-0.064706133) q[0];
x q[1];
rz(2.0080645) q[2];
sx q[2];
rz(-1.184706) q[2];
sx q[2];
rz(-1.3646477) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5408052) q[1];
sx q[1];
rz(-1.729319) q[1];
sx q[1];
rz(1.7608789) q[1];
x q[2];
rz(-1.2536212) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(-1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(0.3953735) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.19668002) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5644154) q[0];
sx q[0];
rz(-2.5860997) q[0];
sx q[0];
rz(-0.76340686) q[0];
rz(-pi) q[1];
rz(1.0388971) q[2];
sx q[2];
rz(-2.7173018) q[2];
sx q[2];
rz(-1.588856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45914868) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(3.0241806) q[1];
x q[2];
rz(1.907903) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(-0.043134886) q[3];
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
rz(-1.1941236) q[2];
rz(1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(-0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.4032723) q[0];
rz(1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-0.9202252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.614914) q[0];
sx q[0];
rz(-0.81256142) q[0];
sx q[0];
rz(2.148669) q[0];
rz(-1.4843066) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(2.6038225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3355545) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(0.6946509) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7761049) q[3];
sx q[3];
rz(-1.218154) q[3];
sx q[3];
rz(0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(0.81930339) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.4047594) q[0];
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
rz(-2.4015409) q[0];
sx q[0];
rz(-1.4972367) q[0];
sx q[0];
rz(0.037365035) q[0];
x q[1];
rz(-3.1066936) q[2];
sx q[2];
rz(-1.4279832) q[2];
sx q[2];
rz(-2.3226483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52757712) q[1];
sx q[1];
rz(-1.3506538) q[1];
sx q[1];
rz(0.084053587) q[1];
rz(-pi) q[2];
rz(2.1987678) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(0.065746106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(-0.72367469) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-2.2912912) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1171922) q[0];
sx q[0];
rz(-1.8937366) q[0];
sx q[0];
rz(-0.49844235) q[0];
x q[1];
rz(1.1210576) q[2];
sx q[2];
rz(-1.5151086) q[2];
sx q[2];
rz(2.4186717) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7442419) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(-0.97126295) q[1];
rz(-2.7558277) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.05904077) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(0.13051662) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966973) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(2.0969735) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0824532) q[2];
sx q[2];
rz(-1.5416317) q[2];
sx q[2];
rz(1.189122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.051441593) q[1];
sx q[1];
rz(-1.4370343) q[1];
sx q[1];
rz(0.83530463) q[1];
rz(-0.39077057) q[3];
sx q[3];
rz(-0.6676538) q[3];
sx q[3];
rz(2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54661575) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(2.0440408) q[0];
rz(-pi) q[1];
rz(2.0975935) q[2];
sx q[2];
rz(-1.5287405) q[2];
sx q[2];
rz(-2.2823357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(1.8712908) q[1];
rz(-pi) q[2];
rz(0.86352591) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(-1.2442148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8490303) q[2];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(2.0064328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508182) q[0];
sx q[0];
rz(-2.7047815) q[0];
sx q[0];
rz(2.6831021) q[0];
rz(-pi) q[1];
rz(1.4831545) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(1.9784387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2444297) q[1];
sx q[1];
rz(-2.0532554) q[1];
sx q[1];
rz(0.06628118) q[1];
x q[2];
rz(2.1938571) q[3];
sx q[3];
rz(-1.590168) q[3];
sx q[3];
rz(-1.7957578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(0.35183516) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(-2.9454254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575532) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(-2.2474225) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.050373366) q[2];
sx q[2];
rz(-1.435624) q[2];
sx q[2];
rz(0.86063517) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4539459) q[1];
sx q[1];
rz(-1.5861662) q[1];
sx q[1];
rz(2.3049298) q[1];
rz(-0.65532834) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(-0.67656803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(-1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(2.0666163) q[3];
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
rz(-pi) q[0];
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
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(-2.9076004) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(1.3410939) q[3];
sx q[3];
rz(-1.615974) q[3];
sx q[3];
rz(1.6254197) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];