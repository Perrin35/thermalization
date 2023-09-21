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
rz(0.35559911) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(-1.0441138) q[1];
sx q[1];
rz(1.27966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0227538) q[0];
sx q[0];
rz(-1.7099329) q[0];
sx q[0];
rz(1.820178) q[0];
rz(-2.6718164) q[2];
sx q[2];
rz(-2.6758286) q[2];
sx q[2];
rz(-2.2848406) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(-2.5024662) q[1];
rz(1.7765331) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(-0.3068026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(2.4025829) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-0.48746902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4741164) q[0];
sx q[0];
rz(-1.6352909) q[0];
sx q[0];
rz(1.6518031) q[0];
rz(2.7198118) q[2];
sx q[2];
rz(-1.9739208) q[2];
sx q[2];
rz(-0.031907206) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65709719) q[1];
sx q[1];
rz(-2.8946981) q[1];
sx q[1];
rz(-2.2730278) q[1];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.653169) q[3];
sx q[3];
rz(1.9272643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-0.3953735) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(-3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19668002) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(-3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-2.9188459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4091464) q[0];
sx q[0];
rz(-1.961686) q[0];
sx q[0];
rz(1.976165) q[0];
x q[1];
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
x q[0];
rz(-2.1258771) q[1];
sx q[1];
rz(-1.6385957) q[1];
sx q[1];
rz(-0.61342872) q[1];
rz(-pi) q[2];
rz(-1.2336897) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2591851) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(-2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(2.2213675) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.675925) q[0];
sx q[0];
rz(-1.9786069) q[0];
sx q[0];
rz(-0.84665914) q[0];
x q[1];
rz(-1.4843066) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(-0.53777018) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0238266) q[1];
sx q[1];
rz(-2.1823366) q[1];
sx q[1];
rz(1.0002506) q[1];
rz(-pi) q[2];
rz(-2.6357189) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(1.7368332) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3081007) q[0];
sx q[0];
rz(-1.5335324) q[0];
sx q[0];
rz(1.6444071) q[0];
rz(-1.808851) q[2];
sx q[2];
rz(-0.14698725) q[2];
sx q[2];
rz(-2.5631995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0799775) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(-1.7916937) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-2.3036028) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0673163) q[0];
sx q[0];
rz(-0.58642504) q[0];
sx q[0];
rz(-0.61074722) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.020535) q[2];
sx q[2];
rz(-1.626484) q[2];
sx q[2];
rz(-2.4186717) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7442419) q[1];
sx q[1];
rz(-2.2480818) q[1];
sx q[1];
rz(2.1703297) q[1];
x q[2];
rz(-0.38576491) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(-0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.05904077) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1259574) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(-1.9546753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163527) q[0];
sx q[0];
rz(-1.851474) q[0];
sx q[0];
rz(-1.0513845) q[0];
rz(-pi) q[1];
rz(0.45849623) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(-0.83930783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3992607) q[1];
sx q[1];
rz(-2.2982344) q[1];
sx q[1];
rz(-2.9620693) q[1];
rz(-1.8625453) q[3];
sx q[3];
rz(-2.1803133) q[3];
sx q[3];
rz(-2.9677344) q[3];
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
rz(-1.7857893) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.81925201) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(2.1112679) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54661575) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(1.0975518) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0439992) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(0.85925697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(1.8712908) q[1];
rz(-pi) q[2];
rz(2.2780667) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(1.2442148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.29256233) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(-0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(2.0064328) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90056706) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(0.39649773) q[0];
rz(-1.8912002) q[2];
sx q[2];
rz(-1.5984629) q[2];
sx q[2];
rz(-0.49081341) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3865349) q[1];
sx q[1];
rz(-2.6549576) q[1];
sx q[1];
rz(1.6965894) q[1];
rz(1.5376066) q[3];
sx q[3];
rz(-0.6233218) q[3];
sx q[3];
rz(2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5967782) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(0.19616729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0553186) q[0];
sx q[0];
rz(-2.1325169) q[0];
sx q[0];
rz(0.67027153) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0912193) q[2];
sx q[2];
rz(-1.435624) q[2];
sx q[2];
rz(0.86063517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6876467) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(0.83666283) q[1];
x q[2];
rz(2.2979126) q[3];
sx q[3];
rz(-2.0920678) q[3];
sx q[3];
rz(-2.664444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(-1.1817415) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.9076004) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(1.7668016) q[3];
sx q[3];
rz(-0.23402611) q[3];
sx q[3];
rz(0.2454161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
