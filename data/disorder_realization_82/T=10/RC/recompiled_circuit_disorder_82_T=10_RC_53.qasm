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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1188388) q[0];
sx q[0];
rz(-1.4316598) q[0];
sx q[0];
rz(-1.3214146) q[0];
rz(-pi) q[1];
rz(-0.46977629) q[2];
sx q[2];
rz(-0.46576408) q[2];
sx q[2];
rz(-2.2848406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3763334) q[1];
sx q[1];
rz(-1.9783124) q[1];
sx q[1];
rz(-2.5024662) q[1];
x q[2];
rz(-0.7198556) q[3];
sx q[3];
rz(-0.30656439) q[3];
sx q[3];
rz(-0.43678624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(2.546229) q[0];
rz(-3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4741164) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(1.6518031) q[0];
rz(-0.42178085) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(-3.1096854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5408052) q[1];
sx q[1];
rz(-1.4122737) q[1];
sx q[1];
rz(-1.7608789) q[1];
x q[2];
rz(3.0549166) q[3];
sx q[3];
rz(-1.254734) q[3];
sx q[3];
rz(-2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(0.3953735) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-0.22274676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4653141) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(2.7200384) q[0];
x q[1];
rz(1.0388971) q[2];
sx q[2];
rz(-2.7173018) q[2];
sx q[2];
rz(1.5527366) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.682444) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(-0.1174121) q[1];
rz(2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2591851) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(-2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(0.96364337) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(-1.7383204) q[0];
rz(-1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738149) q[0];
sx q[0];
rz(-0.91705634) q[0];
sx q[0];
rz(0.52315229) q[0];
rz(-pi) q[1];
rz(-0.36107365) q[2];
sx q[2];
rz(-1.6517342) q[2];
sx q[2];
rz(1.0635478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.72318075) q[1];
sx q[1];
rz(-2.3310117) q[1];
sx q[1];
rz(-0.65631887) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50587378) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(-0.94569262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9437287) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(-1.7839446) q[3];
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
x q[0];
rz(-pi/2) q[0];
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
rz(1.7368332) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-2.6884902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74005175) q[0];
sx q[0];
rz(-1.4972367) q[0];
sx q[0];
rz(-3.1042276) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1066936) q[2];
sx q[2];
rz(-1.4279832) q[2];
sx q[2];
rz(-0.81894433) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6140155) q[1];
sx q[1];
rz(-1.7909389) q[1];
sx q[1];
rz(3.0575391) q[1];
rz(-pi) q[2];
rz(-0.94282486) q[3];
sx q[3];
rz(-0.70126611) q[3];
sx q[3];
rz(-0.065746106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(1.6298693) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.5135117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1171922) q[0];
sx q[0];
rz(-1.2478561) q[0];
sx q[0];
rz(-0.49844235) q[0];
rz(-pi) q[1];
rz(-1.6983301) q[2];
sx q[2];
rz(-0.45293929) q[2];
sx q[2];
rz(-2.4085101) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7442419) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(-2.1703297) q[1];
rz(0.38576491) q[3];
sx q[3];
rz(-1.1499975) q[3];
sx q[3];
rz(0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0156353) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(-1.827084) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.9546753) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163527) q[0];
sx q[0];
rz(-1.851474) q[0];
sx q[0];
rz(-2.0902082) q[0];
x q[1];
rz(0.45849623) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(-0.83930783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3992607) q[1];
sx q[1];
rz(-2.2982344) q[1];
sx q[1];
rz(2.9620693) q[1];
x q[2];
rz(0.62992923) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.94901) q[0];
sx q[0];
rz(-0.858925) q[0];
sx q[0];
rz(-2.683995) q[0];
x q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.7221754) q[1];
sx q[1];
rz(1.8712908) q[1];
rz(-pi) q[2];
rz(-0.86352591) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(1.2442148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(0.070177468) q[2];
rz(-0.82693806) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(-0.31627396) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(1.1351599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410256) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(-2.7450949) q[0];
x q[1];
rz(-1.6584381) q[2];
sx q[2];
rz(-2.8200375) q[2];
sx q[2];
rz(-1.9784387) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89716298) q[1];
sx q[1];
rz(-2.0532554) q[1];
sx q[1];
rz(-0.06628118) q[1];
rz(-pi) q[2];
x q[2];
rz(0.023852392) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6241374) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
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
rz(-2.7668178) q[1];
sx q[1];
rz(-2.9454254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(0.79191533) q[0];
rz(-pi) q[1];
rz(-1.4354544) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(-0.70336715) q[2];
rz(-pi) q[3];
x q[3];
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
rz(1.5937362) q[1];
rz(0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(0.47714864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-3.0604559) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(0.20028533) q[2];
sx q[2];
rz(-2.9030048) q[2];
sx q[2];
rz(1.5422274) q[2];
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
