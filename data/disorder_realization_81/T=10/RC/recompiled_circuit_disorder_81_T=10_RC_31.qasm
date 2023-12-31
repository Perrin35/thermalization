OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(4.9217304) q[0];
sx q[0];
rz(11.187727) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987797) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(0.16422693) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3583463) q[2];
sx q[2];
rz(-0.48848402) q[2];
sx q[2];
rz(-2.8993895) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3528459) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-0.84428865) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-0.60602337) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79016722) q[0];
sx q[0];
rz(-1.5520436) q[0];
sx q[0];
rz(3.0856531) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8430014) q[2];
sx q[2];
rz(-2.9332187) q[2];
sx q[2];
rz(-1.5339472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9895049) q[1];
sx q[1];
rz(-0.67645914) q[1];
sx q[1];
rz(1.771404) q[1];
rz(-pi) q[2];
rz(-2.9119592) q[3];
sx q[3];
rz(-2.4335055) q[3];
sx q[3];
rz(1.1585483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.9821232) q[2];
rz(0.37108478) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(-2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7263111) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(-2.2583654) q[0];
rz(-2.9224612) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(-2.6205274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1097475) q[1];
sx q[1];
rz(-2.577563) q[1];
sx q[1];
rz(2.2721223) q[1];
x q[2];
rz(0.1180325) q[3];
sx q[3];
rz(-1.1088088) q[3];
sx q[3];
rz(0.0052009728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(-1.8211676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5399649) q[0];
sx q[0];
rz(-0.40973445) q[0];
sx q[0];
rz(-1.5967303) q[0];
rz(2.690372) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(-2.1860683) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79975407) q[1];
sx q[1];
rz(-0.8952039) q[1];
sx q[1];
rz(-1.9104596) q[1];
rz(-pi) q[2];
rz(0.60453316) q[3];
sx q[3];
rz(-2.2570838) q[3];
sx q[3];
rz(-1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(0.17677447) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(1.140973) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-0.86529055) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(1.3006166) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039935) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(0.07261891) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6426716) q[2];
sx q[2];
rz(-1.2051029) q[2];
sx q[2];
rz(0.14771151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3805441) q[1];
sx q[1];
rz(-0.64861464) q[1];
sx q[1];
rz(2.5785239) q[1];
x q[2];
rz(0.50118581) q[3];
sx q[3];
rz(-0.29705829) q[3];
sx q[3];
rz(-2.3230769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(0.011750301) q[0];
rz(-2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(-1.5531497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1588622) q[0];
sx q[0];
rz(-1.2356865) q[0];
sx q[0];
rz(1.9431252) q[0];
x q[1];
rz(1.9366829) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(-2.6423955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8223871) q[1];
sx q[1];
rz(-2.2777595) q[1];
sx q[1];
rz(-1.6824526) q[1];
rz(-2.2736069) q[3];
sx q[3];
rz(-1.6276974) q[3];
sx q[3];
rz(-1.9753319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-0.67725956) q[0];
rz(2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-2.1645434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517294) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(-0.59318869) q[0];
rz(-pi) q[1];
rz(-1.5708959) q[2];
sx q[2];
rz(-1.4387555) q[2];
sx q[2];
rz(0.11167234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6923185) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(1.3432137) q[1];
x q[2];
rz(-1.2491751) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(-1.4332989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3892422) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4417393) q[0];
sx q[0];
rz(-0.2158567) q[0];
sx q[0];
rz(-0.21780832) q[0];
rz(-2.6685733) q[2];
sx q[2];
rz(-2.0447391) q[2];
sx q[2];
rz(-0.75616403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3932712) q[1];
sx q[1];
rz(-1.5763348) q[1];
sx q[1];
rz(0.022151532) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9551679) q[3];
sx q[3];
rz(-0.66925183) q[3];
sx q[3];
rz(-0.8347019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32968783) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(-1.1784941) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(-1.8485803) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(2.5440149) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4353838) q[0];
sx q[0];
rz(-1.4357114) q[0];
sx q[0];
rz(-0.038890966) q[0];
rz(-pi) q[1];
rz(-2.3915646) q[2];
sx q[2];
rz(-2.4366597) q[2];
sx q[2];
rz(0.36545576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55652009) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(2.3162566) q[1];
rz(-pi) q[2];
rz(-0.58595539) q[3];
sx q[3];
rz(-1.5065985) q[3];
sx q[3];
rz(2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(0.67217174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32165747) q[0];
sx q[0];
rz(-1.5655087) q[0];
sx q[0];
rz(1.5810285) q[0];
x q[1];
rz(-1.5615084) q[2];
sx q[2];
rz(-1.0592959) q[2];
sx q[2];
rz(1.4052504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8476665) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(0.49077175) q[1];
rz(-0.32398128) q[3];
sx q[3];
rz(-0.45352648) q[3];
sx q[3];
rz(1.6739664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(-1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-0.72369408) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(0.15974113) q[2];
sx q[2];
rz(-1.1764871) q[2];
sx q[2];
rz(2.7574678) q[2];
rz(1.1643812) q[3];
sx q[3];
rz(-1.2953399) q[3];
sx q[3];
rz(-3.0084707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
