OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(-0.52283302) q[0];
sx q[0];
rz(-0.62358207) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5335124) q[0];
sx q[0];
rz(-0.97097662) q[0];
sx q[0];
rz(-1.5661201) q[0];
rz(3.0481553) q[2];
sx q[2];
rz(-1.7200025) q[2];
sx q[2];
rz(2.0051533) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1653633) q[1];
sx q[1];
rz(-2.2803218) q[1];
sx q[1];
rz(2.707259) q[1];
rz(-1.1327098) q[3];
sx q[3];
rz(-1.4287018) q[3];
sx q[3];
rz(-0.70664584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7913251) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(-1.4482927) q[3];
sx q[3];
rz(-0.99213123) q[3];
sx q[3];
rz(-2.1712415) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(1.0789385) q[0];
rz(1.3868015) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(2.4761377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5011713) q[0];
sx q[0];
rz(-1.8542395) q[0];
sx q[0];
rz(-0.47407504) q[0];
rz(-1.3999248) q[2];
sx q[2];
rz(-1.8881646) q[2];
sx q[2];
rz(-2.3651809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.39311073) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(1.3039116) q[1];
x q[2];
rz(-1.5311702) q[3];
sx q[3];
rz(-2.431776) q[3];
sx q[3];
rz(-0.94330793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60454303) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(3.0664505) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(2.4501734) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103257) q[0];
sx q[0];
rz(-1.9570436) q[0];
sx q[0];
rz(1.8589742) q[0];
x q[1];
rz(2.5163469) q[2];
sx q[2];
rz(-2.9125104) q[2];
sx q[2];
rz(0.89876995) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5249467) q[1];
sx q[1];
rz(-2.4541306) q[1];
sx q[1];
rz(2.0928659) q[1];
x q[2];
rz(-1.0812976) q[3];
sx q[3];
rz(-2.7989417) q[3];
sx q[3];
rz(-3.0150068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65163461) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(1.970132) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30663438) q[0];
sx q[0];
rz(-1.8638896) q[0];
sx q[0];
rz(-0.65960633) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9910827) q[2];
sx q[2];
rz(-1.2005271) q[2];
sx q[2];
rz(-2.9589257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0749803) q[1];
sx q[1];
rz(-2.1790824) q[1];
sx q[1];
rz(3.0327256) q[1];
rz(-pi) q[2];
rz(0.69339852) q[3];
sx q[3];
rz(-0.20892538) q[3];
sx q[3];
rz(0.25017504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7248914) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-2.7187738) q[2];
rz(-2.4041798) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(-2.6111531) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.8431429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4209375) q[0];
sx q[0];
rz(-1.4268095) q[0];
sx q[0];
rz(-1.7341341) q[0];
rz(1.8680044) q[2];
sx q[2];
rz(-2.1564335) q[2];
sx q[2];
rz(-0.98779087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0188705) q[1];
sx q[1];
rz(-1.4713305) q[1];
sx q[1];
rz(0.36032569) q[1];
x q[2];
rz(-0.88921806) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(1.0970955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-0.47362622) q[2];
rz(-0.099362699) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6376003) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(0.9978869) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(2.6748437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186819) q[0];
sx q[0];
rz(-2.4372299) q[0];
sx q[0];
rz(-2.8055311) q[0];
x q[1];
rz(-2.7141063) q[2];
sx q[2];
rz(-1.3689405) q[2];
sx q[2];
rz(2.0919378) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0135865) q[1];
sx q[1];
rz(-2.9832573) q[1];
sx q[1];
rz(-2.8360785) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0454037) q[3];
sx q[3];
rz(-1.4078762) q[3];
sx q[3];
rz(-0.067226203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85990396) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(1.5647282) q[2];
rz(0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307813) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(2.7217641) q[0];
rz(-2.9201674) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5931169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80207434) q[0];
sx q[0];
rz(-0.087856494) q[0];
sx q[0];
rz(-0.095803424) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19303796) q[2];
sx q[2];
rz(-1.8503354) q[2];
sx q[2];
rz(0.9370196) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.78649) q[1];
sx q[1];
rz(-3.0349602) q[1];
sx q[1];
rz(-2.998888) q[1];
rz(-pi) q[2];
rz(2.1920131) q[3];
sx q[3];
rz(-0.55286828) q[3];
sx q[3];
rz(2.051193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5856813) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(-1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(-2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.8274868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24459141) q[0];
sx q[0];
rz(-0.76652157) q[0];
sx q[0];
rz(-2.9826829) q[0];
x q[1];
rz(-1.2411225) q[2];
sx q[2];
rz(-2.2909819) q[2];
sx q[2];
rz(0.61851293) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96771679) q[1];
sx q[1];
rz(-1.516725) q[1];
sx q[1];
rz(1.6915583) q[1];
rz(-pi) q[2];
rz(1.2460327) q[3];
sx q[3];
rz(-2.357558) q[3];
sx q[3];
rz(-1.3336381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(-1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6190417) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(1.2458941) q[0];
rz(-0.11101162) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(2.5949809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36771691) q[0];
sx q[0];
rz(-1.9426632) q[0];
sx q[0];
rz(2.0433776) q[0];
rz(2.7466752) q[2];
sx q[2];
rz(-1.9153708) q[2];
sx q[2];
rz(-0.049023703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48601549) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(-1.9678712) q[1];
x q[2];
rz(-0.48874493) q[3];
sx q[3];
rz(-0.97526032) q[3];
sx q[3];
rz(0.83317703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.6938422) q[2];
rz(1.1374121) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.85957134) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-0.71722537) q[0];
rz(1.9316797) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(0.70770121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74117888) q[0];
sx q[0];
rz(-0.78010633) q[0];
sx q[0];
rz(-2.0692503) q[0];
rz(-pi) q[1];
x q[1];
rz(2.131358) q[2];
sx q[2];
rz(-0.26294225) q[2];
sx q[2];
rz(-1.5901142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.087079436) q[1];
sx q[1];
rz(-2.2594249) q[1];
sx q[1];
rz(2.7282532) q[1];
x q[2];
rz(-0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(0.64630634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.306504) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(2.3564561) q[2];
sx q[2];
rz(-0.61500906) q[2];
sx q[2];
rz(-1.014819) q[2];
rz(-0.054374183) q[3];
sx q[3];
rz(-1.5297223) q[3];
sx q[3];
rz(-2.3755304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];