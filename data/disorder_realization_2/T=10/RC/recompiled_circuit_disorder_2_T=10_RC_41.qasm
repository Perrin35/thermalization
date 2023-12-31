OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(-2.5180106) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5335124) q[0];
sx q[0];
rz(-0.97097662) q[0];
sx q[0];
rz(-1.5754726) q[0];
rz(2.1262769) q[2];
sx q[2];
rz(-0.17586389) q[2];
sx q[2];
rz(1.6989087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9762293) q[1];
sx q[1];
rz(-0.86127087) q[1];
sx q[1];
rz(2.707259) q[1];
x q[2];
rz(2.9849103) q[3];
sx q[3];
rz(-2.0041668) q[3];
sx q[3];
rz(2.3436848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3502675) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(-1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7704849) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(1.0789385) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(2.4761377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6404214) q[0];
sx q[0];
rz(-1.8542395) q[0];
sx q[0];
rz(2.6675176) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3999248) q[2];
sx q[2];
rz(-1.8881646) q[2];
sx q[2];
rz(-2.3651809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7484819) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(1.837681) q[1];
rz(-pi) q[2];
rz(3.1075675) q[3];
sx q[3];
rz(-2.2799387) q[3];
sx q[3];
rz(2.2505086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5370496) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(-1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(2.4501734) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.8483775) q[1];
sx q[1];
rz(-1.4000777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83126691) q[0];
sx q[0];
rz(-1.184549) q[0];
sx q[0];
rz(-1.2826184) q[0];
rz(2.954735) q[2];
sx q[2];
rz(-1.4374905) q[2];
sx q[2];
rz(1.2847628) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6058265) q[1];
sx q[1];
rz(-1.2488135) q[1];
sx q[1];
rz(0.95226007) q[1];
rz(-pi) q[2];
rz(-1.8758043) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(-2.4839694) q[2];
rz(1.1714606) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-0.79384971) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.4452274) q[0];
rz(1.4472648) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(0.34805527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570243) q[0];
sx q[0];
rz(-0.94385249) q[0];
sx q[0];
rz(1.2059962) q[0];
x q[1];
rz(1.9449171) q[2];
sx q[2];
rz(-1.4305563) q[2];
sx q[2];
rz(1.8082878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.12236774) q[1];
sx q[1];
rz(-2.5248563) q[1];
sx q[1];
rz(-1.7255746) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69339852) q[3];
sx q[3];
rz(-2.9326673) q[3];
sx q[3];
rz(-0.25017504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(2.7187738) q[2];
rz(2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(-2.6111531) q[0];
rz(2.2166705) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(-1.2984498) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7206551) q[0];
sx q[0];
rz(-1.4268095) q[0];
sx q[0];
rz(1.7341341) q[0];
rz(0.41579397) q[2];
sx q[2];
rz(-0.64877629) q[2];
sx q[2];
rz(2.6598038) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4106771) q[1];
sx q[1];
rz(-1.2123322) q[1];
sx q[1];
rz(1.6770384) q[1];
rz(-pi) q[2];
rz(-0.88921806) q[3];
sx q[3];
rz(-0.71435706) q[3];
sx q[3];
rz(2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9412781) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-0.47362622) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(2.30106) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6376003) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(0.9978869) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(2.6748437) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40690639) q[0];
sx q[0];
rz(-0.91306251) q[0];
sx q[0];
rz(1.8440194) q[0];
x q[1];
rz(0.45852197) q[2];
sx q[2];
rz(-0.47007559) q[2];
sx q[2];
rz(0.93570659) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0135865) q[1];
sx q[1];
rz(-0.15833536) q[1];
sx q[1];
rz(2.8360785) q[1];
rz(2.09957) q[3];
sx q[3];
rz(-2.9526132) q[3];
sx q[3];
rz(-0.46940645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(-1.5647282) q[2];
rz(2.5148897) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3395183) q[0];
sx q[0];
rz(-3.0537362) q[0];
sx q[0];
rz(0.095803424) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19303796) q[2];
sx q[2];
rz(-1.2912573) q[2];
sx q[2];
rz(-0.9370196) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0737887) q[1];
sx q[1];
rz(-1.5556591) q[1];
sx q[1];
rz(3.036036) q[1];
rz(2.0358622) q[3];
sx q[3];
rz(-1.8814058) q[3];
sx q[3];
rz(0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5856813) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(2.3445271) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(-1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(-2.8523493) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(1.3141059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1159191) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(-1.7220108) q[0];
rz(0.35347519) q[2];
sx q[2];
rz(-0.77958737) q[2];
sx q[2];
rz(1.0970864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1738759) q[1];
sx q[1];
rz(-1.516725) q[1];
sx q[1];
rz(-1.6915583) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81359158) q[3];
sx q[3];
rz(-1.7980669) q[3];
sx q[3];
rz(0.0031301216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(1.7129664) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6190417) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(-1.8956986) q[0];
rz(3.030581) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(-0.54661173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0194191) q[0];
sx q[0];
rz(-2.008736) q[0];
sx q[0];
rz(-0.41284783) q[0];
x q[1];
rz(-1.941628) q[2];
sx q[2];
rz(-1.2002581) q[2];
sx q[2];
rz(-1.759699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48601549) q[1];
sx q[1];
rz(-1.3994819) q[1];
sx q[1];
rz(-1.9678712) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9162174) q[3];
sx q[3];
rz(-1.1715874) q[3];
sx q[3];
rz(-2.6938714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0299915) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.4477504) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(0.70770121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74117888) q[0];
sx q[0];
rz(-0.78010633) q[0];
sx q[0];
rz(1.0723423) q[0];
rz(1.3466481) q[2];
sx q[2];
rz(-1.4321616) q[2];
sx q[2];
rz(-2.5773406) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0545132) q[1];
sx q[1];
rz(-0.88216773) q[1];
sx q[1];
rz(0.41333945) q[1];
rz(-pi) q[2];
rz(2.1901312) q[3];
sx q[3];
rz(-1.8639495) q[3];
sx q[3];
rz(-0.64630634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(-1.5385657) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(0.8159591) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(-2.3564561) q[2];
sx q[2];
rz(-2.5265836) q[2];
sx q[2];
rz(2.1267736) q[2];
rz(3.0872185) q[3];
sx q[3];
rz(-1.5297223) q[3];
sx q[3];
rz(-2.3755304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
