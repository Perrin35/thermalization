OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24902046) q[0];
sx q[0];
rz(5.3244642) q[0];
sx q[0];
rz(9.1992314) q[0];
rz(-1.072999) q[1];
sx q[1];
rz(-2.6838611) q[1];
sx q[1];
rz(-0.20888858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07889121) q[0];
sx q[0];
rz(-1.3008504) q[0];
sx q[0];
rz(0.45227082) q[0];
x q[1];
rz(-1.8525847) q[2];
sx q[2];
rz(-2.1827841) q[2];
sx q[2];
rz(-0.91152292) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1465197) q[1];
sx q[1];
rz(-2.0834298) q[1];
sx q[1];
rz(-1.9216838) q[1];
x q[2];
rz(-2.8257583) q[3];
sx q[3];
rz(-1.8973992) q[3];
sx q[3];
rz(2.0706319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0245584) q[2];
sx q[2];
rz(-2.890675) q[2];
sx q[2];
rz(-2.2158902) q[2];
rz(0.77445817) q[3];
sx q[3];
rz(-1.9175994) q[3];
sx q[3];
rz(1.2831877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9722612) q[0];
sx q[0];
rz(-0.65003482) q[0];
sx q[0];
rz(-1.4647123) q[0];
rz(-0.88893923) q[1];
sx q[1];
rz(-2.2496532) q[1];
sx q[1];
rz(1.3533786) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037925) q[0];
sx q[0];
rz(-0.97957459) q[0];
sx q[0];
rz(-1.2205489) q[0];
rz(1.5487461) q[2];
sx q[2];
rz(-1.2047909) q[2];
sx q[2];
rz(-0.61701202) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6135237) q[1];
sx q[1];
rz(-0.98586833) q[1];
sx q[1];
rz(1.8208353) q[1];
rz(-pi) q[2];
rz(-1.8436271) q[3];
sx q[3];
rz(-2.3300814) q[3];
sx q[3];
rz(-1.2173513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9951524) q[2];
sx q[2];
rz(-2.1285987) q[2];
sx q[2];
rz(-2.1688482) q[2];
rz(1.3257596) q[3];
sx q[3];
rz(-1.9917859) q[3];
sx q[3];
rz(0.61746922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4165118) q[0];
sx q[0];
rz(-1.4840115) q[0];
sx q[0];
rz(-3.1120279) q[0];
rz(1.0211396) q[1];
sx q[1];
rz(-0.28426668) q[1];
sx q[1];
rz(-0.31612843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164564) q[0];
sx q[0];
rz(-3.1119149) q[0];
sx q[0];
rz(-0.82334955) q[0];
rz(-pi) q[1];
rz(-1.3695551) q[2];
sx q[2];
rz(-1.6627208) q[2];
sx q[2];
rz(2.6215009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8384605) q[1];
sx q[1];
rz(-2.7632398) q[1];
sx q[1];
rz(1.0835453) q[1];
x q[2];
rz(1.6484327) q[3];
sx q[3];
rz(-2.2901553) q[3];
sx q[3];
rz(-0.84859728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11008392) q[2];
sx q[2];
rz(-1.2981334) q[2];
sx q[2];
rz(0.52424866) q[2];
rz(-0.60018572) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(2.0704827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780739) q[0];
sx q[0];
rz(-1.9215895) q[0];
sx q[0];
rz(-2.779261) q[0];
rz(2.433297) q[1];
sx q[1];
rz(-2.7999122) q[1];
sx q[1];
rz(-1.1452311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6714685) q[0];
sx q[0];
rz(-0.13397476) q[0];
sx q[0];
rz(2.2132316) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87605642) q[2];
sx q[2];
rz(-1.183702) q[2];
sx q[2];
rz(0.53384483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.143024) q[1];
sx q[1];
rz(-2.2557229) q[1];
sx q[1];
rz(-2.8110503) q[1];
rz(-1.6207375) q[3];
sx q[3];
rz(-1.7347214) q[3];
sx q[3];
rz(-2.0633179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6142673) q[2];
sx q[2];
rz(-0.82794398) q[2];
sx q[2];
rz(3.1329727) q[2];
rz(2.4873867) q[3];
sx q[3];
rz(-0.71627408) q[3];
sx q[3];
rz(-2.8692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24097405) q[0];
sx q[0];
rz(-2.8130377) q[0];
sx q[0];
rz(-0.7884489) q[0];
rz(1.0139326) q[1];
sx q[1];
rz(-1.3194111) q[1];
sx q[1];
rz(-3.0533155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006314) q[0];
sx q[0];
rz(-2.2128803) q[0];
sx q[0];
rz(-0.93812801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9806251) q[2];
sx q[2];
rz(-2.5720398) q[2];
sx q[2];
rz(-0.19191027) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8535117) q[1];
sx q[1];
rz(-1.2895023) q[1];
sx q[1];
rz(-2.0916677) q[1];
rz(-pi) q[2];
rz(1.5203736) q[3];
sx q[3];
rz(-2.6459624) q[3];
sx q[3];
rz(-0.69678604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91904116) q[2];
sx q[2];
rz(-1.4107979) q[2];
sx q[2];
rz(0.50773531) q[2];
rz(-0.12568411) q[3];
sx q[3];
rz(-1.9583743) q[3];
sx q[3];
rz(-0.79508933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3313726) q[0];
sx q[0];
rz(-0.75646821) q[0];
sx q[0];
rz(3.0561225) q[0];
rz(2.5147009) q[1];
sx q[1];
rz(-1.9848928) q[1];
sx q[1];
rz(1.289182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9753646) q[0];
sx q[0];
rz(-1.3718107) q[0];
sx q[0];
rz(0.10597845) q[0];
x q[1];
rz(-1.7250546) q[2];
sx q[2];
rz(-2.3084894) q[2];
sx q[2];
rz(1.2546033) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8823008) q[1];
sx q[1];
rz(-1.8083824) q[1];
sx q[1];
rz(-0.17452328) q[1];
rz(-pi) q[2];
rz(1.3335532) q[3];
sx q[3];
rz(-1.5565775) q[3];
sx q[3];
rz(2.1923101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36279303) q[2];
sx q[2];
rz(-0.94393602) q[2];
sx q[2];
rz(1.8865406) q[2];
rz(-3.0826027) q[3];
sx q[3];
rz(-1.2041644) q[3];
sx q[3];
rz(0.98439938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4208218) q[0];
sx q[0];
rz(-1.2728007) q[0];
sx q[0];
rz(2.3765748) q[0];
rz(1.4014686) q[1];
sx q[1];
rz(-0.88686371) q[1];
sx q[1];
rz(2.9041362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260587) q[0];
sx q[0];
rz(-2.1977673) q[0];
sx q[0];
rz(-3.0978908) q[0];
rz(-pi) q[1];
rz(1.6638905) q[2];
sx q[2];
rz(-2.4590465) q[2];
sx q[2];
rz(-2.1381706) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.69370338) q[1];
sx q[1];
rz(-0.2940601) q[1];
sx q[1];
rz(-2.990688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6066437) q[3];
sx q[3];
rz(-0.69679835) q[3];
sx q[3];
rz(0.076059503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6852297) q[2];
sx q[2];
rz(-0.80358973) q[2];
sx q[2];
rz(-2.2881499) q[2];
rz(1.338909) q[3];
sx q[3];
rz(-1.5693376) q[3];
sx q[3];
rz(-2.9760823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6859739) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(0.18534216) q[0];
rz(-2.7475884) q[1];
sx q[1];
rz(-1.3961671) q[1];
sx q[1];
rz(-1.0967163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7506517) q[0];
sx q[0];
rz(-2.0454198) q[0];
sx q[0];
rz(0.43584521) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11285891) q[2];
sx q[2];
rz(-0.90870171) q[2];
sx q[2];
rz(2.2619132) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85307676) q[1];
sx q[1];
rz(-0.6986874) q[1];
sx q[1];
rz(-0.53650155) q[1];
x q[2];
rz(0.59456749) q[3];
sx q[3];
rz(-1.2487186) q[3];
sx q[3];
rz(-0.26628029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77157053) q[2];
sx q[2];
rz(-3.0574419) q[2];
sx q[2];
rz(2.8328075) q[2];
rz(1.2540981) q[3];
sx q[3];
rz(-1.2718688) q[3];
sx q[3];
rz(2.0103644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9517188) q[0];
sx q[0];
rz(-1.254344) q[0];
sx q[0];
rz(0.21981123) q[0];
rz(-1.1687357) q[1];
sx q[1];
rz(-1.7857779) q[1];
sx q[1];
rz(1.8692325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74368661) q[0];
sx q[0];
rz(-2.1748161) q[0];
sx q[0];
rz(-0.029725909) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4859305) q[2];
sx q[2];
rz(-2.0800635) q[2];
sx q[2];
rz(-2.9743663) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90779623) q[1];
sx q[1];
rz(-1.4005473) q[1];
sx q[1];
rz(0.40791224) q[1];
rz(-pi) q[2];
rz(2.6401889) q[3];
sx q[3];
rz(-1.1165757) q[3];
sx q[3];
rz(2.8387031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9610338) q[2];
sx q[2];
rz(-0.84332931) q[2];
sx q[2];
rz(-0.45832222) q[2];
rz(0.35442963) q[3];
sx q[3];
rz(-1.3684045) q[3];
sx q[3];
rz(1.9663158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4273222) q[0];
sx q[0];
rz(-1.2027807) q[0];
sx q[0];
rz(0.74139968) q[0];
rz(1.7085913) q[1];
sx q[1];
rz(-0.42867908) q[1];
sx q[1];
rz(-0.0892078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36486711) q[0];
sx q[0];
rz(-1.3898384) q[0];
sx q[0];
rz(-2.0511766) q[0];
x q[1];
rz(1.9544425) q[2];
sx q[2];
rz(-2.128278) q[2];
sx q[2];
rz(1.5889744) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20232378) q[1];
sx q[1];
rz(-1.9096194) q[1];
sx q[1];
rz(-3.0975049) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44762917) q[3];
sx q[3];
rz(-2.3530686) q[3];
sx q[3];
rz(2.727689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0006813) q[2];
sx q[2];
rz(-0.050364308) q[2];
sx q[2];
rz(0.64765206) q[2];
rz(0.40376136) q[3];
sx q[3];
rz(-1.253456) q[3];
sx q[3];
rz(-3.0103179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972926) q[0];
sx q[0];
rz(-0.97351749) q[0];
sx q[0];
rz(1.0607006) q[0];
rz(-2.3464959) q[1];
sx q[1];
rz(-0.92082321) q[1];
sx q[1];
rz(-0.77650741) q[1];
rz(-3.0349229) q[2];
sx q[2];
rz(-2.223458) q[2];
sx q[2];
rz(-1.3724422) q[2];
rz(2.285801) q[3];
sx q[3];
rz(-0.11920155) q[3];
sx q[3];
rz(-3.0536065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
