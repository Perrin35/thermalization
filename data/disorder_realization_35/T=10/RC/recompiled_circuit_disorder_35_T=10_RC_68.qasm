OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(1.024363) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3496075) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(2.4277359) q[0];
x q[1];
rz(0.77180441) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(2.8231205) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6572666) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(2.3372997) q[1];
rz(1.4952881) q[3];
sx q[3];
rz(-1.2592053) q[3];
sx q[3];
rz(-2.989245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99012016) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(2.5734076) q[0];
rz(-pi) q[1];
rz(0.77387626) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(0.80640031) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2270826) q[1];
sx q[1];
rz(-2.7343379) q[1];
sx q[1];
rz(1.6362908) q[1];
x q[2];
rz(-2.5492937) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(-2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(3.0569055) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.9975196) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(-2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(0.57317615) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5432376) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(1.7217365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2531883) q[2];
sx q[2];
rz(-0.8237969) q[2];
sx q[2];
rz(-2.2044646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14904505) q[1];
sx q[1];
rz(-0.48679513) q[1];
sx q[1];
rz(-2.4942314) q[1];
x q[2];
rz(2.2872541) q[3];
sx q[3];
rz(-2.5594098) q[3];
sx q[3];
rz(0.69609387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1406143) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(-2.5909246) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013214839) q[0];
sx q[0];
rz(-0.88071874) q[0];
sx q[0];
rz(2.8280558) q[0];
x q[1];
rz(1.5393125) q[2];
sx q[2];
rz(-1.5523947) q[2];
sx q[2];
rz(3.1061663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4094761) q[1];
sx q[1];
rz(-0.64142694) q[1];
sx q[1];
rz(-2.9684121) q[1];
rz(-3.0691931) q[3];
sx q[3];
rz(-2.1941333) q[3];
sx q[3];
rz(2.6423531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.682122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0844903) q[0];
sx q[0];
rz(-1.717289) q[0];
sx q[0];
rz(1.811972) q[0];
x q[1];
rz(2.1826477) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(-0.81629717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0046878) q[1];
sx q[1];
rz(-1.3740174) q[1];
sx q[1];
rz(-2.911527) q[1];
rz(-pi) q[2];
rz(2.8206283) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(-0.62063673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-0.84189502) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(-1.9110514) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16620557) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(0.020394527) q[0];
rz(2.9406592) q[2];
sx q[2];
rz(-1.4954508) q[2];
sx q[2];
rz(-0.88694015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87673346) q[1];
sx q[1];
rz(-0.85163341) q[1];
sx q[1];
rz(0.77244669) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3994282) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(-0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(0.53652525) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1034575) q[0];
sx q[0];
rz(-2.3893444) q[0];
sx q[0];
rz(-0.30970807) q[0];
rz(-pi) q[1];
rz(2.1937218) q[2];
sx q[2];
rz(-1.9588753) q[2];
sx q[2];
rz(-0.58448234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.045947) q[1];
sx q[1];
rz(-2.1150757) q[1];
sx q[1];
rz(-2.4023158) q[1];
rz(-pi) q[2];
rz(-0.94654406) q[3];
sx q[3];
rz(-2.4157899) q[3];
sx q[3];
rz(-3.1012227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(0.12891842) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0607973) q[0];
sx q[0];
rz(-1.9658488) q[0];
sx q[0];
rz(0.11443826) q[0];
rz(-pi) q[1];
rz(-2.244197) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(-2.430254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8184549) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(2.6998991) q[1];
x q[2];
rz(1.586732) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(-2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-0.62409419) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(-2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(1.7519959) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1662912) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(1.5750969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8920039) q[1];
sx q[1];
rz(-0.87000404) q[1];
sx q[1];
rz(0.65111098) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8072855) q[3];
sx q[3];
rz(-0.83179501) q[3];
sx q[3];
rz(2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(-2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(-0.84035994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6662153) q[2];
sx q[2];
rz(-1.4350841) q[2];
sx q[2];
rz(-0.00098468653) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1018131) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(1.5323557) q[1];
rz(-pi) q[2];
rz(-0.54639001) q[3];
sx q[3];
rz(-2.8791109) q[3];
sx q[3];
rz(-2.4266092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(0.51013851) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1417086) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(3.1267358) q[2];
sx q[2];
rz(-1.6179634) q[2];
sx q[2];
rz(0.85214324) q[2];
rz(0.32904939) q[3];
sx q[3];
rz(-1.1209189) q[3];
sx q[3];
rz(0.95595595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
