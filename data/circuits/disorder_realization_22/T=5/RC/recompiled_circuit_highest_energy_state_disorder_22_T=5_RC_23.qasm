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
rz(2.4616315) q[0];
sx q[0];
rz(7.1890561) q[0];
sx q[0];
rz(11.472975) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(-1.1727762) q[1];
sx q[1];
rz(-0.42458951) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49219201) q[0];
sx q[0];
rz(-1.7723506) q[0];
sx q[0];
rz(3.1026476) q[0];
rz(-pi) q[1];
x q[1];
rz(1.254804) q[2];
sx q[2];
rz(-2.9796763) q[2];
sx q[2];
rz(0.94053167) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54271419) q[1];
sx q[1];
rz(-1.492036) q[1];
sx q[1];
rz(1.3244737) q[1];
rz(-pi) q[2];
rz(3.0077259) q[3];
sx q[3];
rz(-2.2765719) q[3];
sx q[3];
rz(-5.2701252e-05) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37497416) q[2];
sx q[2];
rz(-1.582229) q[2];
sx q[2];
rz(-1.9826822) q[2];
rz(2.9186987) q[3];
sx q[3];
rz(-1.736172) q[3];
sx q[3];
rz(-2.8515653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463773) q[0];
sx q[0];
rz(-1.615849) q[0];
sx q[0];
rz(-2.0624397) q[0];
rz(-1.4214628) q[1];
sx q[1];
rz(-2.454897) q[1];
sx q[1];
rz(-0.064528331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2384528) q[0];
sx q[0];
rz(-1.6232326) q[0];
sx q[0];
rz(-1.658123) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47625292) q[2];
sx q[2];
rz(-1.5956399) q[2];
sx q[2];
rz(0.78643878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1074984) q[1];
sx q[1];
rz(-2.3281257) q[1];
sx q[1];
rz(-0.33593049) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3702303) q[3];
sx q[3];
rz(-0.75955694) q[3];
sx q[3];
rz(-2.650039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1530389) q[2];
sx q[2];
rz(-1.3554074) q[2];
sx q[2];
rz(-1.1174196) q[2];
rz(2.2828263) q[3];
sx q[3];
rz(-0.78341165) q[3];
sx q[3];
rz(-2.7045238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7607464) q[0];
sx q[0];
rz(-0.28457156) q[0];
sx q[0];
rz(0.17307702) q[0];
rz(-2.1848047) q[1];
sx q[1];
rz(-2.1134977) q[1];
sx q[1];
rz(-1.0622567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5231402) q[0];
sx q[0];
rz(-1.5644531) q[0];
sx q[0];
rz(0.84756084) q[0];
x q[1];
rz(2.144084) q[2];
sx q[2];
rz(-1.0110223) q[2];
sx q[2];
rz(0.56247518) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30327163) q[1];
sx q[1];
rz(-2.1049989) q[1];
sx q[1];
rz(2.3278585) q[1];
rz(-pi) q[2];
rz(0.62543243) q[3];
sx q[3];
rz(-0.98821001) q[3];
sx q[3];
rz(2.7435477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7303077) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(2.3968706) q[2];
rz(2.5005285) q[3];
sx q[3];
rz(-1.3499667) q[3];
sx q[3];
rz(0.49066576) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500126) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(-2.2963754) q[0];
rz(-0.40329626) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(-0.32803112) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83932861) q[0];
sx q[0];
rz(-1.7516881) q[0];
sx q[0];
rz(-0.14206391) q[0];
rz(-pi) q[1];
rz(2.9035527) q[2];
sx q[2];
rz(-1.5063707) q[2];
sx q[2];
rz(0.062459613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.58807948) q[1];
sx q[1];
rz(-2.5633286) q[1];
sx q[1];
rz(0.2167313) q[1];
rz(1.0260887) q[3];
sx q[3];
rz(-1.8063365) q[3];
sx q[3];
rz(2.1229471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24568096) q[2];
sx q[2];
rz(-0.74378219) q[2];
sx q[2];
rz(0.15920676) q[2];
rz(-3.1253452) q[3];
sx q[3];
rz(-2.1135606) q[3];
sx q[3];
rz(-0.73133674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.74723393) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(-2.0514945) q[0];
rz(-0.12995003) q[1];
sx q[1];
rz(-0.81472412) q[1];
sx q[1];
rz(-1.5257588) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9620044) q[0];
sx q[0];
rz(-1.9854913) q[0];
sx q[0];
rz(2.5432822) q[0];
rz(-pi) q[1];
rz(-3.0661929) q[2];
sx q[2];
rz(-1.7345718) q[2];
sx q[2];
rz(2.5679156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.255521) q[1];
sx q[1];
rz(-2.204245) q[1];
sx q[1];
rz(-1.6048163) q[1];
rz(0.035338621) q[3];
sx q[3];
rz(-1.7515079) q[3];
sx q[3];
rz(-2.904195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7502363) q[2];
sx q[2];
rz(-2.3199234) q[2];
sx q[2];
rz(-2.6643122) q[2];
rz(0.20398772) q[3];
sx q[3];
rz(-0.19652772) q[3];
sx q[3];
rz(-2.5175214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9452962) q[0];
sx q[0];
rz(-2.3260703) q[0];
sx q[0];
rz(0.77769172) q[0];
rz(-2.5241191) q[1];
sx q[1];
rz(-1.6505046) q[1];
sx q[1];
rz(2.2106574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27706596) q[0];
sx q[0];
rz(-2.0210938) q[0];
sx q[0];
rz(-0.45876512) q[0];
x q[1];
rz(-2.0495049) q[2];
sx q[2];
rz(-1.1943814) q[2];
sx q[2];
rz(-2.7680754) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27241421) q[1];
sx q[1];
rz(-2.8003516) q[1];
sx q[1];
rz(1.7029087) q[1];
x q[2];
rz(-0.39998885) q[3];
sx q[3];
rz(-2.1533186) q[3];
sx q[3];
rz(-1.3056684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4882539) q[2];
sx q[2];
rz(-1.808017) q[2];
sx q[2];
rz(2.874157) q[2];
rz(-2.2087162) q[3];
sx q[3];
rz(-1.8578015) q[3];
sx q[3];
rz(-2.647184) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6701508) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(0.66584051) q[0];
rz(-0.2039856) q[1];
sx q[1];
rz(-2.1603656) q[1];
sx q[1];
rz(1.1873672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68722938) q[0];
sx q[0];
rz(-0.47894127) q[0];
sx q[0];
rz(-2.9402551) q[0];
x q[1];
rz(2.6186278) q[2];
sx q[2];
rz(-0.93237662) q[2];
sx q[2];
rz(-2.0948727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1923869) q[1];
sx q[1];
rz(-1.306136) q[1];
sx q[1];
rz(-0.70313518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3755619) q[3];
sx q[3];
rz(-0.25849202) q[3];
sx q[3];
rz(-3.0743161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5953956) q[2];
sx q[2];
rz(-2.6232145) q[2];
sx q[2];
rz(-2.4714244) q[2];
rz(-0.3012805) q[3];
sx q[3];
rz(-1.2944841) q[3];
sx q[3];
rz(2.7231351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8381074) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(2.0759034) q[0];
rz(-0.65713716) q[1];
sx q[1];
rz(-2.4333351) q[1];
sx q[1];
rz(-1.7117975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098786548) q[0];
sx q[0];
rz(-1.0142066) q[0];
sx q[0];
rz(-0.3996398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7523319) q[2];
sx q[2];
rz(-1.2444656) q[2];
sx q[2];
rz(1.2943314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7969574) q[1];
sx q[1];
rz(-1.9785787) q[1];
sx q[1];
rz(-1.8090882) q[1];
rz(-pi) q[2];
rz(-2.9518806) q[3];
sx q[3];
rz(-1.4290775) q[3];
sx q[3];
rz(-0.82822463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8810205) q[2];
sx q[2];
rz(-1.842097) q[2];
sx q[2];
rz(-1.0183498) q[2];
rz(1.5492505) q[3];
sx q[3];
rz(-1.6377623) q[3];
sx q[3];
rz(0.75756592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8757979) q[0];
sx q[0];
rz(-2.6261411) q[0];
sx q[0];
rz(0.96762586) q[0];
rz(0.34010092) q[1];
sx q[1];
rz(-0.83931559) q[1];
sx q[1];
rz(2.1077154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5115625) q[0];
sx q[0];
rz(-2.1242737) q[0];
sx q[0];
rz(-0.16279499) q[0];
rz(-pi) q[1];
rz(-0.56615717) q[2];
sx q[2];
rz(-0.70933178) q[2];
sx q[2];
rz(0.69401238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8845773) q[1];
sx q[1];
rz(-2.3126763) q[1];
sx q[1];
rz(-2.1121426) q[1];
rz(-pi) q[2];
rz(0.22397174) q[3];
sx q[3];
rz(-1.4272318) q[3];
sx q[3];
rz(0.75907133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0391482) q[2];
sx q[2];
rz(-1.6951963) q[2];
sx q[2];
rz(0.038912494) q[2];
rz(-0.14032042) q[3];
sx q[3];
rz(-0.084429927) q[3];
sx q[3];
rz(1.3154359) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4843531) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(-1.8835618) q[0];
rz(-2.925442) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(0.36453882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.011019) q[0];
sx q[0];
rz(-0.51968658) q[0];
sx q[0];
rz(-2.0849865) q[0];
rz(-pi) q[1];
rz(0.80354013) q[2];
sx q[2];
rz(-1.2422971) q[2];
sx q[2];
rz(-0.9847509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6072104) q[1];
sx q[1];
rz(-0.057091968) q[1];
sx q[1];
rz(2.9109138) q[1];
rz(-pi) q[2];
rz(-3.0836948) q[3];
sx q[3];
rz(-1.5385043) q[3];
sx q[3];
rz(2.3253289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2898499) q[2];
sx q[2];
rz(-1.2052636) q[2];
sx q[2];
rz(2.5002948) q[2];
rz(1.4801721) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(0.67354584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4079473) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(-2.1570878) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(1.3732688) q[2];
sx q[2];
rz(-1.6850204) q[2];
sx q[2];
rz(-1.9584283) q[2];
rz(-0.52624191) q[3];
sx q[3];
rz(-1.9947697) q[3];
sx q[3];
rz(1.5823732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
