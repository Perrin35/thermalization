OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9094698) q[0];
sx q[0];
rz(-2.1535518) q[0];
sx q[0];
rz(2.7890132) q[0];
rz(2.4352788) q[1];
sx q[1];
rz(-2.3568454) q[1];
sx q[1];
rz(3.1142601) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3685212) q[0];
sx q[0];
rz(-1.7681244) q[0];
sx q[0];
rz(-2.7993566) q[0];
rz(-pi) q[1];
rz(0.66157056) q[2];
sx q[2];
rz(-0.21851608) q[2];
sx q[2];
rz(1.0929581) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3308709) q[1];
sx q[1];
rz(-1.6020944) q[1];
sx q[1];
rz(0.16783952) q[1];
rz(-pi) q[2];
rz(-0.23689278) q[3];
sx q[3];
rz(-1.4290571) q[3];
sx q[3];
rz(-0.6976034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0722384) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(1.5864774) q[2];
rz(-2.41411) q[3];
sx q[3];
rz(-2.1860217) q[3];
sx q[3];
rz(0.47310841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.469406) q[0];
sx q[0];
rz(-1.7962026) q[0];
sx q[0];
rz(2.192002) q[0];
rz(-2.8166215) q[1];
sx q[1];
rz(-1.800671) q[1];
sx q[1];
rz(0.51345888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0810222) q[0];
sx q[0];
rz(-0.96310421) q[0];
sx q[0];
rz(-0.49852235) q[0];
x q[1];
rz(0.32294257) q[2];
sx q[2];
rz(-2.2720798) q[2];
sx q[2];
rz(2.8632134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9565485) q[1];
sx q[1];
rz(-2.0868851) q[1];
sx q[1];
rz(-0.98770492) q[1];
rz(-pi) q[2];
rz(-0.98987119) q[3];
sx q[3];
rz(-1.9912915) q[3];
sx q[3];
rz(0.48244259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.010633858) q[2];
sx q[2];
rz(-1.3291239) q[2];
sx q[2];
rz(-1.0884253) q[2];
rz(-0.66793495) q[3];
sx q[3];
rz(-2.9974944) q[3];
sx q[3];
rz(-1.4151423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3208991) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(1.8841085) q[0];
rz(-3.056774) q[1];
sx q[1];
rz(-3.0501922) q[1];
sx q[1];
rz(2.4841097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456928) q[0];
sx q[0];
rz(-2.6728711) q[0];
sx q[0];
rz(2.8492941) q[0];
x q[1];
rz(-0.74983238) q[2];
sx q[2];
rz(-0.74374357) q[2];
sx q[2];
rz(2.2735571) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0850683) q[1];
sx q[1];
rz(-0.51125647) q[1];
sx q[1];
rz(0.35447094) q[1];
x q[2];
rz(-1.9593997) q[3];
sx q[3];
rz(-0.69903421) q[3];
sx q[3];
rz(-0.70632589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76454863) q[2];
sx q[2];
rz(-0.42082861) q[2];
sx q[2];
rz(-1.2505038) q[2];
rz(1.4416384) q[3];
sx q[3];
rz(-1.7504033) q[3];
sx q[3];
rz(-0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4568951) q[0];
sx q[0];
rz(-0.95624113) q[0];
sx q[0];
rz(0.60234219) q[0];
rz(-2.1777228) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(-0.22183713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93942964) q[0];
sx q[0];
rz(-1.5785281) q[0];
sx q[0];
rz(-0.011223249) q[0];
rz(2.8922563) q[2];
sx q[2];
rz(-2.0866924) q[2];
sx q[2];
rz(-3.0003386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.027699359) q[1];
sx q[1];
rz(-2.4287819) q[1];
sx q[1];
rz(2.5068587) q[1];
rz(-pi) q[2];
rz(0.45000123) q[3];
sx q[3];
rz(-1.5061989) q[3];
sx q[3];
rz(-2.3774176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7092789) q[2];
sx q[2];
rz(-1.9807434) q[2];
sx q[2];
rz(-3.0529037) q[2];
rz(2.6426219) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(0.078908198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0093507) q[0];
sx q[0];
rz(-0.049533822) q[0];
sx q[0];
rz(2.2392739) q[0];
rz(-2.8514476) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(1.1660928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6571044) q[0];
sx q[0];
rz(-2.1323626) q[0];
sx q[0];
rz(-1.7620371) q[0];
rz(1.3275215) q[2];
sx q[2];
rz(-1.8406879) q[2];
sx q[2];
rz(0.44651595) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0048827) q[1];
sx q[1];
rz(-0.45082475) q[1];
sx q[1];
rz(2.2656127) q[1];
rz(-pi) q[2];
rz(-0.10106713) q[3];
sx q[3];
rz(-0.25288452) q[3];
sx q[3];
rz(2.853228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12639283) q[2];
sx q[2];
rz(-0.79605278) q[2];
sx q[2];
rz(1.2882721) q[2];
rz(-2.1081693) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(1.6667295) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.425659) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(-2.9820251) q[0];
rz(-2.5345934) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(2.0807696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55218766) q[0];
sx q[0];
rz(-2.1656143) q[0];
sx q[0];
rz(-2.908559) q[0];
x q[1];
rz(-1.3449613) q[2];
sx q[2];
rz(-1.3094007) q[2];
sx q[2];
rz(0.022646101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6276803) q[1];
sx q[1];
rz(-0.6100756) q[1];
sx q[1];
rz(1.808829) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5626866) q[3];
sx q[3];
rz(-1.488655) q[3];
sx q[3];
rz(0.77041799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0173505) q[2];
sx q[2];
rz(-0.81581798) q[2];
sx q[2];
rz(-2.7537277) q[2];
rz(1.3346437) q[3];
sx q[3];
rz(-0.67238656) q[3];
sx q[3];
rz(-2.1147125) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2514728) q[0];
sx q[0];
rz(-0.9820745) q[0];
sx q[0];
rz(-2.2482596) q[0];
rz(-0.54126254) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(1.062324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4058901) q[0];
sx q[0];
rz(-2.2828809) q[0];
sx q[0];
rz(-1.2232365) q[0];
x q[1];
rz(-0.063060305) q[2];
sx q[2];
rz(-0.82411924) q[2];
sx q[2];
rz(-2.7450739) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.17328611) q[1];
sx q[1];
rz(-1.5472551) q[1];
sx q[1];
rz(2.5858155) q[1];
rz(-pi) q[2];
rz(-2.9422308) q[3];
sx q[3];
rz(-0.81962517) q[3];
sx q[3];
rz(-0.79341753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5541151) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(-1.1678196) q[2];
rz(2.3217412) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(2.4988373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.576936) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(-1.8336953) q[0];
rz(1.4564184) q[1];
sx q[1];
rz(-1.5142199) q[1];
sx q[1];
rz(-3.0110722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042386656) q[0];
sx q[0];
rz(-1.4975647) q[0];
sx q[0];
rz(-1.8797148) q[0];
rz(-2.0575664) q[2];
sx q[2];
rz(-0.83510885) q[2];
sx q[2];
rz(1.6693008) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91192818) q[1];
sx q[1];
rz(-2.4897631) q[1];
sx q[1];
rz(-0.064464557) q[1];
rz(-0.7217312) q[3];
sx q[3];
rz(-2.2829307) q[3];
sx q[3];
rz(-0.16068383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9556344) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(1.7308509) q[2];
rz(-3.020982) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(-1.6284778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738344) q[0];
sx q[0];
rz(-2.5262316) q[0];
sx q[0];
rz(-0.14061418) q[0];
rz(-0.73453844) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(0.75070423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44245369) q[0];
sx q[0];
rz(-1.8635731) q[0];
sx q[0];
rz(-1.8163766) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9042013) q[2];
sx q[2];
rz(-1.2331881) q[2];
sx q[2];
rz(-0.92729502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0486601) q[1];
sx q[1];
rz(-1.5092998) q[1];
sx q[1];
rz(1.7503121) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4846913) q[3];
sx q[3];
rz(-0.70843452) q[3];
sx q[3];
rz(1.751386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1371896) q[2];
sx q[2];
rz(-0.69978324) q[2];
sx q[2];
rz(-0.75263754) q[2];
rz(2.802012) q[3];
sx q[3];
rz(-1.3656253) q[3];
sx q[3];
rz(-0.021421758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.705377) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(2.5035653) q[0];
rz(0.72884196) q[1];
sx q[1];
rz(-1.3325656) q[1];
sx q[1];
rz(-0.80150882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7201231) q[0];
sx q[0];
rz(-2.3178708) q[0];
sx q[0];
rz(-0.27702866) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24157816) q[2];
sx q[2];
rz(-2.2890586) q[2];
sx q[2];
rz(-0.11213839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90097016) q[1];
sx q[1];
rz(-2.3055446) q[1];
sx q[1];
rz(-3.088984) q[1];
rz(-0.16149811) q[3];
sx q[3];
rz(-0.93220369) q[3];
sx q[3];
rz(1.0486918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30404299) q[2];
sx q[2];
rz(-1.3648698) q[2];
sx q[2];
rz(-1.1207646) q[2];
rz(-0.70901999) q[3];
sx q[3];
rz(-2.6778383) q[3];
sx q[3];
rz(-1.2254865) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40502248) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(-0.90514056) q[1];
sx q[1];
rz(-0.97733472) q[1];
sx q[1];
rz(2.6505145) q[1];
rz(-2.2974986) q[2];
sx q[2];
rz(-1.8432968) q[2];
sx q[2];
rz(-2.7805614) q[2];
rz(0.2855339) q[3];
sx q[3];
rz(-1.0665421) q[3];
sx q[3];
rz(1.4530696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
