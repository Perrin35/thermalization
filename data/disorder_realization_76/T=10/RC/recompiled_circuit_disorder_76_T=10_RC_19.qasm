OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(-0.13248086) q[0];
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(-2.4490228) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4953295) q[0];
sx q[0];
rz(-1.6705992) q[0];
sx q[0];
rz(-1.050759) q[0];
rz(-pi) q[1];
rz(0.64451005) q[2];
sx q[2];
rz(-1.1877726) q[2];
sx q[2];
rz(0.884998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.270826) q[1];
sx q[1];
rz(-1.4987136) q[1];
sx q[1];
rz(1.6986398) q[1];
rz(-pi) q[2];
rz(1.0971783) q[3];
sx q[3];
rz(-1.3485104) q[3];
sx q[3];
rz(1.3003295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(1.6202554) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(0.56150395) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(-0.5805648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8037572) q[0];
sx q[0];
rz(-2.8135186) q[0];
sx q[0];
rz(2.8003545) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3524019) q[2];
sx q[2];
rz(-0.94422715) q[2];
sx q[2];
rz(-0.44109694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7932574) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(-2.4130915) q[1];
x q[2];
rz(-0.92506261) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(2.1893082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(0.77366775) q[0];
rz(-3.1402918) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(-0.032827854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87356991) q[0];
sx q[0];
rz(-2.7497254) q[0];
sx q[0];
rz(-0.64143945) q[0];
rz(-pi) q[1];
rz(0.29067729) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(-0.79613396) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.153861) q[1];
sx q[1];
rz(-0.78480936) q[1];
sx q[1];
rz(1.3373109) q[1];
x q[2];
rz(0.34605108) q[3];
sx q[3];
rz(-2.0882574) q[3];
sx q[3];
rz(-0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34439987) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(-2.8642505) q[2];
rz(-2.7456361) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(2.4424281) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-0.26279703) q[0];
rz(0.2335877) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(-2.3707726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.748427) q[0];
sx q[0];
rz(-1.59032) q[0];
sx q[0];
rz(0.016419134) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0371573) q[2];
sx q[2];
rz(-2.0586788) q[2];
sx q[2];
rz(1.3654396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.96344906) q[1];
sx q[1];
rz(-1.2091067) q[1];
sx q[1];
rz(-2.5835035) q[1];
rz(1.7951489) q[3];
sx q[3];
rz(-1.7566578) q[3];
sx q[3];
rz(0.88691521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16584855) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(-0.7129933) q[2];
rz(2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35904303) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(0.094093181) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-2.9715911) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3225587) q[0];
sx q[0];
rz(-1.337262) q[0];
sx q[0];
rz(2.3755431) q[0];
rz(-pi) q[1];
rz(2.0427809) q[2];
sx q[2];
rz(-1.7760279) q[2];
sx q[2];
rz(-0.99265487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1483037) q[1];
sx q[1];
rz(-1.9699886) q[1];
sx q[1];
rz(0.94435512) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65808987) q[3];
sx q[3];
rz(-0.98539017) q[3];
sx q[3];
rz(0.64774367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99469441) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(1.7374932) q[2];
rz(1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(-1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(0.25587747) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(2.4564254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4165669) q[0];
sx q[0];
rz(-1.6402813) q[0];
sx q[0];
rz(2.2216703) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7148758) q[2];
sx q[2];
rz(-2.3428168) q[2];
sx q[2];
rz(-0.069552334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8243858) q[1];
sx q[1];
rz(-2.3134391) q[1];
sx q[1];
rz(0.87888996) q[1];
rz(-1.54322) q[3];
sx q[3];
rz(-0.58610361) q[3];
sx q[3];
rz(0.54857777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7489862) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(1.7123327) q[2];
rz(-2.0424992) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-0.18923047) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012506164) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(-0.72934735) q[0];
rz(-0.29306456) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(-1.1475295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.731819) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(-0.41418196) q[0];
rz(-0.36979923) q[2];
sx q[2];
rz(-2.1714006) q[2];
sx q[2];
rz(-0.72593867) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7123375) q[1];
sx q[1];
rz(-1.0865679) q[1];
sx q[1];
rz(-1.5244563) q[1];
rz(-2.494874) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(-1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78836936) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(2.3925171) q[2];
rz(2.4979112) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(0.914004) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99326837) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(-2.3102982) q[0];
rz(-1.7656901) q[1];
sx q[1];
rz(-2.3283236) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9015394) q[0];
sx q[0];
rz(-2.0365289) q[0];
sx q[0];
rz(-0.17850152) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6767098) q[2];
sx q[2];
rz(-0.97230655) q[2];
sx q[2];
rz(-0.61818365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7290013) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(0.97357915) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3066735) q[3];
sx q[3];
rz(-2.4432126) q[3];
sx q[3];
rz(-1.8861119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9514256) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(1.3893611) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(-2.0432037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7669945) q[0];
sx q[0];
rz(-0.54531389) q[0];
sx q[0];
rz(-2.0314412) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82010014) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(2.1095914) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2213649) q[1];
sx q[1];
rz(-0.5792633) q[1];
sx q[1];
rz(-3.1406162) q[1];
rz(-pi) q[2];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.8212574) q[3];
sx q[3];
rz(0.40303883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8998469) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(-2.6220654) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(-1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7463503) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(-2.675132) q[0];
rz(-0.17164224) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(2.5126273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6105729) q[0];
sx q[0];
rz(-1.6197228) q[0];
sx q[0];
rz(-1.1958836) q[0];
x q[1];
rz(0.92832698) q[2];
sx q[2];
rz(-1.3767585) q[2];
sx q[2];
rz(2.5460668) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.5728587) q[1];
sx q[1];
rz(0.44209977) q[1];
rz(-0.38871308) q[3];
sx q[3];
rz(-2.3832088) q[3];
sx q[3];
rz(0.10314108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(2.0824599) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(-2.2476851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28329904) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(-2.4181096) q[2];
sx q[2];
rz(-1.537848) q[2];
sx q[2];
rz(1.6707735) q[2];
rz(-2.1577088) q[3];
sx q[3];
rz(-1.513653) q[3];
sx q[3];
rz(-0.37099864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
