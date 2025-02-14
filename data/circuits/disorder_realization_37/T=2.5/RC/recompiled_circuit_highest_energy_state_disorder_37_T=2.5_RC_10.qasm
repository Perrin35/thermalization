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
rz(1.2020741) q[0];
sx q[0];
rz(-2.095686) q[0];
sx q[0];
rz(2.3910971) q[0];
rz(0.22076386) q[1];
sx q[1];
rz(-2.0808487) q[1];
sx q[1];
rz(1.9579252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7322217) q[0];
sx q[0];
rz(-2.3187162) q[0];
sx q[0];
rz(-2.5636682) q[0];
rz(-0.92249845) q[2];
sx q[2];
rz(-2.9335409) q[2];
sx q[2];
rz(1.7788943) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4422801) q[1];
sx q[1];
rz(-2.2506425) q[1];
sx q[1];
rz(1.5857497) q[1];
x q[2];
rz(-0.99917717) q[3];
sx q[3];
rz(-2.7795254) q[3];
sx q[3];
rz(1.3130275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3262647) q[2];
sx q[2];
rz(-1.1417737) q[2];
sx q[2];
rz(-3.0495138) q[2];
rz(2.0103256) q[3];
sx q[3];
rz(-2.0720033) q[3];
sx q[3];
rz(2.2856975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6446514) q[0];
sx q[0];
rz(-2.6583789) q[0];
sx q[0];
rz(2.605873) q[0];
rz(3.1247395) q[1];
sx q[1];
rz(-0.87937513) q[1];
sx q[1];
rz(-3.1244315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494228) q[0];
sx q[0];
rz(-2.0606406) q[0];
sx q[0];
rz(-0.36618284) q[0];
rz(-pi) q[1];
rz(2.6748528) q[2];
sx q[2];
rz(-1.6306976) q[2];
sx q[2];
rz(1.8869099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8507599) q[1];
sx q[1];
rz(-0.28181812) q[1];
sx q[1];
rz(2.5435104) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4799177) q[3];
sx q[3];
rz(-1.7044997) q[3];
sx q[3];
rz(0.17482486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0422334) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(-0.24822203) q[2];
rz(-0.92784268) q[3];
sx q[3];
rz(-0.78301269) q[3];
sx q[3];
rz(0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8645653) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(2.9789341) q[0];
rz(-0.76205572) q[1];
sx q[1];
rz(-0.22480741) q[1];
sx q[1];
rz(2.3077097) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91438566) q[0];
sx q[0];
rz(-1.4482131) q[0];
sx q[0];
rz(-0.57182619) q[0];
x q[1];
rz(-0.35807407) q[2];
sx q[2];
rz(-2.798657) q[2];
sx q[2];
rz(-1.5141405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2411338) q[1];
sx q[1];
rz(-0.83158439) q[1];
sx q[1];
rz(2.1421188) q[1];
x q[2];
rz(1.6271851) q[3];
sx q[3];
rz(-2.0109551) q[3];
sx q[3];
rz(1.7007916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0791846) q[2];
sx q[2];
rz(-2.0024313) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(1.0583813) q[3];
sx q[3];
rz(-2.1975785) q[3];
sx q[3];
rz(-0.03037608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3707378) q[0];
sx q[0];
rz(-2.3069032) q[0];
sx q[0];
rz(-1.7133065) q[0];
rz(-0.76438534) q[1];
sx q[1];
rz(-1.9995707) q[1];
sx q[1];
rz(1.8199325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40250889) q[0];
sx q[0];
rz(-0.54259491) q[0];
sx q[0];
rz(-0.71904166) q[0];
rz(-pi) q[1];
rz(0.051647746) q[2];
sx q[2];
rz(-1.2394636) q[2];
sx q[2];
rz(-0.29690642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5668754) q[1];
sx q[1];
rz(-1.8758901) q[1];
sx q[1];
rz(-1.8691236) q[1];
x q[2];
rz(-1.0760078) q[3];
sx q[3];
rz(-2.3632999) q[3];
sx q[3];
rz(2.8552804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4798212) q[2];
sx q[2];
rz(-1.4134553) q[2];
sx q[2];
rz(2.856355) q[2];
rz(1.2626922) q[3];
sx q[3];
rz(-1.9243536) q[3];
sx q[3];
rz(-1.7181905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406463) q[0];
sx q[0];
rz(-0.74746376) q[0];
sx q[0];
rz(0.88229156) q[0];
rz(2.4709985) q[1];
sx q[1];
rz(-1.0958442) q[1];
sx q[1];
rz(0.98463279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66259164) q[0];
sx q[0];
rz(-3.0714095) q[0];
sx q[0];
rz(-0.46371882) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.042858275) q[2];
sx q[2];
rz(-1.0295261) q[2];
sx q[2];
rz(-2.0861911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24422503) q[1];
sx q[1];
rz(-2.3363718) q[1];
sx q[1];
rz(-0.77181384) q[1];
x q[2];
rz(2.601131) q[3];
sx q[3];
rz(-1.1532047) q[3];
sx q[3];
rz(0.77089522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.275445) q[2];
sx q[2];
rz(-0.53469849) q[2];
sx q[2];
rz(0.61694413) q[2];
rz(2.3894737) q[3];
sx q[3];
rz(-1.3277466) q[3];
sx q[3];
rz(-0.45929685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1521456) q[0];
sx q[0];
rz(-0.716827) q[0];
sx q[0];
rz(2.3710807) q[0];
rz(2.41467) q[1];
sx q[1];
rz(-2.9194071) q[1];
sx q[1];
rz(-2.4368584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7677697) q[0];
sx q[0];
rz(-0.85568586) q[0];
sx q[0];
rz(0.82775791) q[0];
x q[1];
rz(-2.0915085) q[2];
sx q[2];
rz(-1.2479223) q[2];
sx q[2];
rz(-0.31393696) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2134691) q[1];
sx q[1];
rz(-1.537461) q[1];
sx q[1];
rz(-2.7916905) q[1];
rz(-0.25382692) q[3];
sx q[3];
rz(-1.0146662) q[3];
sx q[3];
rz(1.1185916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7695693) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(0.43323576) q[2];
rz(-0.20089928) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(2.6065629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3055426) q[0];
sx q[0];
rz(-1.9925646) q[0];
sx q[0];
rz(-2.7591163) q[0];
rz(1.0410694) q[1];
sx q[1];
rz(-1.6761227) q[1];
sx q[1];
rz(-2.6952851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8167717) q[0];
sx q[0];
rz(-1.8292301) q[0];
sx q[0];
rz(1.6245317) q[0];
rz(2.7591797) q[2];
sx q[2];
rz(-1.3578041) q[2];
sx q[2];
rz(-1.2455924) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.81627264) q[1];
sx q[1];
rz(-0.98794663) q[1];
sx q[1];
rz(-2.7144949) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84844037) q[3];
sx q[3];
rz(-2.1506243) q[3];
sx q[3];
rz(2.8715796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81991601) q[2];
sx q[2];
rz(-1.1175464) q[2];
sx q[2];
rz(-2.4134911) q[2];
rz(0.32650945) q[3];
sx q[3];
rz(-1.510334) q[3];
sx q[3];
rz(-2.2842893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9051055) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(1.2976728) q[0];
rz(-2.1406651) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(2.7330858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4057874) q[0];
sx q[0];
rz(-1.5716432) q[0];
sx q[0];
rz(1.5849831) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3199072) q[2];
sx q[2];
rz(-0.61640451) q[2];
sx q[2];
rz(2.6954755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7311915) q[1];
sx q[1];
rz(-2.9653314) q[1];
sx q[1];
rz(0.86430637) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85902416) q[3];
sx q[3];
rz(-0.54104303) q[3];
sx q[3];
rz(-0.67160254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0268176) q[2];
sx q[2];
rz(-1.6208181) q[2];
sx q[2];
rz(2.3180023) q[2];
rz(0.94281998) q[3];
sx q[3];
rz(-1.4225682) q[3];
sx q[3];
rz(-1.6387117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5470062) q[0];
sx q[0];
rz(-0.44438812) q[0];
sx q[0];
rz(-1.8894926) q[0];
rz(-1.7777255) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(-1.3346671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33556834) q[0];
sx q[0];
rz(-2.8360406) q[0];
sx q[0];
rz(-3.0204885) q[0];
rz(-pi) q[1];
rz(-1.7450856) q[2];
sx q[2];
rz(-2.2054923) q[2];
sx q[2];
rz(1.2668934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6393362) q[1];
sx q[1];
rz(-0.4570804) q[1];
sx q[1];
rz(0.16641081) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0638716) q[3];
sx q[3];
rz(-2.6209313) q[3];
sx q[3];
rz(-2.9021104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1198662) q[2];
sx q[2];
rz(-2.848337) q[2];
sx q[2];
rz(-2.8450656) q[2];
rz(-1.4793388) q[3];
sx q[3];
rz(-1.5115503) q[3];
sx q[3];
rz(-0.1571981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.6327561) q[0];
sx q[0];
rz(-0.65082508) q[0];
sx q[0];
rz(-2.346709) q[0];
rz(-0.55566135) q[1];
sx q[1];
rz(-1.8652922) q[1];
sx q[1];
rz(1.4478252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5705161) q[0];
sx q[0];
rz(-0.85502842) q[0];
sx q[0];
rz(2.0809196) q[0];
rz(-pi) q[1];
rz(-2.6437691) q[2];
sx q[2];
rz(-1.6256672) q[2];
sx q[2];
rz(-0.16615088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8965079) q[1];
sx q[1];
rz(-1.0963529) q[1];
sx q[1];
rz(-0.22405476) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9742188) q[3];
sx q[3];
rz(-1.1165301) q[3];
sx q[3];
rz(-0.59859401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8453703) q[2];
sx q[2];
rz(-1.5105379) q[2];
sx q[2];
rz(1.4307865) q[2];
rz(0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(-0.16451612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0059148) q[0];
sx q[0];
rz(-2.4527241) q[0];
sx q[0];
rz(1.2149568) q[0];
rz(1.6750492) q[1];
sx q[1];
rz(-2.2129682) q[1];
sx q[1];
rz(0.60563544) q[1];
rz(0.56500021) q[2];
sx q[2];
rz(-1.0641268) q[2];
sx q[2];
rz(0.63882154) q[2];
rz(2.130065) q[3];
sx q[3];
rz(-0.20023919) q[3];
sx q[3];
rz(-1.7773624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
