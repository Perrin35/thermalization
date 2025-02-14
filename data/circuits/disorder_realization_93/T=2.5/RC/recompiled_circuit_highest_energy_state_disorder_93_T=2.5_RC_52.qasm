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
rz(-2.8060198) q[0];
sx q[0];
rz(-1.6874474) q[0];
sx q[0];
rz(2.1079221) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(2.1318336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5752707) q[0];
sx q[0];
rz(-0.22804865) q[0];
sx q[0];
rz(2.1581677) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0037438914) q[2];
sx q[2];
rz(-1.4883248) q[2];
sx q[2];
rz(-1.0325026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7518471) q[1];
sx q[1];
rz(-1.6865786) q[1];
sx q[1];
rz(2.9751361) q[1];
x q[2];
rz(3.1107424) q[3];
sx q[3];
rz(-2.9992963) q[3];
sx q[3];
rz(-2.9205585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7734163) q[2];
sx q[2];
rz(-2.2645617) q[2];
sx q[2];
rz(-0.82484335) q[2];
rz(-2.8090737) q[3];
sx q[3];
rz(-2.2907292) q[3];
sx q[3];
rz(-0.48377812) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9124209) q[0];
sx q[0];
rz(-3.0569172) q[0];
sx q[0];
rz(-2.605865) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.3923233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24015445) q[0];
sx q[0];
rz(-1.9124228) q[0];
sx q[0];
rz(1.0761976) q[0];
rz(2.8004592) q[2];
sx q[2];
rz(-1.5617512) q[2];
sx q[2];
rz(1.1760528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.630188) q[1];
sx q[1];
rz(-1.5492445) q[1];
sx q[1];
rz(-2.3261916) q[1];
x q[2];
rz(1.4658607) q[3];
sx q[3];
rz(-1.6989048) q[3];
sx q[3];
rz(2.6554633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9048189) q[2];
sx q[2];
rz(-2.3084013) q[2];
sx q[2];
rz(2.0429677) q[2];
rz(-0.83313471) q[3];
sx q[3];
rz(-1.5435217) q[3];
sx q[3];
rz(-2.0875077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5611834) q[0];
sx q[0];
rz(-1.1495178) q[0];
sx q[0];
rz(2.4533601) q[0];
rz(-1.8904103) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(1.9293264) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040846856) q[0];
sx q[0];
rz(-2.0424574) q[0];
sx q[0];
rz(0.77420401) q[0];
rz(2.278944) q[2];
sx q[2];
rz(-0.35637636) q[2];
sx q[2];
rz(-1.4087806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.41238585) q[1];
sx q[1];
rz(-2.4826035) q[1];
sx q[1];
rz(-0.0091730793) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.631725) q[3];
sx q[3];
rz(-0.079515545) q[3];
sx q[3];
rz(-0.36859504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4091829) q[2];
sx q[2];
rz(-2.1467291) q[2];
sx q[2];
rz(-0.44962064) q[2];
rz(0.85363394) q[3];
sx q[3];
rz(-2.4946404) q[3];
sx q[3];
rz(2.4322815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39321008) q[0];
sx q[0];
rz(-2.125183) q[0];
sx q[0];
rz(2.7384695) q[0];
rz(2.368811) q[1];
sx q[1];
rz(-1.9085596) q[1];
sx q[1];
rz(-0.82537878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51096237) q[0];
sx q[0];
rz(-1.2494506) q[0];
sx q[0];
rz(-2.3031631) q[0];
x q[1];
rz(-1.5947444) q[2];
sx q[2];
rz(-1.5439111) q[2];
sx q[2];
rz(-1.8906817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9627208) q[1];
sx q[1];
rz(-2.555518) q[1];
sx q[1];
rz(1.4337711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55565662) q[3];
sx q[3];
rz(-2.1534277) q[3];
sx q[3];
rz(-0.61911303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0394502) q[2];
sx q[2];
rz(-1.8723698) q[2];
sx q[2];
rz(-0.5086745) q[2];
rz(1.2447119) q[3];
sx q[3];
rz(-1.0735984) q[3];
sx q[3];
rz(-1.3332453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73388571) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(2.7852614) q[0];
rz(1.412926) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(-1.2948571) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.750941) q[0];
sx q[0];
rz(-1.7897507) q[0];
sx q[0];
rz(2.697437) q[0];
rz(-pi) q[1];
rz(-2.9239976) q[2];
sx q[2];
rz(-0.67309531) q[2];
sx q[2];
rz(-2.4270647) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8581681) q[1];
sx q[1];
rz(-1.8008261) q[1];
sx q[1];
rz(-2.7344804) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.242022) q[3];
sx q[3];
rz(-2.482467) q[3];
sx q[3];
rz(-0.61263436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20092189) q[2];
sx q[2];
rz(-1.7565497) q[2];
sx q[2];
rz(-2.5015639) q[2];
rz(-0.43404964) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(1.1205589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5459179) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(2.3906999) q[0];
rz(0.76260507) q[1];
sx q[1];
rz(-1.7060988) q[1];
sx q[1];
rz(0.5336175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92350875) q[0];
sx q[0];
rz(-1.8663915) q[0];
sx q[0];
rz(-1.8718998) q[0];
rz(-0.4328753) q[2];
sx q[2];
rz(-0.85202571) q[2];
sx q[2];
rz(-1.5804039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4009559) q[1];
sx q[1];
rz(-1.6888685) q[1];
sx q[1];
rz(-0.90098874) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98762767) q[3];
sx q[3];
rz(-1.0915874) q[3];
sx q[3];
rz(-0.85239172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0642455) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(2.2557491) q[2];
rz(-1.8387851) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(2.6634789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122022) q[0];
sx q[0];
rz(-2.2181856) q[0];
sx q[0];
rz(0.59342629) q[0];
rz(1.5133096) q[1];
sx q[1];
rz(-1.9624886) q[1];
sx q[1];
rz(2.5679307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5705171) q[0];
sx q[0];
rz(-0.62643753) q[0];
sx q[0];
rz(-2.4721739) q[0];
rz(-2.7717765) q[2];
sx q[2];
rz(-2.5313583) q[2];
sx q[2];
rz(3.1265697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5630514) q[1];
sx q[1];
rz(-2.6669901) q[1];
sx q[1];
rz(2.5792349) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0629458) q[3];
sx q[3];
rz(-1.6212731) q[3];
sx q[3];
rz(-0.61928643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4297428) q[2];
sx q[2];
rz(-1.2028799) q[2];
sx q[2];
rz(1.708606) q[2];
rz(-1.9728707) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(-1.5168813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8716005) q[0];
sx q[0];
rz(-2.2347436) q[0];
sx q[0];
rz(0.20371833) q[0];
rz(-1.3153971) q[1];
sx q[1];
rz(-1.306465) q[1];
sx q[1];
rz(1.809583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3888793) q[0];
sx q[0];
rz(-0.95001555) q[0];
sx q[0];
rz(-0.78691532) q[0];
rz(-pi) q[1];
rz(-2.5723348) q[2];
sx q[2];
rz(-0.98460403) q[2];
sx q[2];
rz(-1.7050515) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.2329008) q[1];
sx q[1];
rz(-0.76592839) q[1];
sx q[1];
rz(2.5604064) q[1];
rz(-pi) q[2];
rz(-1.4077987) q[3];
sx q[3];
rz(-1.6987112) q[3];
sx q[3];
rz(0.35234141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3709162) q[2];
sx q[2];
rz(-0.081826536) q[2];
sx q[2];
rz(1.3814242) q[2];
rz(0.55111876) q[3];
sx q[3];
rz(-1.328732) q[3];
sx q[3];
rz(-0.74530017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4296752) q[0];
sx q[0];
rz(-1.4497919) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(1.7970386) q[1];
sx q[1];
rz(-2.5090736) q[1];
sx q[1];
rz(-1.22619) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525402) q[0];
sx q[0];
rz(-1.8551833) q[0];
sx q[0];
rz(-0.14099462) q[0];
rz(-pi) q[1];
rz(2.2974469) q[2];
sx q[2];
rz(-2.591989) q[2];
sx q[2];
rz(-1.5908102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3873432) q[1];
sx q[1];
rz(-0.74973901) q[1];
sx q[1];
rz(0.096258817) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1192083) q[3];
sx q[3];
rz(-0.58956857) q[3];
sx q[3];
rz(-1.1353253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58174497) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(-1.6017412) q[2];
rz(1.3567443) q[3];
sx q[3];
rz(-2.609085) q[3];
sx q[3];
rz(-1.9073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8643879) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(-0.10449115) q[0];
rz(1.744572) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(-1.5207965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6817271) q[0];
sx q[0];
rz(-2.2978373) q[0];
sx q[0];
rz(-1.2894425) q[0];
rz(-1.4112597) q[2];
sx q[2];
rz(-1.9070093) q[2];
sx q[2];
rz(-1.6415063) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7193774) q[1];
sx q[1];
rz(-1.6746192) q[1];
sx q[1];
rz(0.28909282) q[1];
x q[2];
rz(-3.0836068) q[3];
sx q[3];
rz(-2.1515905) q[3];
sx q[3];
rz(-0.58422663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4349159) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.8527156) q[2];
rz(-2.1364818) q[3];
sx q[3];
rz(-2.0821327) q[3];
sx q[3];
rz(0.10678261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2224251) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(0.87986058) q[1];
sx q[1];
rz(-2.6251371) q[1];
sx q[1];
rz(0.62216204) q[1];
rz(2.6420319) q[2];
sx q[2];
rz(-1.6274263) q[2];
sx q[2];
rz(-2.2123663) q[2];
rz(0.60012695) q[3];
sx q[3];
rz(-1.5059581) q[3];
sx q[3];
rz(-3.0653421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
