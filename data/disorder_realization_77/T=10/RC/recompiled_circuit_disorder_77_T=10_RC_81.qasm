OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(-2.9736829) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9553298) q[0];
sx q[0];
rz(-1.7832527) q[0];
sx q[0];
rz(-1.0832548) q[0];
rz(-pi) q[1];
x q[1];
rz(2.216823) q[2];
sx q[2];
rz(-1.7416818) q[2];
sx q[2];
rz(-2.8303763) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15966283) q[1];
sx q[1];
rz(-2.5291981) q[1];
sx q[1];
rz(2.0484522) q[1];
rz(0.25986259) q[3];
sx q[3];
rz(-1.5279603) q[3];
sx q[3];
rz(0.44997893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5713455) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(2.0061357) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0552911) q[2];
sx q[2];
rz(-1.0420348) q[2];
sx q[2];
rz(0.27258401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0484867) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(0.47936819) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6395256) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-3.085014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6037613) q[0];
sx q[0];
rz(-1.0147525) q[0];
sx q[0];
rz(0.53107302) q[0];
rz(-pi) q[1];
rz(2.8851606) q[2];
sx q[2];
rz(-1.7000546) q[2];
sx q[2];
rz(-1.3916707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1287071) q[1];
sx q[1];
rz(-1.7696847) q[1];
sx q[1];
rz(0.30602869) q[1];
x q[2];
rz(1.9583086) q[3];
sx q[3];
rz(-2.1790677) q[3];
sx q[3];
rz(-2.0351978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0915506) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(1.1372304) q[0];
x q[1];
rz(-1.5092588) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(0.57519826) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5220118) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(2.3492858) q[1];
rz(-pi) q[2];
rz(2.6632518) q[3];
sx q[3];
rz(-2.0739177) q[3];
sx q[3];
rz(0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(-3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5807242) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(0.046037721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40446754) q[2];
sx q[2];
rz(-0.81768113) q[2];
sx q[2];
rz(-0.58194619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.797232) q[1];
sx q[1];
rz(-1.6318983) q[1];
sx q[1];
rz(3.1394715) q[1];
rz(-2.6258351) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(-0.26667903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8828316) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(2.8172857) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(3.086673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391101) q[0];
sx q[0];
rz(-1.6533028) q[0];
sx q[0];
rz(1.8399747) q[0];
rz(-pi) q[1];
rz(-1.5149649) q[2];
sx q[2];
rz(-0.32368127) q[2];
sx q[2];
rz(-1.8813546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7683148) q[1];
sx q[1];
rz(-2.1149966) q[1];
sx q[1];
rz(-2.0139704) q[1];
rz(-pi) q[2];
rz(2.218607) q[3];
sx q[3];
rz(-0.69159782) q[3];
sx q[3];
rz(1.637961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383354) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(-2.2708562) q[0];
rz(-pi) q[1];
rz(-2.8001627) q[2];
sx q[2];
rz(-1.9651946) q[2];
sx q[2];
rz(-2.9801126) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0525166) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(1.7757225) q[1];
x q[2];
rz(0.63776871) q[3];
sx q[3];
rz(-0.83003269) q[3];
sx q[3];
rz(-0.72731599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-1.0197619) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(-0.30050373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(2.6595594) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1883165) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(0.88027871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44725542) q[1];
sx q[1];
rz(-1.9112504) q[1];
sx q[1];
rz(-1.7561595) q[1];
rz(2.3281906) q[3];
sx q[3];
rz(-1.8040931) q[3];
sx q[3];
rz(0.63026159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(2.382747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6700867) q[0];
sx q[0];
rz(-2.717088) q[0];
sx q[0];
rz(1.5295117) q[0];
rz(-pi) q[1];
rz(0.27300948) q[2];
sx q[2];
rz(-0.62676478) q[2];
sx q[2];
rz(3.0221456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2821741) q[1];
sx q[1];
rz(-0.71600435) q[1];
sx q[1];
rz(0.23846682) q[1];
rz(-pi) q[2];
rz(-2.2364053) q[3];
sx q[3];
rz(-1.5822516) q[3];
sx q[3];
rz(0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35995099) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2595554) q[0];
sx q[0];
rz(-1.2872739) q[0];
sx q[0];
rz(-0.77772227) q[0];
x q[1];
rz(1.2953193) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(-1.3678577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.016644195) q[1];
sx q[1];
rz(-0.84608191) q[1];
sx q[1];
rz(1.2490586) q[1];
x q[2];
rz(0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-2.2617214) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.7667608) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(-2.0416904) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];