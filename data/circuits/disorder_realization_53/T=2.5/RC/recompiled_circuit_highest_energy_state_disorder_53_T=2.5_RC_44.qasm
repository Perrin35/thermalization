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
rz(1.9617957) q[0];
sx q[0];
rz(-2.0834162) q[0];
sx q[0];
rz(0.41184586) q[0];
rz(0.62077275) q[1];
sx q[1];
rz(-0.70561916) q[1];
sx q[1];
rz(-0.128428) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0000879) q[0];
sx q[0];
rz(-2.1997118) q[0];
sx q[0];
rz(-1.5733243) q[0];
rz(-2.7343636) q[2];
sx q[2];
rz(-1.1127278) q[2];
sx q[2];
rz(1.973733) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8516396) q[1];
sx q[1];
rz(-1.3760097) q[1];
sx q[1];
rz(2.508393) q[1];
rz(-pi) q[2];
rz(-1.1731568) q[3];
sx q[3];
rz(-1.4825882) q[3];
sx q[3];
rz(2.8928066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2423557) q[2];
sx q[2];
rz(-2.5981116) q[2];
sx q[2];
rz(-2.4282954) q[2];
rz(-1.9244309) q[3];
sx q[3];
rz(-1.7141914) q[3];
sx q[3];
rz(-3.0349019) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941968) q[0];
sx q[0];
rz(-1.0260181) q[0];
sx q[0];
rz(-0.062653616) q[0];
rz(0.90473908) q[1];
sx q[1];
rz(-0.86974564) q[1];
sx q[1];
rz(-3.0135801) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5103914) q[0];
sx q[0];
rz(-1.2423853) q[0];
sx q[0];
rz(-1.9241707) q[0];
rz(-pi) q[1];
rz(0.92526154) q[2];
sx q[2];
rz(-1.0224258) q[2];
sx q[2];
rz(2.546371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3941895) q[1];
sx q[1];
rz(-0.48626562) q[1];
sx q[1];
rz(-3.1267605) q[1];
rz(-pi) q[2];
rz(-2.9826134) q[3];
sx q[3];
rz(-0.74408722) q[3];
sx q[3];
rz(-1.1733623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5855828) q[2];
sx q[2];
rz(-1.2051008) q[2];
sx q[2];
rz(-2.0466059) q[2];
rz(-2.5321142) q[3];
sx q[3];
rz(-2.1938117) q[3];
sx q[3];
rz(-1.8402717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0037917) q[0];
sx q[0];
rz(-2.4584558) q[0];
sx q[0];
rz(-2.0004499) q[0];
rz(-1.2992651) q[1];
sx q[1];
rz(-1.3971993) q[1];
sx q[1];
rz(-0.13967791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3052854) q[0];
sx q[0];
rz(-1.4023614) q[0];
sx q[0];
rz(1.3091716) q[0];
rz(2.9150117) q[2];
sx q[2];
rz(-1.7936631) q[2];
sx q[2];
rz(1.8755166) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91979474) q[1];
sx q[1];
rz(-1.9842902) q[1];
sx q[1];
rz(1.3516803) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4927718) q[3];
sx q[3];
rz(-2.088097) q[3];
sx q[3];
rz(-2.9981138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7695158) q[2];
sx q[2];
rz(-0.76411) q[2];
sx q[2];
rz(-2.8238943) q[2];
rz(2.5959173) q[3];
sx q[3];
rz(-0.92394865) q[3];
sx q[3];
rz(-1.8359449) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9724156) q[0];
sx q[0];
rz(-1.4117874) q[0];
sx q[0];
rz(-0.2485982) q[0];
rz(-1.2480805) q[1];
sx q[1];
rz(-1.0634407) q[1];
sx q[1];
rz(2.4724821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1951576) q[0];
sx q[0];
rz(-1.4753597) q[0];
sx q[0];
rz(-0.49458953) q[0];
x q[1];
rz(1.8774154) q[2];
sx q[2];
rz(-2.3117522) q[2];
sx q[2];
rz(-1.6774943) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65655639) q[1];
sx q[1];
rz(-1.4244231) q[1];
sx q[1];
rz(-1.3629517) q[1];
rz(-pi) q[2];
rz(2.2150854) q[3];
sx q[3];
rz(-2.4725719) q[3];
sx q[3];
rz(-0.45243087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3847547) q[2];
sx q[2];
rz(-0.33051312) q[2];
sx q[2];
rz(-0.11088863) q[2];
rz(-0.97892654) q[3];
sx q[3];
rz(-1.7612235) q[3];
sx q[3];
rz(0.10262779) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0961618) q[0];
sx q[0];
rz(-1.0935421) q[0];
sx q[0];
rz(-2.9683215) q[0];
rz(2.5896367) q[1];
sx q[1];
rz(-2.5872784) q[1];
sx q[1];
rz(0.90131235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93126585) q[0];
sx q[0];
rz(-1.0823609) q[0];
sx q[0];
rz(-0.68426056) q[0];
x q[1];
rz(-2.5576711) q[2];
sx q[2];
rz(-1.1135227) q[2];
sx q[2];
rz(1.538313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3696987) q[1];
sx q[1];
rz(-1.7189184) q[1];
sx q[1];
rz(-3.1225608) q[1];
rz(-pi) q[2];
rz(-2.3278045) q[3];
sx q[3];
rz(-0.57238676) q[3];
sx q[3];
rz(3.0021964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99894714) q[2];
sx q[2];
rz(-1.877942) q[2];
sx q[2];
rz(-0.20550814) q[2];
rz(2.5213126) q[3];
sx q[3];
rz(-2.5871103) q[3];
sx q[3];
rz(-1.8581355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.810629) q[0];
sx q[0];
rz(-2.0331419) q[0];
sx q[0];
rz(-2.1883645) q[0];
rz(-2.3586121) q[1];
sx q[1];
rz(-0.51376659) q[1];
sx q[1];
rz(-0.89797529) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6655639) q[0];
sx q[0];
rz(-1.7930287) q[0];
sx q[0];
rz(2.1118947) q[0];
rz(-1.1994409) q[2];
sx q[2];
rz(-1.1716649) q[2];
sx q[2];
rz(0.45796203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46875698) q[1];
sx q[1];
rz(-1.3870981) q[1];
sx q[1];
rz(3.0335547) q[1];
rz(-0.19583585) q[3];
sx q[3];
rz(-1.8552498) q[3];
sx q[3];
rz(-1.0468512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39752412) q[2];
sx q[2];
rz(-2.4351138) q[2];
sx q[2];
rz(1.2495329) q[2];
rz(-1.2161829) q[3];
sx q[3];
rz(-1.774923) q[3];
sx q[3];
rz(-1.3528311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-0.80094308) q[0];
sx q[0];
rz(-0.018434374) q[0];
sx q[0];
rz(-1.8311485) q[0];
rz(-2.387923) q[1];
sx q[1];
rz(-0.33652702) q[1];
sx q[1];
rz(0.15187844) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4820712) q[0];
sx q[0];
rz(-0.11733012) q[0];
sx q[0];
rz(2.4767673) q[0];
rz(-pi) q[1];
rz(-0.23252587) q[2];
sx q[2];
rz(-0.33862414) q[2];
sx q[2];
rz(1.2004289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8712403) q[1];
sx q[1];
rz(-2.0899822) q[1];
sx q[1];
rz(-0.097071807) q[1];
rz(-pi) q[2];
rz(-0.88942727) q[3];
sx q[3];
rz(-1.7173036) q[3];
sx q[3];
rz(-1.5832988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14022961) q[2];
sx q[2];
rz(-1.9915853) q[2];
sx q[2];
rz(1.7062194) q[2];
rz(2.5367149) q[3];
sx q[3];
rz(-1.6717654) q[3];
sx q[3];
rz(2.308037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3816933) q[0];
sx q[0];
rz(-0.23867358) q[0];
sx q[0];
rz(2.6651486) q[0];
rz(-2.6051688) q[1];
sx q[1];
rz(-1.8073578) q[1];
sx q[1];
rz(2.776872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8275673) q[0];
sx q[0];
rz(-2.7430581) q[0];
sx q[0];
rz(2.3588595) q[0];
rz(2.6047037) q[2];
sx q[2];
rz(-1.6308954) q[2];
sx q[2];
rz(-2.8837332) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5231001) q[1];
sx q[1];
rz(-0.87311166) q[1];
sx q[1];
rz(-0.24055753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.164806) q[3];
sx q[3];
rz(-1.6315809) q[3];
sx q[3];
rz(-1.2244299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81858188) q[2];
sx q[2];
rz(-1.3280758) q[2];
sx q[2];
rz(-2.2200457) q[2];
rz(-1.796683) q[3];
sx q[3];
rz(-1.4321046) q[3];
sx q[3];
rz(-0.014605453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0607249) q[0];
sx q[0];
rz(-2.9089071) q[0];
sx q[0];
rz(-3.0423855) q[0];
rz(1.5469145) q[1];
sx q[1];
rz(-2.4989936) q[1];
sx q[1];
rz(0.1086455) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0616665) q[0];
sx q[0];
rz(-0.52046973) q[0];
sx q[0];
rz(1.5093278) q[0];
rz(2.7127162) q[2];
sx q[2];
rz(-2.1575973) q[2];
sx q[2];
rz(-1.9029146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2205926) q[1];
sx q[1];
rz(-2.0751168) q[1];
sx q[1];
rz(0.31208529) q[1];
rz(-0.21317706) q[3];
sx q[3];
rz(-1.2090599) q[3];
sx q[3];
rz(-0.26644275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28086153) q[2];
sx q[2];
rz(-2.9555369) q[2];
sx q[2];
rz(0.59530386) q[2];
rz(-0.91296494) q[3];
sx q[3];
rz(-0.76415092) q[3];
sx q[3];
rz(-1.3179774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9106134) q[0];
sx q[0];
rz(-1.9934886) q[0];
sx q[0];
rz(2.8379111) q[0];
rz(0.23823711) q[1];
sx q[1];
rz(-2.1447825) q[1];
sx q[1];
rz(-0.69746596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3179585) q[0];
sx q[0];
rz(-2.1058361) q[0];
sx q[0];
rz(-1.7798502) q[0];
x q[1];
rz(1.6666621) q[2];
sx q[2];
rz(-1.2156475) q[2];
sx q[2];
rz(-2.716316) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1917603) q[1];
sx q[1];
rz(-1.0737541) q[1];
sx q[1];
rz(2.0326896) q[1];
rz(2.8805299) q[3];
sx q[3];
rz(-2.0213642) q[3];
sx q[3];
rz(2.0416259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4136228) q[2];
sx q[2];
rz(-1.1009078) q[2];
sx q[2];
rz(0.093416365) q[2];
rz(1.5591239) q[3];
sx q[3];
rz(-0.070135442) q[3];
sx q[3];
rz(-0.087433405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.14417905) q[0];
sx q[0];
rz(-1.1638673) q[0];
sx q[0];
rz(-1.5735004) q[0];
rz(-0.46918029) q[1];
sx q[1];
rz(-0.73742417) q[1];
sx q[1];
rz(-0.87541568) q[1];
rz(2.1067262) q[2];
sx q[2];
rz(-1.9725261) q[2];
sx q[2];
rz(-0.50289666) q[2];
rz(-1.260626) q[3];
sx q[3];
rz(-2.3295129) q[3];
sx q[3];
rz(-2.2214291) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
