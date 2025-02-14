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
rz(-2.7297468) q[0];
rz(-2.5208199) q[1];
sx q[1];
rz(-2.4359735) q[1];
sx q[1];
rz(0.128428) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0043853) q[0];
sx q[0];
rz(-0.62891987) q[0];
sx q[0];
rz(0.0034751277) q[0];
x q[1];
rz(2.24755) q[2];
sx q[2];
rz(-0.6031361) q[2];
sx q[2];
rz(-1.2007842) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.023019869) q[1];
sx q[1];
rz(-0.65852877) q[1];
sx q[1];
rz(0.3218201) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3462696) q[3];
sx q[3];
rz(-0.40679625) q[3];
sx q[3];
rz(-1.612839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2423557) q[2];
sx q[2];
rz(-0.54348102) q[2];
sx q[2];
rz(-0.71329722) q[2];
rz(1.2171618) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84739581) q[0];
sx q[0];
rz(-1.0260181) q[0];
sx q[0];
rz(-3.078939) q[0];
rz(-2.2368536) q[1];
sx q[1];
rz(-0.86974564) q[1];
sx q[1];
rz(-3.0135801) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0580095) q[0];
sx q[0];
rz(-1.2370652) q[0];
sx q[0];
rz(-0.34837848) q[0];
x q[1];
rz(0.77777783) q[2];
sx q[2];
rz(-0.82077089) q[2];
sx q[2];
rz(0.37015265) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81028013) q[1];
sx q[1];
rz(-1.5638651) q[1];
sx q[1];
rz(-0.48622017) q[1];
rz(-2.4038257) q[3];
sx q[3];
rz(-1.6782266) q[3];
sx q[3];
rz(-0.28001912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55600986) q[2];
sx q[2];
rz(-1.2051008) q[2];
sx q[2];
rz(-1.0949868) q[2];
rz(-0.60947841) q[3];
sx q[3];
rz(-0.947781) q[3];
sx q[3];
rz(1.3013209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13780093) q[0];
sx q[0];
rz(-0.68313685) q[0];
sx q[0];
rz(-2.0004499) q[0];
rz(-1.8423276) q[1];
sx q[1];
rz(-1.7443934) q[1];
sx q[1];
rz(-0.13967791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3052854) q[0];
sx q[0];
rz(-1.7392312) q[0];
sx q[0];
rz(1.832421) q[0];
rz(-pi) q[1];
rz(-1.3422826) q[2];
sx q[2];
rz(-1.7916792) q[2];
sx q[2];
rz(-2.8877838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2217979) q[1];
sx q[1];
rz(-1.9842902) q[1];
sx q[1];
rz(-1.7899124) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3862429) q[3];
sx q[3];
rz(-2.3358029) q[3];
sx q[3];
rz(-0.84945018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37207681) q[2];
sx q[2];
rz(-2.3774827) q[2];
sx q[2];
rz(-0.3176983) q[2];
rz(-0.54567537) q[3];
sx q[3];
rz(-2.217644) q[3];
sx q[3];
rz(-1.3056477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1691771) q[0];
sx q[0];
rz(-1.7298052) q[0];
sx q[0];
rz(0.2485982) q[0];
rz(-1.2480805) q[1];
sx q[1];
rz(-1.0634407) q[1];
sx q[1];
rz(-0.6691106) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1951576) q[0];
sx q[0];
rz(-1.666233) q[0];
sx q[0];
rz(2.6470031) q[0];
x q[1];
rz(-0.76478473) q[2];
sx q[2];
rz(-1.7953821) q[2];
sx q[2];
rz(-3.0377521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8324171) q[1];
sx q[1];
rz(-2.8879893) q[1];
sx q[1];
rz(0.95043851) q[1];
rz(-2.1345244) q[3];
sx q[3];
rz(-1.9525212) q[3];
sx q[3];
rz(-1.6508716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3847547) q[2];
sx q[2];
rz(-0.33051312) q[2];
sx q[2];
rz(0.11088863) q[2];
rz(-2.1626661) q[3];
sx q[3];
rz(-1.7612235) q[3];
sx q[3];
rz(-0.10262779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0454309) q[0];
sx q[0];
rz(-2.0480506) q[0];
sx q[0];
rz(0.1732711) q[0];
rz(2.5896367) q[1];
sx q[1];
rz(-0.55431429) q[1];
sx q[1];
rz(-0.90131235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93126585) q[0];
sx q[0];
rz(-2.0592318) q[0];
sx q[0];
rz(2.4573321) q[0];
x q[1];
rz(2.5576711) q[2];
sx q[2];
rz(-2.02807) q[2];
sx q[2];
rz(1.538313) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9001604) q[1];
sx q[1];
rz(-0.14933085) q[1];
sx q[1];
rz(-1.6976507) q[1];
rz(2.3278045) q[3];
sx q[3];
rz(-0.57238676) q[3];
sx q[3];
rz(0.13939626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1426455) q[2];
sx q[2];
rz(-1.877942) q[2];
sx q[2];
rz(-2.9360845) q[2];
rz(0.62028003) q[3];
sx q[3];
rz(-2.5871103) q[3];
sx q[3];
rz(1.8581355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.810629) q[0];
sx q[0];
rz(-2.0331419) q[0];
sx q[0];
rz(-2.1883645) q[0];
rz(-0.78298059) q[1];
sx q[1];
rz(-0.51376659) q[1];
sx q[1];
rz(-2.2436174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1785068) q[0];
sx q[0];
rz(-2.0971812) q[0];
sx q[0];
rz(0.25776074) q[0];
rz(2.7165604) q[2];
sx q[2];
rz(-1.9117182) q[2];
sx q[2];
rz(1.2630315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46875698) q[1];
sx q[1];
rz(-1.3870981) q[1];
sx q[1];
rz(3.0335547) q[1];
rz(-1.2811019) q[3];
sx q[3];
rz(-1.3829261) q[3];
sx q[3];
rz(-0.46833098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7440685) q[2];
sx q[2];
rz(-0.70647883) q[2];
sx q[2];
rz(1.8920598) q[2];
rz(1.2161829) q[3];
sx q[3];
rz(-1.774923) q[3];
sx q[3];
rz(-1.7887615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80094308) q[0];
sx q[0];
rz(-0.018434374) q[0];
sx q[0];
rz(1.3104441) q[0];
rz(2.387923) q[1];
sx q[1];
rz(-0.33652702) q[1];
sx q[1];
rz(-0.15187844) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5688393) q[0];
sx q[0];
rz(-1.4985159) q[0];
sx q[0];
rz(0.092503017) q[0];
rz(-pi) q[1];
rz(-0.33016684) q[2];
sx q[2];
rz(-1.6474198) q[2];
sx q[2];
rz(-2.9909796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0650349) q[1];
sx q[1];
rz(-0.5273653) q[1];
sx q[1];
rz(-1.7387912) q[1];
rz(-2.9538438) q[3];
sx q[3];
rz(-0.89808447) q[3];
sx q[3];
rz(0.10533145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.001363) q[2];
sx q[2];
rz(-1.1500074) q[2];
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
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3816933) q[0];
sx q[0];
rz(-2.9029191) q[0];
sx q[0];
rz(0.47644404) q[0];
rz(-0.53642383) q[1];
sx q[1];
rz(-1.8073578) q[1];
sx q[1];
rz(0.36472067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8275673) q[0];
sx q[0];
rz(-0.3985346) q[0];
sx q[0];
rz(0.78273313) q[0];
rz(-3.0244855) q[2];
sx q[2];
rz(-2.6016781) q[2];
sx q[2];
rz(-1.929226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9832416) q[1];
sx q[1];
rz(-2.4102306) q[1];
sx q[1];
rz(1.29391) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0682848) q[3];
sx q[3];
rz(-2.1635593) q[3];
sx q[3];
rz(2.7542219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81858188) q[2];
sx q[2];
rz(-1.8135169) q[2];
sx q[2];
rz(-0.92154694) q[2];
rz(1.796683) q[3];
sx q[3];
rz(-1.4321046) q[3];
sx q[3];
rz(-3.1269872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.0607249) q[0];
sx q[0];
rz(-2.9089071) q[0];
sx q[0];
rz(3.0423855) q[0];
rz(-1.5469145) q[1];
sx q[1];
rz(-0.64259905) q[1];
sx q[1];
rz(-3.0329472) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324873) q[0];
sx q[0];
rz(-2.0901839) q[0];
sx q[0];
rz(0.035196134) q[0];
x q[1];
rz(0.4288765) q[2];
sx q[2];
rz(-0.9839954) q[2];
sx q[2];
rz(1.238678) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8092919) q[1];
sx q[1];
rz(-0.58590472) q[1];
sx q[1];
rz(-2.078474) q[1];
rz(-2.9284156) q[3];
sx q[3];
rz(-1.2090599) q[3];
sx q[3];
rz(0.26644275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28086153) q[2];
sx q[2];
rz(-2.9555369) q[2];
sx q[2];
rz(-0.59530386) q[2];
rz(-2.2286277) q[3];
sx q[3];
rz(-0.76415092) q[3];
sx q[3];
rz(1.3179774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9106134) q[0];
sx q[0];
rz(-1.148104) q[0];
sx q[0];
rz(-0.30368152) q[0];
rz(-2.9033555) q[1];
sx q[1];
rz(-2.1447825) q[1];
sx q[1];
rz(2.4441267) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36059067) q[0];
sx q[0];
rz(-1.7502898) q[0];
sx q[0];
rz(0.54474564) q[0];
rz(2.8890144) q[2];
sx q[2];
rz(-0.36732964) q[2];
sx q[2];
rz(0.15548104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38790801) q[1];
sx q[1];
rz(-1.1682967) q[1];
sx q[1];
rz(0.54479568) q[1];
rz(-0.26106278) q[3];
sx q[3];
rz(-2.0213642) q[3];
sx q[3];
rz(-1.0999668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72796983) q[2];
sx q[2];
rz(-1.1009078) q[2];
sx q[2];
rz(-3.0481763) q[2];
rz(-1.5824687) q[3];
sx q[3];
rz(-3.0714572) q[3];
sx q[3];
rz(0.087433405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14417905) q[0];
sx q[0];
rz(-1.9777254) q[0];
sx q[0];
rz(1.5680922) q[0];
rz(-2.6724124) q[1];
sx q[1];
rz(-2.4041685) q[1];
sx q[1];
rz(2.266177) q[1];
rz(-2.1067262) q[2];
sx q[2];
rz(-1.1690665) q[2];
sx q[2];
rz(2.638696) q[2];
rz(-0.78315121) q[3];
sx q[3];
rz(-1.7941536) q[3];
sx q[3];
rz(-0.4335946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
