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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1415048) q[0];
sx q[0];
rz(-0.94188085) q[0];
sx q[0];
rz(1.5733243) q[0];
rz(-pi) q[1];
rz(-2.0635705) q[2];
sx q[2];
rz(-1.2076305) q[2];
sx q[2];
rz(0.21445477) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1185728) q[1];
sx q[1];
rz(-0.65852877) q[1];
sx q[1];
rz(2.8197726) q[1];
rz(-pi) q[2];
rz(1.7953231) q[3];
sx q[3];
rz(-0.40679625) q[3];
sx q[3];
rz(-1.612839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8992369) q[2];
sx q[2];
rz(-2.5981116) q[2];
sx q[2];
rz(-2.4282954) q[2];
rz(-1.9244309) q[3];
sx q[3];
rz(-1.4274012) q[3];
sx q[3];
rz(-0.10669073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.2941968) q[0];
sx q[0];
rz(-1.0260181) q[0];
sx q[0];
rz(0.062653616) q[0];
rz(2.2368536) q[1];
sx q[1];
rz(-2.271847) q[1];
sx q[1];
rz(-3.0135801) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9205039) q[0];
sx q[0];
rz(-2.663923) q[0];
sx q[0];
rz(-2.3484557) q[0];
rz(-pi) q[1];
rz(-0.77777783) q[2];
sx q[2];
rz(-2.3208218) q[2];
sx q[2];
rz(-2.77144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76417962) q[1];
sx q[1];
rz(-1.0845889) q[1];
sx q[1];
rz(1.5629565) q[1];
x q[2];
rz(1.4260728) q[3];
sx q[3];
rz(-2.3033352) q[3];
sx q[3];
rz(-1.7536557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55600986) q[2];
sx q[2];
rz(-1.2051008) q[2];
sx q[2];
rz(1.0949868) q[2];
rz(2.5321142) q[3];
sx q[3];
rz(-2.1938117) q[3];
sx q[3];
rz(-1.3013209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13780093) q[0];
sx q[0];
rz(-2.4584558) q[0];
sx q[0];
rz(-2.0004499) q[0];
rz(-1.2992651) q[1];
sx q[1];
rz(-1.3971993) q[1];
sx q[1];
rz(3.0019147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8477301) q[0];
sx q[0];
rz(-0.31010695) q[0];
sx q[0];
rz(0.98921138) q[0];
rz(1.79931) q[2];
sx q[2];
rz(-1.3499134) q[2];
sx q[2];
rz(2.8877838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91979474) q[1];
sx q[1];
rz(-1.9842902) q[1];
sx q[1];
rz(-1.3516803) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3862429) q[3];
sx q[3];
rz(-0.80578973) q[3];
sx q[3];
rz(-0.84945018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7695158) q[2];
sx q[2];
rz(-2.3774827) q[2];
sx q[2];
rz(0.3176983) q[2];
rz(-0.54567537) q[3];
sx q[3];
rz(-0.92394865) q[3];
sx q[3];
rz(1.3056477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1691771) q[0];
sx q[0];
rz(-1.4117874) q[0];
sx q[0];
rz(0.2485982) q[0];
rz(-1.8935122) q[1];
sx q[1];
rz(-2.0781519) q[1];
sx q[1];
rz(2.4724821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4658858) q[0];
sx q[0];
rz(-1.0786593) q[0];
sx q[0];
rz(-1.6791315) q[0];
x q[1];
rz(-2.3768079) q[2];
sx q[2];
rz(-1.3462105) q[2];
sx q[2];
rz(-3.0377521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30917558) q[1];
sx q[1];
rz(-2.8879893) q[1];
sx q[1];
rz(-2.1911541) q[1];
x q[2];
rz(-2.2150854) q[3];
sx q[3];
rz(-0.66902071) q[3];
sx q[3];
rz(2.6891618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3847547) q[2];
sx q[2];
rz(-0.33051312) q[2];
sx q[2];
rz(0.11088863) q[2];
rz(2.1626661) q[3];
sx q[3];
rz(-1.7612235) q[3];
sx q[3];
rz(0.10262779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0961618) q[0];
sx q[0];
rz(-1.0935421) q[0];
sx q[0];
rz(-0.1732711) q[0];
rz(2.5896367) q[1];
sx q[1];
rz(-0.55431429) q[1];
sx q[1];
rz(-0.90131235) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1364897) q[0];
sx q[0];
rz(-2.1629961) q[0];
sx q[0];
rz(-2.1718958) q[0];
rz(-pi) q[1];
rz(2.5576711) q[2];
sx q[2];
rz(-1.1135227) q[2];
sx q[2];
rz(1.6032796) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9001604) q[1];
sx q[1];
rz(-2.9922618) q[1];
sx q[1];
rz(1.6976507) q[1];
x q[2];
rz(-1.1327733) q[3];
sx q[3];
rz(-1.1896648) q[3];
sx q[3];
rz(-2.1025859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99894714) q[2];
sx q[2];
rz(-1.2636507) q[2];
sx q[2];
rz(-0.20550814) q[2];
rz(0.62028003) q[3];
sx q[3];
rz(-2.5871103) q[3];
sx q[3];
rz(1.8581355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.810629) q[0];
sx q[0];
rz(-1.1084508) q[0];
sx q[0];
rz(-0.95322815) q[0];
rz(-2.3586121) q[1];
sx q[1];
rz(-0.51376659) q[1];
sx q[1];
rz(-0.89797529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.695279) q[0];
sx q[0];
rz(-2.5608664) q[0];
sx q[0];
rz(1.1573791) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9421517) q[2];
sx q[2];
rz(-1.1716649) q[2];
sx q[2];
rz(2.6836306) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0745308) q[1];
sx q[1];
rz(-0.21280405) q[1];
sx q[1];
rz(2.0966543) q[1];
rz(-pi) q[2];
rz(0.19583585) q[3];
sx q[3];
rz(-1.8552498) q[3];
sx q[3];
rz(1.0468512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7440685) q[2];
sx q[2];
rz(-2.4351138) q[2];
sx q[2];
rz(1.2495329) q[2];
rz(1.2161829) q[3];
sx q[3];
rz(-1.3666697) q[3];
sx q[3];
rz(1.7887615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3406496) q[0];
sx q[0];
rz(-0.018434374) q[0];
sx q[0];
rz(-1.3104441) q[0];
rz(0.75366968) q[1];
sx q[1];
rz(-2.8050656) q[1];
sx q[1];
rz(-0.15187844) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150249) q[0];
sx q[0];
rz(-1.6630571) q[0];
sx q[0];
rz(-1.6433861) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33016684) q[2];
sx q[2];
rz(-1.4941729) q[2];
sx q[2];
rz(-0.1506131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.79287) q[1];
sx q[1];
rz(-1.6550437) q[1];
sx q[1];
rz(2.0920175) q[1];
rz(-pi) q[2];
rz(0.88942727) q[3];
sx q[3];
rz(-1.7173036) q[3];
sx q[3];
rz(-1.5582939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.001363) q[2];
sx q[2];
rz(-1.1500074) q[2];
sx q[2];
rz(1.7062194) q[2];
rz(0.6048778) q[3];
sx q[3];
rz(-1.4698272) q[3];
sx q[3];
rz(-0.83355561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3816933) q[0];
sx q[0];
rz(-2.9029191) q[0];
sx q[0];
rz(0.47644404) q[0];
rz(0.53642383) q[1];
sx q[1];
rz(-1.3342349) q[1];
sx q[1];
rz(0.36472067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49052377) q[0];
sx q[0];
rz(-1.2920652) q[0];
sx q[0];
rz(-1.8594478) q[0];
rz(-pi) q[1];
rz(-2.6047037) q[2];
sx q[2];
rz(-1.5106972) q[2];
sx q[2];
rz(-2.8837332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5231001) q[1];
sx q[1];
rz(-0.87311166) q[1];
sx q[1];
rz(0.24055753) q[1];
rz(-pi) q[2];
rz(-2.164806) q[3];
sx q[3];
rz(-1.6315809) q[3];
sx q[3];
rz(1.2244299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3230108) q[2];
sx q[2];
rz(-1.8135169) q[2];
sx q[2];
rz(2.2200457) q[2];
rz(1.796683) q[3];
sx q[3];
rz(-1.4321046) q[3];
sx q[3];
rz(0.014605453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0607249) q[0];
sx q[0];
rz(-2.9089071) q[0];
sx q[0];
rz(3.0423855) q[0];
rz(-1.5946782) q[1];
sx q[1];
rz(-0.64259905) q[1];
sx q[1];
rz(-0.1086455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091054) q[0];
sx q[0];
rz(-2.0901839) q[0];
sx q[0];
rz(-0.035196134) q[0];
rz(-pi) q[1];
rz(-2.1296839) q[2];
sx q[2];
rz(-0.71162738) q[2];
sx q[2];
rz(-2.5932082) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2205926) q[1];
sx q[1];
rz(-1.0664758) q[1];
sx q[1];
rz(0.31208529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0806198) q[3];
sx q[3];
rz(-0.41748306) q[3];
sx q[3];
rz(0.28250697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8607311) q[2];
sx q[2];
rz(-0.18605575) q[2];
sx q[2];
rz(0.59530386) q[2];
rz(2.2286277) q[3];
sx q[3];
rz(-2.3774417) q[3];
sx q[3];
rz(1.3179774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9106134) q[0];
sx q[0];
rz(-1.9934886) q[0];
sx q[0];
rz(0.30368152) q[0];
rz(2.9033555) q[1];
sx q[1];
rz(-0.99681011) q[1];
sx q[1];
rz(2.4441267) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.217933) q[0];
sx q[0];
rz(-0.57070792) q[0];
sx q[0];
rz(-2.804787) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6666621) q[2];
sx q[2];
rz(-1.2156475) q[2];
sx q[2];
rz(0.42527664) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38790801) q[1];
sx q[1];
rz(-1.973296) q[1];
sx q[1];
rz(0.54479568) q[1];
x q[2];
rz(-2.0350214) q[3];
sx q[3];
rz(-1.3363049) q[3];
sx q[3];
rz(0.35500832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4136228) q[2];
sx q[2];
rz(-2.0406849) q[2];
sx q[2];
rz(0.093416365) q[2];
rz(-1.5824687) q[3];
sx q[3];
rz(-0.070135442) q[3];
sx q[3];
rz(3.0541592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9974136) q[0];
sx q[0];
rz(-1.9777254) q[0];
sx q[0];
rz(1.5680922) q[0];
rz(-2.6724124) q[1];
sx q[1];
rz(-2.4041685) q[1];
sx q[1];
rz(2.266177) q[1];
rz(-2.2647247) q[2];
sx q[2];
rz(-0.65779455) q[2];
sx q[2];
rz(0.48566512) q[2];
rz(0.78315121) q[3];
sx q[3];
rz(-1.3474391) q[3];
sx q[3];
rz(2.7079981) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
