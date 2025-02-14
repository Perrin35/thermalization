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
rz(1.5682683) q[0];
rz(-pi) q[1];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-1.8516396) q[1];
sx q[1];
rz(-1.3760097) q[1];
sx q[1];
rz(-2.508393) q[1];
x q[2];
rz(-1.1731568) q[3];
sx q[3];
rz(-1.4825882) q[3];
sx q[3];
rz(2.8928066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8992369) q[2];
sx q[2];
rz(-2.5981116) q[2];
sx q[2];
rz(-0.71329722) q[2];
rz(-1.9244309) q[3];
sx q[3];
rz(-1.7141914) q[3];
sx q[3];
rz(-3.0349019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941968) q[0];
sx q[0];
rz(-1.0260181) q[0];
sx q[0];
rz(3.078939) q[0];
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
rz(-2.9205039) q[0];
sx q[0];
rz(-0.47766968) q[0];
sx q[0];
rz(-2.3484557) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77777783) q[2];
sx q[2];
rz(-0.82077089) q[2];
sx q[2];
rz(-2.77144) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81028013) q[1];
sx q[1];
rz(-1.5638651) q[1];
sx q[1];
rz(2.6553725) q[1];
x q[2];
rz(1.4260728) q[3];
sx q[3];
rz(-0.83825745) q[3];
sx q[3];
rz(1.7536557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.55600986) q[2];
sx q[2];
rz(-1.9364919) q[2];
sx q[2];
rz(-2.0466059) q[2];
rz(-2.5321142) q[3];
sx q[3];
rz(-0.947781) q[3];
sx q[3];
rz(-1.3013209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13780093) q[0];
sx q[0];
rz(-0.68313685) q[0];
sx q[0];
rz(-1.1411427) q[0];
rz(1.8423276) q[1];
sx q[1];
rz(-1.3971993) q[1];
sx q[1];
rz(-0.13967791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6896317) q[0];
sx q[0];
rz(-1.312959) q[0];
sx q[0];
rz(-0.17425107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.79931) q[2];
sx q[2];
rz(-1.7916792) q[2];
sx q[2];
rz(0.2538089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7157212) q[1];
sx q[1];
rz(-0.46498659) q[1];
sx q[1];
rz(-2.6816508) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1909149) q[3];
sx q[3];
rz(-2.1237299) q[3];
sx q[3];
rz(-1.7861136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7695158) q[2];
sx q[2];
rz(-0.76411) q[2];
sx q[2];
rz(0.3176983) q[2];
rz(2.5959173) q[3];
sx q[3];
rz(-2.217644) q[3];
sx q[3];
rz(-1.3056477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9724156) q[0];
sx q[0];
rz(-1.4117874) q[0];
sx q[0];
rz(-0.2485982) q[0];
rz(-1.8935122) q[1];
sx q[1];
rz(-1.0634407) q[1];
sx q[1];
rz(0.6691106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1951576) q[0];
sx q[0];
rz(-1.666233) q[0];
sx q[0];
rz(-0.49458953) q[0];
rz(0.76478473) q[2];
sx q[2];
rz(-1.7953821) q[2];
sx q[2];
rz(-0.10384053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.65655639) q[1];
sx q[1];
rz(-1.7171696) q[1];
sx q[1];
rz(1.778641) q[1];
rz(-pi) q[2];
rz(-2.6982299) q[3];
sx q[3];
rz(-1.0519472) q[3];
sx q[3];
rz(0.31135294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3847547) q[2];
sx q[2];
rz(-0.33051312) q[2];
sx q[2];
rz(0.11088863) q[2];
rz(-0.97892654) q[3];
sx q[3];
rz(-1.3803692) q[3];
sx q[3];
rz(3.0389649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0961618) q[0];
sx q[0];
rz(-2.0480506) q[0];
sx q[0];
rz(-0.1732711) q[0];
rz(0.551956) q[1];
sx q[1];
rz(-0.55431429) q[1];
sx q[1];
rz(-2.2402803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2103268) q[0];
sx q[0];
rz(-2.0592318) q[0];
sx q[0];
rz(2.4573321) q[0];
x q[1];
rz(-0.72867877) q[2];
sx q[2];
rz(-2.4167293) q[2];
sx q[2];
rz(0.55654364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.771894) q[1];
sx q[1];
rz(-1.7189184) q[1];
sx q[1];
rz(-3.1225608) q[1];
rz(-pi) q[2];
rz(-2.0088193) q[3];
sx q[3];
rz(-1.1896648) q[3];
sx q[3];
rz(-1.0390067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1426455) q[2];
sx q[2];
rz(-1.877942) q[2];
sx q[2];
rz(2.9360845) q[2];
rz(2.5213126) q[3];
sx q[3];
rz(-0.55448237) q[3];
sx q[3];
rz(-1.2834572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3309636) q[0];
sx q[0];
rz(-2.0331419) q[0];
sx q[0];
rz(0.95322815) q[0];
rz(0.78298059) q[1];
sx q[1];
rz(-0.51376659) q[1];
sx q[1];
rz(2.2436174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.695279) q[0];
sx q[0];
rz(-0.58072621) q[0];
sx q[0];
rz(1.1573791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71048471) q[2];
sx q[2];
rz(-0.53827751) q[2];
sx q[2];
rz(0.32847095) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46875698) q[1];
sx q[1];
rz(-1.7544946) q[1];
sx q[1];
rz(0.10803799) q[1];
rz(-pi) q[2];
rz(1.2811019) q[3];
sx q[3];
rz(-1.7586666) q[3];
sx q[3];
rz(-0.46833098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39752412) q[2];
sx q[2];
rz(-2.4351138) q[2];
sx q[2];
rz(-1.8920598) q[2];
rz(1.2161829) q[3];
sx q[3];
rz(-1.3666697) q[3];
sx q[3];
rz(1.7887615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3406496) q[0];
sx q[0];
rz(-3.1231583) q[0];
sx q[0];
rz(-1.8311485) q[0];
rz(2.387923) q[1];
sx q[1];
rz(-0.33652702) q[1];
sx q[1];
rz(2.9897142) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5688393) q[0];
sx q[0];
rz(-1.4985159) q[0];
sx q[0];
rz(-3.0490896) q[0];
rz(-pi) q[1];
rz(-2.9090668) q[2];
sx q[2];
rz(-2.8029685) q[2];
sx q[2];
rz(-1.9411638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0765578) q[1];
sx q[1];
rz(-2.6142274) q[1];
sx q[1];
rz(1.4028014) q[1];
rz(-pi) q[2];
rz(1.8009284) q[3];
sx q[3];
rz(-0.69447434) q[3];
sx q[3];
rz(-2.9509773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14022961) q[2];
sx q[2];
rz(-1.9915853) q[2];
sx q[2];
rz(-1.7062194) q[2];
rz(-2.5367149) q[3];
sx q[3];
rz(-1.6717654) q[3];
sx q[3];
rz(0.83355561) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75989938) q[0];
sx q[0];
rz(-0.23867358) q[0];
sx q[0];
rz(2.6651486) q[0];
rz(2.6051688) q[1];
sx q[1];
rz(-1.8073578) q[1];
sx q[1];
rz(0.36472067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428392) q[0];
sx q[0];
rz(-1.2935862) q[0];
sx q[0];
rz(-0.2901092) q[0];
x q[1];
rz(-0.1171072) q[2];
sx q[2];
rz(-0.53991452) q[2];
sx q[2];
rz(-1.929226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9375853) q[1];
sx q[1];
rz(-1.3871925) q[1];
sx q[1];
rz(0.85869243) q[1];
rz(-pi) q[2];
rz(0.97678661) q[3];
sx q[3];
rz(-1.6315809) q[3];
sx q[3];
rz(-1.9171627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3230108) q[2];
sx q[2];
rz(-1.8135169) q[2];
sx q[2];
rz(-2.2200457) q[2];
rz(-1.3449097) q[3];
sx q[3];
rz(-1.709488) q[3];
sx q[3];
rz(-0.014605453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
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
rz(-1.5946782) q[1];
sx q[1];
rz(-0.64259905) q[1];
sx q[1];
rz(3.0329472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324873) q[0];
sx q[0];
rz(-1.0514088) q[0];
sx q[0];
rz(-0.035196134) q[0];
rz(1.0119087) q[2];
sx q[2];
rz(-0.71162738) q[2];
sx q[2];
rz(-2.5932082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33230072) q[1];
sx q[1];
rz(-2.5556879) q[1];
sx q[1];
rz(2.078474) q[1];
x q[2];
rz(-1.2014162) q[3];
sx q[3];
rz(-1.7699852) q[3];
sx q[3];
rz(1.7607813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28086153) q[2];
sx q[2];
rz(-0.18605575) q[2];
sx q[2];
rz(-0.59530386) q[2];
rz(0.91296494) q[3];
sx q[3];
rz(-2.3774417) q[3];
sx q[3];
rz(-1.3179774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2309793) q[0];
sx q[0];
rz(-1.9934886) q[0];
sx q[0];
rz(2.8379111) q[0];
rz(-2.9033555) q[1];
sx q[1];
rz(-0.99681011) q[1];
sx q[1];
rz(-2.4441267) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92365962) q[0];
sx q[0];
rz(-2.5708847) q[0];
sx q[0];
rz(2.804787) q[0];
rz(1.6666621) q[2];
sx q[2];
rz(-1.2156475) q[2];
sx q[2];
rz(-2.716316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38790801) q[1];
sx q[1];
rz(-1.973296) q[1];
sx q[1];
rz(-2.596797) q[1];
x q[2];
rz(-2.0609217) q[3];
sx q[3];
rz(-0.5162067) q[3];
sx q[3];
rz(1.4913477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72796983) q[2];
sx q[2];
rz(-2.0406849) q[2];
sx q[2];
rz(0.093416365) q[2];
rz(1.5591239) q[3];
sx q[3];
rz(-3.0714572) q[3];
sx q[3];
rz(0.087433405) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14417905) q[0];
sx q[0];
rz(-1.9777254) q[0];
sx q[0];
rz(1.5680922) q[0];
rz(2.6724124) q[1];
sx q[1];
rz(-0.73742417) q[1];
sx q[1];
rz(-0.87541568) q[1];
rz(0.45892528) q[2];
sx q[2];
rz(-1.0815207) q[2];
sx q[2];
rz(-1.8455119) q[2];
rz(2.3584414) q[3];
sx q[3];
rz(-1.7941536) q[3];
sx q[3];
rz(-0.4335946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
