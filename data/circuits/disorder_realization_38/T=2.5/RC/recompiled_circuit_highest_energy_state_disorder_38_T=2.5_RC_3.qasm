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
rz(0.83346811) q[0];
sx q[0];
rz(-0.8136189) q[0];
sx q[0];
rz(-1.5669493) q[0];
rz(-1.1235224) q[1];
sx q[1];
rz(-0.82668537) q[1];
sx q[1];
rz(0.5067265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08164601) q[0];
sx q[0];
rz(-1.5525408) q[0];
sx q[0];
rz(0.55299802) q[0];
rz(-pi) q[1];
rz(0.2519744) q[2];
sx q[2];
rz(-1.0720382) q[2];
sx q[2];
rz(0.70372561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6804569) q[1];
sx q[1];
rz(-1.3142921) q[1];
sx q[1];
rz(1.6272904) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2351609) q[3];
sx q[3];
rz(-1.6015823) q[3];
sx q[3];
rz(-2.8828636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7099521) q[2];
sx q[2];
rz(-1.1974572) q[2];
sx q[2];
rz(0.54090995) q[2];
rz(-1.3276118) q[3];
sx q[3];
rz(-1.4541459) q[3];
sx q[3];
rz(-0.66837627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79880181) q[0];
sx q[0];
rz(-1.3725766) q[0];
sx q[0];
rz(-2.0846682) q[0];
rz(2.7150555) q[1];
sx q[1];
rz(-1.6273472) q[1];
sx q[1];
rz(-0.0043409745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26165923) q[0];
sx q[0];
rz(-0.78908217) q[0];
sx q[0];
rz(-2.1270149) q[0];
rz(-1.0356194) q[2];
sx q[2];
rz(-0.7053203) q[2];
sx q[2];
rz(0.83503228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7340034) q[1];
sx q[1];
rz(-0.079691039) q[1];
sx q[1];
rz(-1.4400287) q[1];
x q[2];
rz(0.53496302) q[3];
sx q[3];
rz(-1.5175948) q[3];
sx q[3];
rz(-0.0004079013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9428375) q[2];
sx q[2];
rz(-1.402907) q[2];
sx q[2];
rz(2.7665561) q[2];
rz(0.047920553) q[3];
sx q[3];
rz(-1.5791357) q[3];
sx q[3];
rz(0.16938773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7584943) q[0];
sx q[0];
rz(-2.5653745) q[0];
sx q[0];
rz(-1.6297485) q[0];
rz(2.0747298) q[1];
sx q[1];
rz(-1.8424415) q[1];
sx q[1];
rz(-2.3323434) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95308304) q[0];
sx q[0];
rz(-3.0874757) q[0];
sx q[0];
rz(2.3982766) q[0];
x q[1];
rz(0.079945076) q[2];
sx q[2];
rz(-0.71942753) q[2];
sx q[2];
rz(-1.5768752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36673646) q[1];
sx q[1];
rz(-1.7637758) q[1];
sx q[1];
rz(0.79905135) q[1];
rz(-1.5736846) q[3];
sx q[3];
rz(-1.2089855) q[3];
sx q[3];
rz(2.1199292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92144907) q[2];
sx q[2];
rz(-2.4704832) q[2];
sx q[2];
rz(-1.5044093) q[2];
rz(0.97357059) q[3];
sx q[3];
rz(-0.74876553) q[3];
sx q[3];
rz(1.0934632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022738986) q[0];
sx q[0];
rz(-1.5819419) q[0];
sx q[0];
rz(2.7795025) q[0];
rz(-1.7936961) q[1];
sx q[1];
rz(-1.7575512) q[1];
sx q[1];
rz(2.0787584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2971173) q[0];
sx q[0];
rz(-1.7430291) q[0];
sx q[0];
rz(-0.14396859) q[0];
rz(-1.8824018) q[2];
sx q[2];
rz(-0.62436337) q[2];
sx q[2];
rz(-1.7758689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3058051) q[1];
sx q[1];
rz(-1.2286245) q[1];
sx q[1];
rz(0.31098109) q[1];
x q[2];
rz(-1.7881259) q[3];
sx q[3];
rz(-1.8193075) q[3];
sx q[3];
rz(-2.0686764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41703835) q[2];
sx q[2];
rz(-2.3596256) q[2];
sx q[2];
rz(-0.057272591) q[2];
rz(0.96632424) q[3];
sx q[3];
rz(-0.94760197) q[3];
sx q[3];
rz(1.2053325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22369497) q[0];
sx q[0];
rz(-0.88858336) q[0];
sx q[0];
rz(0.59820557) q[0];
rz(-1.573645) q[1];
sx q[1];
rz(-1.4718098) q[1];
sx q[1];
rz(2.9817458) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2837958) q[0];
sx q[0];
rz(-2.3295053) q[0];
sx q[0];
rz(-1.8079197) q[0];
rz(-pi) q[1];
rz(2.9303126) q[2];
sx q[2];
rz(-1.6307555) q[2];
sx q[2];
rz(2.243239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14434334) q[1];
sx q[1];
rz(-2.1806988) q[1];
sx q[1];
rz(-2.039012) q[1];
x q[2];
rz(0.076301612) q[3];
sx q[3];
rz(-1.5263183) q[3];
sx q[3];
rz(0.22212576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9827031) q[2];
sx q[2];
rz(-0.44124678) q[2];
sx q[2];
rz(-0.59054792) q[2];
rz(-0.96747893) q[3];
sx q[3];
rz(-1.5891821) q[3];
sx q[3];
rz(1.1837186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443053) q[0];
sx q[0];
rz(-0.66216457) q[0];
sx q[0];
rz(2.1024607) q[0];
rz(-2.2100673) q[1];
sx q[1];
rz(-2.4577699) q[1];
sx q[1];
rz(1.1087803) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241348) q[0];
sx q[0];
rz(-2.0724247) q[0];
sx q[0];
rz(-0.62539102) q[0];
x q[1];
rz(-2.4007892) q[2];
sx q[2];
rz(-2.1880297) q[2];
sx q[2];
rz(2.115135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30516827) q[1];
sx q[1];
rz(-0.9988316) q[1];
sx q[1];
rz(-1.9981153) q[1];
rz(-pi) q[2];
rz(-2.8907701) q[3];
sx q[3];
rz(-1.0366813) q[3];
sx q[3];
rz(-3.0297584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80591566) q[2];
sx q[2];
rz(-2.6695965) q[2];
sx q[2];
rz(-0.3042039) q[2];
rz(-0.16600569) q[3];
sx q[3];
rz(-1.3622354) q[3];
sx q[3];
rz(-1.6148875) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45749998) q[0];
sx q[0];
rz(-2.9170211) q[0];
sx q[0];
rz(1.4373047) q[0];
rz(0.2416839) q[1];
sx q[1];
rz(-1.6432523) q[1];
sx q[1];
rz(2.2041352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3110549) q[0];
sx q[0];
rz(-1.8094581) q[0];
sx q[0];
rz(-0.72913362) q[0];
rz(1.8514013) q[2];
sx q[2];
rz(-1.6001424) q[2];
sx q[2];
rz(1.9179232) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73956758) q[1];
sx q[1];
rz(-1.1023542) q[1];
sx q[1];
rz(-0.7206638) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.203905) q[3];
sx q[3];
rz(-2.4613639) q[3];
sx q[3];
rz(2.2474248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23718701) q[2];
sx q[2];
rz(-1.097647) q[2];
sx q[2];
rz(-2.912525) q[2];
rz(1.1711052) q[3];
sx q[3];
rz(-1.7541211) q[3];
sx q[3];
rz(-2.4598725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197297) q[0];
sx q[0];
rz(-2.5010112) q[0];
sx q[0];
rz(1.7280475) q[0];
rz(1.9238663) q[1];
sx q[1];
rz(-1.7037337) q[1];
sx q[1];
rz(0.53122836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41252688) q[0];
sx q[0];
rz(-0.382889) q[0];
sx q[0];
rz(2.114243) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81905535) q[2];
sx q[2];
rz(-1.8764638) q[2];
sx q[2];
rz(1.7421987) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6636017) q[1];
sx q[1];
rz(-0.79959005) q[1];
sx q[1];
rz(-2.5065304) q[1];
x q[2];
rz(1.8645499) q[3];
sx q[3];
rz(-1.2790907) q[3];
sx q[3];
rz(-0.13495378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5765877) q[2];
sx q[2];
rz(-1.7030623) q[2];
sx q[2];
rz(2.3459404) q[2];
rz(-2.3203881) q[3];
sx q[3];
rz(-0.20902769) q[3];
sx q[3];
rz(-2.8388099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32383701) q[0];
sx q[0];
rz(-0.19446401) q[0];
sx q[0];
rz(-1.4870462) q[0];
rz(-1.5215634) q[1];
sx q[1];
rz(-1.227102) q[1];
sx q[1];
rz(1.0618658) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8811765) q[0];
sx q[0];
rz(-1.9435391) q[0];
sx q[0];
rz(0.56074067) q[0];
x q[1];
rz(-1.2765116) q[2];
sx q[2];
rz(-1.6194034) q[2];
sx q[2];
rz(2.693134) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7567021) q[1];
sx q[1];
rz(-1.1497918) q[1];
sx q[1];
rz(-0.16709354) q[1];
rz(-pi) q[2];
rz(2.0031018) q[3];
sx q[3];
rz(-2.6902373) q[3];
sx q[3];
rz(-0.26717523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1294641) q[2];
sx q[2];
rz(-1.1295372) q[2];
sx q[2];
rz(0.10936347) q[2];
rz(0.28956413) q[3];
sx q[3];
rz(-1.3907631) q[3];
sx q[3];
rz(-0.65350986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7107723) q[0];
sx q[0];
rz(-1.9662974) q[0];
sx q[0];
rz(-1.9868504) q[0];
rz(-0.42354241) q[1];
sx q[1];
rz(-1.7083028) q[1];
sx q[1];
rz(-2.4636041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9311614) q[0];
sx q[0];
rz(-0.95936459) q[0];
sx q[0];
rz(1.0125041) q[0];
x q[1];
rz(-1.343849) q[2];
sx q[2];
rz(-1.316202) q[2];
sx q[2];
rz(-0.4190906) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0667693) q[1];
sx q[1];
rz(-0.69157234) q[1];
sx q[1];
rz(-2.2645386) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20991169) q[3];
sx q[3];
rz(-2.5394507) q[3];
sx q[3];
rz(-1.1478375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28422156) q[2];
sx q[2];
rz(-2.325433) q[2];
sx q[2];
rz(2.5428499) q[2];
rz(1.5395928) q[3];
sx q[3];
rz(-1.9339823) q[3];
sx q[3];
rz(-2.1046751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4620062) q[0];
sx q[0];
rz(-2.7988667) q[0];
sx q[0];
rz(-1.3890247) q[0];
rz(3.1287843) q[1];
sx q[1];
rz(-0.3777596) q[1];
sx q[1];
rz(-1.6526745) q[1];
rz(-1.7187865) q[2];
sx q[2];
rz(-0.64633804) q[2];
sx q[2];
rz(2.9895724) q[2];
rz(-0.0010796428) q[3];
sx q[3];
rz(-2.1297364) q[3];
sx q[3];
rz(-0.10441953) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
