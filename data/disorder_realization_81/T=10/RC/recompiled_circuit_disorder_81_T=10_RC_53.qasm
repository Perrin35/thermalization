OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(4.6255914) q[1];
sx q[1];
rz(8.9738823) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8060018) q[0];
sx q[0];
rz(-1.5854892) q[0];
sx q[0];
rz(-3.0528085) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0300006) q[2];
sx q[2];
rz(-2.0473695) q[2];
sx q[2];
rz(-0.0026207844) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7887468) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-pi) q[2];
rz(1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(0.62108921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(2.297304) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(-0.94353765) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77957905) q[0];
sx q[0];
rz(-1.5148666) q[0];
sx q[0];
rz(1.5895784) q[0];
rz(-0.19940168) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(-2.8859438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1520878) q[1];
sx q[1];
rz(-0.67645914) q[1];
sx q[1];
rz(1.771404) q[1];
rz(-pi) q[2];
rz(-0.22963345) q[3];
sx q[3];
rz(-0.7080871) q[3];
sx q[3];
rz(-1.9830444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.9821232) q[2];
rz(0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(-0.82021964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504282) q[0];
sx q[0];
rz(-2.4463852) q[0];
sx q[0];
rz(-1.7466963) q[0];
rz(-2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-0.52106524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1097475) q[1];
sx q[1];
rz(-0.56402962) q[1];
sx q[1];
rz(-2.2721223) q[1];
rz(-pi) q[2];
rz(1.3385653) q[3];
sx q[3];
rz(-0.47577061) q[3];
sx q[3];
rz(-0.25482086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31072581) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(-0.91127515) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.320425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62990084) q[0];
sx q[0];
rz(-1.1612079) q[0];
sx q[0];
rz(3.1303309) q[0];
x q[1];
rz(2.690372) q[2];
sx q[2];
rz(-1.3845452) q[2];
sx q[2];
rz(2.1860683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79975407) q[1];
sx q[1];
rz(-0.8952039) q[1];
sx q[1];
rz(-1.9104596) q[1];
rz(-pi) q[2];
rz(2.5370595) q[3];
sx q[3];
rz(-2.2570838) q[3];
sx q[3];
rz(-1.7741007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(2.9648182) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3115561) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(1.3006166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375992) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(-0.07261891) q[0];
x q[1];
rz(0.36655764) q[2];
sx q[2];
rz(-1.5036811) q[2];
sx q[2];
rz(1.7442489) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0504426) q[1];
sx q[1];
rz(-2.1069408) q[1];
sx q[1];
rz(1.1863143) q[1];
x q[2];
rz(-0.50118581) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(0.81851573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(-0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5531497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1588622) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(1.1984675) q[0];
rz(-pi) q[1];
rz(3.1197238) q[2];
sx q[2];
rz(-1.936603) q[2];
sx q[2];
rz(-2.0621698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1482684) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(-0.12970129) q[1];
x q[2];
rz(3.0670777) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(-2.6889192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-1.1184568) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(2.4643331) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(2.1645434) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517294) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(0.59318869) q[0];
rz(-pi) q[1];
rz(3.1408429) q[2];
sx q[2];
rz(-0.13204083) q[2];
sx q[2];
rz(-0.11242871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44927412) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(1.3432137) q[1];
rz(-1.2491751) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(-1.4332989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(-2.1288669) q[2];
rz(1.9536473) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(-2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(0.69865984) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.8922071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218924) q[0];
sx q[0];
rz(-1.7814753) q[0];
sx q[0];
rz(-1.6181437) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0480568) q[2];
sx q[2];
rz(-1.9881696) q[2];
sx q[2];
rz(-1.0440895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74832143) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(-3.1194411) q[1];
rz(-0.93805712) q[3];
sx q[3];
rz(-1.3360099) q[3];
sx q[3];
rz(0.42890047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(-1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15411988) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(-1.2921635) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0963247) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(1.8906821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9012007) q[1];
sx q[1];
rz(-1.8933834) q[1];
sx q[1];
rz(-0.31369536) q[1];
rz(-0.11573128) q[3];
sx q[3];
rz(-0.58905187) q[3];
sx q[3];
rz(-2.2310886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.2333599) q[2];
rz(0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-0.023660252) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(-0.67217174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32165747) q[0];
sx q[0];
rz(-1.5655087) q[0];
sx q[0];
rz(1.5810285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6300738) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(-0.17009232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2939261) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(-0.49077175) q[1];
rz(-pi) q[2];
rz(0.43283312) q[3];
sx q[3];
rz(-1.4308617) q[3];
sx q[3];
rz(2.7452552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621915) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-0.72369408) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(-1.9359246) q[2];
sx q[2];
rz(-0.42386133) q[2];
sx q[2];
rz(-0.78122666) q[2];
rz(1.1643812) q[3];
sx q[3];
rz(-1.2953399) q[3];
sx q[3];
rz(-3.0084707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];