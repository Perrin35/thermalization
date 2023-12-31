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
rz(4.9217304) q[0];
sx q[0];
rz(11.187727) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987797) q[0];
sx q[0];
rz(-3.0516041) q[0];
sx q[0];
rz(-2.9773657) q[0];
x q[1];
rz(1.3583463) q[2];
sx q[2];
rz(-2.6531086) q[2];
sx q[2];
rz(0.24220315) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1554943) q[1];
sx q[1];
rz(-1.3487779) q[1];
sx q[1];
rz(-0.51198952) q[1];
rz(-0.095590683) q[3];
sx q[3];
rz(-0.88926892) q[3];
sx q[3];
rz(2.3694627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(-0.84428865) q[2];
rz(0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(-2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.9553604) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037558) q[0];
sx q[0];
rz(-3.0825966) q[0];
sx q[0];
rz(0.32365139) q[0];
rz(0.19940168) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(2.8859438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4071719) q[1];
sx q[1];
rz(-0.91033519) q[1];
sx q[1];
rz(0.15863005) q[1];
x q[2];
rz(1.7632742) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(-1.457085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3804669) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(-2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(-2.321373) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263111) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(-0.88322722) q[0];
rz(-pi) q[1];
rz(2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-2.6205274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0318451) q[1];
sx q[1];
rz(-2.577563) q[1];
sx q[1];
rz(2.2721223) q[1];
x q[2];
rz(-1.3385653) q[3];
sx q[3];
rz(-2.665822) q[3];
sx q[3];
rz(2.8867718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(-0.43831929) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.8211676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5116918) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(3.1303309) q[0];
x q[1];
rz(2.690372) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(0.95552432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1530694) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(-2.4371229) q[1];
x q[2];
rz(0.60453316) q[3];
sx q[3];
rz(-0.88450888) q[3];
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
rz(-0.92869174) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(-1.226549) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(1.3006166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9601701) q[0];
sx q[0];
rz(-1.2167861) q[0];
sx q[0];
rz(1.5439073) q[0];
rz(-pi) q[1];
rz(-2.956203) q[2];
sx q[2];
rz(-0.37237793) q[2];
sx q[2];
rz(0.34639726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.865766) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(-2.5715716) q[1];
rz(0.50118581) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(2.3230769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(2.9495083) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-0.011750301) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5884429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4295411) q[0];
sx q[0];
rz(-0.49563227) q[0];
sx q[0];
rz(0.80722157) q[0];
rz(-pi) q[1];
rz(-1.2049098) q[2];
sx q[2];
rz(-1.591218) q[2];
sx q[2];
rz(2.6423955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1482684) q[1];
sx q[1];
rz(-2.4273708) q[1];
sx q[1];
rz(-3.0118914) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-0.67725956) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-2.1645434) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58986321) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(-0.59318869) q[0];
rz(-pi) q[1];
rz(-3.0095519) q[2];
sx q[2];
rz(-1.5706976) q[2];
sx q[2];
rz(1.6824818) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3358826) q[1];
sx q[1];
rz(-1.6478331) q[1];
sx q[1];
rz(-1.9105934) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2491751) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(-1.7082937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-2.1288669) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-0.69865984) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.2493856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.218924) q[0];
sx q[0];
rz(-1.3601174) q[0];
sx q[0];
rz(-1.6181437) q[0];
rz(0.84455372) q[2];
sx q[2];
rz(-2.484998) q[2];
sx q[2];
rz(0.086364634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57751209) q[1];
sx q[1];
rz(-3.1187594) q[1];
sx q[1];
rz(-0.24502416) q[1];
rz(-1.1864248) q[3];
sx q[3];
rz(-0.66925183) q[3];
sx q[3];
rz(2.3068908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.1784941) q[2];
rz(-1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(-1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(0.59757772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4353838) q[0];
sx q[0];
rz(-1.7058813) q[0];
sx q[0];
rz(-0.038890966) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55677982) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(2.5533887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55652009) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(-0.82533605) q[1];
rz(-pi) q[2];
rz(3.0258614) q[3];
sx q[3];
rz(-0.58905187) q[3];
sx q[3];
rz(0.91050402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.2333599) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(-1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(3.1179324) q[0];
rz(0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(2.4694209) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7260872) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(1.0938209) q[0];
x q[1];
rz(-2.6300738) q[2];
sx q[2];
rz(-1.5626972) q[2];
sx q[2];
rz(0.17009232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0553531) q[1];
sx q[1];
rz(-1.2330016) q[1];
sx q[1];
rz(2.4243381) q[1];
x q[2];
rz(2.7087595) q[3];
sx q[3];
rz(-1.7107309) q[3];
sx q[3];
rz(-0.39633745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-2.9818515) q[2];
sx q[2];
rz(-1.1764871) q[2];
sx q[2];
rz(2.7574678) q[2];
rz(-2.1914992) q[3];
sx q[3];
rz(-0.48662574) q[3];
sx q[3];
rz(1.1403198) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
