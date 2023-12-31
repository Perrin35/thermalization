OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(1.3786432) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(0.4508957) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(-0.16422693) q[0];
x q[1];
rz(-3.0300006) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(-0.0026207844) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3528459) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(-2.7098141) q[1];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(-1.9553604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3620136) q[0];
sx q[0];
rz(-1.626726) q[0];
sx q[0];
rz(1.5895784) q[0];
x q[1];
rz(-0.29859121) q[2];
sx q[2];
rz(-0.208374) q[2];
sx q[2];
rz(1.5339472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9895049) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(-1.771404) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69497377) q[3];
sx q[3];
rz(-1.7193828) q[3];
sx q[3];
rz(-2.9050764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(-0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(-2.321373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0778753) q[0];
sx q[0];
rz(-0.88839196) q[0];
sx q[0];
rz(-2.9966485) q[0];
rz(1.0772466) q[2];
sx q[2];
rz(-1.764467) q[2];
sx q[2];
rz(2.1950302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1097475) q[1];
sx q[1];
rz(-0.56402962) q[1];
sx q[1];
rz(-2.2721223) q[1];
x q[2];
rz(3.0235602) q[3];
sx q[3];
rz(-2.0327838) q[3];
sx q[3];
rz(0.0052009728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(2.9860949) q[3];
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
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(-1.320425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5399649) q[0];
sx q[0];
rz(-0.40973445) q[0];
sx q[0];
rz(1.5967303) q[0];
rz(-0.45122066) q[2];
sx q[2];
rz(-1.3845452) q[2];
sx q[2];
rz(2.1860683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1530694) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(0.70446976) q[1];
rz(-pi) q[2];
rz(2.5370595) q[3];
sx q[3];
rz(-2.2570838) q[3];
sx q[3];
rz(1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(-2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(2.2763021) q[0];
rz(1.226549) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.8409761) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9601701) q[0];
sx q[0];
rz(-1.9248065) q[0];
sx q[0];
rz(1.5439073) q[0];
x q[1];
rz(0.36655764) q[2];
sx q[2];
rz(-1.5036811) q[2];
sx q[2];
rz(1.7442489) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2758267) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(-2.5715716) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50118581) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(2.3230769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(2.664393) q[2];
rz(-0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-0.011750301) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6019183) q[0];
sx q[0];
rz(-1.2201021) q[0];
sx q[0];
rz(-2.7838216) q[0];
rz(-1.5137709) q[2];
sx q[2];
rz(-0.36643039) q[2];
sx q[2];
rz(2.1232405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8223871) q[1];
sx q[1];
rz(-0.86383312) q[1];
sx q[1];
rz(-1.6824526) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0670777) q[3];
sx q[3];
rz(-0.86935589) q[3];
sx q[3];
rz(-0.45267347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(1.2825512) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(1.1184568) q[3];
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
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1099694) q[0];
sx q[0];
rz(-0.97951802) q[0];
sx q[0];
rz(-1.4800319) q[0];
rz(-3.1408429) q[2];
sx q[2];
rz(-0.13204083) q[2];
sx q[2];
rz(-3.0291639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44927412) q[1];
sx q[1];
rz(-0.34808967) q[1];
sx q[1];
rz(-1.7983789) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3903923) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(1.8922071) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69985336) q[0];
sx q[0];
rz(-2.925736) q[0];
sx q[0];
rz(0.21780832) q[0];
rz(-pi) q[1];
rz(0.84455372) q[2];
sx q[2];
rz(-0.65659467) q[2];
sx q[2];
rz(-0.086364634) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57751209) q[1];
sx q[1];
rz(-3.1187594) q[1];
sx q[1];
rz(-0.24502416) q[1];
rz(-0.93805712) q[3];
sx q[3];
rz(-1.8055827) q[3];
sx q[3];
rz(-0.42890047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(-1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4998528) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-0.59757772) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874728) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(1.8494291) q[0];
x q[1];
rz(-2.5848128) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(2.5533887) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55652009) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(-2.3162566) q[1];
rz(-pi) q[2];
rz(-3.0258614) q[3];
sx q[3];
rz(-0.58905187) q[3];
sx q[3];
rz(-0.91050402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(-1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(-0.67217174) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7260872) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(-2.0477717) q[0];
rz(-3.1250481) q[2];
sx q[2];
rz(-0.51157727) q[2];
sx q[2];
rz(-1.3862773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8476665) q[1];
sx q[1];
rz(-0.77984174) q[1];
sx q[1];
rz(0.49077175) q[1];
rz(-pi) q[2];
x q[2];
rz(1.724733) q[3];
sx q[3];
rz(-1.9991176) q[3];
sx q[3];
rz(1.1101013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5621915) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.205668) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(2.1914992) q[3];
sx q[3];
rz(-2.6549669) q[3];
sx q[3];
rz(-2.0012729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
