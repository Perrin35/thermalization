OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.05623090267181) q[0];
sx q[0];
rz(3.41688543756539) q[0];
sx q[0];
rz(9.36981535925671) q[0];
rz(0.742052495479584) q[1];
sx q[1];
rz(3.25221780140931) q[1];
sx q[1];
rz(10.1799588560979) q[1];
cx q[1],q[0];
rz(0.395121455192566) q[0];
sx q[0];
rz(3.95693478186662) q[0];
sx q[0];
rz(9.4035042527984) q[0];
rz(0.490253835916519) q[2];
sx q[2];
rz(4.09614917834336) q[2];
sx q[2];
rz(10.1736035108487) q[2];
cx q[2],q[1];
rz(1.28581893444061) q[1];
sx q[1];
rz(3.57276728947694) q[1];
sx q[1];
rz(9.41533177103057) q[1];
rz(1.82014358043671) q[3];
sx q[3];
rz(3.7488450725847) q[3];
sx q[3];
rz(9.65676917730972) q[3];
cx q[3],q[2];
rz(1.55219686031342) q[2];
sx q[2];
rz(4.29751363595063) q[2];
sx q[2];
rz(10.8249400615613) q[2];
rz(0.333411484956741) q[3];
sx q[3];
rz(3.78940752347047) q[3];
sx q[3];
rz(9.92223337887927) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.196315035223961) q[0];
sx q[0];
rz(2.54402867157991) q[0];
sx q[0];
rz(9.51826076804801) q[0];
rz(0.641583323478699) q[1];
sx q[1];
rz(3.58110588987405) q[1];
sx q[1];
rz(9.61693393289253) q[1];
cx q[1],q[0];
rz(1.29633247852325) q[0];
sx q[0];
rz(3.98284527857835) q[0];
sx q[0];
rz(9.45838623716637) q[0];
rz(0.208835229277611) q[2];
sx q[2];
rz(4.50307658513124) q[2];
sx q[2];
rz(9.58766241966888) q[2];
cx q[2],q[1];
rz(-0.836289942264557) q[1];
sx q[1];
rz(3.94149670203263) q[1];
sx q[1];
rz(11.6453635454099) q[1];
rz(-0.00475625786930323) q[3];
sx q[3];
rz(5.03638246853883) q[3];
sx q[3];
rz(9.47106688319846) q[3];
cx q[3],q[2];
rz(0.415003478527069) q[2];
sx q[2];
rz(3.78724613984162) q[2];
sx q[2];
rz(10.4307000398557) q[2];
rz(0.0350602194666862) q[3];
sx q[3];
rz(3.70220032532746) q[3];
sx q[3];
rz(10.184643483154) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.189089179039) q[0];
sx q[0];
rz(4.94223478634889) q[0];
sx q[0];
rz(9.61103882490798) q[0];
rz(0.261789053678513) q[1];
sx q[1];
rz(4.34773114522035) q[1];
sx q[1];
rz(8.9153902888219) q[1];
cx q[1],q[0];
rz(0.875686883926392) q[0];
sx q[0];
rz(3.15799158637459) q[0];
sx q[0];
rz(9.4333351612012) q[0];
rz(1.2234753370285) q[2];
sx q[2];
rz(4.06463042100007) q[2];
sx q[2];
rz(7.5743093252103) q[2];
cx q[2],q[1];
rz(0.898159742355347) q[1];
sx q[1];
rz(3.75603440602357) q[1];
sx q[1];
rz(10.8168488502423) q[1];
rz(0.842864274978638) q[3];
sx q[3];
rz(5.21173659165437) q[3];
sx q[3];
rz(9.55019492506191) q[3];
cx q[3],q[2];
rz(0.190532609820366) q[2];
sx q[2];
rz(3.81579265196855) q[2];
sx q[2];
rz(9.72257489561244) q[2];
rz(0.315362453460693) q[3];
sx q[3];
rz(3.44999242027337) q[3];
sx q[3];
rz(11.1657697915952) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0676588714122772) q[0];
sx q[0];
rz(3.58185523946817) q[0];
sx q[0];
rz(9.42585607957236) q[0];
rz(-0.0623809695243835) q[1];
sx q[1];
rz(4.29310396512086) q[1];
sx q[1];
rz(11.07114813327) q[1];
cx q[1],q[0];
rz(0.304946810007095) q[0];
sx q[0];
rz(4.78157559235627) q[0];
sx q[0];
rz(9.92970982789203) q[0];
rz(4.09286785125732) q[2];
sx q[2];
rz(3.84529516299302) q[2];
sx q[2];
rz(8.32623813151523) q[2];
cx q[2],q[1];
rz(0.298410356044769) q[1];
sx q[1];
rz(3.94552645285661) q[1];
sx q[1];
rz(9.01421949862644) q[1];
rz(1.39950001239777) q[3];
sx q[3];
rz(5.34775153000886) q[3];
sx q[3];
rz(10.664441561691) q[3];
cx q[3],q[2];
rz(1.87947869300842) q[2];
sx q[2];
rz(2.41411903698976) q[2];
sx q[2];
rz(10.4811864852826) q[2];
rz(0.941117167472839) q[3];
sx q[3];
rz(4.42347029049928) q[3];
sx q[3];
rz(9.34697253852292) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.625174045562744) q[0];
sx q[0];
rz(3.84835526545579) q[0];
sx q[0];
rz(10.6112154483716) q[0];
rz(0.93900740146637) q[1];
sx q[1];
rz(3.86986580689485) q[1];
sx q[1];
rz(10.200441813461) q[1];
cx q[1],q[0];
rz(0.0308886580169201) q[0];
sx q[0];
rz(3.49912399251992) q[0];
sx q[0];
rz(9.97414193152591) q[0];
rz(0.553023457527161) q[2];
sx q[2];
rz(4.34081056912477) q[2];
sx q[2];
rz(9.34360404907867) q[2];
cx q[2],q[1];
rz(1.75550019741058) q[1];
sx q[1];
rz(3.73962745268876) q[1];
sx q[1];
rz(8.61382064818546) q[1];
rz(-0.0449599660933018) q[3];
sx q[3];
rz(4.49237695534761) q[3];
sx q[3];
rz(10.2507014632146) q[3];
cx q[3],q[2];
rz(-0.0808072164654732) q[2];
sx q[2];
rz(3.74413362343843) q[2];
sx q[2];
rz(10.9382422923963) q[2];
rz(1.3249409198761) q[3];
sx q[3];
rz(2.50382682879502) q[3];
sx q[3];
rz(9.54824662803813) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0476646460592747) q[0];
sx q[0];
rz(4.14911154110963) q[0];
sx q[0];
rz(10.0602112770002) q[0];
rz(-1.11349976062775) q[1];
sx q[1];
rz(3.57734647591645) q[1];
sx q[1];
rz(9.81916344761058) q[1];
cx q[1],q[0];
rz(-0.522882342338562) q[0];
sx q[0];
rz(3.25636457850272) q[0];
sx q[0];
rz(10.2563249826352) q[0];
rz(1.9362108707428) q[2];
sx q[2];
rz(3.63966042001779) q[2];
sx q[2];
rz(8.65145347117587) q[2];
cx q[2],q[1];
rz(0.043076153844595) q[1];
sx q[1];
rz(4.35788455803926) q[1];
sx q[1];
rz(10.5160860776822) q[1];
rz(0.907690584659576) q[3];
sx q[3];
rz(4.61293974717195) q[3];
sx q[3];
rz(9.58060838877364) q[3];
cx q[3],q[2];
rz(0.340279430150986) q[2];
sx q[2];
rz(4.49584391911561) q[2];
sx q[2];
rz(9.67042761146232) q[2];
rz(0.736860156059265) q[3];
sx q[3];
rz(4.10131749709184) q[3];
sx q[3];
rz(10.0148963689725) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0678294748067856) q[0];
sx q[0];
rz(4.98637917836244) q[0];
sx q[0];
rz(10.2762141585271) q[0];
rz(0.403594523668289) q[1];
sx q[1];
rz(3.74945566256578) q[1];
sx q[1];
rz(9.28458424507781) q[1];
cx q[1],q[0];
rz(0.644634127616882) q[0];
sx q[0];
rz(4.347929032641) q[0];
sx q[0];
rz(10.7551343202512) q[0];
rz(-1.0561101436615) q[2];
sx q[2];
rz(3.36544262071187) q[2];
sx q[2];
rz(11.4374630212705) q[2];
cx q[2],q[1];
rz(1.27800095081329) q[1];
sx q[1];
rz(3.88872710068757) q[1];
sx q[1];
rz(9.87001604437038) q[1];
rz(0.987366557121277) q[3];
sx q[3];
rz(4.32434585888917) q[3];
sx q[3];
rz(11.5311987161557) q[3];
cx q[3],q[2];
rz(0.957895278930664) q[2];
sx q[2];
rz(4.01686731179292) q[2];
sx q[2];
rz(9.62747000753089) q[2];
rz(1.96285879611969) q[3];
sx q[3];
rz(3.4027425964647) q[3];
sx q[3];
rz(8.32835195063754) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.751007676124573) q[0];
sx q[0];
rz(3.17706925992901) q[0];
sx q[0];
rz(10.6246506929319) q[0];
rz(1.22127819061279) q[1];
sx q[1];
rz(2.58836356003816) q[1];
sx q[1];
rz(7.94649229048892) q[1];
cx q[1],q[0];
rz(1.28658986091614) q[0];
sx q[0];
rz(3.15373942454393) q[0];
sx q[0];
rz(10.6323524475019) q[0];
rz(-0.488566845655441) q[2];
sx q[2];
rz(4.37395265896852) q[2];
sx q[2];
rz(10.6158962011258) q[2];
cx q[2],q[1];
rz(0.850466132164001) q[1];
sx q[1];
rz(4.15500298340852) q[1];
sx q[1];
rz(9.20961700975105) q[1];
rz(-1.0411365032196) q[3];
sx q[3];
rz(3.47222620447213) q[3];
sx q[3];
rz(9.66208178400203) q[3];
cx q[3],q[2];
rz(-0.0259098019450903) q[2];
sx q[2];
rz(2.21757105191285) q[2];
sx q[2];
rz(9.68682882785007) q[2];
rz(0.154843643307686) q[3];
sx q[3];
rz(3.34539242287213) q[3];
sx q[3];
rz(9.38504499792262) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.24209751188755) q[0];
sx q[0];
rz(3.84083053668077) q[0];
sx q[0];
rz(9.98748681544467) q[0];
rz(1.57446765899658) q[1];
sx q[1];
rz(4.4270238002115) q[1];
sx q[1];
rz(10.2345289945523) q[1];
cx q[1],q[0];
rz(0.278036773204803) q[0];
sx q[0];
rz(3.56780994136865) q[0];
sx q[0];
rz(9.99300137757465) q[0];
rz(0.631671965122223) q[2];
sx q[2];
rz(4.15847531159455) q[2];
sx q[2];
rz(9.88782543539211) q[2];
cx q[2],q[1];
rz(-0.875314712524414) q[1];
sx q[1];
rz(3.6326728781038) q[1];
sx q[1];
rz(11.8572642564695) q[1];
rz(0.169500201940536) q[3];
sx q[3];
rz(3.79825052817399) q[3];
sx q[3];
rz(9.8781693637292) q[3];
cx q[3],q[2];
rz(0.338538199663162) q[2];
sx q[2];
rz(2.70923167665536) q[2];
sx q[2];
rz(10.836245751373) q[2];
rz(0.15905299782753) q[3];
sx q[3];
rz(4.40109959443147) q[3];
sx q[3];
rz(10.0892310500066) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.148010462522507) q[0];
sx q[0];
rz(2.9382885565334) q[0];
sx q[0];
rz(10.0547263979833) q[0];
rz(0.826076209545135) q[1];
sx q[1];
rz(3.59495470126206) q[1];
sx q[1];
rz(10.425316309921) q[1];
cx q[1],q[0];
rz(0.31622177362442) q[0];
sx q[0];
rz(4.66620686848695) q[0];
sx q[0];
rz(10.4530070781629) q[0];
rz(-1.10648083686829) q[2];
sx q[2];
rz(4.45261135895784) q[2];
sx q[2];
rz(11.8360533475797) q[2];
cx q[2],q[1];
rz(0.368808209896088) q[1];
sx q[1];
rz(3.93064013321931) q[1];
sx q[1];
rz(10.3943575382154) q[1];
rz(0.482021898031235) q[3];
sx q[3];
rz(4.06159737904603) q[3];
sx q[3];
rz(9.61930482684776) q[3];
cx q[3],q[2];
rz(0.916538119316101) q[2];
sx q[2];
rz(3.44480163057382) q[2];
sx q[2];
rz(10.4851104974668) q[2];
rz(0.669250845909119) q[3];
sx q[3];
rz(3.80689415534074) q[3];
sx q[3];
rz(10.1917993783872) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.457220613956451) q[0];
sx q[0];
rz(3.88296380837495) q[0];
sx q[0];
rz(10.1610430836599) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.951238214969635) q[1];
sx q[1];
rz(3.57499754627282) q[1];
sx q[1];
rz(9.08666623233959) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.204516842961311) q[2];
sx q[2];
rz(4.1585295518213) q[2];
sx q[2];
rz(10.1392331480901) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.942237198352814) q[3];
sx q[3];
rz(3.88820293744142) q[3];
sx q[3];
rz(9.14456916450664) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
