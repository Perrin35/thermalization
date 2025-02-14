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
rz(0.451568901538849) q[0];
sx q[0];
rz(1.418192537623) q[0];
sx q[0];
rz(12.1374921560208) q[0];
rz(-3.28473258018494) q[1];
sx q[1];
rz(-1.09094556967681) q[1];
sx q[1];
rz(9.34468909203216) q[1];
cx q[1],q[0];
rz(-1.03196656703949) q[0];
sx q[0];
rz(4.63320759137208) q[0];
sx q[0];
rz(12.5530884027402) q[0];
rz(0.360882490873337) q[2];
sx q[2];
rz(4.61607280571992) q[2];
sx q[2];
rz(10.4173866271894) q[2];
cx q[2],q[1];
rz(6.73495721817017) q[1];
sx q[1];
rz(4.0233965237909) q[1];
sx q[1];
rz(13.7932982206266) q[1];
rz(3.30580043792725) q[3];
sx q[3];
rz(0.999421509104319) q[3];
sx q[3];
rz(8.64561281203433) q[3];
cx q[3],q[2];
rz(3.21637153625488) q[2];
sx q[2];
rz(5.23181596596772) q[2];
sx q[2];
rz(7.70522663592502) q[2];
rz(1.72857964038849) q[3];
sx q[3];
rz(1.54141095479066) q[3];
sx q[3];
rz(11.1181035995404) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.01953518390656) q[0];
sx q[0];
rz(1.39665654500062) q[0];
sx q[0];
rz(9.09120703338786) q[0];
rz(3.30752539634705) q[1];
sx q[1];
rz(6.6270000060373) q[1];
sx q[1];
rz(8.67267451285526) q[1];
cx q[1],q[0];
rz(0.0108757941052318) q[0];
sx q[0];
rz(4.11303308804566) q[0];
sx q[0];
rz(10.4457787036817) q[0];
rz(-1.78954589366913) q[2];
sx q[2];
rz(1.59652117093141) q[2];
sx q[2];
rz(15.778220152847) q[2];
cx q[2],q[1];
rz(-1.18114066123962) q[1];
sx q[1];
rz(4.60327580769593) q[1];
sx q[1];
rz(12.0383341073911) q[1];
rz(2.27769231796265) q[3];
sx q[3];
rz(3.48034733732278) q[3];
sx q[3];
rz(8.21364829539462) q[3];
cx q[3],q[2];
rz(-0.0779452621936798) q[2];
sx q[2];
rz(0.912194879847117) q[2];
sx q[2];
rz(13.223409152023) q[2];
rz(-0.71197235584259) q[3];
sx q[3];
rz(0.651980551081248) q[3];
sx q[3];
rz(10.1696514248769) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.70653581619263) q[0];
sx q[0];
rz(2.5223794897371) q[0];
sx q[0];
rz(10.4662277460019) q[0];
rz(2.76677465438843) q[1];
sx q[1];
rz(4.77980664570863) q[1];
sx q[1];
rz(12.0094218015592) q[1];
cx q[1],q[0];
rz(-0.14675635099411) q[0];
sx q[0];
rz(0.868491562204905) q[0];
sx q[0];
rz(13.0348529577176) q[0];
rz(0.15424732863903) q[2];
sx q[2];
rz(2.68523124058778) q[2];
sx q[2];
rz(10.3452211975972) q[2];
cx q[2],q[1];
rz(-1.53656888008118) q[1];
sx q[1];
rz(5.60615578492219) q[1];
sx q[1];
rz(10.8320278882901) q[1];
rz(-1.87761247158051) q[3];
sx q[3];
rz(0.674248846369334) q[3];
sx q[3];
rz(13.6698360204618) q[3];
cx q[3],q[2];
rz(0.53691428899765) q[2];
sx q[2];
rz(5.35891357262666) q[2];
sx q[2];
rz(12.1311924219052) q[2];
rz(3.66629719734192) q[3];
sx q[3];
rz(3.47514224250848) q[3];
sx q[3];
rz(6.14840528964206) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.101145766675472) q[0];
sx q[0];
rz(4.20596936543519) q[0];
sx q[0];
rz(10.9392077684323) q[0];
rz(2.09284043312073) q[1];
sx q[1];
rz(5.36043396790559) q[1];
sx q[1];
rz(11.1305939912717) q[1];
cx q[1],q[0];
rz(1.74833166599274) q[0];
sx q[0];
rz(2.9611640890413) q[0];
sx q[0];
rz(8.0951092004697) q[0];
rz(1.02748274803162) q[2];
sx q[2];
rz(1.26769533951814) q[2];
sx q[2];
rz(11.1064426660459) q[2];
cx q[2],q[1];
rz(-0.175711199641228) q[1];
sx q[1];
rz(3.66722026665742) q[1];
sx q[1];
rz(6.15532324313327) q[1];
rz(-2.51715612411499) q[3];
sx q[3];
rz(2.58899596531922) q[3];
sx q[3];
rz(13.275156235687) q[3];
cx q[3],q[2];
rz(3.54179573059082) q[2];
sx q[2];
rz(4.92533162434632) q[2];
sx q[2];
rz(12.5915271997373) q[2];
rz(1.48704206943512) q[3];
sx q[3];
rz(2.74367228348786) q[3];
sx q[3];
rz(8.60856959819003) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0954892784357071) q[0];
sx q[0];
rz(5.63366475899751) q[0];
sx q[0];
rz(12.4073920011441) q[0];
rz(-2.03646302223206) q[1];
sx q[1];
rz(2.41598937113816) q[1];
sx q[1];
rz(9.64908484219714) q[1];
cx q[1],q[0];
rz(-1.10625779628754) q[0];
sx q[0];
rz(2.43199369509751) q[0];
sx q[0];
rz(7.2105808019559) q[0];
rz(2.58134293556213) q[2];
sx q[2];
rz(1.5429826100641) q[2];
sx q[2];
rz(7.82913408278629) q[2];
cx q[2],q[1];
rz(-0.0807066559791565) q[1];
sx q[1];
rz(4.59703675110871) q[1];
sx q[1];
rz(8.89875141381427) q[1];
rz(1.37782692909241) q[3];
sx q[3];
rz(5.69747653801972) q[3];
sx q[3];
rz(11.5095293283383) q[3];
cx q[3],q[2];
rz(0.205324560403824) q[2];
sx q[2];
rz(7.64536491234834) q[2];
sx q[2];
rz(8.67617604731723) q[2];
rz(3.25207471847534) q[3];
sx q[3];
rz(4.36910036404664) q[3];
sx q[3];
rz(12.9777953386228) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.16141378879547) q[0];
sx q[0];
rz(5.91465154488618) q[0];
sx q[0];
rz(10.1550990104596) q[0];
rz(1.6397362947464) q[1];
sx q[1];
rz(4.5316917022043) q[1];
sx q[1];
rz(9.65912643670245) q[1];
cx q[1],q[0];
rz(0.653559505939484) q[0];
sx q[0];
rz(3.97598752577836) q[0];
sx q[0];
rz(11.60919044017) q[0];
rz(-2.16970562934875) q[2];
sx q[2];
rz(5.2674094756418) q[2];
sx q[2];
rz(10.3235211729924) q[2];
cx q[2],q[1];
rz(-1.15336561203003) q[1];
sx q[1];
rz(4.54980638821656) q[1];
sx q[1];
rz(6.23500893115207) q[1];
rz(2.73208785057068) q[3];
sx q[3];
rz(6.43610111077363) q[3];
sx q[3];
rz(11.9730717897336) q[3];
cx q[3],q[2];
rz(6.47649097442627) q[2];
sx q[2];
rz(-0.813934652013234) q[2];
sx q[2];
rz(10.706143474571) q[2];
rz(2.53656005859375) q[3];
sx q[3];
rz(2.04219594796235) q[3];
sx q[3];
rz(10.6572928190152) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.1518816947937) q[0];
sx q[0];
rz(1.87937703927095) q[0];
sx q[0];
rz(10.0378189444463) q[0];
rz(0.680606544017792) q[1];
sx q[1];
rz(5.17203346093232) q[1];
sx q[1];
rz(8.0994615316312) q[1];
cx q[1],q[0];
rz(2.31917357444763) q[0];
sx q[0];
rz(-1.02860101859038) q[0];
sx q[0];
rz(9.1117172896783) q[0];
rz(1.20300137996674) q[2];
sx q[2];
rz(6.05524006684358) q[2];
sx q[2];
rz(11.1190334319989) q[2];
cx q[2],q[1];
rz(0.0541067086160183) q[1];
sx q[1];
rz(4.1517985184961) q[1];
sx q[1];
rz(9.27521949111625) q[1];
rz(0.988729000091553) q[3];
sx q[3];
rz(4.37211600144441) q[3];
sx q[3];
rz(13.2162768602292) q[3];
cx q[3],q[2];
rz(0.947205424308777) q[2];
sx q[2];
rz(1.94018653233583) q[2];
sx q[2];
rz(7.09124467372104) q[2];
rz(3.78747797012329) q[3];
sx q[3];
rz(3.40951258142526) q[3];
sx q[3];
rz(6.58655450343295) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.5757669210434) q[0];
sx q[0];
rz(3.86820265849168) q[0];
sx q[0];
rz(6.72595355509921) q[0];
rz(5.69113874435425) q[1];
sx q[1];
rz(5.56180110772187) q[1];
sx q[1];
rz(9.60046645104095) q[1];
cx q[1],q[0];
rz(0.739787757396698) q[0];
sx q[0];
rz(1.29538384278352) q[0];
sx q[0];
rz(11.3642486095349) q[0];
rz(2.02234792709351) q[2];
sx q[2];
rz(4.67651871045167) q[2];
sx q[2];
rz(13.2269394159238) q[2];
cx q[2],q[1];
rz(4.09548139572144) q[1];
sx q[1];
rz(-3.7842172066397) q[1];
sx q[1];
rz(10.9235371112744) q[1];
rz(-2.06667041778564) q[3];
sx q[3];
rz(2.12258330185945) q[3];
sx q[3];
rz(10.8896358966748) q[3];
cx q[3],q[2];
rz(0.597779333591461) q[2];
sx q[2];
rz(2.39763489563996) q[2];
sx q[2];
rz(5.37607667445346) q[2];
rz(-0.907265663146973) q[3];
sx q[3];
rz(4.88884952862794) q[3];
sx q[3];
rz(9.54246917962238) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.18427813053131) q[0];
sx q[0];
rz(4.04342720110948) q[0];
sx q[0];
rz(9.61482313870593) q[0];
rz(4.70050716400146) q[1];
sx q[1];
rz(4.68853667576844) q[1];
sx q[1];
rz(11.2591680049817) q[1];
cx q[1],q[0];
rz(-0.350185394287109) q[0];
sx q[0];
rz(3.53505063255365) q[0];
sx q[0];
rz(12.948821759216) q[0];
rz(1.90938782691956) q[2];
sx q[2];
rz(4.04110011656816) q[2];
sx q[2];
rz(6.29582045077487) q[2];
cx q[2],q[1];
rz(4.29663991928101) q[1];
sx q[1];
rz(1.53801706631715) q[1];
sx q[1];
rz(9.31704025565788) q[1];
rz(-0.115477412939072) q[3];
sx q[3];
rz(4.28040495713288) q[3];
sx q[3];
rz(9.83006504773303) q[3];
cx q[3],q[2];
rz(3.21004581451416) q[2];
sx q[2];
rz(5.19606510003144) q[2];
sx q[2];
rz(12.3652188539426) q[2];
rz(1.02260947227478) q[3];
sx q[3];
rz(5.60066572030122) q[3];
sx q[3];
rz(11.0267796277921) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.967118382453918) q[0];
sx q[0];
rz(5.23723903496797) q[0];
sx q[0];
rz(11.6899046659391) q[0];
rz(3.76409006118774) q[1];
sx q[1];
rz(6.49496236641938) q[1];
sx q[1];
rz(12.2236566305081) q[1];
cx q[1],q[0];
rz(0.087992824614048) q[0];
sx q[0];
rz(2.19810328085954) q[0];
sx q[0];
rz(4.79377934931918) q[0];
rz(2.80439376831055) q[2];
sx q[2];
rz(3.49683720071847) q[2];
sx q[2];
rz(11.5955564737241) q[2];
cx q[2],q[1];
rz(0.491492807865143) q[1];
sx q[1];
rz(2.01801851590211) q[1];
sx q[1];
rz(11.6866948366086) q[1];
rz(1.18311715126038) q[3];
sx q[3];
rz(4.94693103631074) q[3];
sx q[3];
rz(12.78253028392) q[3];
cx q[3],q[2];
rz(-1.18913018703461) q[2];
sx q[2];
rz(4.5834870656305) q[2];
sx q[2];
rz(8.98098898529216) q[2];
rz(5.16443061828613) q[3];
sx q[3];
rz(1.54645624955232) q[3];
sx q[3];
rz(8.87094763516589) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.58862805366516) q[0];
sx q[0];
rz(3.08969670732553) q[0];
sx q[0];
rz(12.8320858240049) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.705950319766998) q[1];
sx q[1];
rz(6.00395050843293) q[1];
sx q[1];
rz(10.3719424962918) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-3.90883898735046) q[2];
sx q[2];
rz(3.01363483269746) q[2];
sx q[2];
rz(12.0899669885556) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.48349797725677) q[3];
sx q[3];
rz(2.56364205678041) q[3];
sx q[3];
rz(10.7156289577405) q[3];
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
