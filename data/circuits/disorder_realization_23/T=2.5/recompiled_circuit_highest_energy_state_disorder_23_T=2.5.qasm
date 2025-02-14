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
rz(1.79534018039703) q[0];
sx q[0];
rz(4.59023085434968) q[0];
sx q[0];
rz(11.4147051334302) q[0];
rz(1.08980739116669) q[1];
sx q[1];
rz(1.66482189496095) q[1];
sx q[1];
rz(9.27381684481307) q[1];
cx q[1],q[0];
rz(3.52447056770325) q[0];
sx q[0];
rz(2.84235137899453) q[0];
sx q[0];
rz(6.10422465800449) q[0];
rz(-0.086920291185379) q[2];
sx q[2];
rz(1.88332322438294) q[2];
sx q[2];
rz(9.45058236121341) q[2];
cx q[2],q[1];
rz(-0.261502951383591) q[1];
sx q[1];
rz(3.12633985572542) q[1];
sx q[1];
rz(9.18061397074863) q[1];
rz(1.47411513328552) q[3];
sx q[3];
rz(1.61143174965913) q[3];
sx q[3];
rz(10.796480870239) q[3];
cx q[3],q[2];
rz(1.02074086666107) q[2];
sx q[2];
rz(3.11817017582292) q[2];
sx q[2];
rz(9.91607097386524) q[2];
rz(1.32257199287415) q[3];
sx q[3];
rz(4.74880019028718) q[3];
sx q[3];
rz(8.09692845343753) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.781119287014008) q[0];
sx q[0];
rz(0.554471166925975) q[0];
sx q[0];
rz(8.72517523764774) q[0];
rz(-1.55050015449524) q[1];
sx q[1];
rz(3.63643041451509) q[1];
sx q[1];
rz(12.7362327337186) q[1];
cx q[1],q[0];
rz(1.24015974998474) q[0];
sx q[0];
rz(3.22831911792094) q[0];
sx q[0];
rz(7.39282414912387) q[0];
rz(3.00699400901794) q[2];
sx q[2];
rz(1.87301174004609) q[2];
sx q[2];
rz(8.8394043803136) q[2];
cx q[2],q[1];
rz(-0.583860337734222) q[1];
sx q[1];
rz(2.52164444525773) q[1];
sx q[1];
rz(11.9220475912015) q[1];
rz(-1.49165320396423) q[3];
sx q[3];
rz(1.11496606667573) q[3];
sx q[3];
rz(10.8223998308103) q[3];
cx q[3],q[2];
rz(0.262131720781326) q[2];
sx q[2];
rz(3.50395822723443) q[2];
sx q[2];
rz(8.00786910056278) q[2];
rz(2.05876851081848) q[3];
sx q[3];
rz(8.53743043740327) q[3];
sx q[3];
rz(8.59757558106586) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.56034922599792) q[0];
sx q[0];
rz(4.06169048150117) q[0];
sx q[0];
rz(11.0071186780851) q[0];
rz(-0.70445591211319) q[1];
sx q[1];
rz(1.99439409573609) q[1];
sx q[1];
rz(6.81998965739414) q[1];
cx q[1],q[0];
rz(3.19422578811646) q[0];
sx q[0];
rz(2.45287022193009) q[0];
sx q[0];
rz(10.9227092027585) q[0];
rz(-3.10498905181885) q[2];
sx q[2];
rz(4.62219444115693) q[2];
sx q[2];
rz(3.62834498881503) q[2];
cx q[2],q[1];
rz(2.06728625297546) q[1];
sx q[1];
rz(3.32664036949212) q[1];
sx q[1];
rz(7.23547885417148) q[1];
rz(-0.290729016065598) q[3];
sx q[3];
rz(-0.717462627095632) q[3];
sx q[3];
rz(7.19943783282443) q[3];
cx q[3],q[2];
rz(-2.39293599128723) q[2];
sx q[2];
rz(4.94586041768128) q[2];
sx q[2];
rz(11.7875249147336) q[2];
rz(-3.05741620063782) q[3];
sx q[3];
rz(4.01096233923966) q[3];
sx q[3];
rz(10.8275933027188) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.244514301419258) q[0];
sx q[0];
rz(4.06264403660829) q[0];
sx q[0];
rz(9.89543659090205) q[0];
rz(1.45875191688538) q[1];
sx q[1];
rz(3.146440735674) q[1];
sx q[1];
rz(10.2895392537038) q[1];
cx q[1],q[0];
rz(-1.70474588871002) q[0];
sx q[0];
rz(3.65946665604646) q[0];
sx q[0];
rz(11.1130803584973) q[0];
rz(-0.132238447666168) q[2];
sx q[2];
rz(1.79004278977449) q[2];
sx q[2];
rz(8.17010018824741) q[2];
cx q[2],q[1];
rz(-2.383385181427) q[1];
sx q[1];
rz(4.64666632016236) q[1];
sx q[1];
rz(14.1882490873258) q[1];
rz(1.0799013376236) q[3];
sx q[3];
rz(3.73821172316606) q[3];
sx q[3];
rz(6.48682305812045) q[3];
cx q[3],q[2];
rz(0.0278157591819763) q[2];
sx q[2];
rz(4.00622645218904) q[2];
sx q[2];
rz(10.5817321300428) q[2];
rz(-0.35465869307518) q[3];
sx q[3];
rz(5.04704001744325) q[3];
sx q[3];
rz(9.53531097470924) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.81195867061615) q[0];
sx q[0];
rz(5.34760800202424) q[0];
sx q[0];
rz(8.89114186762973) q[0];
rz(4.35528326034546) q[1];
sx q[1];
rz(3.11849892896647) q[1];
sx q[1];
rz(10.6060485601346) q[1];
cx q[1],q[0];
rz(-1.26425194740295) q[0];
sx q[0];
rz(2.97216572065885) q[0];
sx q[0];
rz(10.7197599172513) q[0];
rz(-2.22186279296875) q[2];
sx q[2];
rz(5.050929697352) q[2];
sx q[2];
rz(12.3719427347104) q[2];
cx q[2],q[1];
rz(-2.27948403358459) q[1];
sx q[1];
rz(4.06753191550309) q[1];
sx q[1];
rz(10.7996676921765) q[1];
rz(1.4328670501709) q[3];
sx q[3];
rz(4.21188071568544) q[3];
sx q[3];
rz(13.1660172700803) q[3];
cx q[3],q[2];
rz(0.108267620205879) q[2];
sx q[2];
rz(5.76755443413789) q[2];
sx q[2];
rz(7.21151659487888) q[2];
rz(2.86322736740112) q[3];
sx q[3];
rz(1.31855777104432) q[3];
sx q[3];
rz(9.4934251293461) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.53348970413208) q[0];
sx q[0];
rz(2.61874959071214) q[0];
sx q[0];
rz(9.23173703848525) q[0];
rz(2.46423840522766) q[1];
sx q[1];
rz(3.11321279046173) q[1];
sx q[1];
rz(7.79321393965884) q[1];
cx q[1],q[0];
rz(2.22644352912903) q[0];
sx q[0];
rz(6.1409102996164) q[0];
sx q[0];
rz(14.2789053678434) q[0];
rz(4.66743755340576) q[2];
sx q[2];
rz(2.03259018261964) q[2];
sx q[2];
rz(9.39797993785843) q[2];
cx q[2],q[1];
rz(-0.905512988567352) q[1];
sx q[1];
rz(1.81549361546571) q[1];
sx q[1];
rz(10.1602561235349) q[1];
rz(-2.12511420249939) q[3];
sx q[3];
rz(-0.807625619573049) q[3];
sx q[3];
rz(9.22845638393565) q[3];
cx q[3],q[2];
rz(-4.50428533554077) q[2];
sx q[2];
rz(3.85393932660157) q[2];
sx q[2];
rz(7.27029821871921) q[2];
rz(-1.08460640907288) q[3];
sx q[3];
rz(6.82925668557221) q[3];
sx q[3];
rz(11.9454600572507) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.5793594121933) q[0];
sx q[0];
rz(2.14913037617738) q[0];
sx q[0];
rz(7.87672755717441) q[0];
rz(0.690018951892853) q[1];
sx q[1];
rz(0.0236121733957013) q[1];
sx q[1];
rz(14.003007388107) q[1];
cx q[1],q[0];
rz(-0.489559173583984) q[0];
sx q[0];
rz(4.192810805636) q[0];
sx q[0];
rz(6.23295066355869) q[0];
rz(0.021996270865202) q[2];
sx q[2];
rz(5.39699593384797) q[2];
sx q[2];
rz(5.0697026014249) q[2];
cx q[2],q[1];
rz(-0.440139055252075) q[1];
sx q[1];
rz(4.65285709698732) q[1];
sx q[1];
rz(9.9671520948331) q[1];
rz(1.11476838588715) q[3];
sx q[3];
rz(1.83979621727998) q[3];
sx q[3];
rz(7.59166524409457) q[3];
cx q[3],q[2];
rz(0.373166084289551) q[2];
sx q[2];
rz(1.23701408703858) q[2];
sx q[2];
rz(10.0096646904866) q[2];
rz(1.59604263305664) q[3];
sx q[3];
rz(4.86036744912202) q[3];
sx q[3];
rz(13.4808001279752) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.57174134254456) q[0];
sx q[0];
rz(5.36570826371247) q[0];
sx q[0];
rz(10.9917719125669) q[0];
rz(-0.517081081867218) q[1];
sx q[1];
rz(5.34773555596406) q[1];
sx q[1];
rz(5.01659009455844) q[1];
cx q[1],q[0];
rz(2.1684730052948) q[0];
sx q[0];
rz(3.13550973462919) q[0];
sx q[0];
rz(10.3951076626699) q[0];
rz(0.5549556016922) q[2];
sx q[2];
rz(5.43640700181062) q[2];
sx q[2];
rz(11.5727329015653) q[2];
cx q[2],q[1];
rz(3.47149467468262) q[1];
sx q[1];
rz(1.92018297513063) q[1];
sx q[1];
rz(8.84600940941974) q[1];
rz(2.35654902458191) q[3];
sx q[3];
rz(0.172212990122386) q[3];
sx q[3];
rz(5.88689777850314) q[3];
cx q[3],q[2];
rz(2.76855850219727) q[2];
sx q[2];
rz(3.94217506249482) q[2];
sx q[2];
rz(7.5558687210004) q[2];
rz(0.439676195383072) q[3];
sx q[3];
rz(5.12377801735932) q[3];
sx q[3];
rz(9.36018022745057) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.83982968330383) q[0];
sx q[0];
rz(6.29945674737031) q[0];
sx q[0];
rz(9.72104484438106) q[0];
rz(3.07731938362122) q[1];
sx q[1];
rz(4.60628202755982) q[1];
sx q[1];
rz(4.77211902140781) q[1];
cx q[1],q[0];
rz(-0.134162694215775) q[0];
sx q[0];
rz(1.90159705479676) q[0];
sx q[0];
rz(7.99148139952823) q[0];
rz(1.82385540008545) q[2];
sx q[2];
rz(2.45292434294755) q[2];
sx q[2];
rz(10.5873989820401) q[2];
cx q[2],q[1];
rz(1.84781968593597) q[1];
sx q[1];
rz(6.6330248435312) q[1];
sx q[1];
rz(9.81004924177333) q[1];
rz(3.04758143424988) q[3];
sx q[3];
rz(3.26324908633763) q[3];
sx q[3];
rz(10.3554992437284) q[3];
cx q[3],q[2];
rz(5.73714733123779) q[2];
sx q[2];
rz(4.57644978364045) q[2];
sx q[2];
rz(7.40724060534641) q[2];
rz(1.83736395835876) q[3];
sx q[3];
rz(5.9874383529001) q[3];
sx q[3];
rz(6.22502205371066) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.87033629417419) q[0];
sx q[0];
rz(-0.791263429326467) q[0];
sx q[0];
rz(13.7331914663236) q[0];
rz(1.60094904899597) q[1];
sx q[1];
rz(2.84380421240861) q[1];
sx q[1];
rz(10.7513311862867) q[1];
cx q[1],q[0];
rz(-0.0139152826741338) q[0];
sx q[0];
rz(1.81550374825532) q[0];
sx q[0];
rz(10.3266148924749) q[0];
rz(-4.30474758148193) q[2];
sx q[2];
rz(4.37657407124574) q[2];
sx q[2];
rz(7.78282079695865) q[2];
cx q[2],q[1];
rz(-1.57591438293457) q[1];
sx q[1];
rz(4.03131726582582) q[1];
sx q[1];
rz(12.2921607255857) q[1];
rz(-3.54285335540771) q[3];
sx q[3];
rz(4.80235925515229) q[3];
sx q[3];
rz(12.0048141240995) q[3];
cx q[3],q[2];
rz(5.39655351638794) q[2];
sx q[2];
rz(2.97802636225755) q[2];
sx q[2];
rz(1.44760415553256) q[2];
rz(0.126689538359642) q[3];
sx q[3];
rz(1.66197517712648) q[3];
sx q[3];
rz(7.39921948908969) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.666176855564117) q[0];
sx q[0];
rz(3.47256952722604) q[0];
sx q[0];
rz(8.49571505784198) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.01429224014282) q[1];
sx q[1];
rz(1.45667723019654) q[1];
sx q[1];
rz(11.8180560827176) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.0735437050461769) q[2];
sx q[2];
rz(4.86427512963349) q[2];
sx q[2];
rz(10.5430915117185) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.89714753627777) q[3];
sx q[3];
rz(4.66367534001405) q[3];
sx q[3];
rz(11.5933508634488) q[3];
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
