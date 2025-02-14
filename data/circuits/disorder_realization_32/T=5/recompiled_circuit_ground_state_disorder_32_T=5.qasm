OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.13429927825928) q[0];
sx q[0];
rz(8.22598949273164) q[0];
sx q[0];
rz(10.2099570393483) q[0];
rz(2.56853175163269) q[1];
sx q[1];
rz(2.1905724128061) q[1];
sx q[1];
rz(11.934926724426) q[1];
cx q[1],q[0];
rz(1.31086897850037) q[0];
sx q[0];
rz(2.00183025200898) q[0];
sx q[0];
rz(11.8628680467527) q[0];
rz(3.98786330223083) q[2];
sx q[2];
rz(4.0738475044542) q[2];
sx q[2];
rz(7.19061539172336) q[2];
cx q[2],q[1];
rz(2.54413032531738) q[1];
sx q[1];
rz(4.2825507243448) q[1];
sx q[1];
rz(11.6882054567258) q[1];
rz(8.29882717132568) q[3];
sx q[3];
rz(1.66037836869294) q[3];
sx q[3];
rz(11.2057022809903) q[3];
cx q[3],q[2];
rz(-2.63739204406738) q[2];
sx q[2];
rz(4.84554341633851) q[2];
sx q[2];
rz(12.8379015684049) q[2];
rz(-0.0463937520980835) q[3];
sx q[3];
rz(4.55191138585145) q[3];
sx q[3];
rz(9.79240573047801) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.725359976291656) q[0];
sx q[0];
rz(3.16874171060557) q[0];
sx q[0];
rz(9.87038258313342) q[0];
rz(-0.360005706548691) q[1];
sx q[1];
rz(1.58411923249299) q[1];
sx q[1];
rz(7.2602212190549) q[1];
cx q[1],q[0];
rz(-0.925218999385834) q[0];
sx q[0];
rz(2.3686102946573) q[0];
sx q[0];
rz(7.08368036746188) q[0];
rz(-2.29270553588867) q[2];
sx q[2];
rz(4.18861916859681) q[2];
sx q[2];
rz(8.61823043822452) q[2];
cx q[2],q[1];
rz(-0.47525542974472) q[1];
sx q[1];
rz(-2.52863916556304) q[1];
sx q[1];
rz(9.68965358137294) q[1];
rz(0.825047969818115) q[3];
sx q[3];
rz(3.99833378394181) q[3];
sx q[3];
rz(6.03458759783908) q[3];
cx q[3],q[2];
rz(0.526860237121582) q[2];
sx q[2];
rz(4.48561850388581) q[2];
sx q[2];
rz(10.3669999599378) q[2];
rz(-3.63803172111511) q[3];
sx q[3];
rz(1.67832914193208) q[3];
sx q[3];
rz(10.4924382924955) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.08388578891754) q[0];
sx q[0];
rz(4.47133913834626) q[0];
sx q[0];
rz(10.0809913039128) q[0];
rz(3.62197113037109) q[1];
sx q[1];
rz(2.18001327117021) q[1];
sx q[1];
rz(6.89613244532748) q[1];
cx q[1],q[0];
rz(-0.975263237953186) q[0];
sx q[0];
rz(3.56113931735093) q[0];
sx q[0];
rz(6.52320954798862) q[0];
rz(3.13425517082214) q[2];
sx q[2];
rz(3.97632447083528) q[2];
sx q[2];
rz(13.875128722183) q[2];
cx q[2],q[1];
rz(2.61005449295044) q[1];
sx q[1];
rz(6.23894682725007) q[1];
sx q[1];
rz(10.2860371231954) q[1];
rz(-0.0691727250814438) q[3];
sx q[3];
rz(1.07284799416596) q[3];
sx q[3];
rz(8.7205353140752) q[3];
cx q[3],q[2];
rz(1.93566846847534) q[2];
sx q[2];
rz(4.19904032548005) q[2];
sx q[2];
rz(12.0419659376065) q[2];
rz(-2.27624011039734) q[3];
sx q[3];
rz(1.05823460419709) q[3];
sx q[3];
rz(10.0852226376454) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.65372896194458) q[0];
sx q[0];
rz(4.97625103791291) q[0];
sx q[0];
rz(9.52104675620004) q[0];
rz(-0.0138215664774179) q[1];
sx q[1];
rz(3.642178388434) q[1];
sx q[1];
rz(7.3355271577756) q[1];
cx q[1],q[0];
rz(3.83214163780212) q[0];
sx q[0];
rz(4.75936845143373) q[0];
sx q[0];
rz(12.6881813764493) q[0];
rz(0.550192952156067) q[2];
sx q[2];
rz(1.32507041295106) q[2];
sx q[2];
rz(9.46745909228131) q[2];
cx q[2],q[1];
rz(-4.01568508148193) q[1];
sx q[1];
rz(3.51928222377832) q[1];
sx q[1];
rz(7.34021804331943) q[1];
rz(0.525548040866852) q[3];
sx q[3];
rz(1.60176721413667) q[3];
sx q[3];
rz(10.2455746889035) q[3];
cx q[3],q[2];
rz(-2.4200165271759) q[2];
sx q[2];
rz(3.67268422444398) q[2];
sx q[2];
rz(9.94377831219836) q[2];
rz(2.04299330711365) q[3];
sx q[3];
rz(4.82692316372926) q[3];
sx q[3];
rz(7.01321265696689) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.622434318065643) q[0];
sx q[0];
rz(3.33463746507699) q[0];
sx q[0];
rz(8.69053599833652) q[0];
rz(-1.28920674324036) q[1];
sx q[1];
rz(5.23808446724946) q[1];
sx q[1];
rz(11.6578693151395) q[1];
cx q[1],q[0];
rz(-3.4953305721283) q[0];
sx q[0];
rz(4.0578928907686) q[0];
sx q[0];
rz(8.70245698689624) q[0];
rz(-1.3650586605072) q[2];
sx q[2];
rz(5.70218792756135) q[2];
sx q[2];
rz(6.57145569323703) q[2];
cx q[2],q[1];
rz(3.31401205062866) q[1];
sx q[1];
rz(10.1611687262827) q[1];
sx q[1];
rz(12.669951415054) q[1];
rz(1.59518837928772) q[3];
sx q[3];
rz(4.44581869442994) q[3];
sx q[3];
rz(9.74256486295863) q[3];
cx q[3],q[2];
rz(3.18613338470459) q[2];
sx q[2];
rz(3.9278246482187) q[2];
sx q[2];
rz(13.2499892473142) q[2];
rz(-1.09400701522827) q[3];
sx q[3];
rz(1.32355991204316) q[3];
sx q[3];
rz(10.9971421718518) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.08786582946777) q[0];
sx q[0];
rz(1.69310668309266) q[0];
sx q[0];
rz(9.85039407610103) q[0];
rz(-5.67641639709473) q[1];
sx q[1];
rz(-0.687770454091481) q[1];
sx q[1];
rz(16.9105763196866) q[1];
cx q[1],q[0];
rz(1.37891137599945) q[0];
sx q[0];
rz(5.27820840676362) q[0];
sx q[0];
rz(13.5480561017911) q[0];
rz(-1.50496971607208) q[2];
sx q[2];
rz(4.69786122639711) q[2];
sx q[2];
rz(9.6043629258792) q[2];
cx q[2],q[1];
rz(-1.76406705379486) q[1];
sx q[1];
rz(4.21852353413636) q[1];
sx q[1];
rz(11.7371625661771) q[1];
rz(0.789389729499817) q[3];
sx q[3];
rz(3.46297368605668) q[3];
sx q[3];
rz(7.70214853285953) q[3];
cx q[3],q[2];
rz(0.468350291252136) q[2];
sx q[2];
rz(3.6709293444925) q[2];
sx q[2];
rz(5.65444419383212) q[2];
rz(0.337442338466644) q[3];
sx q[3];
rz(4.50135925610597) q[3];
sx q[3];
rz(3.90862081050082) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.59359800815582) q[0];
sx q[0];
rz(9.35473743279512) q[0];
sx q[0];
rz(10.7647536754529) q[0];
rz(0.105466708540916) q[1];
sx q[1];
rz(4.12300077279145) q[1];
sx q[1];
rz(6.13890311717197) q[1];
cx q[1],q[0];
rz(0.540504157543182) q[0];
sx q[0];
rz(1.27990368207032) q[0];
sx q[0];
rz(5.11607692240878) q[0];
rz(2.91530919075012) q[2];
sx q[2];
rz(3.86249783833558) q[2];
sx q[2];
rz(12.5180346727292) q[2];
cx q[2],q[1];
rz(-1.07935857772827) q[1];
sx q[1];
rz(6.69382754166658) q[1];
sx q[1];
rz(10.4525134324948) q[1];
rz(-0.543428361415863) q[3];
sx q[3];
rz(5.30967608292634) q[3];
sx q[3];
rz(10.8187887430112) q[3];
cx q[3],q[2];
rz(-1.10628008842468) q[2];
sx q[2];
rz(-1.99297460715239) q[2];
sx q[2];
rz(10.2627354621808) q[2];
rz(2.64181399345398) q[3];
sx q[3];
rz(5.48983732064302) q[3];
sx q[3];
rz(11.8651852369229) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.42262935638428) q[0];
sx q[0];
rz(2.03276494343812) q[0];
sx q[0];
rz(7.03441903590366) q[0];
rz(-5.02960157394409) q[1];
sx q[1];
rz(1.8389848788553) q[1];
sx q[1];
rz(17.4814529180448) q[1];
cx q[1],q[0];
rz(3.87276911735535) q[0];
sx q[0];
rz(5.50395837624604) q[0];
sx q[0];
rz(8.77394524811908) q[0];
rz(1.25475716590881) q[2];
sx q[2];
rz(4.99228671391542) q[2];
sx q[2];
rz(8.95101303457423) q[2];
cx q[2],q[1];
rz(0.66361278295517) q[1];
sx q[1];
rz(3.91938212712342) q[1];
sx q[1];
rz(9.54549090414449) q[1];
rz(3.71552872657776) q[3];
sx q[3];
rz(4.40484026272828) q[3];
sx q[3];
rz(9.06304833888217) q[3];
cx q[3],q[2];
rz(-0.832968354225159) q[2];
sx q[2];
rz(5.91152778466279) q[2];
sx q[2];
rz(12.5652024507444) q[2];
rz(-0.332311689853668) q[3];
sx q[3];
rz(1.46095910866792) q[3];
sx q[3];
rz(9.90440270899936) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.76976561546326) q[0];
sx q[0];
rz(5.06608656247193) q[0];
sx q[0];
rz(8.32869253157779) q[0];
rz(1.74258840084076) q[1];
sx q[1];
rz(4.64693728287751) q[1];
sx q[1];
rz(8.50491169690295) q[1];
cx q[1],q[0];
rz(-0.771147847175598) q[0];
sx q[0];
rz(-0.439752904576711) q[0];
sx q[0];
rz(9.52475765197679) q[0];
rz(0.418678998947144) q[2];
sx q[2];
rz(4.96953883965547) q[2];
sx q[2];
rz(3.74953696726962) q[2];
cx q[2],q[1];
rz(4.0191535949707) q[1];
sx q[1];
rz(4.56207934220368) q[1];
sx q[1];
rz(7.20433948039218) q[1];
rz(-2.60872673988342) q[3];
sx q[3];
rz(1.69543305237825) q[3];
sx q[3];
rz(10.2232208013456) q[3];
cx q[3],q[2];
rz(-1.21907353401184) q[2];
sx q[2];
rz(3.386842237907) q[2];
sx q[2];
rz(17.218107199661) q[2];
rz(-2.50569295883179) q[3];
sx q[3];
rz(3.84680780966813) q[3];
sx q[3];
rz(12.0517914056699) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.95663273334503) q[0];
sx q[0];
rz(4.0086386521631) q[0];
sx q[0];
rz(7.92910561560794) q[0];
rz(3.02368593215942) q[1];
sx q[1];
rz(4.88881090481813) q[1];
sx q[1];
rz(9.16798741220638) q[1];
cx q[1],q[0];
rz(1.19081556797028) q[0];
sx q[0];
rz(5.84638181527192) q[0];
sx q[0];
rz(9.29245824217006) q[0];
rz(1.13638639450073) q[2];
sx q[2];
rz(6.91396013100679) q[2];
sx q[2];
rz(14.3776616811673) q[2];
cx q[2],q[1];
rz(0.499187111854553) q[1];
sx q[1];
rz(2.74912637670571) q[1];
sx q[1];
rz(10.2521719098012) q[1];
rz(4.15494632720947) q[3];
sx q[3];
rz(8.80189576943452) q[3];
sx q[3];
rz(10.4723561763684) q[3];
cx q[3],q[2];
rz(-0.397664964199066) q[2];
sx q[2];
rz(3.88508346875245) q[2];
sx q[2];
rz(14.7579860448758) q[2];
rz(-0.602036118507385) q[3];
sx q[3];
rz(5.0461941083246) q[3];
sx q[3];
rz(5.94502732752963) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.00297701358795) q[0];
sx q[0];
rz(1.00666323502595) q[0];
sx q[0];
rz(10.7465820074002) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-2.7019636631012) q[1];
sx q[1];
rz(1.72929075558717) q[1];
sx q[1];
rz(6.9412092924039) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-1.17912423610687) q[2];
sx q[2];
rz(3.56561708648736) q[2];
sx q[2];
rz(8.67834988831683) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.90956282615662) q[3];
sx q[3];
rz(5.69199791749055) q[3];
sx q[3];
rz(6.81431052683994) q[3];
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
