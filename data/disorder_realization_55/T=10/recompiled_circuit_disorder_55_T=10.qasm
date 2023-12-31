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
rz(0.546859323978424) q[0];
sx q[0];
rz(4.65805533726747) q[0];
sx q[0];
rz(9.16049945949718) q[0];
rz(-0.973794102668762) q[1];
sx q[1];
rz(5.07305398781831) q[1];
sx q[1];
rz(10.1600253343503) q[1];
cx q[1],q[0];
rz(-0.111569464206696) q[0];
sx q[0];
rz(2.22462472518022) q[0];
sx q[0];
rz(9.23810043036147) q[0];
rz(2.63442826271057) q[2];
sx q[2];
rz(5.3274292071634) q[2];
sx q[2];
rz(10.4787215948026) q[2];
cx q[2],q[1];
rz(-0.689025104045868) q[1];
sx q[1];
rz(4.42635587056214) q[1];
sx q[1];
rz(13.9761857747953) q[1];
rz(1.52797448635101) q[3];
sx q[3];
rz(3.72245833476121) q[3];
sx q[3];
rz(9.04915944337055) q[3];
cx q[3],q[2];
rz(3.8111035823822) q[2];
sx q[2];
rz(1.84102037747438) q[2];
sx q[2];
rz(10.5286366701047) q[2];
rz(-1.87072288990021) q[3];
sx q[3];
rz(4.36936012108857) q[3];
sx q[3];
rz(9.69881329535648) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.02741491794586) q[0];
sx q[0];
rz(2.61296209891374) q[0];
sx q[0];
rz(9.86115277408763) q[0];
rz(0.462886810302734) q[1];
sx q[1];
rz(4.17916158040101) q[1];
sx q[1];
rz(12.832481122009) q[1];
cx q[1],q[0];
rz(1.76040959358215) q[0];
sx q[0];
rz(5.09043696721131) q[0];
sx q[0];
rz(10.1860999822538) q[0];
rz(5.70772075653076) q[2];
sx q[2];
rz(1.43337121804292) q[2];
sx q[2];
rz(10.3753568291585) q[2];
cx q[2],q[1];
rz(3.52922868728638) q[1];
sx q[1];
rz(4.1738670190149) q[1];
sx q[1];
rz(10.1613695382993) q[1];
rz(1.38807213306427) q[3];
sx q[3];
rz(5.31126657326753) q[3];
sx q[3];
rz(9.53651982396051) q[3];
cx q[3],q[2];
rz(0.464885473251343) q[2];
sx q[2];
rz(4.41839090188081) q[2];
sx q[2];
rz(8.91328409909412) q[2];
rz(2.33208465576172) q[3];
sx q[3];
rz(4.67291644414003) q[3];
sx q[3];
rz(12.2787122487943) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.706161022186279) q[0];
sx q[0];
rz(1.54787388642366) q[0];
sx q[0];
rz(10.3535154223363) q[0];
rz(1.73548948764801) q[1];
sx q[1];
rz(3.84166005452211) q[1];
sx q[1];
rz(8.07764909266635) q[1];
cx q[1],q[0];
rz(0.980262100696564) q[0];
sx q[0];
rz(7.59248557885224) q[0];
sx q[0];
rz(8.32786509989902) q[0];
rz(3.42682290077209) q[2];
sx q[2];
rz(4.8642305453592) q[2];
sx q[2];
rz(8.38857851027652) q[2];
cx q[2],q[1];
rz(0.24088628590107) q[1];
sx q[1];
rz(0.979292305307933) q[1];
sx q[1];
rz(7.50331995486423) q[1];
rz(1.41955649852753) q[3];
sx q[3];
rz(2.27217069466645) q[3];
sx q[3];
rz(8.69387934207126) q[3];
cx q[3],q[2];
rz(1.99462628364563) q[2];
sx q[2];
rz(4.8965769131952) q[2];
sx q[2];
rz(10.9468580245893) q[2];
rz(-0.264320313930511) q[3];
sx q[3];
rz(4.17804995377595) q[3];
sx q[3];
rz(9.9207308947961) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.792148232460022) q[0];
sx q[0];
rz(5.13481334050233) q[0];
sx q[0];
rz(10.3904785871427) q[0];
rz(0.722158849239349) q[1];
sx q[1];
rz(4.64581433136994) q[1];
sx q[1];
rz(9.98453429936572) q[1];
cx q[1],q[0];
rz(0.00696654664352536) q[0];
sx q[0];
rz(-0.43381938139861) q[0];
sx q[0];
rz(8.04800949095889) q[0];
rz(2.38440561294556) q[2];
sx q[2];
rz(4.15421238740022) q[2];
sx q[2];
rz(7.60814545153781) q[2];
cx q[2],q[1];
rz(-2.8371114730835) q[1];
sx q[1];
rz(4.75803175766999) q[1];
sx q[1];
rz(10.0405121207158) q[1];
rz(2.35999584197998) q[3];
sx q[3];
rz(4.54410246213014) q[3];
sx q[3];
rz(10.0336753487508) q[3];
cx q[3],q[2];
rz(-0.427977830171585) q[2];
sx q[2];
rz(4.6503862460428) q[2];
sx q[2];
rz(10.8493644952695) q[2];
rz(2.88118410110474) q[3];
sx q[3];
rz(4.89070919354493) q[3];
sx q[3];
rz(8.99373075961276) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.9907306432724) q[0];
sx q[0];
rz(4.78514376481111) q[0];
sx q[0];
rz(9.79114029406711) q[0];
rz(1.54619610309601) q[1];
sx q[1];
rz(3.69551089604432) q[1];
sx q[1];
rz(9.76845721005603) q[1];
cx q[1],q[0];
rz(0.300037980079651) q[0];
sx q[0];
rz(2.24047878583009) q[0];
sx q[0];
rz(14.8511543035428) q[0];
rz(-0.992210268974304) q[2];
sx q[2];
rz(4.09698769648606) q[2];
sx q[2];
rz(10.8032421827237) q[2];
cx q[2],q[1];
rz(-3.7136287689209) q[1];
sx q[1];
rz(4.215156467753) q[1];
sx q[1];
rz(7.77003917693301) q[1];
rz(-3.16993546485901) q[3];
sx q[3];
rz(3.78564819891984) q[3];
sx q[3];
rz(10.945004439346) q[3];
cx q[3],q[2];
rz(1.39178240299225) q[2];
sx q[2];
rz(5.42386999924714) q[2];
sx q[2];
rz(12.5126256704251) q[2];
rz(1.73714768886566) q[3];
sx q[3];
rz(5.78564730485017) q[3];
sx q[3];
rz(12.2748088598172) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.486972600221634) q[0];
sx q[0];
rz(2.49168613751466) q[0];
sx q[0];
rz(10.1805374979894) q[0];
rz(0.0251563303172588) q[1];
sx q[1];
rz(5.35592928727204) q[1];
sx q[1];
rz(9.68451724051639) q[1];
cx q[1],q[0];
rz(1.32096695899963) q[0];
sx q[0];
rz(3.62988671858842) q[0];
sx q[0];
rz(10.1426121950071) q[0];
rz(0.0565683171153069) q[2];
sx q[2];
rz(0.336773784952708) q[2];
sx q[2];
rz(11.7925283670346) q[2];
cx q[2],q[1];
rz(0.425279200077057) q[1];
sx q[1];
rz(4.48640874226625) q[1];
sx q[1];
rz(12.9037308454435) q[1];
rz(1.55281937122345) q[3];
sx q[3];
rz(4.26308027108247) q[3];
sx q[3];
rz(9.20088583826228) q[3];
cx q[3],q[2];
rz(0.488668948411942) q[2];
sx q[2];
rz(3.94763890107209) q[2];
sx q[2];
rz(9.52538926749631) q[2];
rz(-0.182090222835541) q[3];
sx q[3];
rz(2.26333496172959) q[3];
sx q[3];
rz(11.2186838149945) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.9942193031311) q[0];
sx q[0];
rz(4.56187024910981) q[0];
sx q[0];
rz(9.7349444091241) q[0];
rz(2.63933300971985) q[1];
sx q[1];
rz(6.80719462235505) q[1];
sx q[1];
rz(10.0307366013448) q[1];
cx q[1],q[0];
rz(0.464424073696136) q[0];
sx q[0];
rz(5.62224546273286) q[0];
sx q[0];
rz(10.6234973430555) q[0];
rz(-0.975257635116577) q[2];
sx q[2];
rz(3.48950511415536) q[2];
sx q[2];
rz(13.2351290941159) q[2];
cx q[2],q[1];
rz(0.0529186055064201) q[1];
sx q[1];
rz(0.624019058542796) q[1];
sx q[1];
rz(8.81558135747119) q[1];
rz(3.68796181678772) q[3];
sx q[3];
rz(5.17672053177888) q[3];
sx q[3];
rz(11.1198096036832) q[3];
cx q[3],q[2];
rz(0.240177482366562) q[2];
sx q[2];
rz(1.18379739125306) q[2];
sx q[2];
rz(11.7135159730832) q[2];
rz(1.77157056331635) q[3];
sx q[3];
rz(4.60026160080964) q[3];
sx q[3];
rz(12.361401295654) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.2296427488327) q[0];
sx q[0];
rz(3.73756799300248) q[0];
sx q[0];
rz(10.8861089706342) q[0];
rz(-4.54449939727783) q[1];
sx q[1];
rz(4.11583391030366) q[1];
sx q[1];
rz(12.5023331403653) q[1];
cx q[1],q[0];
rz(-0.562115728855133) q[0];
sx q[0];
rz(2.17702827055986) q[0];
sx q[0];
rz(7.62530348300144) q[0];
rz(-1.26447880268097) q[2];
sx q[2];
rz(4.22644987900788) q[2];
sx q[2];
rz(10.8998099327008) q[2];
cx q[2],q[1];
rz(3.50913429260254) q[1];
sx q[1];
rz(1.60721746285493) q[1];
sx q[1];
rz(7.49879751204654) q[1];
rz(2.82963085174561) q[3];
sx q[3];
rz(7.81639638741548) q[3];
sx q[3];
rz(12.5288219213407) q[3];
cx q[3],q[2];
rz(-3.23280024528503) q[2];
sx q[2];
rz(3.77241531212861) q[2];
sx q[2];
rz(7.83861181735202) q[2];
rz(4.02979373931885) q[3];
sx q[3];
rz(5.04337647755677) q[3];
sx q[3];
rz(10.3541667818944) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.144261941313744) q[0];
sx q[0];
rz(1.06347957451875) q[0];
sx q[0];
rz(7.99286005496188) q[0];
rz(0.568884670734406) q[1];
sx q[1];
rz(2.60619965394075) q[1];
sx q[1];
rz(7.41102430819675) q[1];
cx q[1],q[0];
rz(-0.237524300813675) q[0];
sx q[0];
rz(4.72840360005433) q[0];
sx q[0];
rz(7.70525453089877) q[0];
rz(-0.682386577129364) q[2];
sx q[2];
rz(3.89355489810044) q[2];
sx q[2];
rz(8.38963732718631) q[2];
cx q[2],q[1];
rz(-1.7833263874054) q[1];
sx q[1];
rz(4.95404675801332) q[1];
sx q[1];
rz(4.70044181346103) q[1];
rz(-0.550171196460724) q[3];
sx q[3];
rz(2.50848254759843) q[3];
sx q[3];
rz(10.6844410657804) q[3];
cx q[3],q[2];
rz(8.896897315979) q[2];
sx q[2];
rz(3.99472937186296) q[2];
sx q[2];
rz(9.25112802385494) q[2];
rz(-0.336378335952759) q[3];
sx q[3];
rz(1.91895464261109) q[3];
sx q[3];
rz(5.89167735575839) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.43534517288208) q[0];
sx q[0];
rz(3.72883156140382) q[0];
sx q[0];
rz(11.1008593797605) q[0];
rz(2.31748747825623) q[1];
sx q[1];
rz(4.75442138512666) q[1];
sx q[1];
rz(9.99720475672885) q[1];
cx q[1],q[0];
rz(3.65206336975098) q[0];
sx q[0];
rz(4.23181024392182) q[0];
sx q[0];
rz(10.3710404395978) q[0];
rz(-0.517357349395752) q[2];
sx q[2];
rz(2.62470188935334) q[2];
sx q[2];
rz(9.38243847563072) q[2];
cx q[2],q[1];
rz(0.13092565536499) q[1];
sx q[1];
rz(4.92831686337525) q[1];
sx q[1];
rz(11.2399092674176) q[1];
rz(-0.862138450145721) q[3];
sx q[3];
rz(4.48816064198548) q[3];
sx q[3];
rz(11.0283546209256) q[3];
cx q[3],q[2];
rz(-0.779990136623383) q[2];
sx q[2];
rz(3.61611616809899) q[2];
sx q[2];
rz(8.88755688666507) q[2];
rz(1.05726635456085) q[3];
sx q[3];
rz(4.03310910065705) q[3];
sx q[3];
rz(10.1183933973233) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.82555496692657) q[0];
sx q[0];
rz(3.54719767172868) q[0];
sx q[0];
rz(11.0219340085904) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(5.91592931747437) q[1];
sx q[1];
rz(5.58722511132295) q[1];
sx q[1];
rz(14.2110066175382) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.831074655056) q[2];
sx q[2];
rz(3.88445338805253) q[2];
sx q[2];
rz(9.30960684864923) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.17814779281616) q[3];
sx q[3];
rz(0.52211061318452) q[3];
sx q[3];
rz(12.6763713121335) q[3];
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
