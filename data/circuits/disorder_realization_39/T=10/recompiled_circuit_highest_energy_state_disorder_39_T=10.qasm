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
rz(0.612528741359711) q[0];
sx q[0];
rz(3.77473685343797) q[0];
sx q[0];
rz(9.78821910022899) q[0];
rz(-0.587030470371246) q[1];
sx q[1];
rz(3.67098495562608) q[1];
sx q[1];
rz(10.9324728012006) q[1];
cx q[1],q[0];
rz(0.310012191534042) q[0];
sx q[0];
rz(3.2477846463495) q[0];
sx q[0];
rz(9.91805783509418) q[0];
rz(2.15632390975952) q[2];
sx q[2];
rz(3.94799199898774) q[2];
sx q[2];
rz(7.78091309069797) q[2];
cx q[2],q[1];
rz(0.225442931056023) q[1];
sx q[1];
rz(3.74656418164308) q[1];
sx q[1];
rz(12.1848387479703) q[1];
rz(1.02302193641663) q[3];
sx q[3];
rz(4.0324776490503) q[3];
sx q[3];
rz(9.74469763635799) q[3];
cx q[3],q[2];
rz(-1.86813604831696) q[2];
sx q[2];
rz(3.2324148436361) q[2];
sx q[2];
rz(13.0167503118436) q[2];
rz(-0.069933794438839) q[3];
sx q[3];
rz(4.25692561467225) q[3];
sx q[3];
rz(9.86703676580592) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0623341985046864) q[0];
sx q[0];
rz(3.96919033129747) q[0];
sx q[0];
rz(9.72390273808643) q[0];
rz(0.660807192325592) q[1];
sx q[1];
rz(4.25186899502809) q[1];
sx q[1];
rz(10.3482470273893) q[1];
cx q[1],q[0];
rz(-0.740893006324768) q[0];
sx q[0];
rz(2.77411964734132) q[0];
sx q[0];
rz(8.92273703812763) q[0];
rz(0.820060133934021) q[2];
sx q[2];
rz(2.82446494896943) q[2];
sx q[2];
rz(9.79608381389781) q[2];
cx q[2],q[1];
rz(-0.266580998897552) q[1];
sx q[1];
rz(4.34592667420442) q[1];
sx q[1];
rz(9.49313866942331) q[1];
rz(0.410090357065201) q[3];
sx q[3];
rz(3.76905825932557) q[3];
sx q[3];
rz(10.2347597241323) q[3];
cx q[3],q[2];
rz(-0.4835165143013) q[2];
sx q[2];
rz(2.74077555735643) q[2];
sx q[2];
rz(9.50868035703107) q[2];
rz(0.995332539081573) q[3];
sx q[3];
rz(4.29867282708222) q[3];
sx q[3];
rz(10.9356273174207) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.296265721321106) q[0];
sx q[0];
rz(4.3771255334192) q[0];
sx q[0];
rz(10.3873481511991) q[0];
rz(0.307490974664688) q[1];
sx q[1];
rz(2.25667849381501) q[1];
sx q[1];
rz(9.64371553658649) q[1];
cx q[1],q[0];
rz(0.229029521346092) q[0];
sx q[0];
rz(3.49467259843881) q[0];
sx q[0];
rz(10.7933054923932) q[0];
rz(0.287880152463913) q[2];
sx q[2];
rz(3.84054103692109) q[2];
sx q[2];
rz(10.7538028716962) q[2];
cx q[2],q[1];
rz(2.15561985969543) q[1];
sx q[1];
rz(3.54114479024942) q[1];
sx q[1];
rz(9.41750232958376) q[1];
rz(0.955730974674225) q[3];
sx q[3];
rz(4.13375804026658) q[3];
sx q[3];
rz(9.71075532435581) q[3];
cx q[3],q[2];
rz(1.68300497531891) q[2];
sx q[2];
rz(4.53543630440766) q[2];
sx q[2];
rz(10.9279286622922) q[2];
rz(0.250348687171936) q[3];
sx q[3];
rz(3.86545059283311) q[3];
sx q[3];
rz(9.22631331383392) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.649716258049011) q[0];
sx q[0];
rz(3.46577760775621) q[0];
sx q[0];
rz(9.74157104491397) q[0];
rz(0.0162387955933809) q[1];
sx q[1];
rz(3.93411752780015) q[1];
sx q[1];
rz(9.17646079360648) q[1];
cx q[1],q[0];
rz(0.0288489479571581) q[0];
sx q[0];
rz(4.09738758404786) q[0];
sx q[0];
rz(10.0614951610486) q[0];
rz(1.10359036922455) q[2];
sx q[2];
rz(4.64573696454103) q[2];
sx q[2];
rz(9.20824027656719) q[2];
cx q[2],q[1];
rz(-0.592453062534332) q[1];
sx q[1];
rz(3.49288740952546) q[1];
sx q[1];
rz(11.3681223153989) q[1];
rz(-0.218960717320442) q[3];
sx q[3];
rz(4.90013507206971) q[3];
sx q[3];
rz(10.3565333843152) q[3];
cx q[3],q[2];
rz(1.87796401977539) q[2];
sx q[2];
rz(3.0402277131849) q[2];
sx q[2];
rz(8.2807689666669) q[2];
rz(0.762172400951385) q[3];
sx q[3];
rz(3.6666483600908) q[3];
sx q[3];
rz(10.4163146376531) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.856319725513458) q[0];
sx q[0];
rz(3.05064360250766) q[0];
sx q[0];
rz(10.4772401809613) q[0];
rz(0.490280300378799) q[1];
sx q[1];
rz(4.23043075402314) q[1];
sx q[1];
rz(9.55000219344302) q[1];
cx q[1],q[0];
rz(0.767160475254059) q[0];
sx q[0];
rz(4.31977179844911) q[0];
sx q[0];
rz(10.6292909145276) q[0];
rz(1.13128578662872) q[2];
sx q[2];
rz(3.73488232691819) q[2];
sx q[2];
rz(9.78020725249454) q[2];
cx q[2],q[1];
rz(-1.09088063240051) q[1];
sx q[1];
rz(4.7844250520044) q[1];
sx q[1];
rz(9.17439115642711) q[1];
rz(1.35087418556213) q[3];
sx q[3];
rz(2.27079048951203) q[3];
sx q[3];
rz(9.10784230231448) q[3];
cx q[3],q[2];
rz(-0.483097165822983) q[2];
sx q[2];
rz(2.60396370490129) q[2];
sx q[2];
rz(10.9775013685147) q[2];
rz(0.0565648637712002) q[3];
sx q[3];
rz(4.65571609337861) q[3];
sx q[3];
rz(9.70208627580806) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.561894297599792) q[0];
sx q[0];
rz(3.60598471959169) q[0];
sx q[0];
rz(9.85249597429439) q[0];
rz(2.60042095184326) q[1];
sx q[1];
rz(2.73105570872361) q[1];
sx q[1];
rz(7.64689478873416) q[1];
cx q[1],q[0];
rz(-0.833745241165161) q[0];
sx q[0];
rz(3.7529760320955) q[0];
sx q[0];
rz(9.61672005652591) q[0];
rz(1.14364957809448) q[2];
sx q[2];
rz(4.76480546792085) q[2];
sx q[2];
rz(10.2018775105397) q[2];
cx q[2],q[1];
rz(1.02376711368561) q[1];
sx q[1];
rz(3.90899798472459) q[1];
sx q[1];
rz(8.97895214556857) q[1];
rz(-0.0250442922115326) q[3];
sx q[3];
rz(3.55450985034043) q[3];
sx q[3];
rz(9.45191591269478) q[3];
cx q[3],q[2];
rz(0.0428187288343906) q[2];
sx q[2];
rz(3.66136750777299) q[2];
sx q[2];
rz(9.8341393828313) q[2];
rz(0.219710752367973) q[3];
sx q[3];
rz(4.6909797509485) q[3];
sx q[3];
rz(11.0767693281095) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.701942443847656) q[0];
sx q[0];
rz(2.72034287651116) q[0];
sx q[0];
rz(9.56559123694106) q[0];
rz(0.536786496639252) q[1];
sx q[1];
rz(4.49708464940126) q[1];
sx q[1];
rz(9.80145121215984) q[1];
cx q[1],q[0];
rz(0.464640647172928) q[0];
sx q[0];
rz(3.44021481473977) q[0];
sx q[0];
rz(9.99666265248462) q[0];
rz(0.895869016647339) q[2];
sx q[2];
rz(4.39441195328767) q[2];
sx q[2];
rz(9.03261995910808) q[2];
cx q[2],q[1];
rz(-0.721953272819519) q[1];
sx q[1];
rz(3.93414655526216) q[1];
sx q[1];
rz(11.1780352353971) q[1];
rz(1.22760951519012) q[3];
sx q[3];
rz(3.57954120834405) q[3];
sx q[3];
rz(10.818672990791) q[3];
cx q[3],q[2];
rz(-1.10303997993469) q[2];
sx q[2];
rz(4.30781033833558) q[2];
sx q[2];
rz(10.2260116100232) q[2];
rz(1.07711052894592) q[3];
sx q[3];
rz(2.88476497133309) q[3];
sx q[3];
rz(9.56725516020461) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.242566630244255) q[0];
sx q[0];
rz(3.30225973029668) q[0];
sx q[0];
rz(10.3658627033155) q[0];
rz(0.896489918231964) q[1];
sx q[1];
rz(4.56003049214418) q[1];
sx q[1];
rz(10.2394998431127) q[1];
cx q[1],q[0];
rz(0.77234947681427) q[0];
sx q[0];
rz(2.60084781249101) q[0];
sx q[0];
rz(10.2289914846341) q[0];
rz(1.28456711769104) q[2];
sx q[2];
rz(2.96364179451997) q[2];
sx q[2];
rz(9.62309493719741) q[2];
cx q[2],q[1];
rz(-0.434896975755692) q[1];
sx q[1];
rz(3.86623826821382) q[1];
sx q[1];
rz(9.9503109216611) q[1];
rz(0.548050463199615) q[3];
sx q[3];
rz(3.78324678738649) q[3];
sx q[3];
rz(8.28172085284396) q[3];
cx q[3],q[2];
rz(0.616463780403137) q[2];
sx q[2];
rz(2.50696709950502) q[2];
sx q[2];
rz(10.0015998840253) q[2];
rz(-0.473546266555786) q[3];
sx q[3];
rz(1.99279311497743) q[3];
sx q[3];
rz(9.46811961232826) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.438054949045181) q[0];
sx q[0];
rz(3.89710846741731) q[0];
sx q[0];
rz(9.33687142132922) q[0];
rz(1.58806228637695) q[1];
sx q[1];
rz(3.72595885594422) q[1];
sx q[1];
rz(9.16715270876094) q[1];
cx q[1],q[0];
rz(-0.157302975654602) q[0];
sx q[0];
rz(3.57330578764016) q[0];
sx q[0];
rz(9.71233273147746) q[0];
rz(0.0974606052041054) q[2];
sx q[2];
rz(3.77123603423173) q[2];
sx q[2];
rz(9.82277608513042) q[2];
cx q[2],q[1];
rz(-0.963660717010498) q[1];
sx q[1];
rz(3.63832193811471) q[1];
sx q[1];
rz(11.5050921201627) q[1];
rz(0.576453983783722) q[3];
sx q[3];
rz(4.21070841153199) q[3];
sx q[3];
rz(9.69493002294704) q[3];
cx q[3],q[2];
rz(-0.0755001753568649) q[2];
sx q[2];
rz(4.82176414330537) q[2];
sx q[2];
rz(9.50256564318343) q[2];
rz(1.05456876754761) q[3];
sx q[3];
rz(3.881795616942) q[3];
sx q[3];
rz(10.0757004380147) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.24423778057098) q[0];
sx q[0];
rz(3.17503447656567) q[0];
sx q[0];
rz(9.86053205131694) q[0];
rz(-0.0325953662395477) q[1];
sx q[1];
rz(3.87899413903291) q[1];
sx q[1];
rz(11.2184294223706) q[1];
cx q[1],q[0];
rz(1.45300579071045) q[0];
sx q[0];
rz(3.96191379626329) q[0];
sx q[0];
rz(9.41508268228873) q[0];
rz(1.16849195957184) q[2];
sx q[2];
rz(3.70554468234117) q[2];
sx q[2];
rz(10.0590908884923) q[2];
cx q[2],q[1];
rz(2.09501981735229) q[1];
sx q[1];
rz(3.1802675147825) q[1];
sx q[1];
rz(9.05564347504779) q[1];
rz(1.25601422786713) q[3];
sx q[3];
rz(4.48988583882386) q[3];
sx q[3];
rz(9.5137078076522) q[3];
cx q[3],q[2];
rz(1.21604478359222) q[2];
sx q[2];
rz(4.17570844491059) q[2];
sx q[2];
rz(10.1810837745587) q[2];
rz(0.0821527242660522) q[3];
sx q[3];
rz(3.53393975098664) q[3];
sx q[3];
rz(10.1160121321599) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.49741804599762) q[0];
sx q[0];
rz(4.17876270611817) q[0];
sx q[0];
rz(8.39338729380771) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.299401044845581) q[1];
sx q[1];
rz(4.3871615250879) q[1];
sx q[1];
rz(9.88927564620181) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.13819909095764) q[2];
sx q[2];
rz(3.90220943291719) q[2];
sx q[2];
rz(9.67379598914787) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.321728318929672) q[3];
sx q[3];
rz(4.46321299870545) q[3];
sx q[3];
rz(9.95285311936542) q[3];
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
