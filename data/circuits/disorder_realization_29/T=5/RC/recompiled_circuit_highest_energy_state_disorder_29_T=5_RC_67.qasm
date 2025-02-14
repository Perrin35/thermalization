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
rz(-0.40773243) q[0];
sx q[0];
rz(4.3507504) q[0];
sx q[0];
rz(11.103295) q[0];
rz(0.2585803) q[1];
sx q[1];
rz(-1.9039896) q[1];
sx q[1];
rz(3.0078476) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8834849) q[0];
sx q[0];
rz(-2.9843669) q[0];
sx q[0];
rz(-1.5327461) q[0];
rz(0.57853477) q[2];
sx q[2];
rz(-1.0867439) q[2];
sx q[2];
rz(0.59147385) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3773777) q[1];
sx q[1];
rz(-1.0395607) q[1];
sx q[1];
rz(1.2197391) q[1];
rz(-pi) q[2];
rz(1.5720346) q[3];
sx q[3];
rz(-0.15286013) q[3];
sx q[3];
rz(0.78939702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75970834) q[2];
sx q[2];
rz(-2.3298161) q[2];
sx q[2];
rz(1.831057) q[2];
rz(-1.367502) q[3];
sx q[3];
rz(-1.1279227) q[3];
sx q[3];
rz(0.016157063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1644208) q[0];
sx q[0];
rz(-1.9855969) q[0];
sx q[0];
rz(1.0154065) q[0];
rz(2.7482391) q[1];
sx q[1];
rz(-1.669408) q[1];
sx q[1];
rz(-2.4400087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3053646) q[0];
sx q[0];
rz(-0.71963632) q[0];
sx q[0];
rz(1.4426065) q[0];
rz(-pi) q[1];
rz(-3.1308514) q[2];
sx q[2];
rz(-1.2678384) q[2];
sx q[2];
rz(-3.0634653) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7483731) q[1];
sx q[1];
rz(-1.4956359) q[1];
sx q[1];
rz(1.2542226) q[1];
x q[2];
rz(-0.75922482) q[3];
sx q[3];
rz(-0.99731481) q[3];
sx q[3];
rz(-2.218117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1141438) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(2.0360428) q[2];
rz(0.28087273) q[3];
sx q[3];
rz(-0.61087817) q[3];
sx q[3];
rz(0.15225473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0356782) q[0];
sx q[0];
rz(-0.67258251) q[0];
sx q[0];
rz(2.0689082) q[0];
rz(-3.1283123) q[1];
sx q[1];
rz(-1.599879) q[1];
sx q[1];
rz(-0.57659155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.428513) q[0];
sx q[0];
rz(-1.2831076) q[0];
sx q[0];
rz(1.6384533) q[0];
rz(-pi) q[1];
rz(1.0112314) q[2];
sx q[2];
rz(-0.15655993) q[2];
sx q[2];
rz(2.8383534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2999163) q[1];
sx q[1];
rz(-0.51989976) q[1];
sx q[1];
rz(1.9701881) q[1];
x q[2];
rz(-2.1284038) q[3];
sx q[3];
rz(-1.3609145) q[3];
sx q[3];
rz(-0.55279532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85792929) q[2];
sx q[2];
rz(-0.89202213) q[2];
sx q[2];
rz(-2.499495) q[2];
rz(2.5765431) q[3];
sx q[3];
rz(-0.27500209) q[3];
sx q[3];
rz(2.9716085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2791486) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(2.8992262) q[0];
rz(-0.53388059) q[1];
sx q[1];
rz(-0.68478525) q[1];
sx q[1];
rz(-1.6530316) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2676754) q[0];
sx q[0];
rz(-2.5416964) q[0];
sx q[0];
rz(-0.12592648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62612957) q[2];
sx q[2];
rz(-0.51373791) q[2];
sx q[2];
rz(2.0425678) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10746126) q[1];
sx q[1];
rz(-2.3337769) q[1];
sx q[1];
rz(0.40523578) q[1];
x q[2];
rz(3.0128128) q[3];
sx q[3];
rz(-0.30211651) q[3];
sx q[3];
rz(-2.0649892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7786467) q[2];
sx q[2];
rz(-1.8774418) q[2];
sx q[2];
rz(-1.470559) q[2];
rz(3.0342024) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(-1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11112467) q[0];
sx q[0];
rz(-2.680439) q[0];
sx q[0];
rz(-1.9386468) q[0];
rz(0.89490926) q[1];
sx q[1];
rz(-1.1824965) q[1];
sx q[1];
rz(0.21200171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0184191) q[0];
sx q[0];
rz(-2.1668262) q[0];
sx q[0];
rz(-0.53431781) q[0];
x q[1];
rz(-0.52258018) q[2];
sx q[2];
rz(-1.3088639) q[2];
sx q[2];
rz(-0.22522989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4350076) q[1];
sx q[1];
rz(-0.75453484) q[1];
sx q[1];
rz(2.2526787) q[1];
rz(-pi) q[2];
rz(2.1240953) q[3];
sx q[3];
rz(-1.5626255) q[3];
sx q[3];
rz(-3.0752575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2682858) q[2];
sx q[2];
rz(-2.3640859) q[2];
sx q[2];
rz(3.0262465) q[2];
rz(0.98571944) q[3];
sx q[3];
rz(-1.4109572) q[3];
sx q[3];
rz(2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2301521) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(2.4440515) q[0];
rz(2.7279834) q[1];
sx q[1];
rz(-1.6802843) q[1];
sx q[1];
rz(-2.4381309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7005806) q[0];
sx q[0];
rz(-1.8412151) q[0];
sx q[0];
rz(-1.7187814) q[0];
rz(-0.38833924) q[2];
sx q[2];
rz(-1.0000129) q[2];
sx q[2];
rz(2.6663189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20774923) q[1];
sx q[1];
rz(-1.2022756) q[1];
sx q[1];
rz(-0.50362419) q[1];
rz(-pi) q[2];
rz(-1.2715152) q[3];
sx q[3];
rz(-1.8095152) q[3];
sx q[3];
rz(-1.4786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0387663) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(1.212567) q[2];
rz(-0.28193998) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(0.50659242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3459699) q[0];
sx q[0];
rz(-0.75241929) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(0.94547358) q[1];
sx q[1];
rz(-1.5676326) q[1];
sx q[1];
rz(-2.4623154) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51157975) q[0];
sx q[0];
rz(-1.094745) q[0];
sx q[0];
rz(-1.6569754) q[0];
rz(2.7921259) q[2];
sx q[2];
rz(-2.4878056) q[2];
sx q[2];
rz(-0.9330627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5838562) q[1];
sx q[1];
rz(-2.5917555) q[1];
sx q[1];
rz(1.932895) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6369989) q[3];
sx q[3];
rz(-0.28488628) q[3];
sx q[3];
rz(-1.6411622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88963738) q[2];
sx q[2];
rz(-2.1869982) q[2];
sx q[2];
rz(-1.6471222) q[2];
rz(0.55388081) q[3];
sx q[3];
rz(-1.6051822) q[3];
sx q[3];
rz(1.7249829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81903356) q[0];
sx q[0];
rz(-2.366134) q[0];
sx q[0];
rz(2.0490647) q[0];
rz(0.93387261) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(-2.3563103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0781495) q[0];
sx q[0];
rz(-2.5012996) q[0];
sx q[0];
rz(0.75115738) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44194989) q[2];
sx q[2];
rz(-2.0734678) q[2];
sx q[2];
rz(-1.82077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8173302) q[1];
sx q[1];
rz(-1.3249389) q[1];
sx q[1];
rz(-2.861633) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1340895) q[3];
sx q[3];
rz(-0.99815449) q[3];
sx q[3];
rz(2.9331895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87172047) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(-2.5093057) q[2];
rz(-2.2271683) q[3];
sx q[3];
rz(-1.2605896) q[3];
sx q[3];
rz(0.15973346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81894994) q[0];
sx q[0];
rz(-2.2845848) q[0];
sx q[0];
rz(2.6981165) q[0];
rz(-0.30287287) q[1];
sx q[1];
rz(-0.58285204) q[1];
sx q[1];
rz(-1.833896) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46638495) q[0];
sx q[0];
rz(-1.2699915) q[0];
sx q[0];
rz(1.3928901) q[0];
x q[1];
rz(1.2595749) q[2];
sx q[2];
rz(-1.0838795) q[2];
sx q[2];
rz(0.68070179) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.16567) q[1];
sx q[1];
rz(-0.56708401) q[1];
sx q[1];
rz(2.0815064) q[1];
rz(2.0409706) q[3];
sx q[3];
rz(-1.8368145) q[3];
sx q[3];
rz(3.0373552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12094721) q[2];
sx q[2];
rz(-1.0409313) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(-0.81765085) q[3];
sx q[3];
rz(-2.652467) q[3];
sx q[3];
rz(0.31942719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15449512) q[0];
sx q[0];
rz(-1.3767865) q[0];
sx q[0];
rz(-0.35926551) q[0];
rz(-0.66419762) q[1];
sx q[1];
rz(-1.3281053) q[1];
sx q[1];
rz(2.7204303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4818648) q[0];
sx q[0];
rz(-2.7057428) q[0];
sx q[0];
rz(1.0208561) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4975411) q[2];
sx q[2];
rz(-3.0689291) q[2];
sx q[2];
rz(-0.49468985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2639271) q[1];
sx q[1];
rz(-1.3093728) q[1];
sx q[1];
rz(1.1007453) q[1];
rz(-pi) q[2];
rz(-2.5362016) q[3];
sx q[3];
rz(-1.0773813) q[3];
sx q[3];
rz(2.8132954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8757214) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(1.9270012) q[2];
rz(2.752979) q[3];
sx q[3];
rz(-2.2190084) q[3];
sx q[3];
rz(-2.0428366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.494396) q[0];
sx q[0];
rz(-1.5460486) q[0];
sx q[0];
rz(-1.8478951) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(1.1399798) q[2];
sx q[2];
rz(-0.93968334) q[2];
sx q[2];
rz(0.46776007) q[2];
rz(-1.8372336) q[3];
sx q[3];
rz(-2.022701) q[3];
sx q[3];
rz(-1.4823784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
