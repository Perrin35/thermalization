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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(0.32455197) q[0];
rz(2.1895154) q[1];
sx q[1];
rz(-2.9157186) q[1];
sx q[1];
rz(1.8332551) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150779) q[0];
sx q[0];
rz(-0.9008207) q[0];
sx q[0];
rz(-1.7448533) q[0];
rz(-2.8988774) q[2];
sx q[2];
rz(-2.6449892) q[2];
sx q[2];
rz(-0.14641031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1349012) q[1];
sx q[1];
rz(-0.54062245) q[1];
sx q[1];
rz(-2.0903793) q[1];
rz(-pi) q[2];
rz(-0.69921391) q[3];
sx q[3];
rz(-1.4627856) q[3];
sx q[3];
rz(2.8349769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41210458) q[2];
sx q[2];
rz(-1.2594014) q[2];
sx q[2];
rz(2.7903902) q[2];
rz(1.2900194) q[3];
sx q[3];
rz(-2.0100287) q[3];
sx q[3];
rz(1.2235519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(2.2112041) q[0];
rz(1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-2.5321541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7284645) q[0];
sx q[0];
rz(-0.43914686) q[0];
sx q[0];
rz(1.3261611) q[0];
rz(1.8571505) q[2];
sx q[2];
rz(-0.30977466) q[2];
sx q[2];
rz(1.3510493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6516782) q[1];
sx q[1];
rz(-2.6082572) q[1];
sx q[1];
rz(-0.79483219) q[1];
x q[2];
rz(-0.16709631) q[3];
sx q[3];
rz(-1.0067938) q[3];
sx q[3];
rz(0.43844863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(-3.0668104) q[2];
rz(2.6212202) q[3];
sx q[3];
rz(-0.64295355) q[3];
sx q[3];
rz(-2.5274932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-0.30211788) q[0];
rz(-2.865454) q[1];
sx q[1];
rz(-1.8243676) q[1];
sx q[1];
rz(1.4322697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8373568) q[0];
sx q[0];
rz(-0.68528995) q[0];
sx q[0];
rz(2.1987112) q[0];
x q[1];
rz(1.309607) q[2];
sx q[2];
rz(-2.1495594) q[2];
sx q[2];
rz(0.28489339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70967445) q[1];
sx q[1];
rz(-2.1248105) q[1];
sx q[1];
rz(2.9565503) q[1];
rz(-pi) q[2];
rz(1.4301705) q[3];
sx q[3];
rz(-0.46984497) q[3];
sx q[3];
rz(1.8610561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8351195) q[2];
sx q[2];
rz(-2.6758631) q[2];
sx q[2];
rz(-1.845537) q[2];
rz(2.6895788) q[3];
sx q[3];
rz(-1.1072423) q[3];
sx q[3];
rz(2.7590416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52962676) q[0];
sx q[0];
rz(-1.9744248) q[0];
sx q[0];
rz(1.3324598) q[0];
rz(-1.0491071) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(-1.5528991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72399607) q[0];
sx q[0];
rz(-1.9105457) q[0];
sx q[0];
rz(1.8906361) q[0];
rz(-2.0715782) q[2];
sx q[2];
rz(-1.5827279) q[2];
sx q[2];
rz(2.3119761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5415948) q[1];
sx q[1];
rz(-1.3728956) q[1];
sx q[1];
rz(2.2210414) q[1];
rz(-pi) q[2];
rz(-0.96943198) q[3];
sx q[3];
rz(-0.80925377) q[3];
sx q[3];
rz(-1.7084029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4730452) q[2];
sx q[2];
rz(-1.0205525) q[2];
sx q[2];
rz(-0.2505396) q[2];
rz(2.1446832) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(2.7610049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32960358) q[0];
sx q[0];
rz(-0.488509) q[0];
sx q[0];
rz(1.3478152) q[0];
rz(-0.40428058) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(2.0124729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92496678) q[0];
sx q[0];
rz(-0.89168404) q[0];
sx q[0];
rz(-2.62519) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9086228) q[2];
sx q[2];
rz(-2.4475636) q[2];
sx q[2];
rz(2.239486) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68539219) q[1];
sx q[1];
rz(-0.90185748) q[1];
sx q[1];
rz(-3.0078956) q[1];
rz(-pi) q[2];
rz(1.4604767) q[3];
sx q[3];
rz(-1.2577783) q[3];
sx q[3];
rz(0.68389326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79444844) q[2];
sx q[2];
rz(-0.930154) q[2];
sx q[2];
rz(1.244119) q[2];
rz(1.0602903) q[3];
sx q[3];
rz(-0.62842193) q[3];
sx q[3];
rz(2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619549) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(2.7144077) q[0];
rz(3.1071013) q[1];
sx q[1];
rz(-1.5692312) q[1];
sx q[1];
rz(0.010206612) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50960582) q[0];
sx q[0];
rz(-1.5731531) q[0];
sx q[0];
rz(1.5725333) q[0];
rz(-1.353327) q[2];
sx q[2];
rz(-0.64369338) q[2];
sx q[2];
rz(0.74483192) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2040703) q[1];
sx q[1];
rz(-0.88562219) q[1];
sx q[1];
rz(3.0584945) q[1];
rz(-pi) q[2];
rz(-2.0766076) q[3];
sx q[3];
rz(-2.0772604) q[3];
sx q[3];
rz(2.9714366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5222142) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(-2.348032) q[2];
rz(-2.2788952) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(2.4376552) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70235395) q[0];
sx q[0];
rz(-2.2965501) q[0];
sx q[0];
rz(1.1335565) q[0];
rz(2.0639065) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(-2.8394707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7493499) q[0];
sx q[0];
rz(-1.8729405) q[0];
sx q[0];
rz(-2.7878615) q[0];
rz(-pi) q[1];
rz(2.5100915) q[2];
sx q[2];
rz(-0.69984791) q[2];
sx q[2];
rz(-1.2597142) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97000098) q[1];
sx q[1];
rz(-2.0920894) q[1];
sx q[1];
rz(1.820351) q[1];
rz(-0.98432912) q[3];
sx q[3];
rz(-2.1696089) q[3];
sx q[3];
rz(2.5564155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9083531) q[2];
sx q[2];
rz(-2.2519604) q[2];
sx q[2];
rz(-1.7944149) q[2];
rz(1.2498648) q[3];
sx q[3];
rz(-1.1976539) q[3];
sx q[3];
rz(3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5938479) q[0];
sx q[0];
rz(-0.43314728) q[0];
sx q[0];
rz(-3.0331842) q[0];
rz(-2.0938865) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(1.0188867) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8150367) q[0];
sx q[0];
rz(-2.360965) q[0];
sx q[0];
rz(-2.3579979) q[0];
rz(-2.3802929) q[2];
sx q[2];
rz(-1.5660962) q[2];
sx q[2];
rz(2.0705303) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0483886) q[1];
sx q[1];
rz(-2.0282413) q[1];
sx q[1];
rz(-1.1253352) q[1];
rz(-pi) q[2];
rz(-1.4298444) q[3];
sx q[3];
rz(-2.3059855) q[3];
sx q[3];
rz(2.6648389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5041647) q[2];
sx q[2];
rz(-2.1647858) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(-0.49312433) q[3];
sx q[3];
rz(-0.92032856) q[3];
sx q[3];
rz(-1.5776618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2757932) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(-1.5419143) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(-0.42617282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754643) q[0];
sx q[0];
rz(-1.5541561) q[0];
sx q[0];
rz(0.017701935) q[0];
rz(-0.8502281) q[2];
sx q[2];
rz(-2.2022708) q[2];
sx q[2];
rz(-0.15635083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.328769) q[1];
sx q[1];
rz(-1.6881094) q[1];
sx q[1];
rz(-1.3552865) q[1];
rz(-pi) q[2];
rz(-0.79094633) q[3];
sx q[3];
rz(-2.2344347) q[3];
sx q[3];
rz(-2.5928465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.49491945) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(-2.7827175) q[3];
sx q[3];
rz(-1.6152629) q[3];
sx q[3];
rz(-2.0487962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36295715) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(-0.92078513) q[0];
rz(-0.50865632) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(2.7210534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2131491) q[0];
sx q[0];
rz(-0.9125114) q[0];
sx q[0];
rz(0.64853213) q[0];
rz(-1.2447717) q[2];
sx q[2];
rz(-2.4851126) q[2];
sx q[2];
rz(-0.028353779) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1433761) q[1];
sx q[1];
rz(-2.3399379) q[1];
sx q[1];
rz(0.90513913) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.693142) q[3];
sx q[3];
rz(-2.4469355) q[3];
sx q[3];
rz(-0.26315022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4113808) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(2.8249557) q[2];
rz(0.5591048) q[3];
sx q[3];
rz(-0.63566256) q[3];
sx q[3];
rz(-2.3234698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.485514) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(-2.6651233) q[1];
sx q[1];
rz(-1.0404027) q[1];
sx q[1];
rz(-1.7146005) q[1];
rz(1.6638577) q[2];
sx q[2];
rz(-0.75329594) q[2];
sx q[2];
rz(2.6483173) q[2];
rz(0.093802916) q[3];
sx q[3];
rz(-1.1891014) q[3];
sx q[3];
rz(0.56567241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
