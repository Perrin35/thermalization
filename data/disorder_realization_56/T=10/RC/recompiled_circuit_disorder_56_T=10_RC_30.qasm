OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(-0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(2.1138432) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42067895) q[0];
sx q[0];
rz(-1.5936216) q[0];
sx q[0];
rz(0.37585293) q[0];
rz(-1.1711636) q[2];
sx q[2];
rz(-2.6174449) q[2];
sx q[2];
rz(-1.580796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7386908) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(0.72523592) q[1];
rz(2.5146033) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(-2.2771319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1296922) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(0.29775277) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(-2.0334977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28716921) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(-1.9186583) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1531395) q[2];
sx q[2];
rz(-0.69500143) q[2];
sx q[2];
rz(-2.7671791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83924343) q[1];
sx q[1];
rz(-1.1855372) q[1];
sx q[1];
rz(1.1883931) q[1];
x q[2];
rz(-1.8403345) q[3];
sx q[3];
rz(-0.88629913) q[3];
sx q[3];
rz(-2.136363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(-2.1767445) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7643825) q[0];
sx q[0];
rz(-1.2493734) q[0];
sx q[0];
rz(1.2755323) q[0];
x q[1];
rz(-0.13871128) q[2];
sx q[2];
rz(-0.29964667) q[2];
sx q[2];
rz(1.2219929) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.112903) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.4494004) q[1];
rz(1.996702) q[3];
sx q[3];
rz(-0.45255462) q[3];
sx q[3];
rz(-1.9354613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.09459153) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(-2.0084521) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(-0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621181) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(-1.3765155) q[0];
rz(0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(2.8994697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0644558) q[0];
sx q[0];
rz(-0.69364871) q[0];
sx q[0];
rz(-2.8855188) q[0];
rz(1.5577199) q[2];
sx q[2];
rz(-1.1411975) q[2];
sx q[2];
rz(2.122288) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33430797) q[1];
sx q[1];
rz(-2.0172999) q[1];
sx q[1];
rz(1.7590894) q[1];
x q[2];
rz(1.4297156) q[3];
sx q[3];
rz(-1.9851079) q[3];
sx q[3];
rz(2.8431818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(-1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(1.1345908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15123385) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(-2.4609341) q[0];
rz(-0.35933944) q[2];
sx q[2];
rz(-1.1655072) q[2];
sx q[2];
rz(-2.5422424) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6280917) q[1];
sx q[1];
rz(-2.6673632) q[1];
sx q[1];
rz(-2.7952816) q[1];
x q[2];
rz(-2.9633425) q[3];
sx q[3];
rz(-0.44294391) q[3];
sx q[3];
rz(-2.6435542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(-1.9893601) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1359065) q[0];
sx q[0];
rz(-1.583477) q[0];
sx q[0];
rz(-1.4655051) q[0];
x q[1];
rz(-1.9828898) q[2];
sx q[2];
rz(-1.1043784) q[2];
sx q[2];
rz(-2.7127624) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17105477) q[1];
sx q[1];
rz(-1.6754524) q[1];
sx q[1];
rz(2.9403951) q[1];
rz(2.4310962) q[3];
sx q[3];
rz(-2.3264255) q[3];
sx q[3];
rz(-1.9424155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(-0.1858055) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018935238) q[0];
sx q[0];
rz(-1.1766953) q[0];
sx q[0];
rz(1.3144486) q[0];
rz(-pi) q[1];
rz(-0.37624069) q[2];
sx q[2];
rz(-1.2708086) q[2];
sx q[2];
rz(0.30181995) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9638622) q[1];
sx q[1];
rz(-1.1451671) q[1];
sx q[1];
rz(-1.298794) q[1];
rz(-pi) q[2];
rz(0.15345927) q[3];
sx q[3];
rz(-1.6073265) q[3];
sx q[3];
rz(-3.1094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(-0.31204143) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7438695) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(-1.2851108) q[0];
rz(1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.3407019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2368603) q[0];
sx q[0];
rz(-1.245984) q[0];
sx q[0];
rz(-1.9843285) q[0];
rz(-2.4457744) q[2];
sx q[2];
rz(-1.6098445) q[2];
sx q[2];
rz(-2.2829636) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33412877) q[1];
sx q[1];
rz(-2.6263116) q[1];
sx q[1];
rz(0.79828243) q[1];
x q[2];
rz(2.1741381) q[3];
sx q[3];
rz(-0.52166044) q[3];
sx q[3];
rz(-1.1286917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1061219) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0969365) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(0.83713371) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1820113) q[0];
sx q[0];
rz(-2.2416229) q[0];
sx q[0];
rz(-1.8295893) q[0];
rz(-2.7021072) q[2];
sx q[2];
rz(-2.5186933) q[2];
sx q[2];
rz(-2.5650131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7652119) q[1];
sx q[1];
rz(-1.7227255) q[1];
sx q[1];
rz(-3.0357749) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8505441) q[3];
sx q[3];
rz(-2.0488727) q[3];
sx q[3];
rz(0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(-2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76374617) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996795) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(0.051254) q[0];
rz(-0.51771848) q[2];
sx q[2];
rz(-2.0132408) q[2];
sx q[2];
rz(2.8874318) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53303888) q[1];
sx q[1];
rz(-1.3320919) q[1];
sx q[1];
rz(-1.1281611) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2782709) q[3];
sx q[3];
rz(-2.8046126) q[3];
sx q[3];
rz(1.8388336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-2.3095619) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(-0.17765799) q[2];
sx q[2];
rz(-1.0189609) q[2];
sx q[2];
rz(-0.56448274) q[2];
rz(1.1888614) q[3];
sx q[3];
rz(-1.100913) q[3];
sx q[3];
rz(2.037896) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
