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
rz(-2.8170407) q[0];
rz(-0.95207721) q[1];
sx q[1];
rz(-0.2258741) q[1];
sx q[1];
rz(1.3083375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4885534) q[0];
sx q[0];
rz(-1.4346315) q[0];
sx q[0];
rz(-0.67739886) q[0];
rz(-pi) q[1];
rz(-2.657349) q[2];
sx q[2];
rz(-1.6855557) q[2];
sx q[2];
rz(1.6387303) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5950299) q[1];
sx q[1];
rz(-1.1076704) q[1];
sx q[1];
rz(-0.28966499) q[1];
x q[2];
rz(0.16690688) q[3];
sx q[3];
rz(-0.70611533) q[3];
sx q[3];
rz(-1.1366858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7294881) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(0.35120249) q[2];
rz(1.8515733) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(-1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38607645) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(-2.2112041) q[0];
rz(1.8366086) q[1];
sx q[1];
rz(-2.3001859) q[1];
sx q[1];
rz(0.60943857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6822426) q[0];
sx q[0];
rz(-1.1455904) q[0];
sx q[0];
rz(-0.11328477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8687042) q[2];
sx q[2];
rz(-1.6570083) q[2];
sx q[2];
rz(-2.6484368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7623065) q[1];
sx q[1];
rz(-1.2067144) q[1];
sx q[1];
rz(-1.9696139) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16709631) q[3];
sx q[3];
rz(-1.0067938) q[3];
sx q[3];
rz(-0.43844863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(-2.6212202) q[3];
sx q[3];
rz(-0.64295355) q[3];
sx q[3];
rz(2.5274932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-0.30211788) q[0];
rz(2.865454) q[1];
sx q[1];
rz(-1.3172251) q[1];
sx q[1];
rz(1.4322697) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30423588) q[0];
sx q[0];
rz(-2.4563027) q[0];
sx q[0];
rz(0.94288148) q[0];
rz(-1.8319857) q[2];
sx q[2];
rz(-2.1495594) q[2];
sx q[2];
rz(0.28489339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.773743) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(1.2817205) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4301705) q[3];
sx q[3];
rz(-0.46984497) q[3];
sx q[3];
rz(1.8610561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3064731) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(-1.845537) q[2];
rz(2.6895788) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(-2.7590416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6119659) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(1.3324598) q[0];
rz(1.0491071) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(-1.5886935) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5063254) q[0];
sx q[0];
rz(-0.46231368) q[0];
sx q[0];
rz(-0.72700951) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1279911) q[2];
sx q[2];
rz(-1.0700534) q[2];
sx q[2];
rz(-2.406943) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2237593) q[1];
sx q[1];
rz(-2.4661015) q[1];
sx q[1];
rz(-1.89066) q[1];
rz(-pi) q[2];
rz(-0.96943198) q[3];
sx q[3];
rz(-0.80925377) q[3];
sx q[3];
rz(1.4331897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4730452) q[2];
sx q[2];
rz(-1.0205525) q[2];
sx q[2];
rz(2.891053) q[2];
rz(-0.9969095) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(-0.38058773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8119891) q[0];
sx q[0];
rz(-0.488509) q[0];
sx q[0];
rz(1.3478152) q[0];
rz(-2.7373121) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(-2.0124729) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92496678) q[0];
sx q[0];
rz(-2.2499086) q[0];
sx q[0];
rz(-2.62519) q[0];
rz(-1.9086228) q[2];
sx q[2];
rz(-2.4475636) q[2];
sx q[2];
rz(-0.90210669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89896527) q[1];
sx q[1];
rz(-2.4614442) q[1];
sx q[1];
rz(1.7378452) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4604767) q[3];
sx q[3];
rz(-1.2577783) q[3];
sx q[3];
rz(-2.4576994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3471442) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(1.244119) q[2];
rz(-2.0813023) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(0.30317831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4853972) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(0.42718497) q[0];
rz(-3.1071013) q[1];
sx q[1];
rz(-1.5692312) q[1];
sx q[1];
rz(-0.010206612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1447574) q[0];
sx q[0];
rz(-0.0029276927) q[0];
sx q[0];
rz(-2.5064431) q[0];
x q[1];
rz(1.353327) q[2];
sx q[2];
rz(-2.4978993) q[2];
sx q[2];
rz(0.74483192) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.806655) q[1];
sx q[1];
rz(-2.4522071) q[1];
sx q[1];
rz(1.4695808) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4236492) q[3];
sx q[3];
rz(-2.4417851) q[3];
sx q[3];
rz(1.0216624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5222142) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(-2.348032) q[2];
rz(2.2788952) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(0.70393744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4392387) q[0];
sx q[0];
rz(-2.2965501) q[0];
sx q[0];
rz(-2.0080361) q[0];
rz(2.0639065) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(-2.8394707) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0724807) q[0];
sx q[0];
rz(-1.9078507) q[0];
sx q[0];
rz(1.2500136) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0321192) q[2];
sx q[2];
rz(-1.0240842) q[2];
sx q[2];
rz(-0.49671587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6985642) q[1];
sx q[1];
rz(-2.5686567) q[1];
sx q[1];
rz(0.40614508) q[1];
x q[2];
rz(0.68639836) q[3];
sx q[3];
rz(-1.0960311) q[3];
sx q[3];
rz(-1.7976185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9083531) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(-1.7944149) q[2];
rz(1.2498648) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(0.10716042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5477448) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(-3.0331842) q[0];
rz(-1.0477061) q[1];
sx q[1];
rz(-2.3240418) q[1];
sx q[1];
rz(-2.122706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86025809) q[0];
sx q[0];
rz(-2.0906013) q[0];
sx q[0];
rz(0.61183527) q[0];
rz(-pi) q[1];
rz(2.3802929) q[2];
sx q[2];
rz(-1.5754964) q[2];
sx q[2];
rz(-1.0710623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68542504) q[1];
sx q[1];
rz(-1.1738832) q[1];
sx q[1];
rz(-2.6422068) q[1];
rz(-pi) q[2];
rz(2.4014428) q[3];
sx q[3];
rz(-1.6751846) q[3];
sx q[3];
rz(1.1889282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5041647) q[2];
sx q[2];
rz(-2.1647858) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(-0.49312433) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(-1.5639308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86579943) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(-2.8637874) q[0];
rz(1.5419143) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(0.42617282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3329553) q[0];
sx q[0];
rz(-1.5884958) q[0];
sx q[0];
rz(1.5541535) q[0];
rz(-pi) q[1];
rz(-2.2913646) q[2];
sx q[2];
rz(-2.2022708) q[2];
sx q[2];
rz(-2.9852418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.26697576) q[1];
sx q[1];
rz(-0.24493453) q[1];
sx q[1];
rz(1.0670948) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83287333) q[3];
sx q[3];
rz(-0.98372059) q[3];
sx q[3];
rz(0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6466732) q[2];
sx q[2];
rz(-2.9743331) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(-0.35887512) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36295715) q[0];
sx q[0];
rz(-2.9243922) q[0];
sx q[0];
rz(-0.92078513) q[0];
rz(-0.50865632) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(2.7210534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32146183) q[0];
sx q[0];
rz(-0.88867868) q[0];
sx q[0];
rz(0.90773037) q[0];
x q[1];
rz(1.2447717) q[2];
sx q[2];
rz(-0.65648001) q[2];
sx q[2];
rz(-0.028353779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1433761) q[1];
sx q[1];
rz(-2.3399379) q[1];
sx q[1];
rz(2.2364535) q[1];
rz(3.0402571) q[3];
sx q[3];
rz(-2.2592415) q[3];
sx q[3];
rz(-2.7197414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4113808) q[2];
sx q[2];
rz(-0.79019848) q[2];
sx q[2];
rz(-2.8249557) q[2];
rz(-0.5591048) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(-2.3234698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.086924788) q[2];
sx q[2];
rz(-0.82155052) q[2];
sx q[2];
rz(-0.62053298) q[2];
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
