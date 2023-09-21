OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8310175) q[0];
sx q[0];
rz(-2.0895045) q[0];
sx q[0];
rz(-1.6488099) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(2.066943) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91905347) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(0.28199621) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4104112) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(0.10057848) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5684811) q[1];
sx q[1];
rz(-2.9955578) q[1];
sx q[1];
rz(-2.5558429) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0243458) q[3];
sx q[3];
rz(-1.082886) q[3];
sx q[3];
rz(2.818056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(1.6281698) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.2233541) q[0];
rz(0.16076316) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.5011903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6340856) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(-3.0188942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1378527) q[2];
sx q[2];
rz(-0.54076414) q[2];
sx q[2];
rz(-0.31976779) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7565631) q[1];
sx q[1];
rz(-2.2379413) q[1];
sx q[1];
rz(1.3515616) q[1];
rz(-pi) q[2];
rz(3.1355751) q[3];
sx q[3];
rz(-1.9618417) q[3];
sx q[3];
rz(1.0139549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.5037781) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2730946) q[0];
sx q[0];
rz(-1.4833741) q[0];
sx q[0];
rz(-3.0762061) q[0];
x q[1];
rz(-2.2321646) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(0.94240377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40239247) q[1];
sx q[1];
rz(-1.8715197) q[1];
sx q[1];
rz(-1.5261569) q[1];
rz(-pi) q[2];
rz(-3.010473) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.1733615) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(1.051735) q[0];
rz(-0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(-1.1490885) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3455032) q[0];
sx q[0];
rz(-2.8755099) q[0];
sx q[0];
rz(-2.7059113) q[0];
rz(-pi) q[1];
rz(2.1707702) q[2];
sx q[2];
rz(-0.68471013) q[2];
sx q[2];
rz(-0.1240571) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9149575) q[1];
sx q[1];
rz(-1.8640222) q[1];
sx q[1];
rz(0.83421591) q[1];
rz(1.9150919) q[3];
sx q[3];
rz(-0.85840423) q[3];
sx q[3];
rz(2.815849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(-1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-0.88422424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54908862) q[0];
sx q[0];
rz(-1.3929402) q[0];
sx q[0];
rz(1.3772208) q[0];
rz(-pi) q[1];
rz(-2.7344633) q[2];
sx q[2];
rz(-1.1659157) q[2];
sx q[2];
rz(-1.7083573) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3430749) q[1];
sx q[1];
rz(-1.6371562) q[1];
sx q[1];
rz(-0.14454486) q[1];
rz(-2.1518425) q[3];
sx q[3];
rz(-2.5269066) q[3];
sx q[3];
rz(0.97782545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4051751) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58020619) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(-1.6798621) q[0];
rz(-pi) q[1];
rz(0.078209608) q[2];
sx q[2];
rz(-2.4259896) q[2];
sx q[2];
rz(2.2323334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5507257) q[1];
sx q[1];
rz(-2.496965) q[1];
sx q[1];
rz(-2.4025737) q[1];
x q[2];
rz(-8*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(-1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51450729) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095734) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.3597885) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.70599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5933701) q[0];
sx q[0];
rz(-0.75836327) q[0];
sx q[0];
rz(2.3286657) q[0];
rz(-1.5905459) q[2];
sx q[2];
rz(-1.5891468) q[2];
sx q[2];
rz(-1.8576647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3064733) q[1];
sx q[1];
rz(-1.2700348) q[1];
sx q[1];
rz(1.9496099) q[1];
x q[2];
rz(-2.2613951) q[3];
sx q[3];
rz(-2.5185563) q[3];
sx q[3];
rz(-2.8043889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(-1.0081572) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(0.15667668) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(-2.6538387) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93850219) q[0];
sx q[0];
rz(-1.4346968) q[0];
sx q[0];
rz(0.076119856) q[0];
rz(-pi) q[1];
rz(1.6985967) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(-0.70805659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8147162) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(0.53175064) q[1];
x q[2];
rz(2.7435999) q[3];
sx q[3];
rz(-2.6905305) q[3];
sx q[3];
rz(-0.034702452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(-2.7189642) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(-0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(-0.92393595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7565219) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(-2.4247652) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48859476) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(-2.7985364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82211557) q[1];
sx q[1];
rz(-1.176782) q[1];
sx q[1];
rz(1.0072717) q[1];
rz(-pi) q[2];
rz(0.68848916) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(-1.1779932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.7473934) q[2];
rz(-1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(3.1034234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105135) q[0];
sx q[0];
rz(-1.5974701) q[0];
sx q[0];
rz(-2.1795576) q[0];
rz(-pi) q[1];
rz(-1.1881234) q[2];
sx q[2];
rz(-1.4854991) q[2];
sx q[2];
rz(-1.6558937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0109509) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(-1.6846659) q[1];
rz(-0.60572259) q[3];
sx q[3];
rz(-2.1401792) q[3];
sx q[3];
rz(0.13525087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.9434628) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(0.11001982) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
rz(-0.31717832) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];