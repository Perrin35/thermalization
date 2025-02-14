OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(1.2989651) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(4.4823449) q[1];
sx q[1];
rz(9.525099) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.318905) q[0];
sx q[0];
rz(-1.5189369) q[0];
sx q[0];
rz(2.479631) q[0];
rz(-pi) q[1];
rz(2.5415129) q[2];
sx q[2];
rz(-2.588495) q[2];
sx q[2];
rz(-1.2002522) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8954017) q[1];
sx q[1];
rz(-0.34798393) q[1];
sx q[1];
rz(2.4892508) q[1];
rz(-1.8842949) q[3];
sx q[3];
rz(-2.1380205) q[3];
sx q[3];
rz(2.6424266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4002865) q[2];
sx q[2];
rz(-0.15179673) q[2];
sx q[2];
rz(-0.8068223) q[2];
rz(-0.72871366) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(1.9369283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(-0.30701315) q[0];
rz(1.864805) q[1];
sx q[1];
rz(-2.0025608) q[1];
sx q[1];
rz(1.101864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3392838) q[0];
sx q[0];
rz(-1.3972939) q[0];
sx q[0];
rz(-2.1306778) q[0];
x q[1];
rz(-1.6346143) q[2];
sx q[2];
rz(-1.7498657) q[2];
sx q[2];
rz(3.0916391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9872322) q[1];
sx q[1];
rz(-1.7061966) q[1];
sx q[1];
rz(-1.3920648) q[1];
x q[2];
rz(2.3086433) q[3];
sx q[3];
rz(-2.2112339) q[3];
sx q[3];
rz(-2.4082977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.028194204) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(-0.096435189) q[2];
rz(-2.9891369) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(0.69795394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0518799) q[0];
sx q[0];
rz(-1.1002325) q[0];
sx q[0];
rz(2.7413947) q[0];
rz(1.6290889) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(-2.8834744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5982957) q[0];
sx q[0];
rz(-0.27320293) q[0];
sx q[0];
rz(1.0323204) q[0];
rz(3.0417241) q[2];
sx q[2];
rz(-1.4912919) q[2];
sx q[2];
rz(0.85209633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3052747) q[1];
sx q[1];
rz(-1.5757676) q[1];
sx q[1];
rz(3.1286376) q[1];
rz(-pi) q[2];
rz(1.7430923) q[3];
sx q[3];
rz(-1.3329643) q[3];
sx q[3];
rz(-2.3034277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(-3.1094587) q[2];
rz(-0.26767996) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(1.5816429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716008) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-3.0573523) q[0];
rz(0.035942297) q[1];
sx q[1];
rz(-0.033999559) q[1];
sx q[1];
rz(0.34119225) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81927903) q[0];
sx q[0];
rz(-1.3094762) q[0];
sx q[0];
rz(2.4909291) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.555055) q[2];
sx q[2];
rz(-1.3549202) q[2];
sx q[2];
rz(-3.0160273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1028776) q[1];
sx q[1];
rz(-0.85188085) q[1];
sx q[1];
rz(2.170156) q[1];
rz(-pi) q[2];
rz(1.1197508) q[3];
sx q[3];
rz(-2.1220792) q[3];
sx q[3];
rz(-2.0752843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(-1.6061456) q[2];
rz(-2.8793907) q[3];
sx q[3];
rz(-1.5235135) q[3];
sx q[3];
rz(-0.95482993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3839805) q[0];
sx q[0];
rz(-0.42344991) q[0];
sx q[0];
rz(1.0773995) q[0];
rz(-2.7491838) q[1];
sx q[1];
rz(-0.078374021) q[1];
sx q[1];
rz(2.1108625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1218613) q[0];
sx q[0];
rz(-1.1543589) q[0];
sx q[0];
rz(2.5520578) q[0];
rz(-pi) q[1];
rz(2.0692798) q[2];
sx q[2];
rz(-1.7297812) q[2];
sx q[2];
rz(2.0778401) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8414611) q[1];
sx q[1];
rz(-2.8634954) q[1];
sx q[1];
rz(2.2113731) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0822273) q[3];
sx q[3];
rz(-1.2687195) q[3];
sx q[3];
rz(-1.6405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3173759) q[2];
sx q[2];
rz(-0.65418303) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(0.9134891) q[3];
sx q[3];
rz(-1.8279816) q[3];
sx q[3];
rz(0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4588673) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(-1.6656026) q[0];
rz(0.39235517) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(0.59142339) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77350835) q[0];
sx q[0];
rz(-0.71660935) q[0];
sx q[0];
rz(2.1518097) q[0];
rz(-pi) q[1];
rz(2.5673836) q[2];
sx q[2];
rz(-2.817135) q[2];
sx q[2];
rz(-2.7221219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5352262) q[1];
sx q[1];
rz(-2.6312345) q[1];
sx q[1];
rz(1.4354857) q[1];
rz(2.2675681) q[3];
sx q[3];
rz(-1.2748147) q[3];
sx q[3];
rz(0.10914739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.300294) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(0.57233125) q[2];
rz(2.9193997) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(-0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543024) q[0];
sx q[0];
rz(-0.13986762) q[0];
sx q[0];
rz(2.733316) q[0];
rz(-2.4093742) q[1];
sx q[1];
rz(-0.1258985) q[1];
sx q[1];
rz(2.8439723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0746411) q[0];
sx q[0];
rz(-2.5805001) q[0];
sx q[0];
rz(2.448161) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0760667) q[2];
sx q[2];
rz(-1.552993) q[2];
sx q[2];
rz(2.7731563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0875594) q[1];
sx q[1];
rz(-0.39759025) q[1];
sx q[1];
rz(-1.041574) q[1];
rz(-2.4378889) q[3];
sx q[3];
rz(-1.4559485) q[3];
sx q[3];
rz(-2.5435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1709661) q[2];
sx q[2];
rz(-1.8793224) q[2];
sx q[2];
rz(-0.69620281) q[2];
rz(0.89037406) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625229) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(0.21275511) q[0];
rz(-2.6720324) q[1];
sx q[1];
rz(-0.9649562) q[1];
sx q[1];
rz(0.75417095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5306204) q[0];
sx q[0];
rz(-1.3830796) q[0];
sx q[0];
rz(-1.8646445) q[0];
rz(-pi) q[1];
rz(2.5905072) q[2];
sx q[2];
rz(-2.7646825) q[2];
sx q[2];
rz(2.2152679) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6282916) q[1];
sx q[1];
rz(-2.4835212) q[1];
sx q[1];
rz(-2.9546896) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2175104) q[3];
sx q[3];
rz(-2.1322741) q[3];
sx q[3];
rz(-1.9677066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0065877) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(-2.3593486) q[2];
rz(1.6953281) q[3];
sx q[3];
rz(-0.54636121) q[3];
sx q[3];
rz(2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0924454) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(-2.1771722) q[0];
rz(-1.8655221) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.6395578) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14012155) q[0];
sx q[0];
rz(-1.3706511) q[0];
sx q[0];
rz(-1.3339304) q[0];
rz(0.48351863) q[2];
sx q[2];
rz(-2.3502825) q[2];
sx q[2];
rz(0.4936337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.561475) q[1];
sx q[1];
rz(-2.0881542) q[1];
sx q[1];
rz(-0.026419445) q[1];
rz(-pi) q[2];
rz(1.1427059) q[3];
sx q[3];
rz(-0.16166281) q[3];
sx q[3];
rz(-0.77199329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2532578) q[2];
sx q[2];
rz(-1.8586681) q[2];
sx q[2];
rz(2.1962732) q[2];
rz(0.82593289) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0655521) q[0];
sx q[0];
rz(-1.9726418) q[0];
sx q[0];
rz(0.8031351) q[0];
rz(1.5777292) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(2.8520083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11293791) q[0];
sx q[0];
rz(-1.5585596) q[0];
sx q[0];
rz(1.6763681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31676964) q[2];
sx q[2];
rz(-1.2578336) q[2];
sx q[2];
rz(-2.6027158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3037422) q[1];
sx q[1];
rz(-1.9281862) q[1];
sx q[1];
rz(2.9010495) q[1];
x q[2];
rz(1.2373459) q[3];
sx q[3];
rz(-1.2996009) q[3];
sx q[3];
rz(1.1734133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3760066) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(-0.96013367) q[2];
rz(0.58297408) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(-2.2768903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4509907) q[0];
sx q[0];
rz(-1.4324181) q[0];
sx q[0];
rz(1.6313534) q[0];
rz(-0.040738978) q[1];
sx q[1];
rz(-2.4650885) q[1];
sx q[1];
rz(-3.0104641) q[1];
rz(-2.1292674) q[2];
sx q[2];
rz(-0.53999117) q[2];
sx q[2];
rz(-2.7593031) q[2];
rz(2.1914239) q[3];
sx q[3];
rz(-0.087407268) q[3];
sx q[3];
rz(2.0711318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
