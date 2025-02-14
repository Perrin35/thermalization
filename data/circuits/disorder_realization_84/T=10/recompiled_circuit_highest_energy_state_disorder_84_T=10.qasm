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
rz(0.52047210931778) q[0];
sx q[0];
rz(4.17767241795594) q[0];
sx q[0];
rz(11.9254803419034) q[0];
rz(-0.237352922558784) q[1];
sx q[1];
rz(3.78323981364305) q[1];
sx q[1];
rz(6.19828150271579) q[1];
cx q[1],q[0];
rz(-1.90834021568298) q[0];
sx q[0];
rz(4.80968609650666) q[0];
sx q[0];
rz(10.1301090478818) q[0];
rz(5.20214319229126) q[2];
sx q[2];
rz(3.20784317900474) q[2];
sx q[2];
rz(10.2150172948758) q[2];
cx q[2],q[1];
rz(-0.816048741340637) q[1];
sx q[1];
rz(3.4545688350969) q[1];
sx q[1];
rz(11.3744982242505) q[1];
rz(-0.555198431015015) q[3];
sx q[3];
rz(5.61484304268891) q[3];
sx q[3];
rz(13.2461490392606) q[3];
cx q[3],q[2];
rz(0.251959651708603) q[2];
sx q[2];
rz(4.34428647358949) q[2];
sx q[2];
rz(11.0579654931943) q[2];
rz(2.30973219871521) q[3];
sx q[3];
rz(5.95714035828645) q[3];
sx q[3];
rz(10.7693505048673) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.21075963973999) q[0];
sx q[0];
rz(3.77686390479142) q[0];
sx q[0];
rz(6.0181762933652) q[0];
rz(2.34590721130371) q[1];
sx q[1];
rz(0.345147522287913) q[1];
sx q[1];
rz(8.00365374087497) q[1];
cx q[1],q[0];
rz(4.59992218017578) q[0];
sx q[0];
rz(3.84911814530427) q[0];
sx q[0];
rz(9.83825874923869) q[0];
rz(0.616837978363037) q[2];
sx q[2];
rz(4.32132378418977) q[2];
sx q[2];
rz(13.3727195024411) q[2];
cx q[2],q[1];
rz(2.98275232315063) q[1];
sx q[1];
rz(4.66877916653688) q[1];
sx q[1];
rz(9.67662886380359) q[1];
rz(0.175403401255608) q[3];
sx q[3];
rz(9.59499612649018) q[3];
sx q[3];
rz(8.29797706603214) q[3];
cx q[3],q[2];
rz(5.51100206375122) q[2];
sx q[2];
rz(4.13676092227037) q[2];
sx q[2];
rz(8.76560554503604) q[2];
rz(-2.21994209289551) q[3];
sx q[3];
rz(2.09893360932405) q[3];
sx q[3];
rz(7.8422034740369) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.71218609809875) q[0];
sx q[0];
rz(0.615391882258006) q[0];
sx q[0];
rz(8.14129540919467) q[0];
rz(3.5684928894043) q[1];
sx q[1];
rz(-0.593952981633596) q[1];
sx q[1];
rz(7.45073602198764) q[1];
cx q[1],q[0];
rz(2.14388632774353) q[0];
sx q[0];
rz(0.501811178522654) q[0];
sx q[0];
rz(10.1643995404164) q[0];
rz(-3.64655351638794) q[2];
sx q[2];
rz(5.1866034587198) q[2];
sx q[2];
rz(8.2780492067258) q[2];
cx q[2],q[1];
rz(-1.1192272901535) q[1];
sx q[1];
rz(4.61834731896455) q[1];
sx q[1];
rz(10.7047623157422) q[1];
rz(-1.83774411678314) q[3];
sx q[3];
rz(8.02398935158784) q[3];
sx q[3];
rz(8.50304273366138) q[3];
cx q[3],q[2];
rz(-0.20649641752243) q[2];
sx q[2];
rz(4.79721370537812) q[2];
sx q[2];
rz(9.30698072015449) q[2];
rz(-1.30550336837769) q[3];
sx q[3];
rz(5.42362299759919) q[3];
sx q[3];
rz(8.43390551804706) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(5.97855949401855) q[0];
sx q[0];
rz(7.33494821389253) q[0];
sx q[0];
rz(9.45730970650121) q[0];
rz(0.983719825744629) q[1];
sx q[1];
rz(2.33012208540971) q[1];
sx q[1];
rz(9.96766910552188) q[1];
cx q[1],q[0];
rz(0.98292338848114) q[0];
sx q[0];
rz(0.741159351664134) q[0];
sx q[0];
rz(10.5851902723233) q[0];
rz(-4.81562185287476) q[2];
sx q[2];
rz(4.40156796773011) q[2];
sx q[2];
rz(9.18801205455467) q[2];
cx q[2],q[1];
rz(-3.02631521224976) q[1];
sx q[1];
rz(4.4584780057245) q[1];
sx q[1];
rz(10.1890784859578) q[1];
rz(-0.164435848593712) q[3];
sx q[3];
rz(7.02845135529573) q[3];
sx q[3];
rz(8.28769061564609) q[3];
cx q[3],q[2];
rz(4.72135353088379) q[2];
sx q[2];
rz(4.27158811886842) q[2];
sx q[2];
rz(5.32258365153476) q[2];
rz(-4.11710977554321) q[3];
sx q[3];
rz(5.41169348557527) q[3];
sx q[3];
rz(6.80517385005161) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.31127107143402) q[0];
sx q[0];
rz(5.6011268218332) q[0];
sx q[0];
rz(9.75731605886623) q[0];
rz(2.50348162651062) q[1];
sx q[1];
rz(5.65758410294587) q[1];
sx q[1];
rz(13.0225081205289) q[1];
cx q[1],q[0];
rz(3.92118692398071) q[0];
sx q[0];
rz(4.64962318738038) q[0];
sx q[0];
rz(6.1507241487424) q[0];
rz(2.66253709793091) q[2];
sx q[2];
rz(3.98390397627885) q[2];
sx q[2];
rz(6.23906705378696) q[2];
cx q[2],q[1];
rz(-1.78624534606934) q[1];
sx q[1];
rz(7.94305625756318) q[1];
sx q[1];
rz(8.84013876914188) q[1];
rz(-2.4941520690918) q[3];
sx q[3];
rz(11.2891766150766) q[3];
sx q[3];
rz(7.64094636439487) q[3];
cx q[3],q[2];
rz(-2.07863235473633) q[2];
sx q[2];
rz(2.380206437903) q[2];
sx q[2];
rz(8.97258461116954) q[2];
rz(0.264429867267609) q[3];
sx q[3];
rz(6.13703647454316) q[3];
sx q[3];
rz(10.6615860223691) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.07531571388245) q[0];
sx q[0];
rz(0.831065805750438) q[0];
sx q[0];
rz(8.70691475867435) q[0];
rz(4.73049163818359) q[1];
sx q[1];
rz(7.82387128670747) q[1];
sx q[1];
rz(9.63203447162315) q[1];
cx q[1],q[0];
rz(1.85331809520721) q[0];
sx q[0];
rz(2.69876000483567) q[0];
sx q[0];
rz(14.0733299016873) q[0];
rz(0.422803282737732) q[2];
sx q[2];
rz(7.91885057290132) q[2];
sx q[2];
rz(4.14655492304965) q[2];
cx q[2],q[1];
rz(7.3476710319519) q[1];
sx q[1];
rz(4.23587635357911) q[1];
sx q[1];
rz(11.4795477151792) q[1];
rz(3.40257740020752) q[3];
sx q[3];
rz(4.6063530762964) q[3];
sx q[3];
rz(10.4124452233235) q[3];
cx q[3],q[2];
rz(2.38057160377502) q[2];
sx q[2];
rz(4.19686940510804) q[2];
sx q[2];
rz(5.72201750277683) q[2];
rz(1.05939722061157) q[3];
sx q[3];
rz(4.2582335789972) q[3];
sx q[3];
rz(11.1082750320356) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.21374320983887) q[0];
sx q[0];
rz(2.17084673245484) q[0];
sx q[0];
rz(11.719555592529) q[0];
rz(-2.15325856208801) q[1];
sx q[1];
rz(1.76539316971833) q[1];
sx q[1];
rz(10.257676935188) q[1];
cx q[1],q[0];
rz(3.55254745483398) q[0];
sx q[0];
rz(6.05910030205781) q[0];
sx q[0];
rz(8.33814487456485) q[0];
rz(1.06527352333069) q[2];
sx q[2];
rz(3.30370697577531) q[2];
sx q[2];
rz(5.86866543292209) q[2];
cx q[2],q[1];
rz(-0.786452412605286) q[1];
sx q[1];
rz(5.08140245278413) q[1];
sx q[1];
rz(11.9205353021543) q[1];
rz(1.61991691589355) q[3];
sx q[3];
rz(4.82590809662873) q[3];
sx q[3];
rz(12.6035501718442) q[3];
cx q[3],q[2];
rz(2.67584681510925) q[2];
sx q[2];
rz(3.84444937308366) q[2];
sx q[2];
rz(7.05327532290622) q[2];
rz(-0.13104023039341) q[3];
sx q[3];
rz(1.48210564454133) q[3];
sx q[3];
rz(10.6716942548673) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.347626209259033) q[0];
sx q[0];
rz(2.20632311900193) q[0];
sx q[0];
rz(3.6332144498746) q[0];
rz(2.82880163192749) q[1];
sx q[1];
rz(2.85206949909265) q[1];
sx q[1];
rz(9.34237836896583) q[1];
cx q[1],q[0];
rz(-4.33699750900269) q[0];
sx q[0];
rz(3.2320363466912) q[0];
sx q[0];
rz(9.56611480414077) q[0];
rz(-4.14561700820923) q[2];
sx q[2];
rz(5.32679620583589) q[2];
sx q[2];
rz(5.63592360018894) q[2];
cx q[2],q[1];
rz(5.71620798110962) q[1];
sx q[1];
rz(-1.73213705222075) q[1];
sx q[1];
rz(12.7433247327726) q[1];
rz(1.67284452915192) q[3];
sx q[3];
rz(4.40133610566194) q[3];
sx q[3];
rz(8.21155211924716) q[3];
cx q[3],q[2];
rz(4.58564186096191) q[2];
sx q[2];
rz(5.07515087922151) q[2];
sx q[2];
rz(9.40676508321568) q[2];
rz(-0.87120121717453) q[3];
sx q[3];
rz(2.80499911506707) q[3];
sx q[3];
rz(13.2339322328488) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.19106006622314) q[0];
sx q[0];
rz(3.69611534674699) q[0];
sx q[0];
rz(9.88925830125018) q[0];
rz(2.65060997009277) q[1];
sx q[1];
rz(-1.669848767919) q[1];
sx q[1];
rz(4.81571910380527) q[1];
cx q[1],q[0];
rz(5.75020885467529) q[0];
sx q[0];
rz(-2.51668438116973) q[0];
sx q[0];
rz(10.9956289291303) q[0];
rz(-0.305683165788651) q[2];
sx q[2];
rz(4.22623554070527) q[2];
sx q[2];
rz(6.72355411051914) q[2];
cx q[2],q[1];
rz(2.55821442604065) q[1];
sx q[1];
rz(1.51028564770753) q[1];
sx q[1];
rz(12.317948794357) q[1];
rz(5.58635902404785) q[3];
sx q[3];
rz(7.22369638283784) q[3];
sx q[3];
rz(11.4057458400647) q[3];
cx q[3],q[2];
rz(-3.15624284744263) q[2];
sx q[2];
rz(7.77458206017549) q[2];
sx q[2];
rz(17.5181817769925) q[2];
rz(0.939754486083984) q[3];
sx q[3];
rz(2.08461144764955) q[3];
sx q[3];
rz(10.6120707750241) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.23476552963257) q[0];
sx q[0];
rz(2.65637877781922) q[0];
sx q[0];
rz(5.2454724073331) q[0];
rz(4.19597196578979) q[1];
sx q[1];
rz(4.96759322484071) q[1];
sx q[1];
rz(10.474277830116) q[1];
cx q[1],q[0];
rz(0.188055500388145) q[0];
sx q[0];
rz(-0.0497272888845721) q[0];
sx q[0];
rz(12.570059990875) q[0];
rz(-2.13121891021729) q[2];
sx q[2];
rz(5.13514474232728) q[2];
sx q[2];
rz(7.26809237002536) q[2];
cx q[2],q[1];
rz(4.45815849304199) q[1];
sx q[1];
rz(2.53996429045732) q[1];
sx q[1];
rz(12.511773800842) q[1];
rz(-0.783202707767487) q[3];
sx q[3];
rz(3.90011313756044) q[3];
sx q[3];
rz(8.43057814835712) q[3];
cx q[3],q[2];
rz(4.13186836242676) q[2];
sx q[2];
rz(7.59025779564912) q[2];
sx q[2];
rz(10.8711509466092) q[2];
rz(-0.158987864851952) q[3];
sx q[3];
rz(5.12445274193818) q[3];
sx q[3];
rz(11.8163542509) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.04793000221252) q[0];
sx q[0];
rz(2.20853737195069) q[0];
sx q[0];
rz(8.65439871548816) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.75353074073792) q[1];
sx q[1];
rz(4.63794842560823) q[1];
sx q[1];
rz(7.9807281255643) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.407197147607803) q[2];
sx q[2];
rz(1.25158563454682) q[2];
sx q[2];
rz(11.9162664174955) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.2354748249054) q[3];
sx q[3];
rz(3.64759770234162) q[3];
sx q[3];
rz(10.7312707662503) q[3];
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
